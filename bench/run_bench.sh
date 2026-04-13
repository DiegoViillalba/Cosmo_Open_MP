#!/usr/bin/env bash
# =============================================================================
# run_bench.sh
# Benchmarking automático del simulador PM.
#
# Corre pm_serial y pm_parallel con 1, 2, 4, 8 y 16 hilos, extrae los
# tiempos por etapa desde stdout y genera results.csv.
#
# Uso:
#   cd build
#   bash ../bench/run_bench.sh [n_steps] [output_csv]
#
# Argumentos:
#   n_steps    : pasos de simulación (default: 20, suficiente para medir)
#   output_csv : archivo de salida (default: ../bench/results.csv)
#
# Requisito: pm_serial y pm_parallel compilados en el directorio actual (build/).
# =============================================================================

set -euo pipefail

N_STEPS=${1:-20}
OUTPUT_CSV=${2:-../bench/results.csv}
IC_FLAGS="-g"     # generar IC sintéticas, sin archivo externo

# Hilos a probar
THREAD_COUNTS=(1 2 4 8 16)

# Detectar número máximo de hilos disponibles
MAX_THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)
echo "Núcleos detectados: $MAX_THREADS"

# Filtrar thread counts que no excedan los disponibles
VALID_THREADS=()
for t in "${THREAD_COUNTS[@]}"; do
    if [ "$t" -le "$MAX_THREADS" ]; then
        VALID_THREADS+=("$t")
    fi
done
echo "Hilos a probar: ${VALID_THREADS[*]}"
echo ""

# Inicializar CSV
echo "mode,threads,total_s,cic_s,fft_s,grad_s,interp_s,integr_s" > "$OUTPUT_CSV"

# ─────────────────────────────────────────────────────────────────────────────
# Función para extraer tiempo de una etapa del output de la simulación
# Busca la línea que contiene el nombre de la etapa en la tabla de timers.
# ─────────────────────────────────────────────────────────────────────────────
extract_time() {
    local log="$1"
    local stage="$2"
    grep "$stage" "$log" | grep -oP '\d+\.\d+(?= s)' | head -1 || echo "0"
}

# ─────────────────────────────────────────────────────────────────────────────
# Run un ejecutable y extrae métricas
# ─────────────────────────────────────────────────────────────────────────────
run_and_record() {
    local exe="$1"
    local mode="$2"
    local threads="$3"
    local logfile="/tmp/pm_bench_${mode}_t${threads}.log"

    echo -n "  Corriendo $mode con $threads hilo(s)... "

    # Medir tiempo de pared con 'time' del shell
    local t_start=$(date +%s%N)

    ./"$exe" $IC_FLAGS -t "$threads" -s "$N_STEPS" -o /tmp/pm_bench_out \
        > "$logfile" 2>&1

    local t_end=$(date +%s%N)
    local total_ms=$(( (t_end - t_start) / 1000000 ))
    local total_s=$(echo "scale=3; $total_ms / 1000" | bc)

    echo "${total_s}s"

    # Extraer tiempos por etapa desde la tabla de timers en el log
    local cic=$(extract_time "$logfile" "CIC deposit")
    local fft=$(extract_time "$logfile" "Poisson FFT")
    local grad=$(extract_time "$logfile" "Gradiente")
    local interp=$(extract_time "$logfile" "Force interp")
    local integr=$(extract_time "$logfile" "Integrador")

    # Agregar fila al CSV
    echo "$mode,$threads,$total_s,$cic,$fft,$grad,$interp,$integr" >> "$OUTPUT_CSV"
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Baseline secuencial (pm_serial)
# ─────────────────────────────────────────────────────────────────────────────
if [ -f "./pm_serial" ]; then
    echo "=== Baseline secuencial ==="
    run_and_record "pm_serial" "serial" 1
    echo ""
else
    echo "[WARN] pm_serial no encontrado, saltando baseline"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 2. Versión paralela con distintos números de hilos
# ─────────────────────────────────────────────────────────────────────────────
if [ -f "./pm_parallel" ]; then
    echo "=== Versión paralela ==="
    for t in "${VALID_THREADS[@]}"; do
        run_and_record "pm_parallel" "parallel" "$t"
    done
else
    echo "[ERROR] pm_parallel no encontrado. Compila primero con make."
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# 3. Calcular y mostrar tabla de speedup
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== Resultados guardados en: $OUTPUT_CSV ==="
echo ""

# Leer tiempo serial para calcular speedup
T_SERIAL=$(grep "^serial" "$OUTPUT_CSV" | cut -d',' -f3)

if [ -n "$T_SERIAL" ] && [ "$T_SERIAL" != "0" ]; then
    echo "Tabla de speedup (t_serial = ${T_SERIAL}s):"
    echo ""
    printf "%-12s %-10s %-12s %-12s\n" "modo" "hilos" "tiempo(s)" "speedup"
    printf "%-12s %-10s %-12s %-12s\n" "----" "-----" "---------" "-------"

    while IFS=',' read -r mode threads total rest; do
        [ "$mode" = "mode" ] && continue   # saltar cabecera
        speedup=$(echo "scale=2; $T_SERIAL / $total" | bc 2>/dev/null || echo "N/A")
        printf "%-12s %-10s %-12s %-12s\n" "$mode" "$threads" "$total" "$speedup"
    done < "$OUTPUT_CSV"
fi

echo ""
echo "Para graficar los resultados:"
echo "  python3 ../bench/parse_results.py $OUTPUT_CSV"
