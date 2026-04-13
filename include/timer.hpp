#pragma once

// =============================================================================
// timer.hpp
// Utilidad de medición de tiempos para el análisis de desempeño.
//
// Usa omp_get_wtime() como reloj de pared (wall-clock), que es la función
// recomendada por OpenMP para medir tiempo real de ejecución paralela.
// Esto es correcto: clock() mide tiempo de CPU total (suma de todos los hilos),
// no tiempo transcurrido real, lo que daría speedup = 1 siempre.
//
// StageTimer acumula tiempos por etapa (CIC, FFT, gradiente, etc.) durante
// toda la simulación para reportar el desglose al final.
// =============================================================================

#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

// omp_get_wtime() viene de <omp.h>, pero si compilamos sin OpenMP
// la reemplazamos por clock() / CLOCKS_PER_SEC
#ifdef _OPENMP
  #include <omp.h>
  inline double wall_time() { return omp_get_wtime(); }
#else
  #include <ctime>
  inline double wall_time() {
      return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  }
#endif

namespace pm {

// =============================================================================
// ScopedTimer
// Timer RAII: comienza en construcción, imprime al destruirse.
// Útil para medir bloques de código de forma limpia:
//
//   {
//     ScopedTimer t("CIC deposit");
//     cic_deposit(state);
//   }  // <-- aquí imprime el tiempo
// =============================================================================
class ScopedTimer {
public:
    explicit ScopedTimer(const std::string& label, bool print_on_destroy = true)
        : label_(label), t0_(wall_time()), print_(print_on_destroy) {}

    ~ScopedTimer() {
        if (print_) {
            double elapsed = wall_time() - t0_;
            std::cout << "[timer] " << std::left << std::setw(24) << label_
                      << " : " << std::fixed << std::setprecision(4)
                      << elapsed << " s\n";
        }
    }

    // Devuelve el tiempo transcurrido hasta ahora (sin destruir)
    double elapsed() const { return wall_time() - t0_; }

private:
    std::string label_;
    double t0_;
    bool print_;
};

// =============================================================================
// StageTimer
// Acumula tiempos por etapa a lo largo de múltiples pasos de simulación.
// Al final genera un reporte de desglose con porcentajes y speedup.
//
// Uso:
//   StageTimer st;
//   for (int step = 0; step < N; ++step) {
//     st.start("CIC");   cic_deposit(s);   st.stop("CIC");
//     st.start("FFT");   solve_poisson(s); st.stop("FFT");
//     ...
//   }
//   st.report(n_threads);
// =============================================================================
class StageTimer {
public:
    // Registra el inicio del cronómetro para una etapa
    void start(const std::string& name) {
        starts_[name] = wall_time();
    }

    // Detiene el cronómetro y acumula el tiempo
    void stop(const std::string& name) {
        double dt = wall_time() - starts_.at(name);
        totals_[name] += dt;
        if (std::find(order_.begin(), order_.end(), name) == order_.end())
            order_.push_back(name);  // conserva el orden de inserción
    }

    // Devuelve el tiempo total acumulado de una etapa
    double total(const std::string& name) const {
        auto it = totals_.find(name);
        return (it != totals_.end()) ? it->second : 0.0;
    }

    // Suma de todos los tiempos registrados
    double grand_total() const {
        double s = 0;
        for (auto& kv : totals_) s += kv.second;
        return s;
    }

    // ──────────────────────────────────────────────────────────
    // report()
    // Imprime tabla de desglose:
    //   - tiempo total por etapa
    //   - porcentaje del total
    //   - si se pasa serial_total, calcula speedup por etapa
    //
    // n_threads: número de hilos usados (para la cabecera)
    // serial_total: tiempos del run secuencial para comparar
    // ──────────────────────────────────────────────────────────
    void report(int n_threads = 1,
                const StageTimer* serial = nullptr) const {
        double tot = grand_total();
        int w = 22;

        std::cout << "\n";
        std::cout << "╔══════════════════════════════════════════════════════╗\n";
        std::cout << "║     Desglose de tiempos  |  hilos = "
                  << std::setw(2) << n_threads << "              ║\n";
        std::cout << "╠══════════════════════════════════════════════════════╣\n";
        std::cout << std::left
                  << "║ " << std::setw(w) << "Etapa"
                  << std::right << std::setw(8)  << "t (s)"
                  << std::setw(7)  << "%"
                  << std::setw(10) << "speedup"
                  << " ║\n";
        std::cout << "╠══════════════════════════════════════════════════════╣\n";

        for (auto& name : order_) {
            double t = totals_.at(name);
            double pct = (tot > 0) ? 100.0 * t / tot : 0.0;
            double sp  = 1.0;
            if (serial) {
                double ts = serial->total(name);
                sp = (t > 0 && ts > 0) ? ts / t : 1.0;
            }
            std::cout << "║ "
                      << std::left  << std::setw(w)   << name
                      << std::right << std::fixed << std::setprecision(4)
                      << std::setw(8)  << t
                      << std::setw(6)  << std::setprecision(1) << pct << "%"
                      << std::setw(9)  << std::setprecision(2) << sp  << "x"
                      << " ║\n";
        }

        std::cout << "╠══════════════════════════════════════════════════════╣\n";

        double global_sp = 1.0;
        if (serial) {
            double st = serial->grand_total();
            global_sp = (tot > 0 && st > 0) ? st / tot : 1.0;
        }
        double efficiency = global_sp / std::max(n_threads, 1);

        std::cout << "║ "
                  << std::left  << std::setw(w)   << "TOTAL"
                  << std::right << std::fixed << std::setprecision(4)
                  << std::setw(8)  << tot
                  << std::setw(7)  << "100%"
                  << std::setw(9)  << std::setprecision(2) << global_sp << "x"
                  << " ║\n";
        std::cout << "║ "
                  << std::left  << std::setw(w+14) << "Eficiencia"
                  << std::right << std::setprecision(1) << efficiency * 100.0
                  << "%"
                  << "          ║\n";
        std::cout << "╚══════════════════════════════════════════════════════╝\n\n";
    }

private:
    std::unordered_map<std::string, double> starts_;
    std::unordered_map<std::string, double> totals_;
    std::vector<std::string> order_;   // para preservar orden de impresión
};

} // namespace pm
