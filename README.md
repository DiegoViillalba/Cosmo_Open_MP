# Ilhuca cosmo — Simulador Particle-Mesh Cosmológico con OpenMP

**TSCAD 2026-2 · Física Biomédica · UNAM · Equipo 4**

Simulador de N-cuerpos gravitacional en C++17 usando el método Particle-Mesh (PM).  
Paralelizado con OpenMP. Incluye versiones serial y paralela con medición de speedup por etapa.

---

## Dependencias

| Paquete | Ubuntu/Debian |
|---|---|
| Compilador C++17 | `build-essential` |
| OpenMP | incluido en GCC |
| FFTW3 | `libfftw3-dev` |
| CMake ≥ 3.16 | `cmake` |

```bash
sudo apt install build-essential cmake libfftw3-dev pkg-config
```

---

## Compilación rápida

```bash
# Clonar / descomprimir el proyecto
cd pm_cosmo

# Modo prueba — Ng=128, ~2M partículas (~16 MB de malla)
mkdir build && cd build
cmake .. -DPM_GRID_SIZE=128 -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Modo producción — Ng=256, ~16M partículas (~128 MB de malla)
cmake .. -DPM_GRID_SIZE=256 -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Esto genera tres ejecutables en `build/`:

| Ejecutable | Descripción |
|---|---|
| `pm_serial` | Compilado **sin** OpenMP — línea base de comparación |
| `pm_parallel` | Compilado **con** OpenMP — versión paralela |
| `run_tests` | Suite de 17 tests unitarios |

---

## Uso

```bash
# Prueba rápida con IC generadas internamente
./pm_serial   -g -s 20 -v
./pm_parallel -g -s 20 -t 8 -v

# Con archivo de condiciones iniciales externas
./pm_parallel ../ic/ic_128.dat -t 8 -s 100 -o ../data

# Opciones
#   -t <N>   Número de hilos OpenMP
#   -s <N>   Número de pasos de tiempo
#   -o <dir> Directorio de salida (snapshots + diagnostics.csv)
#   -g       Generar IC sintéticas (grilla regular perturbada)
#   -b       Snapshots en formato binario (más rápido que ASCII)
#   -v       Verbose: tabla de energía en cada paso
#   -h       Ayuda
```

---

## Ejecutar los tests

```bash
cd build
./run_tests              # todos los tests
./run_tests cic          # solo tests del módulo CIC
./run_tests poisson      # solo tests de Poisson
./run_tests integrador   # solo tests del integrador
./run_tests diag         # solo tests de diagnósticos
```

Salida esperada: `17/17 TODOS LOS TESTS PASARON.`

---

## Benchmark de speedup

```bash
cd build
bash ../bench/run_bench.sh 20 ../bench/results.csv
```

Corre `pm_serial` y `pm_parallel` con 1, 2, 4, 8, 16 hilos e imprime la tabla de speedup por etapa.

---

## Parámetros de la simulación

Todos en `include/config.hpp`. Los más importantes:

| Parámetro | Default | Descripción |
|---|---|---|
| `Ng` | 256 (o 128) | Puntos de malla por dimensión |
| `N_PARTICLES` | Ng³ | Total de partículas |
| `L_BOX` | 100.0 Mpc/h | Tamaño de la caja |
| `DT` | 0.01 | Paso de tiempo |
| `N_STEPS` | 100 | Pasos totales |
| `SNAP_INTERVAL` | 10 | Cada cuántos pasos guardar snapshot |

---

## Estructura de archivos

```
pm_cosmo/
├── include/          Headers públicos de cada módulo
│   ├── config.hpp    Parámetros globales (Ng, DT, N_STEPS…)
│   ├── types.hpp     Estructuras Particle, Grid3D, SimState
│   ├── timer.hpp     ScopedTimer y StageTimer (omp_get_wtime)
│   ├── ic_reader.hpp Lector de condiciones iniciales
│   ├── cic.hpp       Cloud-in-Cell deposit
│   ├── poisson_fft.hpp  Solucionador ∇²φ = 4πGρ con FFTW3
│   ├── gradient.hpp  Campo gravitacional g = −∇φ
│   ├── force_interp.hpp Interpolación fuerza malla→partícula
│   ├── integrator.hpp   Integrador Leap-Frog
│   ├── diagnostics.hpp  Energías y momento (reduction OMP)
│   ├── output.hpp    Snapshots ASCII/binario + CSV
│   └── simulation.hpp   Orquestador del loop principal
├── src/              Implementaciones .cpp
├── tests/            Suite de tests unitarios
│   ├── test_framework.hpp  Macros REGISTER_TEST, ASSERT_*
│   ├── test_main.cpp       Runner de tests
│   ├── test_cic.cpp        4 tests del depósito CIC
│   ├── test_poisson.cpp    4 tests del solucionador Poisson
│   ├── test_integrator.cpp 4 tests del integrador Leap-Frog
│   └── test_diagnostics.cpp 5 tests de energías y reduction
├── ic/               Archivos de condiciones iniciales
│   ├── ic_128.dat    IC para Ng=128 (generado externamente)
│   └── ic_256.dat    IC para Ng=256 (generado externamente)
├── data/             Snapshots de salida (generado al correr)
├── bench/
│   ├── run_bench.sh  Script de benchmarking automático
│   └── results.csv   Tabla de speedup (generado por el script)
└── CMakeLists.txt    Build system con targets serial/paralelo/tests
```

---

## Flujo de llamadas y dependencias entre módulos

El diagrama siguiente muestra qué función llama a cuál, con los argumentos exactos que se pasan en cada invocación dentro de una ejecución completa.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  main.cpp                                                                   │
│  argv → RunConfig{ic_file, n_threads, n_steps, output_dir, use_ascii, verb} │
└──────────────────────────┬──────────────────────────────────────────────────┘
                           │ Simulation sim(cfg)
                           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  simulation.cpp :: Simulation                                               │
│                                                                             │
│  ┌── init() ──────────────────────────────────────────────────────────────┐ │
│  │                                                                        │ │
│  │  ┌─ IC ──────────────────────────────────────────────────────────────┐ │ │
│  │  │  read_ic(state_, cfg_.ic_file)                                    │ │ │
│  │  │    └─ ic_reader.cpp: abre archivo, llena state_.particles[]       │ │ │
│  │  │  — o —                                                            │ │ │
│  │  │  generate_uniform_ic(state_, sigma_pos=0.5, sigma_vel=0.1,       │ │ │
│  │  │                       seed=42)                                    │ │ │
│  │  │    └─ ic_reader.cpp: grilla regular + perturbación gaussiana      │ │ │
│  │  └───────────────────────────────────────────────────────────────────┘ │ │
│  │                                                                        │ │
│  │  ┌─ Fuerza inicial (t=0) ────────────────────────────────────────────┐ │ │
│  │  │  cic_deposit(state_)                                              │ │ │
│  │  │    → escribe state_.density[]                                     │ │ │
│  │  │  solve_poisson(state_)                                            │ │ │
│  │  │    → lee state_.density[], escribe state_.potential[]             │ │ │
│  │  │  compute_gradient(state_)                                         │ │ │
│  │  │    → lee state_.potential[], escribe state_.force_x/y/z[]        │ │ │
│  │  │  interpolate_force(state_, accel[N_PART])                        │ │ │
│  │  │    → lee state_.force_x/y/z[], escribe accel[i]={ax,ay,az}       │ │ │
│  │  │  leapfrog_half_kick(state_, accel)                                │ │ │
│  │  │    → modifica state_.particles[i].vel (retrocede dt/2)           │ │ │
│  │  └───────────────────────────────────────────────────────────────────┘ │ │
│  │                                                                        │ │
│  │  ┌─ Diagnóstico t=0 ─────────────────────────────────────────────────┐ │ │
│  │  │  compute_diagnostics(state_) → Diagnostics{Ek,Ep,Etot,P,virial}  │ │ │
│  │  │  print_diagnostics(0, 0.0, diag)          [si verbose]            │ │ │
│  │  │  append_diagnostics_csv(output_dir+"/diagnostics.csv",           │ │ │
│  │  │                         step=0, time=0.0, Ek, Ep, Etot)          │ │ │
│  │  │  write_snapshot_ascii/binary(state_,                             │ │ │
│  │  │                              output_dir+"/snap_0000.dat/.bin")    │ │ │
│  │  └───────────────────────────────────────────────────────────────────┘ │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
│                                                                             │
│  ┌── run() ───────────────────────────────────────────────────────────────┐ │
│  │  for s in 1..N_STEPS:                                                  │ │
│  │    single_step()  ◄────────────── pipeline PM (ver abajo)              │ │
│  │    if s % 10 == 0:                                                     │ │
│  │      compute_diagnostics(state_) → diag                               │ │
│  │      print_diagnostics(s, state_.time, diag)    [si verbose]           │ │
│  │      append_diagnostics_csv(..., s, state_.time, Ek, Ep, Etot)        │ │
│  │    if s % SNAP_INTERVAL == 0:                                          │ │
│  │      write_snapshot_ascii/binary(state_, "snap_NNNN.dat/.bin")        │ │
│  │  timer_.report(n_threads)                                              │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  single_step() — pipeline PM, ejecutado N_STEPS veces                      │
│                                                                             │
│   state_.clear_grids()      limpia density, potential, force_x/y/z         │
│          │                                                                  │
│          ▼                                                                  │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │ 1. cic_deposit(state_)                               [OMP atomic]   │  │
│   │    IN : state_.particles[i].{pos, mass}                             │  │
│   │    OUT: state_.density[Ng³]                                         │  │
│   │    Paralelismo: loop sobre N_PART, atomic en 8 escrituras/partícula │  │
│   └──────────────────────────┬──────────────────────────────────────────┘  │
│                              │ state_.density                               │
│                              ▼                                              │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │ 2. solve_poisson(state_)                         [OMP + FFTW3]     │  │
│   │    IN : state_.density[Ng³]                                         │  │
│   │    OUT: state_.potential[Ng³]                                       │  │
│   │    Pasos: FFT(ρ) → ×kernel(−4πG/k²) → IFFT → normalizar           │  │
│   │    Paralelismo: loop de multiplicación con #pragma omp parallel for │  │
│   └──────────────────────────┬──────────────────────────────────────────┘  │
│                              │ state_.potential                             │
│                              ▼                                              │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │ 3. compute_gradient(state_)                     [OMP collapse(3)]  │  │
│   │    IN : state_.potential[Ng³]                                       │  │
│   │    OUT: state_.force_x/y/z[Ng³]                                    │  │
│   │    Fórmula: g_x = −(φ[i+1]−φ[i−1])/(2·Δx)  (dif. centradas)      │  │
│   │    Paralelismo: triple loop sobre malla, sin condiciones de carrera │  │
│   └──────────────────────────┬──────────────────────────────────────────┘  │
│                              │ state_.force_x/y/z                          │
│                              ▼                                              │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │ 4. interpolate_force(state_, accel[N_PART])     [OMP schedule(static)]│ │
│   │    IN : state_.force_x/y/z[Ng³]                                    │  │
│   │         state_.particles[i].pos                                     │  │
│   │    OUT: accel[i] = {ax, ay, az}  (pesos CIC inverso)               │  │
│   │    Paralelismo: loop sobre N_PART, solo lectura en mallas           │  │
│   └──────────────────────────┬──────────────────────────────────────────┘  │
│                              │ accel[N_PART]                                │
│                              ▼                                              │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │ 5. leapfrog_step(state_, accel)                 [OMP schedule(static)]│ │
│   │    IN : accel[i]                                                    │  │
│   │         state_.particles[i].{pos, vel}                              │  │
│   │    OUT: state_.particles[i].{pos, vel}  (actualizados)              │  │
│   │    Kick:  vel += accel·dt                                           │  │
│   │    Drift: pos += vel·dt  (con wraparound periódico)                 │  │
│   │    Paralelismo: loop sobre N_PART, sin dependencias entre partículas│  │
│   └─────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│   Tiempo medido por etapa con StageTimer (omp_get_wtime)                   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  run_tests — Suite de tests unitarios                                       │
│                                                                             │
│  test_cic.cpp        → cic_deposit(state)                                  │
│    cic_conservacion_masa  · cic_split_50_50                                │
│    cic_no_negativos       · cic_grilla_uniforme                            │
│                                                                             │
│  test_poisson.cpp    → solve_poisson(state)                                │
│    poisson_modo_dc_cero   · poisson_forma_coseno                           │
│    poisson_no_nan_inf     · poisson_linealidad                             │
│                                                                             │
│  test_integrator.cpp → leapfrog_step / leapfrog_half_kick                  │
│    integrador_mru         · integrador_periodicidad                        │
│    integrador_vel_constante · integrador_mua                               │
│                                                                             │
│  test_diagnostics.cpp → compute_diagnostics(state)                         │
│    diag_ek_particula_unica · diag_momento_cero                             │
│    diag_reposo             · diag_energia_total_consistente                │
│    diag_ek_escala_masa  ← verifica reduction OMP                           │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Estrategias de paralelismo por módulo

| Módulo | Directiva OMP | Patrón | Race condition |
|---|---|---|---|
| `cic.cpp` | `parallel for` + `atomic` | N_PART iteraciones, 8 escrituras/iter en malla compartida | **Sí** → resuelta con `atomic` |
| `poisson_fft.cpp` | `parallel for collapse(2)` | Loop sobre modos de Fourier Ng²×(Ng/2+1) | No (cada modo es independiente) |
| `gradient.cpp` | `parallel for collapse(3)` | Loop sobre todas las celdas Ng³ | No (cada celda escribe en posición única) |
| `force_interp.cpp` | `parallel for schedule(static)` | N_PART iters, solo lectura en mallas | No (escritura en `accel[i]` privado) |
| `integrator.cpp` | `parallel for schedule(static)` | N_PART iters, lectura `accel[i]`, escritura `particles[i]` | No (cada partícula es independiente) |
| `diagnostics.cpp` | `parallel for reduction(+:Ek,Px,Py,Pz)` | Suma acumulada sobre N_PART | **Sí** → resuelta con `reduction` |

---

## Datos de salida

```
data/
├── snap_0000.dat     Snapshot inicial  (x y z vx vy vz por línea)
├── snap_0010.dat     Snapshot paso 10
├── snap_NNNN.dat     ...
└── diagnostics.csv   step,time,Ek,Ep,Etot,dE_frac
```

Para visualizar en Python:
```python
import numpy as np
import matplotlib.pyplot as plt

# Leer un snapshot
data = np.loadtxt("data/snap_0100.dat", skiprows=1)
pos = data[:, :3]   # x, y, z  en unidades de celda

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(pos[::100, 0], pos[::100, 1], pos[::100, 2], s=0.1)
plt.show()

# Evolución de energía
import pandas as pd
diag = pd.read_csv("data/diagnostics.csv")
plt.plot(diag.time, diag.Etot, label="E total")
plt.plot(diag.time, diag.Ek,   label="E cinética")
plt.plot(diag.time, diag.Ep,   label="E potencial")
plt.legend(); plt.xlabel("Tiempo"); plt.show()
```

---

## Ejemplo de salida del benchmark

```
modo         hilos      tiempo(s)    speedup
----         -----      ---------    -------
serial       1          18.34        1.00x
parallel     1          19.10        0.96x   ← overhead de OpenMP con 1 hilo
parallel     2          10.21        1.80x
parallel     4           5.48        3.35x
parallel     8           3.12        5.88x
parallel     16          2.19        8.38x   ← escalabilidad limitada por FFT
```

> El solucionador de Poisson (FFT secuencial) es el cuello de botella para > 8 hilos.  
> Los módulos CIC, gradiente e integrador escalan casi linealmente hasta 8 hilos.

---

## Créditos

Proyecto desarrollado para la materia **Temas Selectos de Cómputo de Alto Desempeño 2026-2**  
Facultad de Ciencias, UNAM · Física Biomédica · Equipo 4
