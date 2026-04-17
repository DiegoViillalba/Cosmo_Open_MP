#pragma once

// =============================================================================
// simulation.hpp
// Orquestador: contiene el loop principal de la simulación PM.
//
// El pipeline completo de un paso de tiempo es:
//
//   ┌─────────────────────────────────────────────────┐
//   │  Step n                                         │
//   │                                                 │
//   │  1. cic_deposit()      partículas → ρ(x)        │
//   │  2. solve_poisson()    ρ(x) → φ(k) → φ(x)       │
//   │  3. compute_gradient() φ(x) → g(x)              │
//   │  4. interpolate_force() g(x) → a_i por partícula│
//   │  5. leapfrog_step()    (x,v) → (x',v')          │
//   │  6. [diagnostics]      si step % DIAG_INTERVAL  │
//   │  7. [output]           si step % SNAP_INTERVAL  │
//   └─────────────────────────────────────────────────┘
//
// La clase Simulation encapsula el estado y los timers para que main.cpp
// sea lo más limpio posible.
// =============================================================================

#include <string>
#include "types.hpp"
#include "timer.hpp"

namespace pm {

// Parámetros de configuración pasados desde la línea de comandos
struct RunConfig {
    std::string ic_file;       // archivo de condiciones iniciales
    int n_threads;             // número de hilos OpenMP (0 = automático)
    int n_steps;               // pasos a correr (0 = usa N_STEPS de config.hpp)
    std::string output_dir;    // carpeta de salida para snapshots
    bool use_ascii;            // true = ASCII, false = binario
    bool verbose;              // imprime diagnósticos en cada paso
                               // TODO: Implement a frequancy of saving

    RunConfig()
        : ic_file(""), n_threads(0), n_steps(0),
          output_dir("data"), use_ascii(true), verbose(false) {}
};

class Simulation {
public:
    explicit Simulation(const RunConfig& cfg);

    // Inicializa el estado: lee IC o genera grilla uniforme
    void init();

    // Ejecuta todos los pasos de tiempo
    // Devuelve referencia a los timers para reportar desempeño al final
    const StageTimer& run();

    // Acceso de solo lectura al estado (para los tests)
    const SimState& state() const { return state_; }

private:
    RunConfig  cfg_;
    SimState   state_;
    StageTimer timer_;

    // Ejecuta un solo paso del pipeline PM
    void single_step();

    // Escribe snapshot con el formato y nombre adecuados
    void write_snapshot(int step);
};

} // namespace pm
