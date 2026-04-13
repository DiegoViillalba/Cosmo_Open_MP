// =============================================================================
// simulation.cpp
// Loop principal del simulador Particle-Mesh.
//
// Orquesta el pipeline completo en cada paso de tiempo:
//   CIC → Poisson → Gradiente → Interpolación → Leap-Frog
//
// También gestiona los timers para el reporte de desempeño y la escritura
// de snapshots en los intervalos configurados.
// =============================================================================

#include "simulation.hpp"
#include "ic_reader.hpp"
#include "cic.hpp"
#include "poisson_fft.hpp"
#include "gradient.hpp"
#include "force_interp.hpp"
#include "integrator.hpp"
#include "diagnostics.hpp"
#include "output.hpp"
#include "config.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <filesystem>   // std::filesystem::create_directories (C++17)
#include <stdexcept>
#include <array>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace pm {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Simulation::Simulation(const RunConfig& cfg) : cfg_(cfg) {

    // Configurar número de hilos OpenMP si se especificó
    #ifdef _OPENMP
    if (cfg_.n_threads > 0) {
        omp_set_num_threads(cfg_.n_threads);
        std::cout << "[sim] Hilos OpenMP configurados: " << cfg_.n_threads << "\n";
    } else {
        std::cout << "[sim] Hilos OpenMP (automático): "
                  << omp_get_max_threads() << "\n";
    }
    #else
    std::cout << "[sim] Compilado SIN OpenMP (modo secuencial)\n";
    #endif

    // Crear directorio de salida si no existe
    std::filesystem::create_directories(cfg_.output_dir);
}

// -----------------------------------------------------------------------------
// init
// Carga o genera las condiciones iniciales y prepara el estado inicial.
// -----------------------------------------------------------------------------
void Simulation::init() {

    std::cout << "\n[sim] Inicializando simulación...\n";
    std::cout << "[sim] Malla: " << Ng << "^3 = " << N_CELLS << " celdas\n";
    std::cout << "[sim] Partículas: " << N_PARTICLES << "\n";
    std::cout << "[sim] Pasos: "
              << (cfg_.n_steps > 0 ? cfg_.n_steps : N_STEPS) << "\n\n";

    // ── Leer o generar condiciones iniciales ───────────────────────────
    if (!cfg_.ic_file.empty()) {
        read_ic(state_, cfg_.ic_file);
    } else {
        std::cout << "[sim] No se proporcionó IC → generando grilla uniforme\n";
        generate_uniform_ic(state_);
    }

    // ── Primer cálculo de fuerza (necesario para el half-kick inicial) ─
    // El Leap-Frog requiere las fuerzas en t=0 para retroceder las
    // velocidades al tiempo t = -dt/2.
    std::cout << "[sim] Calculando fuerza inicial...\n";

    state_.clear_grids();
    cic_deposit(state_);
    solve_poisson(state_);
    compute_gradient(state_);

    std::vector<std::array<double,3>> accel(N_PARTICLES);
    interpolate_force(state_, accel);

    // Retroceder velocidades medio paso para el esquema Leap-Frog
    leapfrog_half_kick(state_, accel);

    // Escribir diagnósticos y snapshot inicial
    auto diag = compute_diagnostics(state_);
    if (cfg_.verbose)
        print_diagnostics(0, state_.time, diag);

    const std::string diag_file = cfg_.output_dir + "/diagnostics.csv";
    append_diagnostics_csv(diag_file, 0, 0.0,
                           diag.kinetic_energy,
                           diag.potential_energy,
                           diag.total_energy);

    // Snapshot inicial
    write_snapshot(0);
    std::cout << "[sim] Inicialización completa.\n\n";
}

// -----------------------------------------------------------------------------
// run
// Loop principal: ejecuta todos los pasos de tiempo y reporta desempeño.
// -----------------------------------------------------------------------------
const StageTimer& Simulation::run() {

    const int n_steps = (cfg_.n_steps > 0) ? cfg_.n_steps : N_STEPS;
    const std::string diag_file = cfg_.output_dir + "/diagnostics.csv";

    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║         Iniciando loop de simulación         ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n\n";

    // Tiempo total del run
    const double t_run_start = wall_time();

    for (int s = 1; s <= n_steps; ++s) {
        state_.step = s;
        single_step();

        // ── Diagnósticos (cada 10 pasos o si verbose) ─────────────────
        if (cfg_.verbose || s % 10 == 0) {
            auto diag = compute_diagnostics(state_);
            if (cfg_.verbose)
                print_diagnostics(s, state_.time, diag);

            append_diagnostics_csv(diag_file, s, state_.time,
                                   diag.kinetic_energy,
                                   diag.potential_energy,
                                   diag.total_energy);
        }

        // ── Snapshot ──────────────────────────────────────────────────
        if (s % SNAP_INTERVAL == 0)
            write_snapshot(s);

        // ── Barra de progreso simple ───────────────────────────────────
        if (s % (n_steps / 10 + 1) == 0 || s == n_steps) {
            int pct = static_cast<int>(100.0 * s / n_steps);
            double t_elapsed = wall_time() - t_run_start;
            double t_est = (s > 0) ? t_elapsed * (n_steps - s) / s : 0.0;
            std::cout << "  [" << std::setw(3) << pct << "%] "
                      << "paso " << std::setw(5) << s << "/" << n_steps
                      << "  transcurrido: " << std::fixed
                      << std::setprecision(1) << t_elapsed << "s"
                      << "  restante: ~" << t_est << "s\n";
        }
    }

    // ── Reporte final de tiempos ───────────────────────────────────────
    #ifdef _OPENMP
    const int nt = omp_get_max_threads();
    #else
    const int nt = 1;
    #endif
    timer_.report(nt);

    return timer_;
}

// -----------------------------------------------------------------------------
// single_step (privado)
// Ejecuta una iteración completa del pipeline PM con medición de tiempos.
// -----------------------------------------------------------------------------
void Simulation::single_step() {

    // Limpiar mallas del paso anterior
    state_.clear_grids();

    // ── 1. CIC deposit: partículas → ρ(x) ─────────────────────────────
    timer_.start("1. CIC deposit");
    cic_deposit(state_);
    timer_.stop("1. CIC deposit");

    // ── 2. Poisson FFT: ρ(x) → φ(x) ──────────────────────────────────
    timer_.start("2. Poisson FFT");
    solve_poisson(state_);
    timer_.stop("2. Poisson FFT");

    // ── 3. Gradiente: φ(x) → g(x) ─────────────────────────────────────
    timer_.start("3. Gradiente");
    compute_gradient(state_);
    timer_.stop("3. Gradiente");

    // ── 4. Interpolación: g(x) → a_i ──────────────────────────────────
    timer_.start("4. Force interp");
    std::vector<std::array<double,3>> accel;
    interpolate_force(state_, accel);
    timer_.stop("4. Force interp");

    // ── 5. Integrador: (x,v) → (x',v') ───────────────────────────────
    timer_.start("5. Integrador");
    leapfrog_step(state_, accel);
    timer_.stop("5. Integrador");
}

// -----------------------------------------------------------------------------
// write_snapshot (privado)
// Decide formato y construye el nombre del archivo.
// -----------------------------------------------------------------------------
void Simulation::write_snapshot(int step) {
    std::ostringstream name;
    name << cfg_.output_dir << "/snap_"
         << std::setw(4) << std::setfill('0') << step;

    if (cfg_.use_ascii) {
        write_snapshot_ascii(state_, name.str() + ".dat");
    } else {
        write_snapshot_binary(state_, name.str() + ".bin");
    }
}

} // namespace pm
