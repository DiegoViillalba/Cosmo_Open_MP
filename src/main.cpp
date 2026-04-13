// =============================================================================
// main.cpp
// Punto de entrada del simulador PM cosmológico.
//
// Uso:
//   ./pm_parallel [opciones] [ic_file]
//
// Opciones:
//   -t <N>      Número de hilos OpenMP (default: automático)
//   -s <N>      Número de pasos (default: N_STEPS de config.hpp)
//   -o <dir>    Directorio de salida (default: "data")
//   -b          Usar formato binario para snapshots (default: ASCII)
//   -v          Verbose: imprimir diagnósticos en cada paso
//   -g          Generar IC sintéticas (ignorar ic_file)
//   -h          Mostrar ayuda
//
// Ejemplos:
//   ./pm_parallel ic/ic_128.dat -t 4 -s 50 -v
//   ./pm_serial   -g -s 20 -o resultados/serial
//   ./pm_parallel -g -t 8  -s 20 -o resultados/p8
// =============================================================================

#include "simulation.hpp"
#include "config.hpp"

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstdlib>   // std::atoi

#ifdef _OPENMP
  #include <omp.h>
#endif

// Imprime ayuda y termina
static void print_usage(const char* argv0) {
    std::cout << "\nUso: " << argv0 << " [opciones] [ic_file]\n\n"
              << "Opciones:\n"
              << "  -t <N>   Hilos OpenMP (default: automático)\n"
              << "  -s <N>   Número de pasos (default: " << pm::N_STEPS << ")\n"
              << "  -o <dir> Directorio de salida (default: data)\n"
              << "  -b       Snapshots binarios en vez de ASCII\n"
              << "  -v       Verbose (diagnósticos en cada paso)\n"
              << "  -g       Generar IC sintéticas (no requiere archivo)\n"
              << "  -h       Esta ayuda\n\n"
              << "Configuración actual (config.hpp):\n"
              << "  Ng = " << pm::Ng << "  (tamaño de malla)\n"
              << "  Np = " << pm::Np << "  (partículas por dimensión)\n"
              << "  N_PARTICLES = " << pm::N_PARTICLES << "\n"
              << "  DT = " << pm::DT << "\n\n"
              << "Ejemplos:\n"
              << "  " << argv0 << " ic/ic_128.dat -t 4 -s 100\n"
              << "  " << argv0 << " -g -t 8 -s 50 -v\n\n";
}

int main(int argc, char* argv[]) {

    pm::RunConfig cfg;
    bool force_generate = false;

    // ── Parseo de argumentos ───────────────────────────────────────────
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        else if (arg == "-t" && i + 1 < argc) {
            cfg.n_threads = std::atoi(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc) {
            cfg.n_steps = std::atoi(argv[++i]);
        }
        else if (arg == "-o" && i + 1 < argc) {
            cfg.output_dir = argv[++i];
        }
        else if (arg == "-b") {
            cfg.use_ascii = false;
        }
        else if (arg == "-v") {
            cfg.verbose = true;
        }
        else if (arg == "-g") {
            force_generate = true;
        }
        else if (arg[0] != '-') {
            cfg.ic_file = arg;   // primer argumento no-flag = archivo de IC
        }
        else {
            std::cerr << "Opción desconocida: " << arg
                      << "  (usa -h para ayuda)\n";
            return 1;
        }
    }

    if (force_generate) cfg.ic_file = "";  // forzar generación interna

    // ── Banner inicial ─────────────────────────────────────────────────
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════╗\n";
    std::cout << "║    Simulador Particle-Mesh Cosmológico — UNAM      ║\n";
    std::cout << "║    TSCAD 2026-2  |  Física Biomédica               ║\n";
    std::cout << "╚════════════════════════════════════════════════════╝\n";

    #ifdef PM_SERIAL
    std::cout << "  Modo: SECUENCIAL (pm_serial)\n";
    #else
    std::cout << "  Modo: PARALELO   (pm_parallel)\n";
    #ifdef _OPENMP
    std::cout << "  OpenMP versión:  " << _OPENMP << "\n";
    #endif
    #endif
    std::cout << "\n";

    // ── Ejecutar simulación ────────────────────────────────────────────
    try {
        pm::Simulation sim(cfg);
        sim.init();
        sim.run();
    }
    catch (const std::exception& e) {
        std::cerr << "\n[ERROR] " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n[main] Simulación completada. Resultados en: "
              << cfg.output_dir << "/\n\n";
    return 0;
}
