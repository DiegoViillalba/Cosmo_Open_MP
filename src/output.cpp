// =============================================================================
// output.cpp
// Escritura de snapshots de partículas y archivo de diagnósticos.
//
// La escritura a disco es secuencial por diseño: el cuello de botella es el
// I/O, no el CPU, así que paralelizar no ayuda y complica la implementación.
// Esta es una decisión técnica deliberada que vale la pena mencionar en la
// exposición al discutir qué partes vale la pena paralelizar.
// =============================================================================

#include "output.hpp"
#include "config.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cstdio>    // std::remove para verificar existencia
#include <cstdint>   // uint64_t

namespace pm {

// -----------------------------------------------------------------------------
// write_snapshot_ascii
// Escribe todas las partículas en formato texto legible.
// Cada línea: x y z vx vy vz
// La primera línea contiene el número de partículas (para read_ic).
// -----------------------------------------------------------------------------
void write_snapshot_ascii(const SimState& state, const std::string& filename) {

    std::ofstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("output: no se puede escribir '" + filename + "'");

    const auto& parts = state.particles;
    f << parts.size() << "\n";
    f << std::fixed << std::setprecision(8);

    for (const auto& p : parts) {
        f << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
          << p.vel[0] << " " << p.vel[1] << " " << p.vel[2] << "\n";
    }

    std::cout << "[output] Snapshot ASCII → " << filename
              << "  (" << parts.size() << " partículas)\n";
}

// -----------------------------------------------------------------------------
// write_snapshot_binary
// Escribe las partículas como bytes crudos (struct Particle completo).
//
// Para leer en Python:
//   import numpy as np
//   raw = np.fromfile("snap_000.bin", dtype=np.float64)
//   # Cada partícula tiene 7 doubles: x,y,z, vx,vy,vz, mass
//   data = raw.reshape(-1, 7)
//   pos = data[:, :3]
//   vel = data[:, 3:6]
// -----------------------------------------------------------------------------
void write_snapshot_binary(const SimState& state, const std::string& filename) {

    std::ofstream f(filename, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("output: no se puede escribir '" + filename + "'");

    // Cabecera: número de partículas como uint64
    const uint64_t n = state.particles.size();
    f.write(reinterpret_cast<const char*>(&n), sizeof(n));

    // Datos: volcado directo del vector de Particle
    // Nota: esto asume que Particle no tiene padding inesperado.
    // En la práctica con doubles alineados no hay problema.
    f.write(reinterpret_cast<const char*>(state.particles.data()),
            n * sizeof(Particle));

    std::cout << "[output] Snapshot binario → " << filename
              << "  (" << n << " partículas, "
              << (n * sizeof(Particle)) / (1024*1024) << " MB)\n";
}

// -----------------------------------------------------------------------------
// append_diagnostics_csv
// Agrega una fila al CSV de evolución de energías.
// Crea el archivo con cabecera si no existe todavía.
// -----------------------------------------------------------------------------
void append_diagnostics_csv(const std::string& diag_file,
                            int step, double time,
                            double Ek, double Ep, double Etot) {

    // Verificar si el archivo ya existe para saber si hay que poner cabecera
    bool needs_header = false;
    {
        std::ifstream test(diag_file);
        needs_header = !test.good();
    }

    std::ofstream f(diag_file, std::ios::app);  // append mode
    if (!f.is_open())
        throw std::runtime_error("output: no se puede escribir '" + diag_file + "'");

    if (needs_header)
        f << "step,time,Ek,Ep,Etot,dE_frac\n";

    // Fracción de cambio de energía respecto al paso 0 (se llena externamente)
    const double dE = (std::abs(Etot) > 1e-30)
                      ? (Etot - Ek - Ep) / std::abs(Etot)
                      : 0.0;

    f << step << ","
      << std::fixed << std::setprecision(6) << time << ","
      << std::scientific << std::setprecision(8)
      << Ek   << ","
      << Ep   << ","
      << Etot << ","
      << dE   << "\n";
}

} // namespace pm
