#pragma once

// =============================================================================
// output.hpp
// Escritura de snapshots de la simulación a disco.
//
// Se proveen dos formatos:
//
//   ASCII (write_snapshot_ascii):
//     Legible por humanos y por scripts de Python/MATLAB directamente.
//     Formato: una partícula por línea → "x y z vx vy vz"
//     Lento para 256³ (16M líneas), pero ideal para 128³ y debugging.
//
//   Binario (write_snapshot_binary):
//     Escribe los structs Particle como bytes crudos.
//     ~6x más rápido y ~3x más compacto que ASCII.
//     Para leerlo en Python: np.fromfile(f, dtype=np.float64).reshape(-1,6)
//
// También se escribe un archivo de diagnósticos (energía vs tiempo)
// en formato CSV para graficar la conservación de energía.
// =============================================================================

#include <string>
#include "types.hpp"

namespace pm {

// Escribe posiciones y velocidades en formato ASCII.
// filename: ej. "data/snap_005.dat"
void write_snapshot_ascii(const SimState& state, const std::string& filename);

// Escribe en formato binario (más rápido).
void write_snapshot_binary(const SimState& state, const std::string& filename);

// Agrega una línea al archivo de diagnósticos CSV.
// Crea el archivo con cabecera si no existe.
// diag_file: ej. "data/diagnostics.csv"
void append_diagnostics_csv(const std::string& diag_file,
                            int step, double time,
                            double Ek, double Ep, double Etot);

} // namespace pm
