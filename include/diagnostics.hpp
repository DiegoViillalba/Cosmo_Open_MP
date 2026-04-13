#pragma once

// =============================================================================
// diagnostics.hpp
// Cálculo de energías y diagnósticos físicos de la simulación.
//
// Los diagnósticos sirven para verificar que la simulación es correcta:
//   - La energía total E = Ek + Ep se conserva (aproximadamente) si el
//     paso de tiempo DT es suficientemente pequeño.
//   - El momento total P debe conservarse exactamente (simetría traslacional).
//
// CONEXIÓN CON OPENMP:
// La energía cinética y el momento son sumas sobre todas las partículas.
// Este es el caso paradigmático de la cláusula reduction:
//
//   double Ek = 0;
//   #pragma omp parallel for reduction(+:Ek)
//   for (int i = 0; i < N; i++)
//       Ek += 0.5 * m * (vx²+vy²+vz²);
//
// Sin reduction: múltiples hilos escribirían en Ek simultáneamente → carrera.
// Con reduction: cada hilo tiene su copia privada de Ek, y al final OpenMP
//               las suma automáticamente. Correcto Y eficiente.
// =============================================================================

#include "types.hpp"

namespace pm {

// Resultados de los diagnósticos en un paso dado
struct Diagnostics {
    double kinetic_energy;      // Ek = Σ ½ m v²
    double potential_energy;    // Ep = Σ m φ(x_i)  (interpolada)
    double total_energy;        // E  = Ek + Ep
    double momentum[3];         // P  = Σ m v   (debe conservarse)
    double virial_ratio;        // |Ep| / (2 Ek)  (= 1 en equilibrio virial)
};

// Calcula todos los diagnósticos en el paso actual.
// Usa reduction de OpenMP para las sumas.
Diagnostics compute_diagnostics(const SimState& state);

// Imprime una línea de diagnóstico formateada
void print_diagnostics(int step, double time, const Diagnostics& d);

} // namespace pm
