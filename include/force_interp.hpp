#pragma once

// =============================================================================
// force_interp.hpp
// Interpolación de la fuerza desde la malla a cada partícula (CIC inverso).
//
// Esta es la operación inversa al depósito CIC: en lugar de distribuir masa
// de partícula → malla, ahora distribuimos fuerza de malla → partícula,
// usando exactamente los mismos pesos trilineales.
//
// Para una partícula en posición (px, py, pz):
//   1. Se calcula la celda base (i0,j0,k0) = floor(px,py,pz)
//   2. Se calculan los pesos dx=px-i0, dy=py-j0, dz=pz-k0
//   3. La fuerza se interpola de las 8 celdas vecinas:
//      f = Σ_{δx,δy,δz ∈ {0,1}} w_x^δx · w_y^δy · w_z^δz · force(i0+δx, j0+δy, k0+δz)
//      donde w_x^0 = (1-dx), w_x^1 = dx, etc.
//
// Este paso es puro de lectura en las mallas: cada hilo lee posiciones
// distintas de force_x/y/z (sin escritura) y escribe en su propia partícula.
// No hay condición de carrera → paralelización trivial.
// =============================================================================

#include "types.hpp"

namespace pm {

// Interpola fuerza de malla a partículas usando pesos CIC.
// Lee state.force_x/y/z, escribe aceleración temporal (no persistente).
// La aceleración se aplica directamente al integrador en integrator.cpp.
void interpolate_force(SimState& state,
                       std::vector<std::array<double,3>>& accel);

} // namespace pm
