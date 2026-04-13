#pragma once

// =============================================================================
// gradient.hpp
// Calcula el campo gravitacional g = -∇φ por diferencias finitas centradas.
//
// Una vez que tenemos el potencial φ en cada celda de la malla, la fuerza
// gravitacional por unidad de masa en cada celda es:
//
//   g_x(i,j,k) = -(φ(i+1,j,k) - φ(i-1,j,k)) / (2·Δx)
//   g_y(i,j,k) = -(φ(i,j+1,k) - φ(i,j-1,k)) / (2·Δy)
//   g_z(i,j,k) = -(φ(i,j,k+1) - φ(i,j,k-1)) / (2·Δz)
//
// Las diferencias centradas dan orden de error O(Δx²), mientras que las
// diferencias hacia adelante solo dan O(Δx). Vale la pena el costo extra
// de acceder dos celdas vecinas en lugar de una.
//
// Las condiciones de frontera periódicas se manejan automáticamente con
// Grid3D::at(), que aplica módulo a los índices.
//
// Complejidad: O(N_CELLS) → embarrasingly parallel, el loop se distribuye
// directamente con #pragma omp parallel for sin ninguna condición de carrera
// (cada celda escribe en una posición diferente de force_x/y/z).
// =============================================================================

#include "types.hpp"

namespace pm {

// Calcula g = -∇φ y almacena en state.force_x/y/z.
// Usa diferencias finitas centradas y condiciones de frontera periódicas.
void compute_gradient(SimState& state);

} // namespace pm
