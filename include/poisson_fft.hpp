#pragma once

// =============================================================================
// poisson_fft.hpp
// Resolución de la ecuación de Poisson gravitacional en espacio de Fourier.
//
// La ecuación a resolver es:
//   ∇²φ(x) = 4πG ρ(x)      (espacio real)
//
// En espacio de Fourier esto se convierte en una multiplicación algebraica:
//   φ̂(k) = -4πG ρ̂(k) / k²
//
// donde k² = kx² + ky² + kz² es el módulo cuadrado del vector de onda.
// Este es el núcleo O(N log N) del método PM.
//
// Algoritmo:
//   1. FFT directa:    ρ(x) → ρ̂(k)           [FFTW3: r2c]
//   2. Multiplicar:    φ̂(k) = -4πG ρ̂(k) / k²  [loop paralelo con OMP]
//   3. FFT inversa:    φ̂(k) → φ(x)            [FFTW3: c2r]
//   4. Normalizar:     φ(x) /= N_CELLS
//
// Nota sobre k=0: el modo DC (k=0) corresponde a la densidad media y se
// pone a cero (φ̂(0) = 0), lo que equivale a sustraer la densidad de fondo.
// =============================================================================

#include "types.hpp"

namespace pm {

// Resuelve ∇²φ = 4πGρ usando FFTW3.
// Lee state.density, escribe state.potential.
void solve_poisson(SimState& state);

} // namespace pm
