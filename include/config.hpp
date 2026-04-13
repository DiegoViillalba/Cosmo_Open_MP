#pragma once

// =============================================================================
// config.hpp
// Parámetros globales de la simulación Particle-Mesh.
//
// Todos los valores centrales se concentran aquí para facilitar el cambio
// entre el caso de prueba (128³) y el caso de producción (256³) con una
// sola línea de modificación: #define PM_GRID_SIZE.
//
// Convención de unidades: unidades cosmológicas internas (box = 1.0)
//   - Longitud : fracción de L_box  (posiciones en [0, 1))
//   - Velocidad: km/s (comóviles)
//   - Masa     : fracción de M_total
// =============================================================================

#include <cstddef>   // size_t
#include <cmath>     // M_PI

namespace pm {

// ─────────────────────────────────────────────────────────────
// Tamaño de malla y número de partículas
//
// PM_GRID_SIZE controla el modo de operación:
//   128 → prueba rápida (~2M partículas, ~16 MB de malla)
//   256 → producción   (~16M partículas, ~128 MB de malla)
//
// Para cambiar de modo: recompila pasando -DPM_GRID_SIZE=128
// o modifica la línea de abajo.
// ─────────────────────────────────────────────────────────────
#ifndef PM_GRID_SIZE
#define PM_GRID_SIZE 256
#endif

constexpr std::size_t Ng = PM_GRID_SIZE;      // puntos de malla por dimensión
constexpr std::size_t Np = Ng;                // partículas por dimensión
constexpr std::size_t N_PARTICLES = Np*Np*Np; // total de partículas
constexpr std::size_t N_CELLS     = Ng*Ng*Ng; // total de celdas en la malla

// ─────────────────────────────────────────────────────────────
// Geometría de la caja de simulación
// L_box: tamaño en Mpc/h (solo relevante para unidades físicas;
//        internamente todo trabaja en [0, 1))
// ─────────────────────────────────────────────────────────────
constexpr double L_BOX   = 100.0;             // Mpc/h
constexpr double CELL_SIZE = L_BOX / Ng;      // tamaño de celda en Mpc/h

// ─────────────────────────────────────────────────────────────
// Parámetros de integración temporal
// DT    : paso de tiempo (en unidades internas, ≈ factor de escala)
// N_STEPS: número de pasos totales a correr
// ─────────────────────────────────────────────────────────────
constexpr double DT      = 0.01;
constexpr int    N_STEPS = 100;

// ─────────────────────────────────────────────────────────────
// Constante gravitacional G en unidades internas
// En unidades cosmológicas convenientes, G absorbe factores de
// H_0 y rho_crit. Aquí usamos G=1 para simplificar la demo.
// ─────────────────────────────────────────────────────────────
constexpr double G_GRAV  = 1.0;

// ─────────────────────────────────────────────────────────────
// Salida de datos
// SNAP_INTERVAL: cada cuántos pasos se escribe un snapshot
// ─────────────────────────────────────────────────────────────
constexpr int SNAP_INTERVAL = 10;

// ─────────────────────────────────────────────────────────────
// Softening gravitacional
// Evita divergencias cuando dos partículas quedan muy cerca.
// Típicamente ~ 1/50 del tamaño de celda.
// ─────────────────────────────────────────────────────────────
constexpr double SOFTENING = CELL_SIZE / 50.0;

// ─────────────────────────────────────────────────────────────
// Macro de paralelismo condicional
// Cuando se compila pm_serial (PM_SERIAL definido), las regiones
// OpenMP se convierten en código secuencial sin necesidad de
// tener dos versiones del algoritmo.
// ─────────────────────────────────────────────────────────────
#ifdef PM_SERIAL
  #define OMP_PARALLEL_FOR
  #define OMP_PARALLEL_FOR_REDUCTION(op, var)
  #define OMP_ATOMIC
  #define OMP_CRITICAL(name)
#else
  #define OMP_PARALLEL_FOR        _Pragma("omp parallel for schedule(static)")
  #define OMP_PARALLEL_FOR_REDUCTION(op, var) \
          _Pragma("omp parallel for reduction(" #op ":" #var ")")
  #define OMP_ATOMIC              _Pragma("omp atomic")
  #define OMP_CRITICAL(name)      _Pragma("omp critical(" #name ")")
#endif

} // namespace pm
