#pragma once

// =============================================================================
// cic.hpp
// Cloud-in-Cell (CIC): deposita la masa de las partículas en la malla.
//
// En el método PM la malla discreta de densidad ρ_{ijk} se construye
// asignando la masa de cada partícula a las 8 celdas vecinas según su
// posición fraccionaria (trilinear weighting). Esto equivale a tratar cada
// partícula como un cubo de lado Δx centrado en su posición.
//
// Complejidad: O(N_PART) → perfectamente paralelizable.
//
// El desafío de paralelismo aquí es la condición de carrera: múltiples hilos
// pueden intentar actualizar la misma celda de la malla simultáneamente.
// Se resuelve con acumuladores locales por hilo y una reducción final.
// =============================================================================

#include "types.hpp"

namespace pm {

// Deposita masa de todas las partículas en state.density.
// Internamente usa OMP con acumuladores locales para evitar race conditions.
// Normaliza al final: suma total de density = N_PARTICLES (densidad media = 1)
void cic_deposit(SimState& state);

} // namespace pm
