#pragma once

// =============================================================================
// integrator.hpp
// Integración de las ecuaciones de movimiento con el método Leap-Frog.
//
// El integrador Leap-Frog es el estándar en simulaciones de N-cuerpos por
// su conservación exacta del volumen en el espacio de fase (simpléctico) y
// su segunda orden de precisión con solo un cálculo de fuerza por paso.
//
// Esquema (versión "kick-drift-kick" o simplemente "leapfrog"):
//
//   v(t + dt/2) = v(t - dt/2) + a(t) · dt     [kick: actualiza velocidad]
//   x(t + dt)   = x(t) + v(t + dt/2) · dt     [drift: actualiza posición]
//
// Las velocidades viven a dt/2 desfasadas de las posiciones. Esto es una
// característica del método, no un error.
//
// Condiciones de frontera periódicas: si una partícula sale de [0, Ng)
// se reinserta por el lado opuesto (wrap).
//
// Complejidad: O(N_PART) → embarrasingly parallel por partícula.
// No hay dependencias entre partículas en este paso (solo leen accel[i]
// que ya fue calculado, y escriben en su propia partícula).
// =============================================================================

#include <array>
#include "types.hpp"

namespace pm {

// Aplica un paso Leap-Frog a todas las partículas.
// accel[i] = {ax, ay, az} para la partícula i (calculado por interpolate_force)
// DT se lee de config.hpp
void leapfrog_step(SimState& state,
                   const std::vector<std::array<double,3>>& accel);

// Inicializa las velocidades al primer medio-paso (solo se llama una vez
// al inicio de la simulación, antes del primer step completo).
void leapfrog_half_kick(SimState& state,
                        const std::vector<std::array<double,3>>& accel);

} // namespace pm
