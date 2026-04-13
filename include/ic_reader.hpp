#pragma once

// =============================================================================
// ic_reader.hpp
// Lectura de condiciones iniciales desde archivo.
//
// Formato esperado del archivo ASCII:
//   - Primera línea: N_PART  (número de partículas, entero)
//   - Líneas siguientes: x y z vx vy vz  (6 doubles por línea)
//     donde x,y,z ∈ [0, Ng)  y  vx,vy,vz en unidades internas
//
// Si el archivo no existe o el número de partículas no coincide con
// N_PARTICLES en config.hpp, se lanza std::runtime_error.
//
// También se provee generate_uniform_ic() para generar condiciones iniciales
// simples en memoria (grilla regular + perturbaciones) cuando no se tiene
// un archivo externo disponible, útil para las pruebas de 128³.
// =============================================================================

#include <string>
#include "types.hpp"

namespace pm {

// Lee condiciones iniciales desde un archivo ASCII.
// Llena state.particles con posiciones y velocidades.
void read_ic(SimState& state, const std::string& filename);

// Genera condiciones iniciales en una grilla regular con perturbaciones
// gaussianas pequeñas. No requiere archivo externo.
// sigma_pos : desviación estándar del desplazamiento (en unidades de celda)
// sigma_vel : desviación estándar de la velocidad inicial
// seed      : semilla del generador aleatorio (reproducibilidad)
void generate_uniform_ic(SimState& state,
                         double sigma_pos = 0.5,
                         double sigma_vel = 0.1,
                         unsigned long seed = 42);

} // namespace pm
