// =============================================================================
// ic_reader.cpp
// Implementación de la lectura y generación de condiciones iniciales.
// =============================================================================

#include "ic_reader.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <random>
#include <cmath>
#include <iostream>

namespace pm {

// -----------------------------------------------------------------------------
// read_ic
// Lee condiciones iniciales desde archivo ASCII.
//
// Formato esperado:
//   Línea 1  : N_PART  (entero)
//   Líneas 2…: x y z vx vy vz  (en unidades de celda y vel. interna)
// -----------------------------------------------------------------------------
void read_ic(SimState& state, const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("ic_reader: no se puede abrir '" + filename + "'");

    std::size_t n_part;
    f >> n_part;

    if (n_part != N_PARTICLES) {
        std::ostringstream oss;
        oss << "ic_reader: el archivo tiene " << n_part
            << " partículas, pero la simulación espera " << N_PARTICLES
            << " (Ng=" << Ng << ").\n"
            << "  Solución: recompila con -DPM_GRID_SIZE=N adecuado, "
            << "o genera el IC con el Ng correcto.";
        throw std::runtime_error(oss.str());
    }

    state.particles.resize(n_part);
    const double mass_per_part = 1.0 / static_cast<double>(n_part);

    for (std::size_t i = 0; i < n_part; ++i) {
        auto& p = state.particles[i];
        f >> p.pos[0] >> p.pos[1] >> p.pos[2]
          >> p.vel[0] >> p.vel[1] >> p.vel[2];
        p.mass = mass_per_part;

        if (f.fail())
            throw std::runtime_error(
                "ic_reader: error de lectura en partícula " + std::to_string(i));
    }

    std::cout << "[ic_reader] Leídas " << n_part
              << " partículas desde '" << filename << "'\n";
}

// -----------------------------------------------------------------------------
// generate_uniform_ic
// Genera condiciones iniciales como grilla regular perturbada.
//
// Coloca las partículas en una grilla cúbica de Np^3 puntos, con
// desplazamientos gaussianos y velocidades gaussianas pequeñas.
// Útil para pruebas sin necesidad de archivo externo.
// -----------------------------------------------------------------------------
void generate_uniform_ic(SimState& state,
                         double sigma_pos,
                         double sigma_vel,
                         unsigned long seed) {

    state.particles.resize(N_PARTICLES);
    const double mass_per_part = 1.0 / static_cast<double>(N_PARTICLES);
    const double spacing = static_cast<double>(Ng) / static_cast<double>(Np);

    // Implementación de ICs a partir de desplazamiento
    // double redshift_initial = 50.0;
    // double H0 = 67.4; // Planck collaboration 2018
    
    // double a_initial = 1.0 / (1.0 + redshift_initial);

    // Generador de números aleatorios con distribución normal
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> dist_pos(0.0, sigma_pos);
    std::normal_distribution<double> dist_vel(0.0, sigma_vel);

    std::size_t idx = 0;
    for (std::size_t ix = 0; ix < Np; ++ix) {
        for (std::size_t iy = 0; iy < Np; ++iy) {
            for (std::size_t iz = 0; iz < Np; ++iz) {

                auto& p = state.particles[idx++];

                // Posición inicial en la grilla
                double x = (ix + 0.5) * spacing + dist_pos(rng);
                double y = (iy + 0.5) * spacing + dist_pos(rng);
                double z = (iz + 0.5) * spacing + dist_pos(rng);

                // Aplicar condiciones periódicas
                double ng = static_cast<double>(Ng);
                p.pos[0] = std::fmod(x + ng * 10.0, ng);
                p.pos[1] = std::fmod(y + ng * 10.0, ng);
                p.pos[2] = std::fmod(z + ng * 10.0, ng);

                // Velocidades gaussianas pequeñas
                p.vel[0] = dist_vel(rng);
                p.vel[1] = dist_vel(rng);
                p.vel[2] = dist_vel(rng);

                p.mass = mass_per_part;
            }
        }
    }

    std::cout << "[ic_reader] Generadas " << N_PARTICLES
              << " partículas en grilla uniforme perturbada"
              << " (Ng=" << Ng << ", seed=" << seed << ")\n";
}

} // namespace pm
