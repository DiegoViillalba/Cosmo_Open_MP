// =============================================================================
// test_cic.cpp — Tests para el módulo Cloud-in-Cell
// =============================================================================
#include "test_framework.hpp"
#include "types.hpp"
#include "cic.hpp"
#include "config.hpp"
#include <cmath>
#include <algorithm>

static pm::SimState single_particle(double x, double y, double z) {
    pm::SimState s;
    s.particles.resize(1);
    s.particles[0].pos[0] = x; s.particles[0].pos[1] = y; s.particles[0].pos[2] = z;
    s.particles[0].vel[0] = s.particles[0].vel[1] = s.particles[0].vel[2] = 0;
    s.particles[0].mass   = 1.0;
    return s;
}

// Test 1: suma(density) == N_CELLS después de normalizar
REGISTER_TEST(cic_conservacion_masa) {
    pm::SimState state;
    state.particles.resize(pm::N_PARTICLES);
    const double m = 1.0 / pm::N_PARTICLES;
    for (auto& p : state.particles) {
        p.mass = m; p.pos[0] = p.pos[1] = p.pos[2] = 1.0;
        p.vel[0] = p.vel[1] = p.vel[2] = 0;
    }
    pm::cic_deposit(state);

    double suma = 0;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i) suma += state.density[i];
    ASSERT_NEAR(suma, static_cast<double>(pm::N_CELLS), 1e-6 * pm::N_CELLS);
    return true;
}

// Test 2: partícula en borde → densidad igual en las dos celdas adyacentes
REGISTER_TEST(cic_split_50_50) {
    pm::SimState state = single_particle(5.0, 1.5, 1.5);
    pm::cic_deposit(state);

    const int N = static_cast<int>(pm::Ng);
    double d0 = state.density[4 * N * N + 1 * N + 1];
    double d1 = state.density[5 * N * N + 1 * N + 1];
    ASSERT_GT(d0, 0.0);
    ASSERT_NEAR(d0, d1, 1e-10);
    return true;
}

// Test 3: ningún valor de densidad es negativo
REGISTER_TEST(cic_no_negativos) {
    pm::SimState state;
    state.particles.resize(200);
    const double m = 1.0 / 200.0;
    for (std::size_t i = 0; i < 200; ++i) {
        state.particles[i].mass   = m;
        state.particles[i].pos[0] = std::fmod(i * 7.31,  static_cast<double>(pm::Ng));
        state.particles[i].pos[1] = std::fmod(i * 3.71,  static_cast<double>(pm::Ng));
        state.particles[i].pos[2] = std::fmod(i * 11.13, static_cast<double>(pm::Ng));
        state.particles[i].vel[0] = state.particles[i].vel[1]
                                  = state.particles[i].vel[2] = 0;
    }
    pm::cic_deposit(state);
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        ASSERT_TRUE(state.density[i] >= -1e-12);
    return true;
}

// Test 4: grilla regular → densidad uniforme ≈ 1
REGISTER_TEST(cic_grilla_uniforme) {
    pm::SimState state;
    state.particles.resize(pm::N_PARTICLES);
    const double spacing = static_cast<double>(pm::Ng) / pm::Np;
    const double m = 1.0 / pm::N_PARTICLES;
    std::size_t idx = 0;
    for (std::size_t ix = 0; ix < pm::Np; ++ix)
        for (std::size_t iy = 0; iy < pm::Np; ++iy)
            for (std::size_t iz = 0; iz < pm::Np; ++iz) {
                auto& p = state.particles[idx++];
                p.pos[0] = (ix + 0.5) * spacing;
                p.pos[1] = (iy + 0.5) * spacing;
                p.pos[2] = (iz + 0.5) * spacing;
                p.vel[0] = p.vel[1] = p.vel[2] = 0;
                p.mass   = m;
            }
    pm::cic_deposit(state);

    double min_d = 1e30, max_d = -1e30;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i) {
        min_d = std::min(min_d, state.density[i]);
        max_d = std::max(max_d, state.density[i]);
    }
    ASSERT_NEAR(min_d, 1.0, 0.01);
    ASSERT_NEAR(max_d, 1.0, 0.01);
    return true;
}
