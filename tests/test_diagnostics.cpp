// =============================================================================
// test_diagnostics.cpp — Tests para el módulo de diagnósticos
// =============================================================================
#include "test_framework.hpp"
#include "types.hpp"
#include "diagnostics.hpp"
#include "config.hpp"
#include <cmath>
#include <vector>

// Test 1: Ek = ½·m·v² para partícula única
REGISTER_TEST(diag_ek_particula_unica) {
    pm::SimState state;
    state.particles.resize(1);
    state.particles[0].vel[0] = 3.0; state.particles[0].vel[1] = 4.0; state.particles[0].vel[2] = 0.0;
    state.particles[0].pos[0] = state.particles[0].pos[1] = state.particles[0].pos[2] = 1.0;
    state.particles[0].mass   = 2.0;
    // Ek = ½ × 2 × (9 + 16 + 0) = 25.0
    auto d = pm::compute_diagnostics(state);
    ASSERT_NEAR(d.kinetic_energy, 25.0, 1e-10);
    return true;
}

// Test 2: momento cero con velocidades opuestas iguales
REGISTER_TEST(diag_momento_cero) {
    pm::SimState state;
    state.particles.resize(2);
    const double m = 0.5;
    state.particles[0].mass = m;
    state.particles[0].vel[0] =  2.0; state.particles[0].vel[1] = -1.0; state.particles[0].vel[2] =  3.0;
    state.particles[0].pos[0] = state.particles[0].pos[1] = state.particles[0].pos[2] = 5.0;
    state.particles[1].mass = m;
    state.particles[1].vel[0] = -2.0; state.particles[1].vel[1] =  1.0; state.particles[1].vel[2] = -3.0;
    state.particles[1].pos[0] = state.particles[1].pos[1] = state.particles[1].pos[2] = 5.0;
    auto d = pm::compute_diagnostics(state);
    ASSERT_NEAR(d.momentum[0], 0.0, 1e-12);
    ASSERT_NEAR(d.momentum[1], 0.0, 1e-12);
    ASSERT_NEAR(d.momentum[2], 0.0, 1e-12);
    return true;
}

// Test 3: todas en reposo → Ek = 0, |P| = 0
REGISTER_TEST(diag_reposo) {
    pm::SimState state;
    state.particles.resize(pm::N_PARTICLES);
    const double m = 1.0 / pm::N_PARTICLES;
    for (auto& p : state.particles) {
        p.mass = m;
        p.vel[0] = p.vel[1] = p.vel[2] = 0.0;
        p.pos[0] = p.pos[1] = p.pos[2] = 1.0;
    }
    auto d = pm::compute_diagnostics(state);
    ASSERT_NEAR(d.kinetic_energy, 0.0, 1e-12);
    ASSERT_NEAR(d.momentum[0],   0.0, 1e-12);
    return true;
}

// Test 4: Etot == Ek + Ep
REGISTER_TEST(diag_energia_total_consistente) {
    pm::SimState state;
    state.particles.resize(20);
    for (std::size_t i = 0; i < 20; ++i) {
        state.particles[i].mass   = 0.05;
        state.particles[i].vel[0] = static_cast<double>(i) * 0.1;
        state.particles[i].vel[1] = state.particles[i].vel[2] = 0.0;
        state.particles[i].pos[0] = static_cast<double>(i % pm::Ng) + 0.5;
        state.particles[i].pos[1] = state.particles[i].pos[2] = 1.0;
    }
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        state.potential[i] = 0.01 * std::sin(static_cast<double>(i) * 0.1);
    auto d = pm::compute_diagnostics(state);
    ASSERT_NEAR(d.total_energy, d.kinetic_energy + d.potential_energy, 1e-12);
    return true;
}

// Test 5: Ek escala linealmente con masa (verifica reduction OMP)
REGISTER_TEST(diag_ek_escala_masa) {
    const std::size_t N = 500;
    pm::SimState s1, s2;
    s1.particles.resize(N); s2.particles.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
        const double m = 1.0 / N;
        s1.particles[i].mass   = m;
        s1.particles[i].vel[0] = std::sin(static_cast<double>(i) * 0.1);
        s1.particles[i].vel[1] = std::cos(static_cast<double>(i) * 0.17);
        s1.particles[i].vel[2] = std::sin(static_cast<double>(i) * 0.23 + 1.0);
        s1.particles[i].pos[0] = std::fmod(static_cast<double>(i), static_cast<double>(pm::Ng)) + 0.5;
        s1.particles[i].pos[1] = s1.particles[i].pos[2] = 1.0;
        s2.particles[i]        = s1.particles[i];
        s2.particles[i].mass   = 2.0 * m;
    }
    auto d1 = pm::compute_diagnostics(s1);
    auto d2 = pm::compute_diagnostics(s2);
    ASSERT_NEAR(d2.kinetic_energy, 2.0 * d1.kinetic_energy, 1e-10);
    ASSERT_NEAR(d2.momentum[0],   2.0 * d1.momentum[0],    1e-12);
    return true;
}
