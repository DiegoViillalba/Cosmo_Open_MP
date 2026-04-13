// =============================================================================
// test_integrator.cpp — Tests para el integrador Leap-Frog
// =============================================================================
#include "test_framework.hpp"
#include "types.hpp"
#include "integrator.hpp"
#include "config.hpp"
#include <cmath>
#include <array>
#include <vector>

// Test 1: MRU — sin fuerza x(t) = x0 + v0·N·dt (módulo Ng)
REGISTER_TEST(integrador_mru) {
    pm::SimState state;
    state.particles.resize(1);
    state.particles[0].pos[0] = 10.0; state.particles[0].pos[1] = 20.0; state.particles[0].pos[2] = 30.0;
    state.particles[0].vel[0] =  1.0; state.particles[0].vel[1] =  0.5; state.particles[0].vel[2] = -0.3;
    state.particles[0].mass   = 1.0 / pm::N_PARTICLES;

    std::vector<std::array<double,3>> accel(1);
    accel[0] = {0.0, 0.0, 0.0};

    const double ng = static_cast<double>(pm::Ng);
    const int    NS = 10;
    const double dt = pm::DT;
    const double x0 = 10.0, y0 = 20.0, z0 = 30.0;
    const double vx  =  1.0, vy  =  0.5, vz  = -0.3;

    for (int s = 0; s < NS; ++s) pm::leapfrog_step(state, accel);

    auto wrap = [&](double x, double v) {
        return std::fmod(x + v*dt*NS + ng*100.0, ng);
    };
    ASSERT_NEAR(state.particles[0].pos[0], wrap(x0, vx), 1e-10);
    ASSERT_NEAR(state.particles[0].pos[1], wrap(y0, vy), 1e-10);
    ASSERT_NEAR(state.particles[0].pos[2], wrap(z0, vz), 1e-10);
    return true;
}

// Test 2: periodicidad — partícula cruza borde derecho y reaparece a la izquierda
REGISTER_TEST(integrador_periodicidad) {
    pm::SimState state;
    state.particles.resize(1);
    const double ng = static_cast<double>(pm::Ng);
    const double dt = pm::DT;
    state.particles[0].pos[0] = ng - dt * 0.5;
    state.particles[0].pos[1] = state.particles[0].pos[2] = ng / 2.0;
    state.particles[0].vel[0] = 1.0;
    state.particles[0].vel[1] = state.particles[0].vel[2] = 0.0;
    state.particles[0].mass   = 1.0 / pm::N_PARTICLES;

    std::vector<std::array<double,3>> accel(1);
    accel[0] = {0.0, 0.0, 0.0};
    pm::leapfrog_step(state, accel);

    const double x = state.particles[0].pos[0];
    ASSERT_TRUE(x >= 0.0);
    ASSERT_TRUE(x < ng);
    ASSERT_NEAR(x, dt * 0.5, 1e-10);
    return true;
}

// Test 3: velocidad constante con a = 0 tras 50 pasos
REGISTER_TEST(integrador_vel_constante) {
    pm::SimState state;
    state.particles.resize(1);
    state.particles[0].vel[0] =  2.5; state.particles[0].vel[1] = -1.3; state.particles[0].vel[2] = 0.7;
    state.particles[0].pos[0] = state.particles[0].pos[1] = state.particles[0].pos[2] = 10.0;
    state.particles[0].mass   = 1.0 / pm::N_PARTICLES;

    std::vector<std::array<double,3>> accel(1);
    accel[0] = {0.0, 0.0, 0.0};
    for (int s = 0; s < 50; ++s) pm::leapfrog_step(state, accel);

    ASSERT_NEAR(state.particles[0].vel[0],  2.5, 1e-12);
    ASSERT_NEAR(state.particles[0].vel[1], -1.3, 1e-12);
    ASSERT_NEAR(state.particles[0].vel[2],  0.7, 1e-12);
    return true;
}

// Test 4: MUA — x = ½·a·t² exacto con Leap-Frog y half-kick inicial
REGISTER_TEST(integrador_mua) {
    pm::SimState state;
    state.particles.resize(1);
    state.particles[0].pos[0] = state.particles[0].pos[1] = state.particles[0].pos[2] = 0.0;
    state.particles[0].vel[0] = state.particles[0].vel[1] = state.particles[0].vel[2] = 0.0;
    state.particles[0].mass   = 1.0 / pm::N_PARTICLES;

    const double ax = 0.5;
    const int NS    = 10;
    const double dt = pm::DT;

    std::vector<std::array<double,3>> accel(1);
    accel[0] = {ax, 0.0, 0.0};

    pm::leapfrog_half_kick(state, accel);
    for (int s = 0; s < NS; ++s) pm::leapfrog_step(state, accel);

    const double t = NS * dt;
    ASSERT_NEAR(state.particles[0].pos[0], 0.5 * ax * t * t, 1e-10);
    return true;
}
