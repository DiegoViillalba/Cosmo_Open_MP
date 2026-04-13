// =============================================================================
// test_poisson.cpp — Tests para el solucionador de Poisson con FFT
// =============================================================================
#include "test_framework.hpp"
#include "types.hpp"
#include "poisson_fft.hpp"
#include "config.hpp"
#include <cmath>
#include <vector>

// Test 1: densidad constante → potencial ≈ 0 (modo DC = 0)
REGISTER_TEST(poisson_modo_dc_cero) {
    pm::SimState state;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i) state.density[i] = 1.0;
    pm::solve_poisson(state);
    double max_phi = 0;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        max_phi = std::max(max_phi, std::abs(state.potential[i]));
    ASSERT_NEAR(max_phi, 0.0, 1e-10);
    return true;
}

// Test 2: φ ∝ cos(2πx/N) cuando ρ ∝ cos(2πx/N)
REGISTER_TEST(poisson_forma_coseno) {
    pm::SimState state;
    const int N = static_cast<int>(pm::Ng);
    for (int ix = 0; ix < N; ++ix)
        for (int iy = 0; iy < N; ++iy)
            for (int iz = 0; iz < N; ++iz)
                state.density[ix*N*N + iy*N + iz] = 1.0 + std::cos(2.0*M_PI*ix/N);
    pm::solve_poisson(state);

    std::vector<double> prof(N, 0.0);
    for (int ix = 0; ix < N; ++ix)
        for (int iy = 0; iy < N; ++iy)
            for (int iz = 0; iz < N; ++iz)
                prof[ix] += state.potential[ix*N*N + iy*N + iz];
    for (int ix = 0; ix < N; ++ix) prof[ix] /= (N * N);

    double num = 0, den = 0;
    for (int ix = 0; ix < N; ++ix) {
        double c = std::cos(2.0*M_PI*ix/N);
        num += prof[ix] * c;
        den += prof[ix] * prof[ix];
    }
    double corr = std::abs(num) / (std::sqrt(den) * std::sqrt(N / 2.0) + 1e-30);
    ASSERT_GT(corr, 0.99);
    return true;
}

// Test 3: no NaN / Inf para densidad arbitraria
REGISTER_TEST(poisson_no_nan_inf) {
    pm::SimState state;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        state.density[i] = 1.0 + 0.3 * std::sin(static_cast<double>(i) * 0.007);
    pm::solve_poisson(state);
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        ASSERT_TRUE(std::isfinite(state.potential[i]));
    return true;
}

// Test 4: linealidad φ(α·ρ) = α·φ(ρ)
REGISTER_TEST(poisson_linealidad) {
    const double alpha = 3.7;
    pm::SimState s1;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        s1.density[i] = 1.0 + 0.5 * std::cos(static_cast<double>(i) * 0.05);
    pm::solve_poisson(s1);

    pm::SimState s2;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        s2.density[i] = alpha * s1.density[i];
    pm::solve_poisson(s2);

    double max_err = 0;
    for (std::size_t i = 0; i < pm::N_CELLS; ++i)
        max_err = std::max(max_err, std::abs(s2.potential[i] - alpha * s1.potential[i]));
    ASSERT_NEAR(max_err, 0.0, 1e-8);
    return true;
}
