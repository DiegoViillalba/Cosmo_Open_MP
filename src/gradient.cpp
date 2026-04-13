// =============================================================================
// gradient.cpp
// Gradiente del potencial gravitacional por diferencias finitas centradas.
//
// Calcula g = -∇φ en cada celda de la malla usando el esquema de 2do orden:
//   g_x(i,j,k) = -(φ(i+1,j,k) - φ(i-1,j,k)) / (2·Δx)
//
// Paralelismo: loop 3D sobre la malla → embarrassingly parallel.
// Cada iteración escribe en (force_x[idx], force_y[idx], force_z[idx])
// con idx único → cero condiciones de carrera.
//
// El collapse(3) en el pragma OMP colapsa los tres loops anidados en uno
// solo antes de distribuir, lo que da más trabajo por hilo y reduce el
// overhead de scheduling para mallas grandes.
// =============================================================================

#include "gradient.hpp"
#include "config.hpp"

namespace pm {

void compute_gradient(SimState& state) {

    const int  N    = static_cast<int>(Ng);
    const double inv2dx = 1.0 / (2.0 * CELL_SIZE);  // 1/(2Δx)

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static) collapse(3)
    #endif
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < N; ++iz) {

                // Índice plano de la celda actual
                const std::size_t idx = static_cast<std::size_t>(
                    ix * N * N + iy * N + iz);

                // Vecinos periódicos en cada dirección
                // (N + x - 1) % N    → vecino negativo
                // (x + 1) % N        → vecino positivo
                const int ixm = (N + ix - 1) % N;
                const int ixp = (ix + 1) % N;
                const int iym = (N + iy - 1) % N;
                const int iyp = (iy + 1) % N;
                const int izm = (N + iz - 1) % N;
                const int izp = (iz + 1) % N;

                // Potencial en los 6 vecinos
                const double phi_xm = state.potential[ixm * N * N + iy  * N + iz];
                const double phi_xp = state.potential[ixp * N * N + iy  * N + iz];
                const double phi_ym = state.potential[ix  * N * N + iym * N + iz];
                const double phi_yp = state.potential[ix  * N * N + iyp * N + iz];
                const double phi_zm = state.potential[ix  * N * N + iy  * N + izm];
                const double phi_zp = state.potential[ix  * N * N + iy  * N + izp];

                // Diferencias centradas: g = -∇φ → signo negativo
                state.force_x[idx] = -(phi_xp - phi_xm) * inv2dx;
                state.force_y[idx] = -(phi_yp - phi_ym) * inv2dx;
                state.force_z[idx] = -(phi_zp - phi_zm) * inv2dx;
            }
        }
    }
}

} // namespace pm
