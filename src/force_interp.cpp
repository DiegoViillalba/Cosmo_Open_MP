// =============================================================================
// force_interp.cpp
// Interpolación trilineal (CIC inverso) de la fuerza a cada partícula.
//
// Esta operación es el adjunto exacto del depósito CIC:
//   - En CIC deposit: cada partícula distribuye masa a 8 celdas.
//   - En force_interp: cada partícula recoge fuerza de las mismas 8 celdas
//     con los mismos pesos trilineales.
//
// Paralelismo: solo lectura en las mallas de fuerza + escritura en
// accel[i] que es privado de cada iteración → zero race conditions.
// La paralelización con schedule(static) es óptima aquí porque todas
// las iteraciones tienen exactamente el mismo costo.
// =============================================================================

#include "force_interp.hpp"
#include "config.hpp"

#include <cmath>

namespace pm {

void interpolate_force(SimState& state,
                       std::vector<std::array<double,3>>& accel) {

    const std::size_t N_PART = state.particles.size();
    const int N = static_cast<int>(Ng);

    // Asegurar que el arreglo de aceleración tenga el tamaño correcto
    accel.resize(N_PART);

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t p = 0; p < N_PART; ++p) {

        const auto& part = state.particles[p];

        // Misma convención de posición que en cic_deposit (desplazamiento 0.5)
        const double px = part.pos[0] - 0.5;
        const double py = part.pos[1] - 0.5;
        const double pz = part.pos[2] - 0.5;

        // Celda base
        const int i0 = static_cast<int>(std::floor(px));
        const int j0 = static_cast<int>(std::floor(py));
        const int k0 = static_cast<int>(std::floor(pz));

        // Pesos trilineales (idénticos a los del depósito)
        const double dx = px - i0;
        const double dy = py - j0;
        const double dz = pz - k0;

        const double tx = 1.0 - dx;
        const double ty = 1.0 - dy;
        const double tz = 1.0 - dz;

        // Índices periódicos de las 8 celdas vecinas
        const int i0w = ((i0     % N) + N) % N;
        const int j0w = ((j0     % N) + N) % N;
        const int k0w = ((k0     % N) + N) % N;
        const int i1  = ((i0 + 1) % N + N) % N;
        const int j1  = ((j0 + 1) % N + N) % N;
        const int k1  = ((k0 + 1) % N + N) % N;

        // Macro para índice plano
        #define IDX(a,b,c) (static_cast<std::size_t>((a)*N*N + (b)*N + (c)))

        // Suma ponderada de las 8 celdas vecinas (una por componente de fuerza)
        // Los pesos son exactamente los mismos que en el depósito CIC.
        const double w000 = tx * ty * tz;
        const double w100 = dx * ty * tz;
        const double w010 = tx * dy * tz;
        const double w110 = dx * dy * tz;
        const double w001 = tx * ty * dz;
        const double w101 = dx * ty * dz;
        const double w011 = tx * dy * dz;
        const double w111 = dx * dy * dz;

        accel[p][0] =
            state.force_x[IDX(i0w,j0w,k0w)] * w000 +
            state.force_x[IDX(i1, j0w,k0w)] * w100 +
            state.force_x[IDX(i0w,j1, k0w)] * w010 +
            state.force_x[IDX(i1, j1, k0w)] * w110 +
            state.force_x[IDX(i0w,j0w,k1 )] * w001 +
            state.force_x[IDX(i1, j0w,k1 )] * w101 +
            state.force_x[IDX(i0w,j1, k1 )] * w011 +
            state.force_x[IDX(i1, j1, k1 )] * w111;

        accel[p][1] =
            state.force_y[IDX(i0w,j0w,k0w)] * w000 +
            state.force_y[IDX(i1, j0w,k0w)] * w100 +
            state.force_y[IDX(i0w,j1, k0w)] * w010 +
            state.force_y[IDX(i1, j1, k0w)] * w110 +
            state.force_y[IDX(i0w,j0w,k1 )] * w001 +
            state.force_y[IDX(i1, j0w,k1 )] * w101 +
            state.force_y[IDX(i0w,j1, k1 )] * w011 +
            state.force_y[IDX(i1, j1, k1 )] * w111;

        accel[p][2] =
            state.force_z[IDX(i0w,j0w,k0w)] * w000 +
            state.force_z[IDX(i1, j0w,k0w)] * w100 +
            state.force_z[IDX(i0w,j1, k0w)] * w010 +
            state.force_z[IDX(i1, j1, k0w)] * w110 +
            state.force_z[IDX(i0w,j0w,k1 )] * w001 +
            state.force_z[IDX(i1, j0w,k1 )] * w101 +
            state.force_z[IDX(i0w,j1, k1 )] * w011 +
            state.force_z[IDX(i1, j1, k1 )] * w111;

        #undef IDX
    }
}

} // namespace pm
