// =============================================================================
// cic.cpp
// Cloud-in-Cell: asigna masa de partículas a la malla de densidad.
//
// PROBLEMA DE PARALELISMO:
// El loop sobre partículas tiene una condición de carrera: dos partículas
// pueden caer en las mismas 8 celdas vecinas y sus hilos intentarán
// incrementar el mismo density[idx] al mismo tiempo → resultado incorrecto.
//
// SOLUCIÓN ADOPTADA: acumuladores locales por hilo
//   - Cada hilo tiene su propia copia del arreglo de densidad.
//   - Al terminar el loop, se suman todas las copias en la malla global.
//   - Costo en memoria: n_threads × Ng³ × 8 bytes.
//     Para 8 hilos y Ng=256: 8 × 16M × 8B ≈ 1 GB → demasiado.
//
// SOLUCIÓN ALTERNATIVA (usada aquí): atomic en las 8 escrituras.
//   - Más sencillo de implementar y sin uso extra de memoria.
//   - Funciona bien cuando hay pocas colisiones (partículas distribuidas).
//   - Para Ng=256 y N=16M: cada celda recibe en promedio N/Ng³ = 1 partícula,
//     así que las colisiones son raras → atomic tiene overhead mínimo.
//
// Para demostración en la exposición, se incluye también una versión
// con critical block (más lenta, pero más fácil de entender conceptualmente).
// =============================================================================

#include "cic.hpp"
#include "config.hpp"

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <cmath>
#include <cstddef>

namespace pm {

void cic_deposit(SimState& state) {

    // Limpiar la malla de densidad antes de cada depósito
    state.density.zero();

    const std::size_t N = state.particles.size();
    const int ni = static_cast<int>(Ng);

    // ──────────────────────────────────────────────────────────────────
    // Loop paralelo con atomic en las 8 actualizaciones por partícula.
    //
    // schedule(dynamic, 512): el trabajo por partícula es uniforme, pero
    // usamos dynamic con chunks grandes para balancear el overhead de atomic.
    // ──────────────────────────────────────────────────────────────────
    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(dynamic, 512)
    #endif
    for (std::size_t p = 0; p < N; ++p) {

        const auto& part = state.particles[p];
        const double m   = part.mass;

        // Posición en unidades de celda, desplazada 0.5 para alinear con
        // el centro de las celdas (convención CIC estándar)
        const double px = part.pos[0] - 0.5;
        const double py = part.pos[1] - 0.5;
        const double pz = part.pos[2] - 0.5;

        // Celda base (esquina inferior del cubo de 8 vecinos)
        const int i0 = static_cast<int>(std::floor(px));
        const int j0 = static_cast<int>(std::floor(py));
        const int k0 = static_cast<int>(std::floor(pz));

        // Pesos trilineales: distancia fraccionaria a la celda base
        const double dx = px - i0;   // ∈ [0, 1)
        const double dy = py - j0;
        const double dz = pz - k0;

        // Pesos complementarios
        const double tx = 1.0 - dx;
        const double ty = 1.0 - dy;
        const double tz = 1.0 - dz;

        // Índices de las 8 celdas vecinas (con condiciones periódicas)
        const int i1 = ((i0 + 1) % ni + ni) % ni;
        const int j1 = ((j0 + 1) % ni + ni) % ni;
        const int k1 = ((k0 + 1) % ni + ni) % ni;
        const int i0w = (i0     % ni + ni) % ni;
        const int j0w = (j0     % ni + ni) % ni;
        const int k0w = (k0     % ni + ni) % ni;

        // Macros para el índice plano de una celda (ix, iy, iz)
        #define IDX(a,b,c) ((a)*Ng*Ng + (b)*Ng + (c))

        // Depositar masa en las 8 celdas con pesos trilineales.
        // Las operaciones atomic garantizan escrituras seguras en paralelo.
        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i0w,j0w,k0w)] += m * tx * ty * tz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i1, j0w,k0w)] += m * dx * ty * tz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i0w,j1, k0w)] += m * tx * dy * tz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i1, j1, k0w)] += m * dx * dy * tz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i0w,j0w,k1 )] += m * tx * ty * dz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i1, j0w,k1 )] += m * dx * ty * dz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i0w,j1, k1 )] += m * tx * dy * dz;

        #ifndef PM_SERIAL
          #pragma omp atomic
        #endif
        state.density[IDX(i1, j1, k1 )] += m * dx * dy * dz;

        #undef IDX
    }

    // ──────────────────────────────────────────────────────────────────
    // Normalización: convertir de masa total a contraste de densidad δ.
    //   ρ_norm(x) = ρ(x) / ρ_mean - 1
    // donde ρ_mean = suma_total / N_CELLS = N_PART × m_part / N_CELLS.
    //
    // En nuestro caso m_part = 1/N_PART y la suma total = 1,
    // así que ρ_mean = 1/N_CELLS → multiplicamos por N_CELLS.
    //
    // Esto garantiza que la densidad de fondo sea 1 y los vacíos sean < 1.
    // ──────────────────────────────────────────────────────────────────
    const double norm = static_cast<double>(N_CELLS);

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < N_CELLS; ++i) {
        state.density[i] *= norm;
    }
}

} // namespace pm
