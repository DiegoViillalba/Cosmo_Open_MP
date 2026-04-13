// =============================================================================
// poisson_fft.cpp
// Resolución de ∇²φ = 4πGρ en espacio de Fourier usando FFTW3.
//
// Detalles de implementación con FFTW3:
//
//   r2c (real-to-complex): la FFT de un arreglo real de tamaño N tiene
//   simetría hermítica, así que FFTW solo almacena N/2+1 modos complejos
//   en la última dimensión. Para una malla Ng³:
//     - Entrada:  Ng × Ng × Ng  doubles
//     - Salida:   Ng × Ng × (Ng/2+1)  complex<double>
//
//   c2r (complex-to-real): operación inversa.
//
// Thread-safety de FFTW:
//   FFTW3 con -lfftw3_omp permite planificación y ejecución multi-hilo.
//   En esta implementación usamos fftw_plan_with_nthreads() para activarlo.
//   El loop de multiplicación por el kernel sí lo paralelizamos con OMP.
// =============================================================================

#include "poisson_fft.hpp"
#include "config.hpp"

#include <fftw3.h>
#include <cmath>
#include <complex>
#include <stdexcept>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace pm {

void solve_poisson(SimState& state) {

    const int N = static_cast<int>(Ng);
    const std::size_t real_size   = N * N * N;
    const std::size_t cplx_size   = N * N * (N/2 + 1);

    // ──────────────────────────────────────────────────────────────────
    // Asignar buffers alineados con fftw_malloc.
    // fftw_malloc garantiza alineación de 16 bytes para SIMD (SSE/AVX).
    // ──────────────────────────────────────────────────────────────────
    double*          in  = fftw_alloc_real(real_size);
    fftw_complex*    out = fftw_alloc_complex(cplx_size);

    if (!in || !out)
        throw std::runtime_error("poisson_fft: fallo en fftw_alloc");

    // ──────────────────────────────────────────────────────────────────
    // Nota sobre FFTW y OpenMP:
    // fftw3_omp (fftw_init_threads / fftw_plan_with_nthreads) requiere
    // el paquete libfftw3-bin compilado con soporte OMP, que no siempre
    // está disponible. En su lugar paralelizamos el loop de multiplicación
    // por el kernel con #pragma omp, que es donde está el trabajo real
    // para mallas grandes. Las propias llamadas fftw_execute son secuenciales
    // pero el bottleneck cosmológico está en el kernel, no en la FFT.
    // Para producción con Ng=512+ se puede linkar -lfftw3_omp y descomentar:
    //   fftw_init_threads(); fftw_plan_with_nthreads(omp_get_max_threads());
    // ──────────────────────────────────────────────────────────────────

    // ──────────────────────────────────────────────────────────────────
    // Crear planes de FFT.
    // FFTW_ESTIMATE: no mide, usa heurística. Más lento en ejecución
    // pero evita el overhead de planificación (importante en demos).
    // Para producción real se usaría FFTW_MEASURE con plan guardado.
    // ──────────────────────────────────────────────────────────────────
    fftw_plan plan_fwd = fftw_plan_dft_r2c_3d(N, N, N, in, out, FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_c2r_3d(N, N, N, out, in, FFTW_ESTIMATE);

    if (!plan_fwd || !plan_inv)
        throw std::runtime_error("poisson_fft: fallo al crear planes FFTW");

    // ──────────────────────────────────────────────────────────────────
    // Copiar densidad al buffer de entrada
    // ──────────────────────────────────────────────────────────────────
    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < real_size; ++i)
        in[i] = state.density[i];

    // ── FFT directa: ρ(x) → ρ̂(k) ──────────────────────────────────
    fftw_execute(plan_fwd);

    // ──────────────────────────────────────────────────────────────────
    // Multiplicar por el kernel de Green en espacio de Fourier:
    //   φ̂(k) = -4πG ρ̂(k) / k²
    //
    // Los vectores de onda discretos son:
    //   kx = 2π/N × (ix si ix ≤ N/2, sino ix-N)
    //   ky = 2π/N × (iy si iy ≤ N/2, sino iy-N)
    //   kz = 2π/N × iz   (solo 0..N/2+1 por simetría r2c)
    //
    // La división de Δx (= L_BOX/N) entra en la definición de k:
    //   k_phys = k_discrete / Δx
    // ──────────────────────────────────────────────────────────────────
    const double two_pi_over_N = 2.0 * M_PI / N;
    const double G4pi = 4.0 * M_PI * G_GRAV;
    const double dx_inv2 = 1.0 / (CELL_SIZE * CELL_SIZE);  // 1/Δx²
    const int    Nhalf = N/2 + 1;

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static) collapse(2)
    #endif
    for (int ix = 0; ix < N; ++ix) {
        for (int iy = 0; iy < N; ++iy) {
            for (int iz = 0; iz < Nhalf; ++iz) {

                // Frecuencias discretas (con aliasing para ix,iy > N/2)
                const double kx = two_pi_over_N * (ix <= N/2 ? ix : ix - N);
                const double ky = two_pi_over_N * (iy <= N/2 ? iy : iy - N);
                const double kz = two_pi_over_N * iz;

                // k² en unidades físicas
                const double k2 = (kx*kx + ky*ky + kz*kz) * dx_inv2;

                const std::size_t flat = ix * (N * Nhalf) + iy * Nhalf + iz;

                if (k2 < 1e-30) {
                    // Modo DC (k=0): poner a cero para sustraer densidad media
                    out[flat][0] = 0.0;
                    out[flat][1] = 0.0;
                } else {
                    // φ̂(k) = -4πG ρ̂(k) / k²
                    const double factor = -G4pi / k2;
                    out[flat][0] *= factor;
                    out[flat][1] *= factor;
                }
            }
        }
    }

    // ── FFT inversa: φ̂(k) → φ(x) ──────────────────────────────────
    fftw_execute(plan_inv);

    // ──────────────────────────────────────────────────────────────────
    // Normalizar: FFTW no normaliza la FFT inversa.
    // El factor es 1/N_CELLS para que la IFFT(FFT(f)) = f.
    // ──────────────────────────────────────────────────────────────────
    const double norm = 1.0 / static_cast<double>(real_size);

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < real_size; ++i)
        state.potential[i] = in[i] * norm;

    // Liberar recursos FFTW
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_inv);
    fftw_free(in);
    fftw_free(out);
}

} // namespace pm
