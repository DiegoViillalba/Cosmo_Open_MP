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
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace pm {

namespace {

class FFTWPlanCache {
public:
    FFTWPlanCache()
        : n_(static_cast<int>(Ng)),
          real_size_(static_cast<std::size_t>(n_) * n_ * n_),
          cplx_size_(static_cast<std::size_t>(n_) * n_ * (n_ / 2 + 1)),
          in_(nullptr),
          out_(nullptr),
          plan_fwd_(nullptr),
          plan_inv_(nullptr),
          planned_threads_(0),
          threads_initialized_(false) {
        in_ = fftw_alloc_real(real_size_);
        out_ = fftw_alloc_complex(cplx_size_);
        if (!in_ || !out_)
            throw std::runtime_error("poisson_fft: fallo en fftw_alloc");
    }

    ~FFTWPlanCache() {
        destroy_plans();
        if (in_) fftw_free(in_);
        if (out_) fftw_free(out_);
        if (threads_initialized_)
            fftw_cleanup_threads();
    }

    void ensure_plans_for_threads(int nthreads) {
        const int safe_threads = std::max(1, nthreads);

        if (!threads_initialized_) {
            if (fftw_init_threads() == 0)
                throw std::runtime_error("poisson_fft: fftw_init_threads fallo");
            threads_initialized_ = true;
        }

        if (plan_fwd_ && plan_inv_ && planned_threads_ == safe_threads)
            return;

        destroy_plans();
        fftw_plan_with_nthreads(safe_threads);

        // Reusar planes amortiza el costo de planning en todos los pasos.
        plan_fwd_ = fftw_plan_dft_r2c_3d(n_, n_, n_, in_, out_, FFTW_MEASURE);
        plan_inv_ = fftw_plan_dft_c2r_3d(n_, n_, n_, out_, in_, FFTW_MEASURE);
        if (!plan_fwd_ || !plan_inv_)
            throw std::runtime_error("poisson_fft: fallo al crear planes FFTW");

        planned_threads_ = safe_threads;
    }

    double* in() { return in_; }
    fftw_complex* out() { return out_; }
    fftw_plan forward_plan() const { return plan_fwd_; }
    fftw_plan inverse_plan() const { return plan_inv_; }
    std::size_t real_size() const { return real_size_; }

private:
    void destroy_plans() {
        if (plan_fwd_) {
            fftw_destroy_plan(plan_fwd_);
            plan_fwd_ = nullptr;
        }
        if (plan_inv_) {
            fftw_destroy_plan(plan_inv_);
            plan_inv_ = nullptr;
        }
    }

    int n_;
    std::size_t real_size_;
    std::size_t cplx_size_;
    double* in_;
    fftw_complex* out_;
    fftw_plan plan_fwd_;
    fftw_plan plan_inv_;
    int planned_threads_;
    bool threads_initialized_;
};

FFTWPlanCache& fftw_cache() {
    static FFTWPlanCache cache;
    return cache;
}

int fftw_thread_count() {
#if defined(PM_SERIAL)
    return 1;
#elif defined(_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace

void solve_poisson(SimState& state) {

    const int N = static_cast<int>(Ng);
    auto& cache = fftw_cache();
    cache.ensure_plans_for_threads(fftw_thread_count());
    double* in = cache.in();
    fftw_complex* out = cache.out();
    const std::size_t real_size = cache.real_size();

    // ──────────────────────────────────────────────────────────────────
    // Copiar densidad al buffer de entrada
    // ──────────────────────────────────────────────────────────────────
    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < real_size; ++i)
        in[i] = state.density[i];

    // ── FFT directa: ρ(x) → ρ̂(k) ──────────────────────────────────
    fftw_execute(cache.forward_plan());

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
    fftw_execute(cache.inverse_plan());

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
}

} // namespace pm
