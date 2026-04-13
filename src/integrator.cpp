// =============================================================================
// integrator.cpp
// Integrador Leap-Frog (Verlet de velocidades desfasadas).
//
// El Leap-Frog es un integrador simpléctico de 2do orden:
//   v_{n+1/2} = v_{n-1/2} + a_n · dt     [kick]
//   x_{n+1}   = x_n + v_{n+1/2} · dt     [drift]
//
// "Simpléctico" significa que conserva el volumen en el espacio de fase
// (x, v), lo que en la práctica implica buena conservación de energía
// a largo plazo, sin acumulación de error sistemático.
//
// Condiciones de frontera periódicas:
//   Si x sale de [0, Ng) → se reinserta modularmente.
//   Esto es correcto en cosmología (caja periódica).
//
// Paralelismo: loop over partículas completamente independiente.
//   - Sin lecturas compartidas: cada hilo lee accel[i] y escribe particles[i].
//   - schedule(static) óptimo: costo uniforme por partícula.
// =============================================================================

#include "integrator.hpp"
#include "config.hpp"

#include <cmath>

namespace pm {

// -----------------------------------------------------------------------------
// leapfrog_step
// Un paso completo: kick de velocidad + drift de posición.
// Llamado en cada iteración del loop principal.
// -----------------------------------------------------------------------------
void leapfrog_step(SimState& state,
                   const std::vector<std::array<double,3>>& accel) {

    const std::size_t N_PART = state.particles.size();
    const double dt = DT;
    const double ng = static_cast<double>(Ng);

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < N_PART; ++i) {

        auto& p = state.particles[i];

        // ── Kick: actualizar velocidad con aceleración actual ──────────
        // v_{n+1/2} = v_{n-1/2} + a_n · dt
        p.vel[0] += accel[i][0] * dt;
        p.vel[1] += accel[i][1] * dt;
        p.vel[2] += accel[i][2] * dt;

        // ── Drift: actualizar posición con velocidad nueva ─────────────
        // x_{n+1} = x_n + v_{n+1/2} · dt
        p.pos[0] += p.vel[0] * dt;
        p.pos[1] += p.vel[1] * dt;
        p.pos[2] += p.vel[2] * dt;

        // ── Condiciones de frontera periódicas ─────────────────────────
        // fmod con ajuste para valores negativos:
        //   si pos = -0.1 → fmod(-0.1 + Ng*10, Ng) ≈ Ng - 0.1
        // Multiplicar por 10 asegura que pos + ng*10 > 0 siempre.
        p.pos[0] = std::fmod(p.pos[0] + ng * 10.0, ng);
        p.pos[1] = std::fmod(p.pos[1] + ng * 10.0, ng);
        p.pos[2] = std::fmod(p.pos[2] + ng * 10.0, ng);
    }

    // Avanzar el tiempo global
    state.time += dt;
    state.step  += 1;
}

// -----------------------------------------------------------------------------
// leapfrog_half_kick
// Medio kick inicial: solo se llama UNA VEZ al inicio de la simulación.
//
// El Leap-Frog requiere que las velocidades estén a t = -dt/2 antes del
// primer paso completo. Si las condiciones iniciales dan v en t=0, hay
// que retrasarlas medio paso:
//   v_{-1/2} = v_0 - a_0 · (dt/2)
//
// Esto solo es necesario cuando se parte de IC con v en t=0.
// Si las IC ya dan v_{-1/2} (como las generadas por MUSIC o N-GenIC),
// esta función no debe llamarse.
// -----------------------------------------------------------------------------
void leapfrog_half_kick(SimState& state,
                        const std::vector<std::array<double,3>>& accel) {

    const std::size_t N_PART = state.particles.size();
    const double half_dt = DT * 0.5;

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t i = 0; i < N_PART; ++i) {
        auto& p = state.particles[i];
        // v_{-1/2} = v_0 - a_0 · (dt/2)   (retroceso de medio paso)
        p.vel[0] -= accel[i][0] * half_dt;
        p.vel[1] -= accel[i][1] * half_dt;
        p.vel[2] -= accel[i][2] * half_dt;
    }
}

} // namespace pm
