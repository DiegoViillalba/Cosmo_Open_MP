// =============================================================================
// diagnostics.cpp
// Cálculo de energías cinética, potencial y momento total.
//
// DEMOSTRACIÓN CENTRAL DE LA CLÁUSULA REDUCTION:
//
// Sin reduction, este código sería incorrecto en paralelo:
//
//   double Ek = 0;
//   #pragma omp parallel for          // INCORRECTO: race condition en Ek
//   for (int i = 0; i < N; i++)
//       Ek += 0.5 * m * v2;           // múltiples hilos escriben en Ek
//
// Con reduction, OpenMP crea una copia privada de Ek por hilo, cada hilo
// acumula en su copia local sin conflicto, y al final las suma:
//
//   #pragma omp parallel for reduction(+:Ek)   // CORRECTO
//   for (int i = 0; i < N; i++)
//       Ek += ...;
//
// El resultado es idéntico al secuencial (salvo orden de suma de punto
// flotante, que puede dar diferencias del orden de epsilon de máquina).
// =============================================================================

#include "diagnostics.hpp"
#include "config.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

namespace pm {

Diagnostics compute_diagnostics(const SimState& state) {

    Diagnostics d{};
    const std::size_t N_PART = state.particles.size();
    const int N = static_cast<int>(Ng);

    // ──────────────────────────────────────────────────────────────────
    // Energía cinética: Ek = Σ_i ½ m_i (vx² + vy² + vz²)
    //
    // reduction(+:Ek): cada hilo acumula su porción de la suma en una
    // variable local; OpenMP las reduce (suma) al salir de la región.
    // ──────────────────────────────────────────────────────────────────
    double Ek = 0.0;
    double Px = 0.0, Py = 0.0, Pz = 0.0;

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static) \
        reduction(+:Ek) reduction(+:Px) reduction(+:Py) reduction(+:Pz)
    #endif
    for (std::size_t i = 0; i < N_PART; ++i) {
        const auto& p  = state.particles[i];
        const double m = p.mass;
        const double v2 = p.vel[0]*p.vel[0]
                        + p.vel[1]*p.vel[1]
                        + p.vel[2]*p.vel[2];
        Ek += 0.5 * m * v2;
        Px += m * p.vel[0];
        Py += m * p.vel[1];
        Pz += m * p.vel[2];
    }

    // ──────────────────────────────────────────────────────────────────
    // Energía potencial: Ep = Σ_i m_i · φ(x_i)
    //
    // Interpolamos el potencial en la posición de cada partícula usando
    // el valor de la celda más cercana (NGP: Nearest Grid Point).
    // Para mayor precisión se podría usar interpolación CIC, pero NGP
    // es suficiente para diagnósticos.
    // ──────────────────────────────────────────────────────────────────
    double Ep = 0.0;

    #ifndef PM_SERIAL
    #pragma omp parallel for schedule(static) reduction(+:Ep)
    #endif
    for (std::size_t i = 0; i < N_PART; ++i) {
        const auto& p = state.particles[i];

        // Celda más cercana (NGP)
        const int ix = static_cast<int>(std::round(p.pos[0])) % N;
        const int iy = static_cast<int>(std::round(p.pos[1])) % N;
        const int iz = static_cast<int>(std::round(p.pos[2])) % N;

        const double phi = state.potential[
            static_cast<std::size_t>(ix * N * N + iy * N + iz)];

        Ep += p.mass * phi;
    }
    // El factor ½ evita doble conteo (cada par se cuenta dos veces en
    // la suma sobre partículas individuales)
    Ep *= 0.5;

    d.kinetic_energy   = Ek;
    d.potential_energy = Ep;
    d.total_energy     = Ek + Ep;
    d.momentum[0]      = Px;
    d.momentum[1]      = Py;
    d.momentum[2]      = Pz;
    d.virial_ratio     = (Ek > 1e-30) ? std::abs(Ep) / (2.0 * Ek) : 0.0;

    return d;
}

// -----------------------------------------------------------------------------
// print_diagnostics
// Imprime una línea de diagnóstico legible en stdout.
// El formato CSV en el header facilita redirigir la salida a un archivo.
// -----------------------------------------------------------------------------
void print_diagnostics(int step, double time, const Diagnostics& d) {
    // Header solo en el primer llamado
    static bool header_printed = false;
    if (!header_printed) {
        std::cout << "\n"
                  << std::setw(6)  << "step"  << "  "
                  << std::setw(8)  << "time"   << "  "
                  << std::setw(14) << "Ek"     << "  "
                  << std::setw(14) << "Ep"     << "  "
                  << std::setw(14) << "Etot"   << "  "
                  << std::setw(10) << "|P|"    << "  "
                  << std::setw(8)  << "virial" << "\n"
                  << std::string(80, '-') << "\n";
        header_printed = true;
    }

    const double P_mag = std::sqrt(
        d.momentum[0]*d.momentum[0] +
        d.momentum[1]*d.momentum[1] +
        d.momentum[2]*d.momentum[2]);

    std::cout << std::setw(6)  << step
              << "  " << std::fixed << std::setprecision(4)
              << std::setw(8)  << time
              << "  " << std::scientific << std::setprecision(6)
              << std::setw(14) << d.kinetic_energy
              << "  " << std::setw(14) << d.potential_energy
              << "  " << std::setw(14) << d.total_energy
              << "  " << std::setw(10) << P_mag
              << "  " << std::fixed << std::setprecision(3)
              << std::setw(8)  << d.virial_ratio
              << "\n";
}

} // namespace pm
