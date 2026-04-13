#pragma once

// =============================================================================
// types.hpp
// Estructuras de datos centrales del simulador PM.
//
// Se definen aquí (y no en config.hpp) para separar datos de parámetros,
// siguiendo el principio de responsabilidad única.
// =============================================================================

#include <vector>
#include <array>
#include <cassert>
#include <stdexcept>
#include "config.hpp"

namespace pm {

// =============================================================================
// Particle
// Representa una partícula de materia oscura.
//
// pos[3]: posición en coordenadas de la caja, rango [0, Ng)
//         (es decir, en unidades de celda, no en Mpc/h)
// vel[3]: velocidad en unidades internas
// mass  : masa de la partícula (todas iguales en nuestro caso)
// =============================================================================
struct Particle {
    double pos[3];   // x, y, z en unidades de celda [0, Ng)
    double vel[3];   // vx, vy, vz en unidades internas
    double mass;     // masa (normalizada: suma total = 1)

    Particle() : pos{0,0,0}, vel{0,0,0}, mass(1.0 / N_PARTICLES) {}
};

// =============================================================================
// Grid3D<T>
// Malla tridimensional periódica almacenada en memoria contigua (row-major).
//
// El índice plano para (ix, iy, iz) es:
//   idx = ix * Ng*Ng + iy * Ng + iz
//
// La periodicidad se maneja con wrap(): si un índice sale de [0, Ng),
// se aplica módulo. Esto es correcto para condiciones de frontera periódicas
// estándar en cosmología.
//
// T puede ser double (densidad, potencial, fuerza) o
//   std::complex<double> (transformada de Fourier).
// =============================================================================
template<typename T>
class Grid3D {
public:
    std::vector<T> data;
    const std::size_t n;   // puntos por dimensión (= Ng)

    // Constructor: crea malla cúbica n×n×n inicializada a cero
    explicit Grid3D(std::size_t size = Ng)
        : data(size*size*size, T(0)), n(size) {}

    // Acceso por índice plano (más rápido en loops)
    inline T& operator[](std::size_t i)       { return data[i]; }
    inline const T& operator[](std::size_t i) const { return data[i]; }

    // Acceso 3D con verificación de rango (usa en debug, no en hot loops)
    inline T& at(int ix, int iy, int iz) {
        return data[idx(wrap(ix), wrap(iy), wrap(iz))];
    }
    inline const T& at(int ix, int iy, int iz) const {
        return data[idx(wrap(ix), wrap(iy), wrap(iz))];
    }

    // Convierte (ix,iy,iz) a índice plano SIN verificar rango
    // (ix,iy,iz) deben estar en [0, n)
    inline std::size_t idx(std::size_t ix, std::size_t iy, std::size_t iz) const {
        return ix * n * n + iy * n + iz;
    }

    // Devuelve número total de celdas
    inline std::size_t total() const { return n * n * n; }

    // Pone todos los valores a cero
    void zero() { std::fill(data.begin(), data.end(), T(0)); }

    // Aplica condiciones de frontera periódicas a un índice entero
    inline int wrap(int i) const {
        int ni = static_cast<int>(n);
        return ((i % ni) + ni) % ni;
    }
};

// Alias convenientes para los tipos de malla más usados
using DensityGrid   = Grid3D<double>;   // ρ(x): densidad en espacio real
using PotentialGrid = Grid3D<double>;   // φ(x): potencial gravitacional
using ForceGrid     = Grid3D<double>;   // componente de fuerza (gx, gy, gz)

// =============================================================================
// SimState
// Agrupa todo el estado mutable de la simulación en un solo struct.
// Esto facilita pasar el estado completo entre funciones sin acumulación
// de parámetros.
// =============================================================================
struct SimState {
    std::vector<Particle> particles;   // arreglo de todas las partículas

    DensityGrid   density;             // campo de densidad ρ en la malla
    PotentialGrid potential;           // potencial gravitacional φ
    ForceGrid     force_x;             // fuerza en x: -∂φ/∂x
    ForceGrid     force_y;             // fuerza en y: -∂φ/∂y
    ForceGrid     force_z;             // fuerza en z: -∂φ/∂z

    double time;                       // tiempo actual de la simulación
    int    step;                       // número de paso actual

    SimState()
        : particles(N_PARTICLES),
          density(Ng), potential(Ng),
          force_x(Ng), force_y(Ng), force_z(Ng),
          time(0.0), step(0)
    {}

    // Limpia los campos de malla (se llama al inicio de cada paso)
    void clear_grids() {
        density.zero();
        potential.zero();
        force_x.zero();
        force_y.zero();
        force_z.zero();
    }
};

} // namespace pm
