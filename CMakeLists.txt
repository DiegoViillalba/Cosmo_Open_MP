cmake_minimum_required(VERSION 3.16)
project(pm_cosmo VERSION 1.0 LANGUAGES CXX)

# ─────────────────────────────────────────────
# Estándar C++17 obligatorio en todo el proyecto
# ─────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ─────────────────────────────────────────────
# Localizar OpenMP (requerido para paralelismo)
# ─────────────────────────────────────────────
find_package(OpenMP REQUIRED)

# ─────────────────────────────────────────────
# Localizar FFTW3 (transformada de Fourier)
# En Ubuntu/Debian: sudo apt install libfftw3-dev

# ─────────────────────────────────────────────
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)
pkg_check_modules(FFTW3_OMP QUIET fftw3_omp)

if(FFTW3_OMP_FOUND)
    set(FFTW3_OMP_INCLUDE_PATHS ${FFTW3_OMP_INCLUDE_DIRS})
    set(FFTW3_OMP_LINK_LIBS ${FFTW3_OMP_LIBRARIES})
else()
    find_library(FFTW3_OMP_LIBRARY NAMES fftw3_omp REQUIRED)
    set(FFTW3_OMP_INCLUDE_PATHS "")
    set(FFTW3_OMP_LINK_LIBS ${FFTW3_OMP_LIBRARY})
endif()

# ─────────────────────────────────────────────
# Flags de compilación comunes
#   -O3           optimización máxima
#   -march=native instrucciones del CPU local (AVX, etc.)
#   -Wall -Wextra advertencias útiles de compilador
# ─────────────────────────────────────────────
set(COMMON_FLAGS "-O3 -march=native -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Permitir pasar PM_GRID_SIZE desde la línea de cmake:
#   cmake .. -DPM_GRID_SIZE=128
if(DEFINED PM_GRID_SIZE)
    add_compile_definitions(PM_GRID_SIZE=${PM_GRID_SIZE})
    message(STATUS "  Grid size override: PM_GRID_SIZE=${PM_GRID_SIZE}")
endif()

# ─────────────────────────────────────────────
# Directorio de headers compartidos
# ─────────────────────────────────────────────
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# ─────────────────────────────────────────────
# Fuentes compartidas entre ambos ejecutables
# (todo excepto main.cpp)
# ─────────────────────────────────────────────
set(SHARED_SOURCES
    src/ic_reader.cpp
    src/cic.cpp
    src/poisson_fft.cpp
    src/gradient.cpp
    src/force_interp.cpp
    src/integrator.cpp
    src/diagnostics.cpp
    src/output.cpp
    src/simulation.cpp
)

# ══════════════════════════════════════════════
# TARGET 1: pm_serial
#   Compilado SIN OpenMP para tener la línea base
#   de comparación de desempeño.
# ══════════════════════════════════════════════
add_executable(pm_serial src/main.cpp ${SHARED_SOURCES})
target_include_directories(pm_serial PRIVATE ${INCLUDE_DIR} ${FFTW3_INCLUDE_DIRS} ${FFTW3_OMP_INCLUDE_PATHS})
target_link_libraries(pm_serial PRIVATE ${FFTW3_LIBRARIES} ${FFTW3_OMP_LINK_LIBS})
target_compile_definitions(pm_serial PRIVATE PM_SERIAL)
# Sin OpenMP: no se pasa OpenMP::OpenMP_CXX

# ══════════════════════════════════════════════
# TARGET 2: pm_parallel
#   Compilado CON OpenMP. Las regiones #pragma omp
#   en el código se activan solo cuando PM_SERIAL
#   NO está definido.
# ══════════════════════════════════════════════
add_executable(pm_parallel src/main.cpp ${SHARED_SOURCES})
target_include_directories(pm_parallel PRIVATE ${INCLUDE_DIR} ${FFTW3_INCLUDE_DIRS} ${FFTW3_OMP_INCLUDE_PATHS})
target_link_libraries(pm_parallel PRIVATE OpenMP::OpenMP_CXX ${FFTW3_LIBRARIES} ${FFTW3_OMP_LINK_LIBS})
# OpenMP::OpenMP_CXX ya agrega -fopenmp al compilador y al linker

# ══════════════════════════════════════════════
# TARGET 3: run_tests
#   Ejecutable de pruebas unitarias (sin OpenMP
#   para aislar lógica pura)
# ══════════════════════════════════════════════
add_executable(run_tests
    tests/test_main.cpp
    tests/test_cic.cpp
    tests/test_poisson.cpp
    tests/test_integrator.cpp
    tests/test_diagnostics.cpp
    ${SHARED_SOURCES}
)
target_include_directories(run_tests PRIVATE ${INCLUDE_DIR} ${FFTW3_INCLUDE_DIRS} ${FFTW3_OMP_INCLUDE_PATHS} ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(run_tests PRIVATE OpenMP::OpenMP_CXX ${FFTW3_LIBRARIES} ${FFTW3_OMP_LINK_LIBS})

# ─────────────────────────────────────────────
# Mensaje informativo al configurar con cmake
# ─────────────────────────────────────────────
message(STATUS "")
message(STATUS "=== pm_cosmo build configuration ===")
message(STATUS "  C++ standard  : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type    : ${CMAKE_BUILD_TYPE}")
message(STATUS "  OpenMP found  : ${OpenMP_CXX_FOUND}")
message(STATUS "  FFTW3 found   : ${FFTW3_FOUND}")
if(FFTW3_OMP_FOUND)
    message(STATUS "  FFTW3_OMP found: YES (pkg-config)")
else()
    message(STATUS "  FFTW3_OMP found: YES (find_library)")
endif()
message(STATUS "  Flags         : ${CMAKE_CXX_FLAGS}")
message(STATUS "")
message(STATUS "  Targets disponibles:")
message(STATUS "    pm_serial   -> simulación SIN OpenMP (baseline)")
message(STATUS "    pm_parallel -> simulación CON OpenMP")
message(STATUS "    run_tests   -> suite de pruebas unitarias")
message(STATUS "")
message(STATUS "  Uso rápido:")
message(STATUS "    mkdir build && cd build")
message(STATUS "    cmake ..")
message(STATUS "    make -j\$(nproc)")
message(STATUS "    ./pm_parallel ../ic/ic_128.dat 4 100")
message(STATUS "=======================================")
