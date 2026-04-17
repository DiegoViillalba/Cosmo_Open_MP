# Revision tecnica completa de Cosmo_Open_MP

## 1. Resumen ejecutivo

Este proyecto implementa un simulador cosmologico tipo Particle-Mesh (PM) con pipeline clasico:

1. Deposito de masa en malla (CIC)
2. Resolucion de Poisson por FFT
3. Gradiente del potencial
4. Interpolacion de fuerza a particulas
5. Integracion temporal Leap-Frog

El codigo esta bien estructurado por modulos y la suite de pruebas pasa completa (17/17).

Sobre el cuello de botella:

- En costo numerico por paso, si: la etapa dominante es Poisson FFT.
- En costo total extremo a extremo, el principal bloqueador de speedup no es solo FFT: tambien pesa mucho I/O (snapshot inicial grande) y trabajo fuera del timer por etapa.

## 2. Arquitectura por modulo

### 2.1 Configuracion y tipos base

#### include/config.hpp

- Centraliza parametros globales: tamano de malla, numero de particulas, constantes fisicas y parametros temporales.
- Implementa macros condicionales para activar/desactivar OpenMP segun `PM_SERIAL`.
- Permite sobreescribir `PM_GRID_SIZE` por CMake.

#### include/types.hpp

- Define `Particle`, `Grid3D<T>` y `SimState`.
- `Grid3D<T>` encapsula indice plano, acceso periodico y almacenamiento contiguo.
- `SimState` agrupa particulas + campos de malla + estado temporal de la simulacion.

### 2.2 Orquestacion de simulacion

#### include/simulation.hpp y src/simulation.cpp

- Clase `Simulation` coordina inicializacion, loop principal y snapshots.
- `init()`:
  - Lee IC de archivo o genera IC sinteticas.
  - Ejecuta calculo inicial de fuerzas para `half-kick`.
  - Escribe diagnostico inicial y snapshot de paso 0.
- `run()`:
  - Ejecuta `single_step()` por cada iteracion.
  - Calcula diagnosticos periodicos.
  - Muestra barra de progreso.
  - Reporta desglose de tiempos por etapa.
- `single_step()` implementa el pipeline PM completo.

#### src/main.cpp

- Parsea CLI (`-t`, `-s`, `-o`, `-b`, `-v`, `-g`).
- Construye y ejecuta `Simulation`.
- Muestra banner y modo serial/paralelo.

### 2.3 Modulos numericos del pipeline PM

#### include/cic.hpp y src/cic.cpp

- Deposita masa de particulas en 8 celdas vecinas (trilineal CIC).
- Maneja condiciones periodicas en indices.
- En paralelo usa `omp atomic` en 8 escrituras por particula para evitar race conditions.
- Normaliza densidad al final con factor `N_CELLS`.

#### include/poisson_fft.hpp y src/poisson_fft.cpp

- Resuelve `nabla^2 phi = 4*pi*G*rho` via FFTW3:
  - `r2c` directa
  - multiplicacion por kernel de Green en k-espacio
  - `c2r` inversa
  - normalizacion
- El modo DC (`k=0`) se fuerza a cero.
- Solo el bucle del kernel esta paralelizado con OpenMP.
- Las llamadas `fftw_execute` se ejecutan en modo secuencial en esta implementacion.

#### include/gradient.hpp y src/gradient.cpp

- Calcula `g = -grad(phi)` con diferencias centradas de segundo orden.
- Usa `collapse(3)` para paralelizar barrido 3D de la malla.

#### include/force_interp.hpp y src/force_interp.cpp

- Interpola fuerza de malla a particulas usando CIC inverso (mismos pesos trilineales).
- Solo lectura de malla y escritura por particula: paralelizacion limpia.

#### include/integrator.hpp y src/integrator.cpp

- Integrador Leap-Frog (`kick + drift`) para evolucion temporal.
- Incluye `leapfrog_half_kick()` para arranque correcto en medio paso.
- Aplica periodicidad con `fmod`.

### 2.4 Diagnosticos y salida

#### include/diagnostics.hpp y src/diagnostics.cpp

- Calcula energia cinetica, potencial, total, momento y cociente virial.
- Usa reducciones OpenMP para sumas globales.
- Interpola potencial de forma NGP para estimar `Ep`.

#### include/output.hpp y src/output.cpp

- Escribe snapshots en ASCII o binario.
- ASCII: legible, muy pesado para 256^3.
- Binario: mas rapido, aun asi muy grande en paso 0 (cabecera + vector completo de particulas).
- Agrega fila en CSV de diagnosticos.

### 2.5 Medicion de tiempos

#### include/timer.hpp

- `ScopedTimer` RAII y `StageTimer` acumulado por nombre de etapa.
- `wall_time()` usa `omp_get_wtime()` cuando hay OpenMP.

### 2.6 Pruebas y benchmark

#### tests/*.cpp y tests/test_framework.hpp

- Framework minimalista propio con registro automatico de pruebas.
- Cobertura funcional de CIC, Poisson, integrador y diagnosticos.
- Resultado observado: 17 pruebas, 17 exitos.

#### bench/run_bench.sh

- Ejecuta serial y paralelo para varios hilos.
- Intenta extraer tiempos por etapa de logs y escribir CSV.
- Tiene un fallo de parseo: la expresion regular busca patron con `" s"` que no coincide con el formato actual de `StageTimer`, por eso quedan ceros en columnas por etapa.

## 3. Analisis de cuello de botella

## 3.1 Evidencia de tiempos por etapa (run corto, Ng=256, 2 pasos)

Medicion serial (`pm_serial -g -b -s 2`):

- CIC: 0.5210 s (13.8%)
- Poisson FFT: 1.3938 s (36.9%)
- Gradiente: 0.1190 s (3.1%)
- Force interp: 1.1445 s (30.3%)
- Integrador: 0.6009 s (15.9%)
- Total pipeline: 3.7791 s

Medicion paralela 8 hilos (`pm_parallel -g -b -t 8 -s 2`):

- CIC: 0.5193 s (18.5%)
- Poisson FFT: 1.2839 s (45.8%)
- Gradiente: 0.1255 s (4.5%)
- Force interp: 0.6546 s (23.4%)
- Integrador: 0.2183 s (7.8%)
- Total pipeline: 2.8017 s

## 3.2 Conclusiones del cuello de botella

1. Si, FFT es cuello de botella principal dentro del pipeline numerico.
2. No toda FFT esta aprovechando todos los hilos:
   - El codigo paraleliza el kernel en k-espacio.
   - Pero `fftw_execute()` no esta configurado en modo multihilo (`fftw_init_threads`, `fftw_plan_with_nthreads`, y enlace a `fftw3_omp`).
3. CIC casi no escala por uso intensivo de `omp atomic` en 8 escrituras por particula.
4. `force_interp` e `integrator` si escalan claramente.
5. El rendimiento total extremo a extremo puede quedar plano por costos fuera del timer por etapa:
   - Generacion de IC en `init()` (costosa y no paralelizada).
   - Snapshot inicial muy grande en paso 0 (ASCII peor; binario tambien pesado).

## 4. Por que parece que OpenMP no reparte todo el trabajo

La percepcion es correcta por tres razones:

1. Parte de la FFT permanece secuencial (llamadas FFTW).
2. La etapa CIC usa atomics, que pueden serializar memoria y limitar speedup.
3. El benchmark actual mezcla costos de I/O y setup con el calculo, ocultando mejoras reales de partes paralelizables.

## 5. Hallazgos tecnicos adicionales

1. `bench/run_bench.sh` no captura tiempos por etapa correctamente (columnas en cero).
2. El CSV de diagnosticos calcula `dE_frac` como `(Etot - Ek - Ep)/|Etot|`, que siempre da 0 por identidad algebraica en el mismo paso.
3. `run_tests` se enlaza con OpenMP aunque el comentario dice que es para aislar logica sin OpenMP.

## 6. Recomendaciones priorizadas

1. Activar FFTW multihilo real:
   - Enlazar `fftw3_omp`.
   - Inicializar hilos FFTW una sola vez.
   - Crear y reutilizar planes con `fftw_plan_with_nthreads(nthreads)`.

2. Reducir impacto de CIC:
   - Probar estrategia de buffers privados por bloque/tile y reduccion parcial.
   - Alternativamente, ordenar particulas por celda para bajar colisiones atomicas.

3. Separar benchmark de computo vs I/O:
   - Opcion para desactivar snapshot inicial y CSV durante benchmark.
   - Medir `init`, `compute`, `output` por separado.

4. Corregir parser de `bench/run_bench.sh` para extraer tiempos reales.

5. Revisar medicion de energia relativa:
   - Guardar `E0` del paso inicial y calcular `dE_frac = (Etot - E0)/|E0|`.

## 7. Estado actual de calidad

- Diseno modular: bueno.
- Legibilidad y comentarios: alta.
- Correctitud funcional (tests): buena.
- Escalado OpenMP global: limitado por FFT parcial, atomics en CIC y sobrecosto de I/O/setup.

En resumen: tu sospecha sobre FFT es valida dentro del nucleo numerico, y tambien es cierto que el uso de hilos no esta plenamente distribuido en toda la ruta critica.