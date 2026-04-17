import numpy as np
import matplotlib.pyplot as plt

# FILENAME = "build/data/snap_0000.bin"   
FILENAME = "data/snap_0050.bin"

# 1. Leer el archivo ignorando el primer valor problemático por ahora
# Sabemos que el total de datos es 14,680,064
raw_data = np.fromfile(FILENAME, dtype=np.float64)

# 2. Si el tamaño total es 14,680,065, el primero es el header.
# Si es 14,680,064, significa que NO hay header y son puras partículas.
if len(raw_data) % 7 == 0:
    print("Detectado: Archivo sin header (solo partículas).")
    particles_data = raw_data.reshape(-1, 7)
else:
    print(f"Detectado: Archivo con header. Valor del header: {raw_data[0]}")
    # Forzamos el reshape con los datos restantes
    particles_data = raw_data[1:].reshape(-1, 7)

# 3. Extraer Posiciones (x, y, z son las primeras 3 columnas)
pos = particles_data[:, 0:3]

print(f"Total de partículas procesadas: {len(particles_data)}")

# 4. Histograma de Densidad 2D
plt.figure(figsize=(10, 10), facecolor='black')
# Ajustamos bins a 256 para ver la estructura clara
h = plt.hist2d(pos[:, 0], pos[:, 1], bins=256, cmap='magma')
plt.title("Distribución de Materia - pm_cosmo", color='white')
plt.axis('off')
plt.savefig('resultado.png', dpi=300)
plt.show()