import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import os

# --- CONFIGURACIÓN DE LA ANIMACIÓN ---
DATA_DIR = "data"              # Directorio donde están tus snapshots .bin
BASE_FILENAME = "snap_"       # Prefijo de los archivos
START_SNAP = 0                # Primer snapshot a animar
END_SNAP = 300                 # Último snapshot
BINS = 256                    # Resolución del histograma (ajustar a Ng)
FPS = 30                      # Cuadros por segundo del video final
OUTPUT_VIDEO = "evolucion_cosmica.mp4" # Nombre del archivo de salida

# Configuración visual global
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10, 10), facecolor='black')

# Creamos marcadores vacíos para la barra de color y el histograma
quadmesh = None
cbar = None

# --- FUNCIÓN PARA LEER UN SNAPSHOT INDIVIDUAL (Tu código adaptado) ---
def read_snapshot(snap_id):
    filename = os.path.join(DATA_DIR, f"{BASE_FILENAME}{snap_id:04d}.bin")
    
    if not os.path.exists(filename):
        print(f"Advertencia: Archivo {filename} no encontrado. Saltando.")
        return None

    # 1. Leer el archivo ignorando el primer valor problemático por ahora
    raw_data = np.fromfile(filename, dtype=np.float64)

    # 2. Si el tamaño total es divisible por 7, son puras partículas.
    # Si no, significa que SÍ hay header y son puras partículas.
    if len(raw_data) % 7 == 0:
        particles_data = raw_data.reshape(-1, 7)
    else:
        # Forzamos el reshape con los datos restantes
        particles_data = raw_data[1:].reshape(-1, 7)
    
    # 3. Extraer Posiciones (x, y, z son las primeras 3 columnas)
    return particles_data[:, 0:2] # Solo necesitamos X e Y para el hist2d

# --- FUNCIÓN DE INICIALIZACIÓN DE LA ANIMACIÓN ---
def init():
    ax.set_title("Evolución de Materia - pm_cosmo (Paso 0)", color='white', pad=20)
    ax.axis('off')
    # Definimos los límites basados en tu Ng (ej. 128 o 256)
    # Si usas la convención [0, Ng), ajusta aquí:
    ax.set_xlim(0, BINS)
    ax.set_ylim(0, BINS)
    return []

# --- FUNCIÓN DE ACTUALIZACIÓN DE CUADRO (Frame update) ---
def update(snap_id):
    global quadmesh, cbar
    
    pos_2d = read_snapshot(snap_id)
    if pos_2d is None:
        return []

    # Limpiamos el histograma anterior para no superponer imágenes
    ax.clear()
    ax.axis('off')
    
    # Actualizamos el título con el número de snapshot
    ax.set_title(f"Evolución de Materia - pm_cosmo (Paso {snap_id})", color='white', pad=20)

    # 4. Histograma de Densidad 2D con escala logarítmica
    # Usamos LogNorm para que el colormap 'magma' represente el log de la densidad.
    counts, xedges, yedges, im = ax.hist2d(pos_2d[:, 0], pos_2d[:, 1], 
                                           bins=BINS, cmap='magma')
    
    # Manejo de la barra de color (solo la creamos la primera vez para no duplicarla)
    if cbar is None:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(colors='white')
        cbar.set_label('Densidad Logarítmica', color='white')

    return [im]

# --- GENERAR LA ANIMACIÓN ---
print(f"Iniciando animación de {END_SNAP - START_SNAP + 1} snapshots...")

# Creamos el objeto Animation
ani = animation.FuncAnimation(fig, update, frames=range(START_SNAP, END_SNAP + 1),
                              init_func=init, blit=False)

# Guardar la animación como video .mp4 (requiere ffmpeg instalado)
# Si no tienes ffmpeg, cambia writer='ffmpeg' por writer='pillow' y OUTPUT_VIDEO por .gif
try:
    print(f"Guardando video como {OUTPUT_VIDEO} (esto puede tardar unos minutos)...")
    ani.save(OUTPUT_VIDEO, writer='ffmpeg', fps=FPS, dpi=100, savefig_kwargs={'facecolor':'black'})
    print("¡Animación guardada con éxito!")
except Exception as e:
    print(f"Error al guardar el video: {e}")
    print("Asegúrate de tener 'ffmpeg' instalado o cambia el formato a .gif")

plt.close(fig) # Cierra la figura para liberar memoria