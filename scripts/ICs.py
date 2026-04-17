import numpy as np

def power_spectrum_eh(k, h=0.674, Om=0.315, Ob=0.049, ns=0.965):
    """Espectro de Eisenstein & Hu para cosmología de referencia Planck 2018."""
    theta_2p7 = 2.725 / 2.7 
    keq = 7.46e-2 * Om * h**2 * theta_2p7**-2
    q = k / (keq * h)
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731 / (1 + 62.5 * q)
    T = L0 / (L0 + C0 * q**2)
    return k**ns * T**2

def generate_ics_for_cpp(Ng, Lbox, seed=42):
    np.random.seed(seed)
    
    # 1. Grilla en unidades de celda (0, 1, ..., Ng-1)
    k_lin = np.fft.fftfreq(Ng) * (2 * np.pi)
    kx, ky, kz = np.meshgrid(k_lin, k_lin, k_lin, indexing='ij')
    k_sq = kx**2 + ky**2 + kz**2
    k_sq[0, 0, 0] = 1.0 # Evitar división por cero

    # 2. Espectro de Potencia escalado
    k_phys = np.sqrt(k_sq) * (Ng / Lbox)
    Pk = power_spectrum_eh(k_phys)
    Pk[0,0,0] = 0

    # 3. Campo de desplazamiento de Zel'dovich (PSI)
    white_noise = np.random.normal(size=(Ng, Ng, Ng)) + 1j * np.random.normal(size=(Ng, Ng, Ng))
    delta_k = white_noise * np.sqrt(Pk)

    # El factor -1j * k / k^2 viene de resolver la ec. de Poisson en Fourier
    psi_x = np.fft.ifftn(-1j * kx * delta_k / k_sq).real
    psi_y = np.fft.ifftn(-1j * ky * delta_k / k_sq).real
    psi_z = np.fft.ifftn(-1j * kz * delta_k / k_sq).real

    # 4. Posiciones en unidades de celda [0, Ng)
    # Tu C++ usa std::fmod(pos, ng), así que aquí nos mantenemos en ese rango
    amplitude = 0.1  # Ajuste de escala de la perturbación inicial
    q_lin = np.arange(Ng)
    QX, QY, QZ = np.meshgrid(q_lin, q_lin, q_lin, indexing='ij')
    
    pos = np.stack([QX + amplitude * psi_x, 
                    QY + amplitude * psi_y, 
                    QZ + amplitude * psi_z], axis=-1).reshape(-1, 3)
    
    # 5. Velocidades (Zel'dovich Kick inicial)
    # Estas v0 serán retrasadas medio paso por el leapfrog_half_kick en C++
    v_factor = 0.02 
    vel = v_factor * np.stack([psi_x, psi_y, psi_z], axis=-1).reshape(-1, 3)

    return np.mod(pos, Ng), vel

# --- PARÁMETROS IGUALES A TU CONFIG.HPP ---
NG = 128
L_BOX = 100.0

pos, vel = generate_ics_for_cpp(NG, L_BOX)
data = np.hstack([pos, vel])

# --- GUARDADO EN FORMATO ASCII PARA IC_READER.CPP ---
filename = "ic/ic_128.dat"
with open(filename, "w") as f:
    # Línea 1: N_PART (Tu C++ hace f >> n_part)
    f.write(f"{len(pos)}\n")
    # Líneas 2 en adelante: x y z vx vy vz
    np.savetxt(f, data, fmt='%.8f')

print(f"¡Éxito! Archivo {filename} generado para tu simulador.")
print(f"Recuerda correrlo como: ./pm_parallel -i {filename}")