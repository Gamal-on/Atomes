# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:54:20 2026

@author: damie
"""

# -*- coding: utf-8 -*-
"""
Calcul et affichage 3D de la densité résiduelle Delta_rho
"""
import numpy as np
import matplotlib.pyplot as plt

# === IMPORTS ===
from FS_expérimentaux_dict import *
from FS_théorique_dynamique_dict import *
from FS_théorique_statique_dict import *

# ==============================================================
# PARAMÈTRES DE MAILLE
# ==============================================================
a = 5.578   # Å
b = 5.578
c = 4.686
alpha = 90 * np.pi / 180
beta  = 90 * np.pi / 180
gamma = 90 * np.pi / 180

# Volume réduit
v = np.sqrt(
    1
    - np.cos(alpha)**2
    - np.cos(beta)**2
    - np.cos(gamma)**2
    + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
)

# Volume réduit
v = np.sqrt(
    1
    - np.cos(alpha)**2
    - np.cos(beta)**2
    - np.cos(gamma)**2
    + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
)

# Vecteurs directs (base cartésienne)
a1 = np.array([a, 0, 0])
a2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
a3 = np.array([
    c * np.cos(beta),
    c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
    c * v / np.sin(gamma)
])

# Vecteurs réciproques (base cartésienne)
b1 = np.array([
    1 / a,
    -np.cos(gamma) / (a * np.sin(gamma)),
    (np.cos(gamma) * np.cos(alpha) - np.cos(beta)) / (a * v * np.sin(gamma))
])
b2 = np.array([
    0,
    1 / (b * np.sin(gamma)),
    (np.cos(gamma) * np.cos(beta) - np.cos(alpha)) / (b * v * np.sin(gamma))
])
b3 = np.array([
    0,
    0,
    np.sin(gamma) / (c * v)
])


# ==============================================================
# INDICES ET PAS
# ==============================================================
hmin, hmax = -14, 14
kmin, kmax = -14, 14
lmin, lmax = -14, 14

pas_x = 1.0   # Å
pas_y = 1.0
pas_z = 1.0

dmin1, dmax1 = hmin * a, hmax * a
dmin2, dmax2 = kmin * b, kmax * b
dmin3, dmax3 = lmin * c, lmax * c

# ==============================================================
# PRÉ-CALCUL DES AMPLITUDES ET PHASES (hors grille)
# CORRECTION 3 : dans le code original, amp_th_stat et amp_th étaient
# recalculées à chaque point de la grille (N_points fois), alors qu'elles
# ne dépendent que de (h,k,l). On les sort ici une bonne fois pour toutes.
# ==============================================================

# CORRECTION 4 : vérification que les clés sont communes aux trois dicts
common_keys = (set(FS_théorique_dynamique)
               & set(FS_théorique_statique)
               #& set(FS_expérimentaux)"    
               )

missing = set(FS_théorique_dynamique) - common_keys
if missing:
    print(f"⚠ {len(missing)} réflexions absentes d'un dictionnaire, ignorées.")

reflections = []
for e in common_keys:
    h, k, l = e

    amp_stat = np.sqrt(
        FS_théorique_statique[e][0]**2 +
        FS_théorique_statique[e][1]**2
    )
    amp_dyn = np.sqrt(
        FS_théorique_dynamique[e][0]**2 +
        FS_théorique_dynamique[e][1]**2
    )
    # CORRECTION 5 : dans la version précédente corrigée amp_th_stat venait
    # bien de FS_théorique_statique, mais amp_th et amp_th_stat étaient
    # recalculées à chaque itération de la boucle sur les points.
    # Ici delta_amp est calculé une seule fois par réflexion.
    delta_amp = amp_stat - amp_dyn

    phi = FS_théorique_dynamique[e][2]

    # Coefficients de G·r = k1*(h*b1x) + k2*(h*b1y + k*b2y) + k3*(h*b1z + k*b2z + l*b3z)
    # CORRECTION 6 : les coefficients Gx, Gy, Gz sont précalculés ici
    # plutôt que reconstruits à chaque point de la grille.
    Gx = h * b1[0]
    Gy = h * b1[1] + k * b2[1]
    Gz = h * b1[2] + k * b2[2] + l * b3[2]

    reflections.append((delta_amp, phi, Gx, Gy, Gz))

# Conversion en tableaux NumPy pour vectorisation
delta_amps = np.array([r[0] for r in reflections])   # (M,)
phis       = np.array([r[1] for r in reflections])   # (M,)
Gx_arr     = np.array([r[2] for r in reflections])   # (M,)
Gy_arr     = np.array([r[3] for r in reflections])   # (M,)
Gz_arr     = np.array([r[4] for r in reflections])   # (M,)

# Facteurs d'amplitude complexe : delta_amp * exp(i*phi), shape (M,)
F_hkl = delta_amps * np.exp(1j * phis)

# ==============================================================
# CALCUL DE LA DENSITÉ RÉSIDUELLE — VECTORISÉ
# CORRECTION 7 : remplacement des 3 while imbriqués par meshgrid +
# broadcasting NumPy. La boucle Python sur les points disparaît.
# On conserve une petite boucle sur les réflexions (M itérations)
# plutôt qu'une matrice (N_pts × M) qui pourrait saturer la RAM.
# ==============================================================
x_vals = np.arange(dmin1, dmax1, pas_x)
y_vals = np.arange(dmin2, dmax2, pas_y)
z_vals = np.arange(dmin3, dmax3, pas_z)

Nx, Ny, Nz = len(x_vals), len(y_vals), len(z_vals)
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')  # (Nx, Ny, Nz)

Delta_rho_grid = np.zeros((Nx, Ny, Nz), dtype=complex)

total = len(reflections)
for i, (f, gx, gy, gz) in enumerate(zip(F_hkl, Gx_arr, Gy_arr, Gz_arr)):

    # Barre de progression sur les réflexions (beaucoup plus rapide à afficher)
    pct    = (i + 1) / total
    filled = int(40 * pct)
    bar    = '█' * filled + '░' * (40 - filled)
    print(f'\r[{bar}] {pct*100:5.1f}%  réflexion {i+1}/{total}',
          end='', flush=True)

    G_dot_r = X * gx + Y * gy + Z * gz          # broadcasting (Nx, Ny, Nz)
    Delta_rho_grid += f * np.exp(1j * G_dot_r)  # accumulation

print(f'\nCalcul terminé — {Delta_rho_grid.size} points.')

# CORRECTION 8 : on travaille directement sur le tableau 3D (plus de dict)
values = np.real(Delta_rho_grid)   # (Nx, Ny, Nz)



x_flat = X.ravel()
y_flat = Y.ravel()
z_flat = Z.ravel()
v_flat = values.ravel()

# ==============================================================
# SAUVEGARDE DE x_flat, y_flat, z_flat, v_flat
# ==============================================================
np.savez('delta_rho_flat.npz',
         x_flat=x_flat,
         y_flat=y_flat,
         z_flat=z_flat,
         v_flat=v_flat)
print("Données plates sauvegardées dans delta_rho_flat.npz")


