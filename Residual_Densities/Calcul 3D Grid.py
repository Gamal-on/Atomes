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

# Angle de la maille
alpha = 90 * np.pi / 180
beta  = 90 * np.pi / 180
gamma = 90 * np.pi / 180

# Volume réduit
v = np.sqrt(
    1
    - np.cos(alpha)**2
    - np.cos(beta)**2
    - np.cos(gamma)**2
    + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

# Volume de la maille
V= a * b * c * v

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


# le pas est calculé avec le critère de Shannon pas< 1/(2*hmax/a) avec hmax/a = Gmax 
# Ici il n'est pas respecté car mon ordinnateur n'était pas suffisament puissant pour 
# le faire fonctionner sinon.
pas_x = 1   # Å
pas_y = 1
pas_z = 1

dmin1, dmax1 = hmin * a, hmax * a
dmin2, dmax2 = kmin * b, kmax * b
dmin3, dmax3 = lmin * c, lmax * c


# vérification que les clés sont communes aux trois dicts
common_keys = (set(FS_théorique_dynamique)
               & set(FS_théorique_statique)
               #& set(FS_expérimentaux)"    
               ) - {(0, 0, 0)} # On exclue la réflection directe

missing = set(FS_théorique_dynamique) - common_keys
if missing:
    print(f"⚠ {len(missing)} réflexions absentes d'un dictionnaire, ignorées.")

F_hkl_list = []
Gx_arr     = []
Gy_arr     = []
Gz_arr     = []

for e in common_keys:
    h_r, k_r, l_r = e
    amp_stat  = np.sqrt(FS_théorique_statique[e][0]**2  + FS_théorique_statique[e][1]**2)
    amp_dyn   = np.sqrt(FS_théorique_dynamique[e][0]**2 + FS_théorique_dynamique[e][1]**2)
    phi       = FS_théorique_dynamique[e][2]

    F_hkl_list.append((amp_stat - amp_dyn) * np.exp(1j * phi))
    Gx_arr.append(h_r * b1[0])
    Gy_arr.append(h_r * b1[1] + k_r * b2[1])
    Gz_arr.append(h_r * b1[2] + k_r * b2[2] + l_r * b3[2])

F_hkl  = np.array(F_hkl_list)
Gx_arr = np.array(Gx_arr)
Gy_arr = np.array(Gy_arr)
Gz_arr = np.array(Gz_arr)

x_vals = np.arange(dmin1, dmax1, pas_x)
y_vals = np.arange(dmin2, dmax2, pas_y)
z_vals = np.arange(dmin3, dmax3, pas_z)

Nx, Ny, Nz = len(x_vals), len(y_vals), len(z_vals)
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')  # (Nx, Ny, Nz)

Delta_rho_grid = np.zeros((Nx, Ny, Nz), dtype=complex)

total = len(F_hkl)
for i, (f, gx, gy, gz) in enumerate(zip(F_hkl, Gx_arr, Gy_arr, Gz_arr)):

    # Barre de progression sur les réflexions 
    pct    = (i + 1) / total
    filled = int(40 * pct)
    bar    = '█' * filled + '░' * (40 - filled)
    print(f'\r[{bar}] {pct*100:5.1f}%  réflexion {i+1}/{total}',
          end='', flush=True)

    G_dot_r = - 2 * np.pi *( X * gx + Y * gy + Z * gz )         # broadcasting (Nx, Ny, Nz)
    Delta_rho_grid += f * np.exp(1j * G_dot_r)  # accumulation

print(f'\nCalcul terminé — {Delta_rho_grid.size} points.')


values = np.real(Delta_rho_grid) / V   # (Nx, Ny, Nz)

x_flat = X.ravel()
y_flat = Y.ravel()
z_flat = Z.ravel()
v_flat = values.ravel()

# ==============================================================
# RECHERCHE ET AFFICHAGE DU MAXIMUM
# ==============================================================
# 1. On trouve l'indice "flat" (linéaire) du maximum
idx_flat_max = np.argmax(values)

# 2. On convertit cet indice linéaire en indices 3D (i, j, k)
idx_3d_max = np.unravel_index(idx_flat_max, values.shape)

# 3. On récupère la valeur et les coordonnées physiques correspondantes
max_val = values[idx_3d_max]
x_max = X[idx_3d_max]
y_max = Y[idx_3d_max]
z_max = Z[idx_3d_max]

print("\n" + "="*40)
print("ANALYSE DU MAXIMUM DE DENSITÉ")
print("="*40)
print(f"Valeur maximale de Delta_rho : {max_val:.6f} e/Å³")
print(f"Indices dans la grille (i, j, k) : {idx_3d_max}")
print(f"Coordonnées cartésiennes (Å)   : x = {x_max:.3f}, y = {y_max:.3f}, z = {z_max:.3f}")
print("="*40)

# ==============================================================
# SAUVEGARDE DE x_flat, y_flat, z_flat, v_flat
# ==============================================================
np.savez('delta_rho_flat.npz',
         x_flat=x_flat,
         y_flat=y_flat,
         z_flat=z_flat,
         v_flat=v_flat)
print("Données plates sauvegardées dans delta_rho_flat.npz")


