# -*- coding: utf-8 -*-
"""
Calcul et affichage 2D de la densité résiduelle Delta_rho
dans le plan contenant la molécule sélectionnée.

SÉLECTION DE LA MOLÉCULE :
  - Type de base : 'A' ou 'B'  (les deux molécules de la maille primitive)
  - Triplet (h, k, l)          : indices de la maille → translation h·a1 + k·a2 + l·a3
  La molécule sélectionnée est la copie de la molécule de base translatée
  de h mailles selon a, k mailles selon b, l mailles selon c.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

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

v = np.sqrt(
    1
    - np.cos(alpha)**2
    - np.cos(beta)**2
    - np.cos(gamma)**2
    + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

V = a * b * c * v

# Vecteurs directs (base cartésienne)
a1 = np.array([a, 0, 0])
a2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
a3 = np.array([
    c * np.cos(beta),
    c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
    c * v / np.sin(gamma)
])

# Vecteurs réciproques
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
# CONVERTISSEURS
# ==============================================================

def frac_to_cart(xf, yf, zf):
    """Convertit des coordonnées fractionnaires en coordonnées cartésiennes."""
    return xf * a1 + yf * a2 + zf * a3


# ==============================================================
# DÉFINITION DES DEUX MOLÉCULES DE BASE DE LA MAILLE PRIMITIVE
#
# Molécule A : carbone C1 en (x=0, y=2.789 Å, z=1.538 Å)
# Molécule B : carbone C2 en (x=2.789 Å, y=0, z=-1.538 Å)
#
# Toutes les positions sont données en coordonnées FRACTIONNAIRES
# dans la maille de référence (h=k=l=0).
# ==============================================================

MOLECULES_BASE = {
    'A': {
        'nom': 'Molécule A — Urée (C en x≈0)',
        'atomes': [
            {'label': 'C1',  'symbol': 'C', 'frac': ( 0.000/a,  2.789/b,  1.538/c)},
            {'label': 'O3',  'symbol': 'O', 'frac': ( 0.000/a,  2.789/b, (c - 1.892)/c)},
            {'label': 'N5',  'symbol': 'N', 'frac': ( 0.807/a, (b - 1.982)/b,  0.836/c)},
            {'label': 'N6',  'symbol': 'N', 'frac': (-0.807/a,  1.982/b,  0.836/c)},
            {'label': 'H9',  'symbol': 'H', 'frac': ( 1.426/a, (b - 1.363)/b,  1.331/c)},
            {'label': 'H10', 'symbol': 'H', 'frac': (-1.426/a,  1.363/b,  1.331/c)},
            {'label': 'H13', 'symbol': 'H', 'frac': ( 0.798/a, (b - 1.991)/b, -0.163/c)},
            {'label': 'H14', 'symbol': 'H', 'frac': (-0.798/a,  1.991/b, -0.163/c)},
        ],
    },
    'B': {
        'nom': 'Molécule B — Urée (C en y≈0)',
        'atomes': [
            {'label': 'C2',  'symbol': 'C', 'frac': ( 2.789/a,  0.000/b, (c - 1.538)/c)},
            {'label': 'O4',  'symbol': 'O', 'frac': ( 2.789/a,  0.000/b,  1.892/c)},
            {'label': 'N7',  'symbol': 'N', 'frac': ((a - 1.982)/a, -0.807/b, (c - 0.836)/c)},
            {'label': 'N8',  'symbol': 'N', 'frac': ( 1.982/a,  0.807/b, (c - 0.836)/c)},
            {'label': 'H11', 'symbol': 'H', 'frac': ((a - 1.363)/a, -1.426/b, (c - 1.331)/c)},
            {'label': 'H12', 'symbol': 'H', 'frac': ( 1.363/a,  1.426/b, (c - 1.331)/c)},
            {'label': 'H15', 'symbol': 'H', 'frac': ((a - 1.991)/a, -0.798/b, (0.163 + c)/c)},
            {'label': 'H16', 'symbol': 'H', 'frac': ( 1.991/a,  0.798/b, (0.163 + c)/c)},
        ],
    },
}


def get_molecule(base_type, h, k, l):
    """
    Retourne la molécule de type `base_type` ('A' ou 'B') translatée
    de h mailles selon a1, k mailles selon a2, l mailles selon a3.

    Parameters
    ----------
    base_type : 'A' ou 'B'
    h, k, l   : entiers (indices de la maille — peuvent être négatifs)

    Returns
    -------
    dict avec 'nom', 'description', 'atomes' (coordonnées cartésiennes déjà translatées)
    """
    if base_type not in MOLECULES_BASE:
        raise ValueError(f"Type de base inconnu : '{base_type}'. Choisissez 'A' ou 'B'.")

    base = MOLECULES_BASE[base_type]
    # Vecteur de translation en Å
    translation = h * a1 + k * a2 + l * a3

    atomes_traduits = []
    for at in base['atomes']:
        pos_cart = frac_to_cart(*at['frac']) + translation
        atomes_traduits.append({
            'label':  f"{at['label']}",
            'symbol': at['symbol'],
            'cart':   pos_cart,         # coordonnées cartésiennes directement
        })

    return {
        'nom':         f"{base['nom']}  +  ({h},{k},{l})·maille",
        'description': f"Molécule de base {base_type} translatée de h={h}, k={k}, l={l}",
        'atomes':      atomes_traduits,
    }


# ==============================================================
# SÉLECTION DE LA MOLÉCULE
# ==============================================================
print("\n" + "="*60)
print("  DENSITÉ RÉSIDUELLE 2D — SÉLECTION DE MOLÉCULE")
print("="*60)
print("  Deux types de molécules de base dans la maille primitive :")
print("  [A]  Molécule A — C en (x≈0,  y≈2.79 Å, z≈1.54 Å)")
print("  [B]  Molécule B — C en (x≈2.79 Å, y≈0, z≈-1.54 Å)")
print()
print("  Entrez ensuite un triplet (h, k, l) pour sélectionner")
print("  n'importe quelle copie de cette molécule dans le cristal infini.")
print("  Exemples :  (0,0,0) → maille de référence")
print("              (1,0,0) → translatée d'une maille selon a")
print("              (-1,2,0) → translatée de -1a + 2b")
print("="*60)

# --- Choix du type de base ---
while True:
    base_type = input("\nType de molécule de base (A ou B) : ").strip().upper()
    if base_type in MOLECULES_BASE:
        break
    print("  ⚠ Entrez 'A' ou 'B'.")

# --- Choix du triplet (h, k, l) ---
while True:
    try:
        entree = input("Triplet de maille (h k l), séparés par des espaces : ").strip()
        h, k, l = map(int, entree.split())
        break
    except ValueError:
        print("  ⚠ Entrez trois entiers séparés par des espaces, ex : 1 0 -1")

mol_sel = get_molecule(base_type, h, k, l)
print(f"\n✔ Molécule sélectionnée : {mol_sel['nom']}")
print(f"  {mol_sel['description']}")

# ==============================================================
# CALCUL DU PLAN MOYEN CONTENANT LA MOLÉCULE  (PCA / SVD)
# ==============================================================

# Positions cartésiennes des atomes (déjà calculées dans get_molecule)
positions_cart = np.array([at['cart'] for at in mol_sel['atomes']])
centre = positions_cart.mean(axis=0)
positions_centrees = positions_cart - centre

U, S, Vt = np.linalg.svd(positions_centrees)

normale = Vt[2]
u_axis  = Vt[0]
v_axis  = Vt[1]

normale = normale / np.linalg.norm(normale)
u_axis  = u_axis  / np.linalg.norm(u_axis)
v_axis  = v_axis  / np.linalg.norm(v_axis)

print(f"\n  Centre du plan    : ({centre[0]:.3f}, {centre[1]:.3f}, {centre[2]:.3f}) Å")
print(f"  Axe u (plan)      : ({u_axis[0]:.3f}, {u_axis[1]:.3f}, {u_axis[2]:.3f})")
print(f"  Axe v (plan)      : ({v_axis[0]:.3f}, {v_axis[1]:.3f}, {v_axis[2]:.3f})")
print(f"  Normale au plan   : ({normale[0]:.3f}, {normale[1]:.3f}, {normale[2]:.3f})")

# ==============================================================
# GRILLE 2D DANS LE PLAN DE LA MOLÉCULE
# ==============================================================
marge_u = max(a, b) * 1.5    # Å
marge_v = max(b, c) * 1.5    # Å
pas     = 0.05                # Å

u_vals = np.arange(-marge_u, marge_u + pas, pas)
v_vals = np.arange(-marge_v, marge_v + pas, pas)
Nu, Nv = len(u_vals), len(v_vals)

Ug, Vg = np.meshgrid(u_vals, v_vals, indexing='ij')   # (Nu, Nv)
Points_3D = (centre[np.newaxis, np.newaxis, :]
             + Ug[:, :, np.newaxis] * u_axis[np.newaxis, np.newaxis, :]
             + Vg[:, :, np.newaxis] * v_axis[np.newaxis, np.newaxis, :])

Px = Points_3D[:, :, 0]
Py = Points_3D[:, :, 1]
Pz = Points_3D[:, :, 2]

# ==============================================================
# PRÉ-CALCUL DES AMPLITUDES ET PHASES
# ==============================================================
common_keys = (set(FS_théorique_dynamique)
               & set(FS_théorique_statique)) - {(0, 0, 0)}

missing = set(FS_théorique_dynamique) - common_keys
if missing:
    print(f"\n⚠ {len(missing)} réflexions absentes d'un dictionnaire, ignorées.")

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

# ==============================================================
# CALCUL DE LA DENSITÉ RÉSIDUELLE SUR LA GRILLE 2D
# ==============================================================
print(f"\nCalcul sur {Nu} × {Nv} = {Nu * Nv} points (plan 2D)…")

Delta_rho_2D = np.zeros((Nu, Nv), dtype=complex)
total = len(F_hkl)

for i, (f, gx, gy, gz) in enumerate(zip(F_hkl, Gx_arr, Gy_arr, Gz_arr)):
    pct    = (i + 1) / total
    filled = int(40 * pct)
    bar    = '█' * filled + '░' * (40 - filled)
    print(f'\r[{bar}] {pct*100:5.1f}%  réflexion {i+1}/{total}',
          end='', flush=True)

    G_dot_r = -2 * np.pi * (Px * gx + Py * gy + Pz * gz)
    Delta_rho_2D += f * np.exp(1j * G_dot_r)

print(f'\nCalcul terminé.')

values_2D = np.real(Delta_rho_2D) / V   # e/Å³

# ==============================================================
# COORDONNÉES DES ATOMES PROJETÉES DANS LE PLAN
# ==============================================================
atomes_u, atomes_v, atomes_label, atomes_symbol = [], [], [], []

for at in mol_sel['atomes']:
    diff   = at['cart'] - centre
    u_proj = np.dot(diff, u_axis)
    v_proj = np.dot(diff, v_axis)
    atomes_u.append(u_proj)
    atomes_v.append(v_proj)
    atomes_label.append(at['label'])
    atomes_symbol.append(at['symbol'])

# ==============================================================
# AFFICHAGE 2D — CARTE DE DENSITÉ RÉSIDUELLE
# ==============================================================
vmax = np.max(np.abs(values_2D))
vmin = -vmax

fig, ax = plt.subplots(figsize=(9, 8))

levels_fill = np.linspace(vmin, vmax, 256)
levels_pos  = np.linspace(0, vmax, 12)[1:]
levels_neg  = np.linspace(-vmax, 0, 12)[:-1]

cf = ax.contourf(u_vals, v_vals, values_2D.T,
                 levels=levels_fill, cmap='RdBu', extend='both')

ct_pos = ax.contour(u_vals, v_vals, values_2D.T,
                    levels=levels_pos, colors='#CC2222',
                    linewidths=0.6, alpha=0.8)
ax.clabel(ct_pos, fmt='%.3f', fontsize=6, inline=True)

ct_neg = ax.contour(u_vals, v_vals, values_2D.T,
                    levels=levels_neg, colors='#1144AA',
                    linewidths=0.6, alpha=0.8)
ax.clabel(ct_neg, fmt='%.3f', fontsize=6, inline=True)


cbar = fig.colorbar(cf, ax=ax, pad=0.02, shrink=0.85)
cbar.set_label('Δρ (e/Å³)', fontsize=12)
cbar.ax.tick_params(labelsize=9)
ticks = np.linspace(vmin, vmax, 9)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{t:+.3f}' for t in ticks])

# ==============================================================
# ATOMES : projection 2D + liaisons par seuil de distance
# ==============================================================
COLORS_ATOM = {'C': '#888888', 'O': '#FF4444', 'N': '#44CC44',
               'H': '#DDDDDD', 'Ti': '#4488FF'}
SIZES_ATOM  = {'C': 220, 'O': 260, 'N': 240, 'H': 120, 'Ti': 350}

# Seuils de distance 3D pour tracer une liaison (en Å)
SEUIL_LIAISON = 1.75   # liaisons lourds–lourds
SEUIL_XH      = 1.30   # liaisons X–H

atomes_u, atomes_v, atomes_label, atomes_symbol = [], [], [], []

for at in mol_sel['atomes']:
    diff   = at['cart'] - centre
    u_proj = np.dot(diff, u_axis)
    v_proj = np.dot(diff, v_axis)
    atomes_u.append(u_proj)
    atomes_v.append(v_proj)
    atomes_label.append(at['label'])
    atomes_symbol.append(at['symbol'])

# --- Liaisons (seuil sur distance 3D cartésienne) ---
n_at = len(mol_sel['atomes'])
for i in range(n_at):
    for j in range(i + 1, n_at):
        sym_i = mol_sel['atomes'][i]['symbol']
        sym_j = mol_sel['atomes'][j]['symbol']
        d3d   = np.linalg.norm(
                    mol_sel['atomes'][i]['cart']
                    - mol_sel['atomes'][j]['cart'])
        seuil = SEUIL_XH if 'H' in {sym_i, sym_j} else SEUIL_LIAISON
        if d3d < seuil:
            ax.plot([atomes_u[i], atomes_u[j]],
                    [atomes_v[i], atomes_v[j]],
                    color='white', lw=2.0, zorder=4, alpha=0.85,
                    solid_capstyle='round')
            # Contour sombre pour lisibilité sur la colormap
            ax.plot([atomes_u[i], atomes_u[j]],
                    [atomes_v[i], atomes_v[j]],
                    color='#333333', lw=3.2, zorder=3, alpha=0.5,
                    solid_capstyle='round')

# --- Atomes ---
for u_at, v_at, lbl, sym in zip(atomes_u, atomes_v, atomes_label, atomes_symbol):
    col  = COLORS_ATOM.get(sym, '#FFAA00')
    size = SIZES_ATOM.get(sym, 200)
    # Halo sombre pour contraste
    ax.scatter(u_at, v_at, s=size * 1.5, c='#111111',
               zorder=5, linewidths=0)
    ax.scatter(u_at, v_at, s=size, c=col,
               zorder=6, edgecolors='white', linewidths=1.0)
    ax.annotate(lbl, (u_at, v_at),
                textcoords='offset points', xytext=(7, 5),
                fontsize=8, fontweight='bold', color='white',
                zorder=7,
                path_effects=[pe.withStroke(linewidth=2.0,
                                            foreground='#111111')])

# --- Légende des espèces présentes ---
from matplotlib.lines import Line2D
symboles_presents = sorted({at['symbol'] for at in mol_sel['atomes']})
legend_atoms = [
    Line2D([0], [0], marker='o', color='w', label=sym,
           markerfacecolor=COLORS_ATOM.get(sym, '#FFAA00'),
           markeredgecolor='white', markersize=9)
    for sym in symboles_presents
]
ax.legend(handles=legend_atoms, loc='upper right',
          fontsize=9, framealpha=0.6, title='Atomes', title_fontsize=9)


ax.set_title(
    f'Densité résiduelle Δρ — Mol. {base_type}  +  ({h},{k},{l})·maille\n',
    fontsize=12)
ax.set_aspect('equal')
ax.grid(False)
plt.tight_layout()
plt.show()

# ==============================================================
# STATISTIQUES
# ==============================================================
idx_max = np.unravel_index(np.argmax(values_2D), values_2D.shape)
idx_min = np.unravel_index(np.argmin(values_2D), values_2D.shape)

print("\n" + "="*55)
print(f"  STATISTIQUES — {mol_sel['nom']}")
print("="*55)
print(f"  Maximum : {values_2D[idx_max]:+.5f} e/Å³"
      f"  @ u={u_vals[idx_max[0]]:.2f} Å, v={v_vals[idx_max[1]]:.2f} Å")
print(f"  Minimum : {values_2D[idx_min]:+.5f} e/Å³"
      f"  @ u={u_vals[idx_min[0]]:.2f} Å, v={v_vals[idx_min[1]]:.2f} Å")
print(f"  Résolution de la grille : {pas} Å")
print(f"  Taille de la grille     : {Nu} × {Nv} points")
print("="*55)

plt.show()

# ==============================================================
# AFFICHAGE 3D DES ATOMES DE LA MOLÉCULE SÉLECTIONNÉE
# ==============================================================
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (requis pour projection='3d')

fig3d = plt.figure(figsize=(9, 8))
ax3d  = fig3d.add_subplot(111, projection='3d')

# Taille et couleur par symbole
COLORS_3D = {'C': '#888888', 'O': '#FF4444', 'N': '#44CC44',
             'H': '#DDDDDD', 'Ti': '#4488FF'}
SIZES_3D  = {'C': 300, 'O': 350, 'N': 300, 'H': 150, 'Ti': 500}

positions_cart = np.array([at['cart'] for at in mol_sel['atomes']])
xs = positions_cart[:, 0]
ys = positions_cart[:, 1]
zs = positions_cart[:, 2]

# --- Liaison simple : chaque paire d'atomes consécutifs ---
# (même logique que la vue 2D ; pour des liaisons réelles,
#  il faudrait un seuil de distance ou une table de connectivité)
for i in range(len(mol_sel['atomes']) - 1):
    p1 = mol_sel['atomes'][i]['cart']
    p2 = mol_sel['atomes'][i + 1]['cart']
    ax3d.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
              'w-', lw=1.5, alpha=0.6, zorder=3)

# --- Atomes ---
for at in mol_sel['atomes']:
    x, y, z = at['cart']
    sym   = at['symbol']
    col   = COLORS_3D.get(sym, '#FFAA00')
    size  = SIZES_3D.get(sym, 250)
    ax3d.scatter(x, y, z, s=size, c=col, edgecolors='k',
                 linewidths=0.8, zorder=5, depthshade=True)
    ax3d.text(x, y, z, f"  {at['label']}", fontsize=8,
              fontweight='bold', color='white',
              path_effects=[pe.withStroke(linewidth=1.5, foreground='k')])

# --- Axes et titre ---
ax3d.set_xlabel('X  [Å]', fontsize=10)
ax3d.set_ylabel('Y  [Å]', fontsize=10)
ax3d.set_zlabel('Z  [Å]', fontsize=10)
ax3d.set_title(
    f'Vue 3D — {mol_sel["nom"]}\n',
    fontsize=11)

# Repère de la maille (vecteurs a1, a2, a3 en pointillés depuis l'origine)
origine = np.zeros(3)
for vec, lbl, col in [(a1, 'a₁', '#FF8800'), (a2, 'a₂', '#00BBFF'), (a3, 'a₃', '#88FF00')]:
    ax3d.quiver(*origine, *vec, length=1.0, normalize=False,
                color=col, linewidth=1.5, linestyle='dashed', alpha=0.5,
                arrow_length_ratio=0.15)
    tip = origine + vec
    ax3d.text(*tip, f' {lbl}', color=col, fontsize=9)

# Légende manuelle des espèces présentes
from matplotlib.lines import Line2D
symboles_presents = sorted({at['symbol'] for at in mol_sel['atomes']})
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label=sym,
           markerfacecolor=COLORS_3D.get(sym, '#FFAA00'),
           markeredgecolor='k', markersize=10)
    for sym in symboles_presents
]
ax3d.legend(handles=legend_elements, loc='upper left', fontsize=9,
            framealpha=0.4, title='Espèces', title_fontsize=9)

ax3d.set_box_aspect([1, 1, 1])   # aspect cubique
plt.tight_layout()
plt.show()

# ==============================================================
# AFFICHAGE 3D DE TOUTES LES MOLÉCULES (grille de mailles)
# ==============================================================
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.widgets as widgets
from itertools import product as iproduct

# --- Couleurs et tailles par espèce ---
COLORS_3D = {'C': '#888888', 'O': '#FF4444', 'N': '#44CC44',
             'H': '#DDDDDD', 'Ti': '#4488FF'}
SIZES_3D  = {'C': 80, 'O': 110, 'N': 90, 'H': 40, 'Ti': 200}

# --- Paramètre : plage de mailles à afficher ---
# On affichera toutes les mailles (h, k, l) avec h,k,l dans [-N_MAILLES, N_MAILLES]
N_MAILLES = 1   # → 2 molécules de base × (2N+1)³ = 2×27 = 54 molécules pour N=1

SEUIL_LIAISON = {'HX': 1.3, 'default': 1.75}  # Å, pour tracer les liaisons

def get_all_molecules(n_mailles):
    """Génère toutes les molécules A et B pour h,k,l dans [-n, n]."""
    molécules = []
    for base_type in ('A', 'B'):
        for h_, k_, l_ in iproduct(range(-n_mailles, n_mailles + 1),
                                   range(-n_mailles, n_mailles + 1),
                                   range(-n_mailles, n_mailles + 1)):
            mol = get_molecule(base_type, h_, k_, l_)
            molécules.append((base_type, h_, k_, l_, mol))
    return molécules

def get_liaisons(atomes_list):
    """Retourne les paires d'indices formant une liaison (seuil de distance)."""
    liaisons = []
    n = len(atomes_list)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(atomes_list[i]['cart'] - atomes_list[j]['cart'])
            seuil = SEUIL_LIAISON['HX'] if 'H' in {atomes_list[i]['symbol'],
                                                    atomes_list[j]['symbol']} \
                    else SEUIL_LIAISON['default']
            if d < seuil:
                liaisons.append((i, j))
    return liaisons

# ==============================================================
# Construction de la figure
# ==============================================================
fig3d = plt.figure(figsize=(12, 9))
plt.subplots_adjust(left=0.05, right=0.88, bottom=0.18, top=0.95)
ax3d = fig3d.add_subplot(111, projection='3d')

# Couleur des molécules A vs B (teinte globale de la liaison)
BOND_COLOR = {'A': '#AAAAFF', 'B': '#FFAAAA'}
ALPHA_BOND = 0.5

scatter_handles = {}   # pour la légende
all_scatter = []       # objets scatter (pour show/hide)
all_lines   = []       # objets Line3D (pour show/hide)

def draw_all_molecules(n_mailles, highlight_sel=True):
    """Dessine toutes les molécules sur ax3d."""
    ax3d.cla()

    molécules = get_all_molecules(n_mailles)
    scatter_handles.clear()

    for base_type, h_, k_, l_, mol in molécules:
        atomes_list = mol['atomes']
        positions   = np.array([at['cart'] for at in atomes_list])

        # Mise en évidence de la molécule sélectionnée initialement
        is_selected = (base_type == base_type and h_ == h and k_ == k and l_ == l)
        edge_col    = 'yellow' if (highlight_sel and is_selected) else 'k'
        edge_lw     = 2.0      if (highlight_sel and is_selected) else 0.5

        # --- Liaisons ---
        liaisons = get_liaisons(atomes_list)
        for i_at, j_at in liaisons:
            p1 = positions[i_at]
            p2 = positions[j_at]
            ax3d.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                      color=BOND_COLOR[base_type], lw=1.2, alpha=ALPHA_BOND, zorder=3)

        # --- Atomes ---
        for at in atomes_list:
            x, y, z = at['cart']
            sym  = at['symbol']
            col  = COLORS_3D.get(sym, '#FFAA00')
            size = SIZES_3D.get(sym, 80)
            sc = ax3d.scatter(x, y, z, s=size, c=col,
                              edgecolors=edge_col, linewidths=edge_lw,
                              depthshade=True, zorder=5)
            if sym not in scatter_handles:
                scatter_handles[sym] = sc

    # --- Vecteurs de maille ---
    origine = np.zeros(3)
    for vec, lbl, col in [(a1, 'a₁', '#FF8800'),
                           (a2, 'a₂', '#00BBFF'),
                           (a3, 'a₃', '#88FF00')]:
        ax3d.quiver(*origine, *vec, color=col, linewidth=1.5,
                    alpha=0.7, arrow_length_ratio=0.15)
        ax3d.text(*(origine + vec * 1.05), f' {lbl}', color=col, fontsize=9)

    # --- Légende espèces ---
    from matplotlib.lines import Line2D
    symboles = sorted(scatter_handles.keys())
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label=sym,
               markerfacecolor=COLORS_3D.get(sym, '#FFAA00'),
               markeredgecolor='k', markersize=9)
        for sym in symboles
    ]

    ax3d.legend(handles=legend_elements, loc='upper left',
                fontsize=8, framealpha=0.4,
                title_fontsize=8)

    ax3d.set_xlabel('X [Å]', fontsize=9)
    ax3d.set_ylabel('Y [Å]', fontsize=9)
    ax3d.set_zlabel('Z [Å]', fontsize=9)
    ax3d.set_title(
        f'Vue 3D des mailles\n',
        fontsize=10)
    ax3d.set_box_aspect([1, 1, 1])
    fig3d.canvas.draw_idle()

# Premier rendu
draw_all_molecules(N_MAILLES)

