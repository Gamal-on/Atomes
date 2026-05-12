# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:13:10 2026

@author: damie
"""

# ==============================================================
# SCRIPT D'AFFICHAGE SEUL — charge les données sans recalculer
# ==============================================================
import numpy as np
import plotly.graph_objects as go

data   = np.load('delta_rho_flat.npz')
x_flat = data['x_flat']
y_flat = data['y_flat']
z_flat = data['z_flat']
v_flat = data['v_flat']

print(f"Données chargées — {len(v_flat)} points")
print(f"  Δρ max : {v_flat.max():.4f} e/Å³")
print(f"  Δρ min : {v_flat.min():.4f} e/Å³")

# ==============================================================
# AFFICHAGE 3D
# ==============================================================

import plotly.graph_objects as go

# --- Sous-échantillonnage ---
n_max = 50_000  # nombre de points affichés — augmente si ton PC le permet
if len(v_flat) > n_max:
    idx = np.random.choice(len(v_flat), size=n_max, replace=False)
    x_plot = x_flat[idx]
    y_plot = y_flat[idx]
    z_plot = z_flat[idx]
    v_plot = v_flat[idx]
    print(f"Sous-échantillonnage : {n_max} points affichés sur {len(v_flat)}")
else:
    x_plot, y_plot, z_plot, v_plot = x_flat, y_flat, z_flat, v_flat

# Échelle symétrique centrée sur 0
v_abs_max = np.max(np.abs(v_flat))   # sur TOUS les points, pas juste l'échantillon

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=x_plot,
    y=y_plot,
    z=z_plot,
    mode='markers',
    marker=dict(
        size=2,
        color=v_plot,
        colorscale=[
            [0.000, '#053061'],
            [0.083, '#1a5b8e'],
            [0.167, '#2166ac'],
            [0.250, '#4393c3'],
            [0.333, '#92c5de'],
            [0.400, '#c6e4f0'],
            [0.450, '#e8f4f8'],
            [0.500, '#ffffff'],   # blanc — Δρ = 0
            [0.550, '#fde8d8'],
            [0.600, '#f4a582'],
            [0.667, '#d6604d'],
            [0.750, '#b2182b'],
            [0.833, '#8b0f1e'],
            [0.917, '#650a14'],
            [1.000, '#3d0008'],
        ],
        cmin=-v_abs_max,
        cmax= v_abs_max,
        opacity=0.8,
        colorbar=dict(
            title='Δρ (e/Å³)',
            tickformat='.3f',
            lenmode='fraction',
            len=0.75,
        ),
    ),
    hovertemplate=(
        'x : %{x:.2f} Å<br>'
        'y : %{y:.2f} Å<br>'
        'z : %{z:.2f} Å<br>'
        'Δρ : %{marker.color:.4f} e/Å³<extra></extra>'
    )
))

fig.update_layout(
    title=dict(
        text='Carte de densité résiduelle Δρ = ρ_stat − ρ_dyn',
        font=dict(size=16)
    ),
    scene=dict(
        aspectmode='cube',
        xaxis=dict(title='x (Å)'),
        yaxis=dict(title='y (Å)'),
        zaxis=dict(title='z (Å)'),
    ),
    height=800,
    paper_bgcolor='#1a1a1a',   # fond sombre pour mieux voir le blanc central
    scene_bgcolor='#1a1a1a',
)

import plotly.io as pio
pio.renderers.default = 'browser'
fig.show()

# --- Résumé console ---
print("\n--- RÉSUMÉ DENSITÉ RÉSIDUELLE ---")
print(f"  Δρ max : {v_flat.max():.4f} e/Å³")
print(f"  Δρ min : {v_flat.min():.4f} e/Å³")
print(f"  Δρ moy : {v_flat.mean():.4f} e/Å³")
print(f"  Points totaux  : {len(v_flat)}")
print(f"  Points affichés: {len(v_plot)}")