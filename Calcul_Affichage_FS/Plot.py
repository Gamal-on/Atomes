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


# Échelle symétrique centrée sur 0 pour que le blanc soit exactement à Δρ = 0
v_abs_max = np.max(np.abs(v_flat))

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=x_flat,
    y=y_flat,
    z=z_flat,
    mode='markers',
    marker=dict(
        size=3,
        color=v_flat,
        colorscale=[
    [0.000, '#053061'],   # bleu très foncé
    [0.083, '#1a5b8e'],
    [0.167, '#2166ac'],   # bleu foncé
    [0.250, '#4393c3'],
    [0.333, '#92c5de'],   # bleu clair
    [0.400, '#c6e4f0'],
    [0.450, '#e8f4f8'],   # bleu très pâle
    [0.500, '#ffffff'],   # blanc — Δρ = 0
    [0.550, '#fde8d8'],   # rouge très pâle
    [0.600, '#f4a582'],
    [0.667, '#d6604d'],   # rouge clair
    [0.750, '#b2182b'],   # rouge foncé
    [0.833, '#8b0f1e'],
    [0.917, '#650a14'],
    [1.000, '#3d0008'],   # rouge très foncé
    ],
        cmin=-v_abs_max,         # ← symétrie autour de 0
        cmax= v_abs_max,
        opacity=0.7,
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
)

import plotly.io as pio
pio.renderers.default = 'browser'
fig.show()

# --- Résumé console ---
print("\n--- RÉSUMÉ DENSITÉ RÉSIDUELLE ---")
print(f"  Δρ max : {v_flat.max():.4f} e/Å³")
print(f"  Δρ min : {v_flat.min():.4f} e/Å³")
print(f"  Δρ moy : {v_flat.mean():.4f} e/Å³")