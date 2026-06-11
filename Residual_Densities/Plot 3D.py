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
import plotly.io as pio

data   = np.load('delta_rho_flat.npz')
x_flat = data['x_flat']
y_flat = data['y_flat']
z_flat = data['z_flat']
v_flat = data['v_flat']
print(f"Données chargées — {len(v_flat)} points")
print(f"  Δρ max : {v_flat.max():.4f} e/Å³")
print(f"  Δρ min : {v_flat.min():.4f} e/Å³")

# ==============================================================
# CONFIGURATION DE L'ÉCHELLE DE COULEUR
# ==============================================================
# Seuil Alpha car la densité résiduelle est très élevé en certains points, 
# ce qui fausse l'échelle pour les plus petites valeurs.
alpha = 10  

# ==============================================================
# AFFICHAGE 3D
# ==============================================================

# --- Sous-échantillonnage --- (On prend moins de points car la fenêtre du navigateur n'arrive 
# pas à charger le graphe sinon.)
n_max = 200_000
if len(v_flat) > n_max:
    idx = np.random.choice(len(v_flat), size=n_max, replace=False) # On choisit aléatoirement les points sélectionnés
    x_plot = x_flat[idx]
    y_plot = y_flat[idx]
    z_plot = z_flat[idx]
    v_plot = v_flat[idx]
    print(f"Sous-échantillonnage : {n_max} points affichés sur {len(v_flat)}")
else:
    x_plot, y_plot, z_plot, v_plot = x_flat, y_flat, z_flat, v_flat

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
            [0.00, '#0000ff'],   # bleu pur  — pour tout Δρ <= -alpha
            [0.50, '#ffffff'],   # blanc     — pour Δρ = 0
            [1.00, '#ff0000'],   # rouge pur — pour tout Δρ >= alpha
        ],
        cmin=-alpha,  # Borne inférieure fixe de l'échelle continue
        cmax=alpha,   # Borne supérieure fixe de l'échelle continue
        opacity=0.8,
        colorbar=dict(
            title=f'Δρ (saturé à ±{alpha})',
            tickformat='.3f',
            lenmode='fraction',
            len=0.75,
            bgcolor='#d0d0d0',
            bordercolor='#aaaaaa',
            borderwidth=1,
        ),
    ),
    hovertemplate=(
        'x : %{x:.2f} Å<br>'
        'y : %{y:.2f} Å<br>'
        'z : %{z:.2f} Å<br>'
        'Δρ réels : %{marker.color:.4f} e/Å³<extra></extra>'
    )
))

fig.update_layout(
    title=dict(
        text=f'Carte de densité résiduelle Δρ (Échelle continue resserrée à ±{alpha})',
        font=dict(size=16, color='black')
    ),
    scene=dict(
        aspectmode='cube',
        xaxis=dict(title='x (Å)', backgroundcolor='#c8c8c8', gridcolor='#aaaaaa', zerolinecolor='#888888'),
        yaxis=dict(title='y (Å)', backgroundcolor='#c8c8c8', gridcolor='#aaaaaa', zerolinecolor='#888888'),
        zaxis=dict(title='z (Å)', backgroundcolor='#c8c8c8', gridcolor='#aaaaaa', zerolinecolor='#888888'),
    ),
    height=800,
    paper_bgcolor='#e8e8e8',   # gris clair — contraste avec bleu, blanc ET rouge
    scene_bgcolor='#e8e8e8',
)

pio.renderers.default = 'browser'
fig.show()

# --- Résumé console ---
print("\n--- RÉSUMÉ DENSITÉ RÉSIDUELLE ---")
print(f"  Seuil alpha choisi : ±{alpha} e/Å³")
print(f"  Δρ max : {v_flat.max():.4f} e/Å³")
print(f"  Δρ min : {v_flat.min():.4f} e/Å³")
print(f"  Δρ moy : {v_flat.mean():.4f} e/Å³")
print(f"  Points totaux  : {len(v_flat)}")
print(f"  Points affichés: {len(v_plot)}")