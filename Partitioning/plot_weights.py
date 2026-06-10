import numpy as np
from plotly import data
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import plotly.express as px


def plot_weight(grid, atom_dict):
    atom_id = int(input("Numéro d'atome : "))

    fig = px.scatter_3d(
    x=grid[:, 0],
    y=grid[:, 1],
    z=grid[:, 2],
    color=atom_dict[atom_id]['hirshfeld_weight'],  # La couleur dépend du poids
    color_continuous_scale="Viridis",  # Échelle de couleur
    labels={"x": "X", "y": "Y", "z": "Z", "color": "Poids"},
    title="Poids de Hirshfeld de l'atome d'oxygène (molécule d'urée)",
    )
    fig.update_traces(marker=dict(size=2))

    fig.show()
