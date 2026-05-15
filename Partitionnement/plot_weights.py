import numpy as np
from plotly import data
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import plotly.express as px

from weights_arrays import grid, weight_carbon_1, weight_oxygen_2
from weights_arrays import weight_azote_3, weight_hydrogen_4, weight_hydrogen_5, weight_azote_6, weight_hydrogen_7, weight_hydrogen_8


def plot_weight(grid, weight):
    fig = px.scatter_3d(
    x=grid[:, 0],
    y=grid[:, 1],
    z=grid[:, 2],
    color=weight,  # La couleur dépend du poids
    color_continuous_scale="Viridis",  # Échelle de couleur
    labels={"x": "X", "y": "Y", "z": "Z", "color": "Poids"},
    title="Poids de Hirshfeld de l'Urée en 3D",)
    fig.update_traces(marker=dict(size=2))

    fig.show()


def main():
    plot_weight(grid, weight_carbon_1)
    plot_weight(grid, weight_oxygen_2)

if __name__ == "__main__":
    main()