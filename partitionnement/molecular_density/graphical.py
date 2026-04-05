import matplotlib
import numpy as np

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def plot_result_1d(position, density):
    """Plot the probility density of each molecular orbital"""
    plt.figure()
    plt.grid(True)
    for i in range(0, len(density)): 
        plt.plot(position, density[i], ".", label=f"orbital {i}")
        plt.xlabel("Position, en Bohr radius")
        plt.ylabel("Densités électroniques")
        i += 1 
    plt.legend()
    plt.show()


def plot_orbital_3d(densities_4d, orbital_index, x_range, y_range, z_range,
                    level=None, title="Orbital-3D"):
    """
    Plot une orbitale 3D (une seule) à partir d'un tableau 4D.
    Parameters  : 
    - densities_4d : ndarray
    - orbital_index : int
    - x_range, y_range, z_range : 1D array
    - title : str or None
    """
    
    orbital = densities_4d[orbital_index]

    X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')

    # Normalisation
    orbital_vis = orbital / np.max(np.abs(orbital))

    if level is None:
        level = 0.05 

    fig = go.Figure()

    # Partie positive
    fig.add_trace(go.Isosurface(
        x=X.ravel(),
        y=Y.ravel(),
        z=Z.ravel(),
        value=orbital_vis.ravel(),
        isomin=level,
        isomax=1.0,
        surface_count=1,
        opacity=0.6,
        caps=dict(x_show=False, y_show=False, z_show=False),
        showscale=False,
    ))

    # Partie négative
    fig.add_trace(go.Isosurface(
        x=X.ravel(),
        y=Y.ravel(),
        z=Z.ravel(),
        value=orbital_vis.ravel(),
        isomin=-1.0,
        isomax=-level,
        surface_count=1,
        opacity=0.6,
        caps=dict(x_show=False, y_show=False, z_show=False),
        showscale=False,
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='x',
            yaxis_title='y',
            zaxis_title='z',
        )
    )

    fig.show()


def main():
    pass


if __name__ == "__main__":
    main()