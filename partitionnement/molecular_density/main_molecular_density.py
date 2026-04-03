import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
from data_read import get_all_mo
from molecular_orbital import molecular_orbitals, molecular_orbitals_3d
import plotly.graph_objects as go

def molecular_orbital_calcul(position, fichier): 
    """
    Compute the probability densities of each molecular orbital,
    for a given molecule, as a function of position.
    Arguments : 
    - position : 1D-ndarray 
    Returns : 
    - probability_denisties : 2D array of the proba densities for each orbital
    - molecular_orbital : 2D array of the OM 
    """
    gaussian_exponents, coeff_molecular_orbital, _, _ = get_all_mo(fichier) 
    molecular_orbital_somme = molecular_orbitals(position, coeff_molecular_orbital, gaussian_exponents)
    probility_densities = np.square(molecular_orbital_somme)
    return probility_densities, molecular_orbital_somme

def molecular_orbital_calcul_3d(x, y, z, fichier): 
    """
    Compute the probability densities of each molecular orbital,
    for a given molecule, as a function of position.
    Arguments : 
    - position : 1D-ndarray 
    Returns : 
    - probability_denisties : 2D array of the proba densities for each orbital
    - molecular_orbital : 2D array of the OM 
    """
    gaussian_exponents, coeff_molecular_orbital, _, _ = get_all_mo(fichier) 
    molecular_orbital_somme = molecular_orbitals_3d(x, y, z, coeff_molecular_orbital, gaussian_exponents)
    probility_densities = np.square(molecular_orbital_somme)
    return probility_densities, molecular_orbital_somme

def plot_result_1d(position, density) :
    """Plot the probility density of each molecular orbital"""
    plt.figure()
    plt.grid(True)
    for i in range(0, len(density)) : 
        plt.plot(position, density[i], ".", label=f"orbital {i}")
        plt.xlabel("Position, en Bohr radius")
        plt.ylabel("Densités électroniques")
        i+=1 
    plt.legend()
    plt.show()


def plot_orbital_3d(density, x_range, y_range, z_range, level=None, title="Orbital 3D"):

    # Grille 3D
    X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')

    # Valeur de seuil
    if level is None:
        level = 0.3 * np.max(np.abs(density))

    fig = go.Figure(data=go.Isosurface(
        x=X.ravel(),
        y=Y.ravel(),
        z=Z.ravel(),
        value=density.ravel(),
        isomin=level,
        isomax=density.max(),
        surface_count=1,
        caps=dict(x_show=False, y_show=False, z_show=False),
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
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str, help="Chemin vers le fichier d'entrée")

    args = parser.parse_args()
    
    x = np.arange(-2, 2, 1e-1) # Attention : valeurs multiples de rayon de Bohr
    y = np.arange(-2, 2, 1e-1)
    z = np.arange(-2, 2, 1e-1)

    densities, mo = molecular_orbital_calcul_3d(x, y, z, args.filepath)
    #plot_result(position, densities)
    print(densities)



if __name__ == "__main__":
    main()