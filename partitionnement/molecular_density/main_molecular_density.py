import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
from data_read import get_all_mo
from molecular_orbital import molecular_orbitals, molecular_orbitals_3d
import plotly.graph_objects as go
from graphical import plot_orbital_3d, plot_result_1d

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


def main():
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str, help="Chemin vers le fichier d'entrée")

    args = parser.parse_args()
    
    x = np.arange(-2,2,1e-1)
    y = np.linspace(-6, 6, 80)
    z = np.linspace(-6, 6, 80)

    densities, mo = molecular_orbital_calcul(x, args.filepath)
    print(mo[0])


if __name__ == "__main__":
    main()