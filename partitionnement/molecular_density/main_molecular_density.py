import matplotlib
import numpy as np

matplotlib.use('TkAgg')
import argparse

from data_read import get_all_mo
from graphical import plot_orbital_3d
from molecular_orbital import calculate_molecular_orbitals, molecular_orbitals


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
    Returns : 
    - orbitals : 4D array of the molecular orbitals (n_orbitals, nx, ny, nz)
    - densities : 4D array of the probability densities (n_or
    """
    gauss_exponents, mo_coefficients, _, powers, centers_geom = get_all_mo(fichier)
    orbitals = calculate_molecular_orbitals(centers_geom, gauss_exponents, powers, mo_coefficients, x, y, z) 
    densities = np.square(orbitals)
    return orbitals, densities


def main():
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str)

    args = parser.parse_args()
    
    x = np.linspace(-1e-1, 2e-1, 80)
    y = np.linspace(-1e-1, 2e-1, 80)
    z = np.linspace(-1e-1, 2e-1, 80)

    mo, densities = molecular_orbital_calcul_3d(x, y, z, args.filepath)
    print(densities)


if __name__ == "__main__":
    main()