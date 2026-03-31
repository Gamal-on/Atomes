import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
from data_read import get_all_mo
from molecular_orbital import molecular_orbitals

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

def plot_result(position, density) :
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

def main():
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str, help="Chemin vers le fichier d'entrée")

    args = parser.parse_args()
    
    position = np.arange(-2, 2, 1e-3) # Attention : valeurs multiples de rayon de Bohr
    
    densities, mo = molecular_orbital_calcul(position, args.filepath)
    plot_result(position, densities)
    #print(score)



if __name__ == "__main__":
    main()