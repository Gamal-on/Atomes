import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
from poids import get_all_mo
from molecular_orbital import molecular_orbitals

def molecular_orbital_calcul(position, fichier): 
    gaussian_exponents, coeff_molecular_orbital, _, _ = get_all_mo(fichier) 
    molecular_orbital_somme = molecular_orbitals(position, coeff_molecular_orbital, gaussian_exponents)
    return molecular_orbital_somme

def plot_result(position, orbital) :
    plt.figure()
    plt.grid(True)
    for i in orbital : 
        plt.plot(position, i, ".", label=i)
        plt.xlabel("Position, en Bohr radius")
        plt.ylabel("Densités électroniques")
        i+=1 
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str, help="Chemin vers le fichier d'entrée")

    args = parser.parse_args()
    
    position = np.arange(0, 3, 1e-3) # Attention : valeurs multiples de rayon de Bohr
    
    score = molecular_orbital_calcul(position, args.filepath)
    plot_result(position, score)
    print(score)



if __name__ == "__main__":
    main()