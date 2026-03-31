import numpy as np
import argparse
from poids import get_all_mo
from molecular_orbital import molecular_orbitals

def molecular_orbital_calcul(position, fichier): 
    gaussian_exponents, coeff_molecular_orbital, _, _ = get_all_mo(fichier) 
    molecular_orbital_somme = molecular_orbitals(position, coeff_molecular_orbital, gaussian_exponents)
    return molecular_orbital_somme

def main():
    parser = argparse.ArgumentParser(description="Calcul des orbitales moléculaires.")
    
    parser.add_argument("filepath", type=str, help="Chemin vers le fichier d'entrée")

    args = parser.parse_args()
    
    position_test = np.arange(0, 4, 1e-1) # Attention : valeurs multiples de rayon de Bohr

    try:
        score = molecular_orbital_calcul(position_test, args.filepath)
        print(f"Résultat pour le fichier '{args.filepath}' :")
        print(score)
    except Exception as e:
        print(f"Une erreur est survenue : {e}")

if __name__ == "__main__":
    main()