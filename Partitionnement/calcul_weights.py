import numpy as np

from read_weights import density_promolecule, density_carbon_1, density_oxygen_2, density_azote_3
from read_weights import density_hydrogen_4, density_hydrogen_5, density_azote_6, density_hydrogen_7, density_hydrogen_8
from read_weights import grid

def hirshfeld_weight(density_atom, density_promolecule):
    poids = np.divide(
    density_atom,
    density_promolecule,
    out=np.zeros_like(density_atom),  
    where=density_promolecule != 0,)

    return poids

def main() : 
    print(len(density_promolecule), len(density_carbon_1))
    weight_carbon = hirshfeld_weight(density_carbon_1, density_promolecule)
    print(weight_carbon)

if __name__ == "__main__":
    main()
