import numpy as np

def calcul_final_densities(atom_dict : dict, electronic_density): 
    """Arguments :
     - nparrray : hirshfeld-weight of the atom
     - nparray : electronic density of the molecule"""

    for atom in atom_dict.items():
        weight = atom["hirshfeld_weight"]
        atom["rho_weighted"] = electronic_density * weight

    return atom_dict

if "name" == "__main__":
    atom_dict = {
        1: {"element": "C", "rho": np.array([0.1, 0.2, 0.3]), "hirshfeld_weight": np.array([0.5, 0.6, 0.7])},
        2: {"element": "H", "rho": np.array([0.4, 0.5, 0.6]), "hirshfeld_weight": np.array([0.8, 0.9, 1.0])},
    }
    electronic_density = np.array([0.2, 0.3, 0.4])

    print(calcul_final_densities(atom_dict, electronic_density))

