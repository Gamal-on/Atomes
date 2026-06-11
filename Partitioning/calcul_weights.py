import numpy as np

def hirshfeld_weight(atom_dict: dict, density_promolecule: np.ndarray) :
    """Calculate the Hirshfeld weight for each atom based on the density of the atom and the promolecule.
    Modifies atom_dict in place by adding a 'hirshfeld_weight' entry for each atom.
    Arguments :
    - atom_dict : dictionnaire {atom_num: {'element': str, 'rho': np.ndarray, ...}}
    - density_promolecule : densité de la promolécule (np.ndarray)"""
    
    for atom_num, atom_data in atom_dict.items():
        rho = atom_data['rho']
        poids = np.divide(
            rho,
            density_promolecule,
            out=np.zeros_like(rho),
            where=density_promolecule != 0,
        )
        atom_dict[atom_num]['hirshfeld_weight'] = poids

    return atom_dict