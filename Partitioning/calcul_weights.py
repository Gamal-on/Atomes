import numpy as np


def hirshfeld_weight(density_atom, density_promolecule):
    """Calculate the Hirshfeld weight for each atom based on the density of the atom and the promolecule.
    Argument : 
    - density_atom : density of the atom
    - density_promolecule : density of the promolecule"""
    poids = np.divide(
    density_atom,
    density_promolecule,
    out=np.zeros_like(density_atom),  
    where=density_promolecule != 0,)

    return poids
