import numpy as np


def hirshfeld_weight(density_atom, density_promolecule):
    poids = np.divide(
    density_atom,
    density_promolecule,
    out=np.zeros_like(density_atom),  
    where=density_promolecule != 0,)

    return poids
