import numpy as np

def load_density_promolecule(filepath: str) -> np.ndarray:
    with open(filepath, 'r') as f:
        lines = f.readlines()

    promolecule_density = np.array([])


    # Data
    data = np.loadtxt(filepath, skiprows=1)
    promolecule_density = np.append(promolecule_density, data[:, 3])

    return promolecule_density