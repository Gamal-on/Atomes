import numpy as np

def load_electronic_density(filepath: str, atom_dict: dict) -> None:
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # En-tête
    header = lines[0].split()
    name, number = header[0], int(header[1])

    # Vérification de cohérence
    assert number in atom_dict and atom_dict[number]['element'] == name, (
        f"Atome ({name}, {number}) incohérent avec atom_dict"
    )

    # Data
    # Only the density column is kept: the grid points are the same for every atom
    atom_dict[number]['rho'] = np.loadtxt(filepath, skiprows=1, usecols=3)

    return atom_dict