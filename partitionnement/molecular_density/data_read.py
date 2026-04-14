import pathlib
import re

import numpy as np


def lire_fichier(filepath: str):
    with pathlib.Path(filepath).open() as f:
        lines = f.readlines()

    gaussienne_vals = []
    type_vals = []  # <--- Ajout pour les types
    centre_vals = []  # <--- Ajout pour les centres
    mo_coeffs = {}

    current_mo = None
    current_vals = []
    atomes_coords = {} 

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        atome_match = re.search(r"CENTRE\s+(\d+)\)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)", line)
        if atome_match:
            idx = int(atome_match.group(1))
            coords = [float(atome_match.group(2)), float(atome_match.group(3)), float(atome_match.group(4))]
            atomes_coords[idx] = coords

        if stripped.upper().startswith("END DATA"):
            break

        if stripped.upper().startswith("EXPONENTS"):
            after = stripped[len("EXPONENTS") :].strip()
            gaussienne_vals.extend(after.split())
            current_mo = None
            continue

        # --- Ajout : Capture des types ---
        if stripped.upper().startswith("TYPE ASSIGNMENTS"):
            after = stripped[len("TYPE ASSIGNMENTS") :].strip()
            type_vals.extend(after.split())
            continue

        # --- Ajout : Capture des centres ---
        if stripped.upper().startswith("CENTRE ASSIGNMENTS"):
            after = stripped[len("CENTRE ASSIGNMENTS") :].strip()
            centre_vals.extend(after.split())
            continue

        mo_match = re.match(r"^MO\s+(\d+)", stripped, re.IGNORECASE)
        if mo_match:
            if current_mo is not None and current_vals:
                mo_coeffs[current_mo] = np.array(
                    [float(v) for v in current_vals]
                )
            current_mo = int(mo_match.group(1))
            current_vals = []
            continue

        if current_mo is not None:
            tokens = stripped.split()
            try:
                [float(t) for t in tokens]
                current_vals.extend(tokens)
            except ValueError:
                pass
            continue

    if current_mo is not None and current_vals:
        mo_coeffs[current_mo] = np.array([float(v) for v in current_vals])

    # --- Traduction des types en puissances l, m, n ---
    TYPE_TO_POWERS = {
        1: [0, 0, 0],  # s
        2: [1, 0, 0],  # x
        3: [0, 1, 0],  # y
        4: [0, 0, 1],  # z
        5: [2, 0, 0],  # xx
        6: [0, 2, 0],  # yy
        7: [0, 0, 2],  # zz
        8: [1, 1, 0],  # xy
        9: [1, 0, 1],  # xz
        10: [0, 1, 1],  # yz
    }

    types = [int(v) for v in type_vals]
    powers = np.array([TYPE_TO_POWERS[t] for t in types])

    centers_indices = [int(v) for v in centre_vals]

    gaussienne = np.array([float(v) for v in gaussienne_vals])
    return gaussienne, mo_coeffs, powers, centers_indices, atomes_coords


def get_all_mo(filepath: str):
    """
    Returns : 
    - gaussian : nd array  (36,) with the exponents
    all_coeffs : np.ndarray with shape (m, 36) of the orbitals coefficients
    mo_numbers : list of int, identifying the orbitals
    powers : np.ndarray 2D with shape (36, 3). For  each gaussian, it gives its quantical numbers [l,n,m]
    centers_geom : np.ndarray (36, 3), giving the center of the gaussian, for each one
    """
    gaussienne, mo_coeffs, powers, centers_indices, atomes_coords = lire_fichier(filepath)
    
    centers_geom = np.array([atomes_coords[int(idx)] for idx in centers_indices])
    
    mo_numbers = sorted(mo_coeffs.keys())
    all_coeffs = np.array([mo_coeffs[i] for i in mo_numbers])
    
    return gaussienne, all_coeffs, mo_numbers, powers, centers_geom


if __name__ == "__main__":
    pass
