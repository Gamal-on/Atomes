import numpy as np
import re
import sys



#lecture du fichier
def lire_fichier(filepath: str):
    with open(filepath, "r") as f:
        lines = f.readlines()

    gaussienne_vals = []
    mo_coeffs = {}

    current_mo = None
    current_vals = []

# On parcourt ligne par ligne pour extraire les données voulues (exposants et coefficients)
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.upper().startswith("END DATA"): #On arrête la lecture dès qu'on rencontre "END DATA" dans le fichier texte
            break

        if stripped.upper().startswith("EXPONENTS"):
            after = stripped[len("EXPONENTS"):].strip() # On récupère tout ce qui suit "EXPONENTS" sur la même ligne, puis on continue à lire les lignes suivantes pour compléter les exposants jusqu'à ce qu'on rencontre une ligne qui ne contient pas de nombre
            gaussienne_vals.extend(after.split())
            current_mo = None
            continue

        mo_match = re.match(r"^MO\s+(\d+)", stripped, re.IGNORECASE) # On cherche les lignes qui commencent par "MO" (molecular orbital)
        if mo_match:
            if current_mo is not None and current_vals:
                mo_coeffs[current_mo] = np.array([float(v) for v in current_vals])
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
        mo_coeffs[current_mo] = np.array([float(v) for v in current_vals]) # On ajoute les coefficients de la dernière MO lue, si elle existe

    if not gaussienne_vals: # Si on n'a pas trouvé de ligne "EXPONENTS", c'est une erreur
        raise ValueError("Mot-clé 'EXPONENTS' introuvable dans le fichier.")
    if not mo_coeffs: # Si on n'a pas trouvé de ligne "MO i", c'est une erreur
        raise ValueError("Aucun bloc 'MO i' trouvé dans le fichier.")

    gaussienne = np.array([float(v) for v in gaussienne_vals]) # On convertit les exposants en tableau numpy de float
    return gaussienne, mo_coeffs


def get_mo(filepath: str, mo_number: int): # On demande à l'utilisateur de spécifier le numéro de la MO qu'il souhaite récupérer, et on vérifie que ce numéro existe bien dans les données extraites du fichier

    """

    filepath  : chemin vers le fichier texte
    mo_number : numéro entier de la MO souhaitée

    """
    gaussienne, mo_coeffs = lire_fichier(filepath)  

    if mo_number not in mo_coeffs: 
        mos_dispo = sorted(mo_coeffs.keys())
        raise KeyError(
            f"MO {mo_number} introuvable. MO disponibles : {mos_dispo}"
        )

    coeff_lin = mo_coeffs[mo_number] 

    if gaussienne.shape != coeff_lin.shape:
        raise ValueError(
            f"Dimensions incohérentes : {len(gaussienne)} exposants "
            f"mais {len(coeff_lin)} coefficients pour MO {mo_number}."
        )

    return gaussienne, coeff_lin


def get_all_mo(filepath: str): # On récupère les exposants et les coefficients de toutes les MO

    """

    filepath : chemin vers le fichier texte

    """
    gaussienne, mo_coeffs = lire_fichier(filepath)

    mo_numbers = sorted(mo_coeffs.keys())
    all_coeffs = np.array([mo_coeffs[i] for i in mo_numbers])

    # tableau final : ligne 0 = exposants, lignes suivantes = coefficients de chaque MO
    # taille : (1 + nb_MO, N)
    tableau_final = np.vstack([gaussienne, all_coeffs])

    return gaussienne, all_coeffs, mo_numbers, tableau_final




if __name__ == "__main__":
    get_all_mo()

