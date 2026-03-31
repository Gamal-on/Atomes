import numpy as np
import re
import sys


def lire_fichier(filepath: str):
    with open(filepath, "r") as f:
        lines = f.readlines()

    gaussienne_vals = []
    mo_coeffs = {}

    current_mo = None
    current_vals = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.upper().startswith("END DATA"):
            break

        if stripped.upper().startswith("EXPONENTS"):
            after = stripped[len("EXPONENTS"):].strip()
            gaussienne_vals.extend(after.split())
            current_mo = None
            continue

        mo_match = re.match(r"^MO\s+(\d+)", stripped, re.IGNORECASE)
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
        mo_coeffs[current_mo] = np.array([float(v) for v in current_vals])

    if not gaussienne_vals:
        raise ValueError("Mot-clé 'EXPONENTS' introuvable dans le fichier.")
    if not mo_coeffs:
        raise ValueError("Aucun bloc 'MO i' trouvé dans le fichier.")

    gaussienne = np.array([float(v) for v in gaussienne_vals])
    return gaussienne, mo_coeffs


def get_mo(filepath: str, mo_number: int):
    """
    Retourne les exposants et les coefficients pour une MO donnée.

    Paramètres
    ----------
    filepath  : chemin vers le fichier texte
    mo_number : numéro entier de la MO souhaitée

    Retourne
    --------
    gaussienne : np.ndarray  -- exposants alpha_j
    coeff_lin  : np.ndarray  -- coefficients a_j pour la MO demandée
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


def get_all_mo(filepath: str):
    """
    Retourne les exposants et les coefficients de toutes les MO.

    Paramètres
    ----------
    filepath : chemin vers le fichier texte

    Retourne
    --------
    gaussienne : np.ndarray de shape (N,)
                 exposants alpha_j communs à toutes les MO

    all_coeffs : np.ndarray de shape (nb_MO, N)
                 all_coeffs[i] contient les coefficients a_j de la MO i+1
                 (les lignes sont ordonnées par numéro de MO croissant)

    mo_numbers : list[int]
                 liste des numéros de MO dans le même ordre que all_coeffs
    """
    gaussienne, mo_coeffs = lire_fichier(filepath)

    mo_numbers = sorted(mo_coeffs.keys())
    all_coeffs = np.array([mo_coeffs[i] for i in mo_numbers])

    # Tableau final : ligne 0 = exposants, lignes suivantes = coefficients de chaque MO
    # shape : (1 + nb_MO, N)
    tableau_final = np.vstack([gaussienne, all_coeffs])

    return gaussienne, all_coeffs, mo_numbers, tableau_final


# ------------------------------------------------------------------
# Point d'entrée CLI
# Usage : python parse_gaussian.py <fichier.txt> [numéro_MO]
#   - sans numéro : affiche toutes les MO
#   - avec numéro : affiche uniquement cette MO
# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage : python parse_gaussian.py <fichier.txt> [numéro_MO]")
        sys.exit(1)

    fichier = sys.argv[1]

    try:
        if len(sys.argv) >= 3:
            # --- Une seule MO ---
            num_mo = int(sys.argv[2])
            gaussienne, coeff_lin = get_mo(fichier, num_mo)
            print(f"\n=== MO {num_mo} ===")
            print(f"Nombre de gaussiennes : {len(gaussienne)}")
            print(f"\nExposants  (alpha_j) :\n{gaussienne}")
            print(f"\nCoefficients (a_j)   :\n{coeff_lin}")

        else:
            # --- Toutes les MO ---
            gaussienne, all_coeffs, mo_numbers, tableau_final = get_all_mo(fichier)
            print(f"\nNombre de gaussiennes : {len(gaussienne)}")
            print(f"Exposants (alpha_j) :\n{gaussienne}\n")
            print(f"{'='*60}")
            for idx, mo_num in enumerate(mo_numbers):
                print(f"\nMO {mo_num} — coefficients (a_j) :")
                print(all_coeffs[idx])
                print(f"{'='*60}")
            print(f"\nTableau final (ligne 0 = exposants, lignes suivantes = coefficients MO) :")
            print(f"Shape : {tableau_final.shape}")
            print(tableau_final)

    except (ValueError, KeyError, FileNotFoundError) as e:
        print(f"Erreur : {e}")
        sys.exit(1)


# renvoie un tableau de tableau avec en 1er les coeffs des exposants des gaussiennes,
# et ensuite les coeffs de chaque MO, dans l'ordre des MO (MO1, MO2, ...)