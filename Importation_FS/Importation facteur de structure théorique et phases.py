"""
Extraction des facteurs de structure depuis un fichier de sortie CRYSTAL14.

Le dictionnaire retourné a pour clés des tuples (h, k, l) et pour valeurs
des tableaux numpy [partie_réelle, partie_imaginaire, phase_en_radians].

Usage :
    python extract_structure_factors.py [fichier.out]
"""

import re
import math
import sys
from collections import defaultdict


def parse_structure_factors(filepath, section="static"):
    """
    Parse les facteurs de structure depuis un fichier de sortie CRYSTAL14.

    Paramètres
    ----------
    filepath : str
        Chemin vers le fichier .out CRYSTAL14.
    section : str
        "static"  → X-RAY STATIC STRUCTURE FACTORS
        "dynamic" → X-RAY DYNAMIC STRUCTURE FACTORS
        "both"    → renvoie un dict avec les deux clés "static" et "dynamic"

    Retourne
    --------
    dict  {(h, k, l): [re, im, phase_rad]}   (ou dict de dicts si section="both")
    """

    # En-tête qui marque le début de chaque section
    SECTION_HEADERS = {
        "static":  "X-RAY STATIC STRUCTURE FACTORS",
        "dynamic": "X-RAY DYNAMIC STRUCTURE FACTORS",
    }

    # Ligne de colonnes : H  K  L  PX  PY  PZ  REAL PART  IMAG. PART  THEORETICAL
    #
    # Deux cas particuliers dans le format CRYSTAL14 :
    #   1. Indices collés : quand deux indices négatifs se suivent, il n'y a pas
    #      d'espace entre eux, ex : " -4-10 10" au lieu de " -4 -10  10"
    #   2. Grand indice H négatif (≤ -10) : la ligne commence sans espace,
    #      ex : "-10  0  0 ..." au lieu de "  -10  0  0 ..."
    #
    # On résout les deux en autorisant \s* en début de ligne et en rendant
    # l'espace entre les indices optionnel (chaque indice peut coller au suivant
    # si ce dernier commence par '-').
    DATA_RE = re.compile(
        r"^\s*(-?\d+)\s*(-?\d+)\s*(-?\d+)"           # H K L (espaces optionnels)
        r"\s+([+-]?\d+\.\d+)"                          # PX
        r"\s+([+-]?\d+\.\d+)"                          # PY
        r"\s+([+-]?\d+\.\d+)"                          # PZ
        r"\s+([+-]?\d+\.\d+E[+-]?\d+)"                # REAL PART
        r"\s+([+-]?\d+\.\d+E[+-]?\d+)"                # IMAG. PART
        r"\s+([+-]?\d+\.\d+E[+-]?\d+)"                # THEORETICAL
    )

    # Lecture du fichier
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    def _extract_section(target_header):
        """Extrait les données d'une section donnée."""
        result = {}
        in_section = False

        for line in lines:
            # Détecter le début de la section
            if target_header in line:
                in_section = True
                result = {}   # réinitialise si la section apparaît plusieurs fois
                continue

            # Fin de section : une nouvelle bannière *** en dehors de la section
            if in_section and "***" in line and target_header not in line:
                # On continue quand même (des séparateurs internes existent)
                pass

            # Détecter début d'une AUTRE section principale → arrêter
            if in_section:
                for key, hdr in SECTION_HEADERS.items():
                    if hdr in line and hdr != target_header:
                        in_section = False
                        break

            if not in_section:
                continue

            m = DATA_RE.match(line)
            if m:
                h, k, l = int(m.group(1)), int(m.group(2)), int(m.group(3))
                re_val = float(m.group(7))
                im_val = float(m.group(8))
                phase  = math.atan2(im_val, re_val)   # en radians
                result[(h, k, l)] = [re_val, im_val, phase]

        return result

    if section == "both":
        return {
            "static":  _extract_section(SECTION_HEADERS["static"]),
            "dynamic": _extract_section(SECTION_HEADERS["dynamic"]),
        }
    else:
        return _extract_section(SECTION_HEADERS[section])


# ── Exemple d'utilisation ──────────────────────────────────────────────────────

if __name__ == "__main__":
    filepath = "C:\\Users\\damie\\Documents\\Projet Wigner Function\\Données facteurs de structures\\FS Théorique + Phase.out"

    print(f"Lecture du fichier : {filepath}\n")

    # Extraire les deux sections
    data = parse_structure_factors(filepath, section="both")

    for section_name, sf_dict in data.items():
        print(f"=== Section : {section_name.upper()} ===")
        print(f"  Nombre de réflexions : {len(sf_dict)}\n")

        # Afficher les 5 premières entrées
        print("  Exemples (h, k, l) → [Re, Im, Phase (rad)]")
        for i, ((h, k, l), vals) in enumerate(sf_dict.items()):
            if i >= 5:
                break
            re_v, im_v, ph = vals
            print(f"    ({h:3d},{k:3d},{l:3d})  Re={re_v:+.4E}  Im={im_v:+.4E}  φ={ph:+.4f} rad  ({math.degrees(ph):+.2f}°)")
        print()

    
# print(data["dynamic"])
# print(data["static"])


