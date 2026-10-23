import os
import sys
import numpy as np

if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

def load_hkl(filepath):
    """
    Lit un fichier .hkl et retourne un dictionnaire indexé par (h,k,l)
    contenant [F_squared_meas, F_squared_sigma].
    """
    hkl_data = {}
    
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    
    # On cherche le début des données (après l'en-tête)
    data_start = False
    for line in lines:
        line = line.strip()
        
        # Détection de la fin de l'en-tête CIF/HKL
        if line.startswith('_refln_F_squared_sigma'):
            data_start = True
            continue
        
        if not data_start:
            continue
            
        # Ligne vide ou commentaire
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 5:
            continue
        
        try:
            h = int(parts[0])
            k = int(parts[1])
            l = int(parts[2])
            f2_meas  = float(parts[3])
            f2_sigma = float(parts[4])
            
            hkl_data[(h, k, l)] = [f2_meas, f2_sigma]
        except ValueError:
            continue
    
    return hkl_data


def query_hkl(hkl_data, h, k, l):
    """
    Interroge le dictionnaire pour un triplet (h,k,l) donné.
    Retourne [F_squared_meas, F_squared_sigma] ou None si absent.
    """
    key = (h, k, l)
    if key in hkl_data:
        return hkl_data[key]
    else:
        print(f"Réflexion ({h},{k},{l}) non trouvée dans les données.")
        return None


# --- Utilisation ---
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "C:\\Users\\damie\\Documents\\Projet Wigner Function\\Données facteurs de structures\\FS expérimentaux.hkl"
        if not os.path.exists(filepath):
            # Chemin alternatif dans le dépôt Atome
            repo_path = os.path.join(os.path.dirname(__file__), "..", "Data_FS", "FS expérimentaux.hkl")
            if os.path.exists(repo_path):
                filepath = os.path.abspath(repo_path)
    
    # Chargement
    hkl_data = load_hkl(filepath)
    print(hkl_data)
    print(f"{len(hkl_data)} réflexions chargées.\n")
    
    # Exemples d'interrogation
    for (h, k, l) in [(0, 0, -6), (1, -1, -1), (-1, -1, 0)]:
        result = query_hkl(hkl_data, h, k, l)
        if result:
            print(f"({h:3d},{k:3d},{l:3d}) -> F²={result[0]:.2f}  σ={result[1]:.2f}")