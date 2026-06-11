import numpy as np

def load_hkl(filepath):
    """
    Lit un fichier .hkl et retourne un dictionnaire indexé par (h,k,l)
    contenant [F_squared_meas, F_squared_sigma].
    """
    hkl_data = {}
    
    with open(filepath, 'r') as f:
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
    filepath = "C:\\Users\\damie\\Documents\\Projet Wigner Function\\Données facteurs de structures\\FS expérimentaux.hkl"
    
    # Chargement
    hkl_data = load_hkl(filepath)
    print(hkl_data)
    print(f"{len(hkl_data)} réflexions chargées.\n")
    
    # Exemples d'interrogation
    for (h, k, l) in [(0, 0, -6), (1, -1, -1), (-1, -1, 0)]:
        result = query_hkl(hkl_data, h, k, l)
        if result:
            print(f"({h:3d},{k:3d},{l:3d}) -> F²={result[0]:.2f}  σ={result[1]:.2f}")