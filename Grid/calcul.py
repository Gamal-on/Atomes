import numpy as np
from pathlib import Path

# Importation de la fonction génératrice
from becke_grid import generate_grid

# =========================================================
# 1. Fonctions de partitionnement de Becke
# =========================================================

def p(mu):
    return (3/2)*mu - (1/2)*(mu**3)

def f_k(mu, k=3):
    resultat = mu 
    for _ in range(k):
        resultat = p(resultat)
    return resultat

def s_k(mu, k=3):
    return (1/2) * (1 - f_k(mu, k))

# =========================================================
# 2. Préparation de la grille et des poids globaux
# =========================================================

def preparer_grille_et_poids(coords_liste, P=100, ordre_choisi=53, r_m=0.5):
    """
    Génère les points de la grille et calcule le poids total pour chaque point
    (intégrant Lebedev, Chebyshev, Jacobien et les poids de Becke).
    """
    N_atomes = len(coords_liste)
    coords_array = np.array(coords_liste)
    
    # Génération de la grille de base
    cart_array = generate_grid(P, coords_liste, ordre_choisi, r_m=r_m, units='angstrom')
    if cart_array is None:
        raise ValueError("Erreur : La grille n'a pas pu être générée.")
        
    N_points_par_grille = cart_array.shape[1]
    
    # --- Poids Radiaux ---
    indices = np.arange(1, P + 1)
    x_i = np.cos(np.pi * indices / (P + 1))
    r_i = r_m * (1 + x_i) / (1 - x_i)
    w_x_cheby = (np.pi / (P + 1)) * np.sin(np.pi * indices / (P + 1))
    jacobien_mapping = (2 * r_m) / (1 - x_i)**2
    w_rad = w_x_cheby * jacobien_mapping * (r_i**2)
    
    # --- Poids Angulaires (Lebedev) ---
    directory_script = Path(__file__).parent
    file_name = f"lebedev_{ordre_choisi}.txt"
    lebedev_data = np.loadtxt(directory_script / file_name)
    w_leb = lebedev_data[:, 2] * 4 * np.pi
    
    w_angulaire_total = np.tile(w_leb, P)
    w_radial_total = np.repeat(w_rad, len(w_leb))
    w_base = w_radial_total * w_angulaire_total # Poids combiné (radial * angulaire)
    
    # --- Poids de Becke (Nécessite coords_liste) ---
    R_ij = np.zeros((N_atomes, N_atomes))
    for i in range(N_atomes):
        for j in range(N_atomes):
            R_ij[i, j] = np.linalg.norm(coords_array[i] - coords_array[j])
            
    points_aplatis = cart_array.reshape(-1, 3)
    M = points_aplatis.shape[0] # Nombre total de points
    r_dist = np.linalg.norm(points_aplatis[:, np.newaxis, :] - coords_array[np.newaxis, :, :], axis=2)
    
    P_vals = np.ones((M, N_atomes))
    for i in range(N_atomes):
        for j in range(N_atomes):
            if i != j:
                mu_ij = (r_dist[:, i] - r_dist[:, j]) / R_ij[i, j]
                P_vals[:, i] *= s_k(mu_ij, k=3)
                
    w_tous_points = P_vals / np.sum(P_vals, axis=1, keepdims=True) 
    w_brut_reshape = w_tous_points.reshape(N_atomes, N_points_par_grille, N_atomes)
    
    # --- Assemblage des vecteurs 1D finaux (x, y, z, poids) ---
    x_tot, y_tot, z_tot, poids_tot = [], [], [], []
    for i in range(N_atomes):
        w_becke_i = w_brut_reshape[i, :, i]
        poids_finaux_i = w_base * w_becke_i
        
        x_tot.extend(cart_array[i, :, 0])
        y_tot.extend(cart_array[i, :, 1])
        z_tot.extend(cart_array[i, :, 2])
        poids_tot.extend(poids_finaux_i)
        
    return np.array(x_tot), np.array(y_tot), np.array(z_tot), np.array(poids_tot)

# =========================================================
# 3. Fonction d'intégration finale (Nettoyée)
# =========================================================

def integrer_fonction(fonction_F, x, y, z, poids):
    """
    Calcule l'intégrale finale.
    Ne dépend plus d'aucune coordonnée atomique ni des paramètres de grille.
    """
    valeurs_F = fonction_F(x, y, z)
    return np.sum(valeurs_F * poids)

# =========================================================
# 4. Exécution et Test
# =========================================================
if __name__ == "__main__":
    
    # Fonction test bidon (pour l'exemple)
    def f(x, y, z):
        return x**2 + y**2 + z**2

    coords_liste = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0])]

    print("1. Préparation de la grille et calcul des poids...")
    x, y, z, poids = preparer_grille_et_poids(coords_liste, P=100, ordre_choisi=53, r_m=0.5)
    
    print("2. Lancement de l'intégration...")
    # Ta fonction finale est maintenant ultra légère !
    resultat = integrer_fonction(f, x, y, z, poids)
    
    print(f"Résultat de l'intégrale : {resultat}")