import numpy as np
from pathlib import Path

# Importation de la fonction génératrice
from becke_grid import generate_grid


#préparation grille et poids 

def preparer_grille_et_poids(coords_liste, P=100, ordre_choisi=53, r_m=0.5, w_hirsh):
    """
    Génère les points de la grille et calcule le poids total pour chaque point
    (intégrant Lebedev, Chebyshev, Jacobien et les poids de Hirshfeld).
    """
    N_atomes = len(coords_liste)
    coords_array = np.array(coords_liste)
    
    # Génération de la grille de base avec le code grid.py
    cart_array = generate_grid(P, coords_liste, ordre_choisi, r_m=r_m, units='angstrom')
    if cart_array is None:
        raise ValueError("Erreur : La grille n'a pas pu être générée.")
        
    N_points_par_grille = cart_array.shape[1]
    
    # Poids et indices radiaux
    indices = np.arange(1, P + 1)
    x_i = np.cos(np.pi * indices / (P + 1))
    r_i = r_m * (1 + x_i) / (1 - x_i)
    w_x_cheby = (np.pi / (P + 1)) * np.sin(np.pi * indices / (P + 1))
    jacobien_mapping = (2 * r_m) / (1 - x_i)**2
    w_rad = w_x_cheby * jacobien_mapping * (r_i**2)
    
    # Poids angulaires
    directory_script = Path(__file__).parent
    file_name = f"lebedev_{ordre_choisi}.txt"
    lebedev_data = np.loadtxt(directory_script / file_name)
    w_leb = lebedev_data[:, 2] * 4 * np.pi
    
    #on cherche à obtenir un tableau [R1*A1,R1*A2,R1*A3,R2*A1...]
    w_angulaire_total = np.tile(w_leb, P)
    w_radial_total = np.repeat(w_rad, len(w_leb))
    w_base = w_radial_total * w_angulaire_total # Poids combiné (radial * angulaire)
    
    # Poids de Becke 
    R_ij = np.zeros((N_atomes, N_atomes))
    for i in range(N_atomes):
        for j in range(N_atomes):
            R_ij[i, j] = np.linalg.norm(coords_array[i] - coords_array[j])
    #on applati les deux premières colonnes        
    points_aplatis = cart_array.reshape(-1, 3)
    M = points_aplatis.shape[0] # Nombre total de points
    #matrice des distances points i_ atomhe j
    r_dist = np.linalg.norm(points_aplatis[:, np.newaxis, :] - coords_array[np.newaxis, :, :], axis=2)
    
    #calcul des corrdonnées cartésiennes et calcul des poids finaux
    x_tot, y_tot, z_tot, poids_tot = [], [], [], []
    for i in range(N_atomes):
        poids_finaux_i = w_base * w_hirsh
        
        x_tot.extend(cart_array[i, :, 0])
        y_tot.extend(cart_array[i, :, 1])
        z_tot.extend(cart_array[i, :, 2])
        poids_tot.extend(poids_finaux_i)
        
    return np.array(x_tot), np.array(y_tot), np.array(z_tot), np.array(poids_tot)

#fonction d'intégration

def integrer_fonction(fonction_F, x, y, z, poids):
    """
    Calcule l'intégrale finale.
    Ne dépend plus d'aucune coordonnée atomique ni des paramètres de grille.
    """
    valeurs_F = fonction_F(x, y, z)
    return np.sum(valeurs_F * poids)

