import time
import numpy as np
import matplotlib.pyplot as plt

# Import des fonctions depuis ton fichier calcul.py
from calcul import preparer_grille_et_poids, integrer_fonction

# =========================================================
# 1. Fonction indépendante pour la complexité
# =========================================================

def tester_complexite(tailles_systeme=range(2, 16), P_val=40, ordre_val=23, afficher_graphique=True):
    """
    Mesure la complexité algorithmique de l'intégration de Becke.
    Retourne les temps de préparation et d'intégration.
    """
    temps_preparation = []
    temps_integration = []

    def f_test(x, y, z):
        return x**2 + y**2

    print(f"Début du test de complexité (Atomes de {min(tailles_systeme)} à {max(tailles_systeme)})...")
    
    for n_atomes in tailles_systeme:
        coords_liste = [np.random.uniform(-2.0, 2.0, 3) for _ in range(n_atomes)]
        
        start_prep = time.time()
        x, y, z, poids = preparer_grille_et_poids(coords_liste, P=P_val, ordre_choisi=ordre_val, r_m=0.5)
        end_prep = time.time()
        temps_preparation.append(end_prep - start_prep)
        
        start_int = time.time()
        integrer_fonction(f_test, x, y, z, poids)
        end_int = time.time()
        temps_integration.append(end_int - start_int)
        
        print(f"Atomes: {n_atomes} | Points: {len(x)} | Temps de préparation: {end_prep - start_prep:.3f} s")

    if afficher_graphique:
        plt.figure(figsize=(10, 6))
        plt.plot(list(tailles_systeme), temps_preparation, marker='o', color='b', label="Préparation (Becke)")
        plt.plot(list(tailles_systeme), temps_integration, marker='s', color='g', label="Intégration vectorisée")
        
        coeff = temps_preparation[-1] / (list(tailles_systeme)[-1]**3)
        ref_cubique = [coeff * (n**3) for n in tailles_systeme]
        plt.plot(list(tailles_systeme), ref_cubique, linestyle=':', color='r', label="Référence $\mathcal{O}(N^3)$")
        
        plt.title("Complexité Algorithmique", fontsize=14)
        plt.xlabel("Nombre d'atomes (N)", fontsize=12)
        plt.ylabel("Temps d'exécution (secondes)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()

    return temps_preparation, temps_integration

# =========================================================
# 2. Fonction indépendante pour les orbitales
# =========================================================

def tester_orbitales(distance_AB=1.0, P_val=70, ordre_val=59):
    """
    Test l'intégration multicentrique sur des orbitales de type Slater.
    Prend en paramètre la distance entre les atomes A et B.
    Retourne un dictionnaire contenant les résultats des intégrales.
    """
    coord_A = np.array([0.0, 0.0, 0.0])
    coord_B = np.array([distance_AB, 0.0, 0.0])
    coords_liste = [coord_A, coord_B]

    # Fonctions d'orbitales (Slater)
    def psi_1s_A(x, y, z):
        r_A = np.sqrt((x - coord_A[0])**2 + (y - coord_A[1])**2 + (z - coord_A[2])**2)
        return (1.0 / np.sqrt(np.pi)) * np.exp(-r_A)

    def psi_1s_B(x, y, z):
        r_B = np.sqrt((x - coord_B[0])**2 + (y - coord_B[1])**2 + (z - coord_B[2])**2)
        return (1.0 / np.sqrt(np.pi)) * np.exp(-r_B)

    def psi_2px_A(x, y, z):
        r_A = np.sqrt((x - coord_A[0])**2 + (y - coord_A[1])**2 + (z - coord_A[2])**2)
        return (x - coord_A[0]) * np.exp(-r_A)

    # Fonctions à intégrer
    def densite_1s_A(x, y, z): return psi_1s_A(x, y, z)**2
    def densite_2px_A(x, y, z): return psi_2px_A(x, y, z)**2
    def recouvrement_AB(x, y, z): return psi_1s_A(x, y, z) * psi_1s_B(x, y, z)

    print(f"Génération de la grille pour test orbitales (P={P_val}, Ordre={ordre_val})...")
    x, y, z, poids = preparer_grille_et_poids(coords_liste, P=P_val, ordre_choisi=ordre_val, r_m=0.5)
    
    # Calculs
    int_1s = integrer_fonction(densite_1s_A, x, y, z, poids)
    int_2px = integrer_fonction(densite_2px_A, x, y, z, poids)
    S_ab = integrer_fonction(recouvrement_AB, x, y, z, poids)
    
    # Valeurs exactes pour comparaison
    S_ab_exact = np.exp(-distance_AB) * (1 + distance_AB + (distance_AB**2)/3.0)
    
    resultats = {
        "norme_1s": {"calcule": int_1s, "exact": 1.0},
        "integrale_2px": {"calcule": int_2px, "exact": np.pi/3},
        "recouvrement_S_ab": {"calcule": S_ab, "exact": S_ab_exact}
    }
    
    return resultats