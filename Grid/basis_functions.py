import time
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# IMPORTATION : On fait uniquement confiance à ton calcul.py
# =========================================================
from calcul import preparer_grille_et_poids, integrer_fonction

# =========================================================
# 1. Fonction indépendante pour la complexité
# =========================================================

def tester_complexite(tailles_systeme=range(2, 8), P_val=40, ordre_val=23, afficher_graphique=True):
    """
    Mesure la complexité algorithmique de l'intégration de Becke.
    """
    temps_preparation = []
    temps_integration = []

    def f_test(x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)
        return (x**2 + y**2) * np.exp(-r)

    print(f"Début du test de complexité (Atomes de {min(tailles_systeme)} à {max(tailles_systeme)})...")
    
    for n_atomes in tailles_systeme:
        # Configuration linéaire stable pour éviter les divisions par zéro
        coords_liste = [np.array([float(i)*2.0, 0.0, 0.0]) for i in range(n_atomes)]
        
        # 1. Temps de calcul des grilles + poids de Becke
        start_prep = time.time()
        x, y, z, poids = preparer_grille_et_poids(coords_liste, P=P_val, ordre_choisi=ordre_val, r_m=0.5)
        end_prep = time.time()
        temps_preparation.append(end_prep - start_prep)
        
        # 2. Temps de l'intégration vectorisée pure
        start_int = time.time()
        integrer_fonction(f_test, x, y, z, poids)
        end_int = time.time()
        temps_integration.append(end_int - start_int)
        
        print(f"Atomes: {n_atomes} | Points totaux: {len(x):,} | Temps préparation: {end_prep - start_prep:.3f} s | Temps calcul: {end_int - start_int:.4f} s")

    if afficher_graphique:
        plt.figure(figsize=(10, 6))
        plt.plot(list(tailles_systeme), temps_preparation, marker='o', color='b', label="Préparation (Grille + Becke)")
        plt.plot(list(tailles_systeme), temps_integration, marker='s', color='g', label="Intégration vectorisée pure")
        
        coeff = temps_preparation[-1] / (list(tailles_systeme)[-1]**3)
        ref_cubique = [coeff * (n**3) for n in tailles_systeme]
        
        plt.plot(list(tailles_systeme), ref_cubique, linestyle=':', color='r', label=r"Référence $\mathcal{O}(N^3)$")
        
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

def tester_orbitales(distance_AB=1.0, P_val=70, ordre_val=53):
    """
    Test l'intégration multicentrique sur des orbitales de type Slater.
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

    # Densités associées
    def densite_1s_A(x, y, z): return psi_1s_A(x, y, z)**2
    def densite_2px_A(x, y, z): return psi_2px_A(x, y, z)**2
    def recouvrement_AB(x, y, z): return psi_1s_A(x, y, z) * psi_1s_B(x, y, z)

    print(f"Génération de la grille via calcul.py (P={P_val}, Ordre={ordre_val})...")
    x, y, z, poids = preparer_grille_et_poids(coords_liste, P=P_val, ordre_choisi=ordre_val, r_m=0.5)
    
    int_1s = integrer_fonction(densite_1s_A, x, y, z, poids)
    int_2px = integrer_fonction(densite_2px_A, x, y, z, poids)
    S_ab = integrer_fonction(recouvrement_AB, x, y, z, poids)
    
    S_ab_exact = np.exp(-distance_AB) * (1 + distance_AB + (distance_AB**2)/3.0)
    
    resultats = {
        "norme_1s": {"calcule": int_1s, "exact": 1.0},
        "integrale_2px": {"calcule": int_2px, "exact": np.pi},
        "recouvrement_S_ab": {"calcule": S_ab, "exact": S_ab_exact}
    }
    
    return resultats


# =========================================================
# 3. Lancement automatique des tests
# =========================================================
if __name__ == "__main__":
    
    # --- Test 1 : Orbitales ---
    print("\n" + "="*50)
    print("TEST NUMÉRIQUE : ORBITALES DE SLATER")
    print("="*50)
    
    resultats_orb = tester_orbitales(distance_AB=1.0, P_val=70, ordre_val=53)
    
    for nom, valeurs in resultats_orb.items():
        calc = valeurs["calcule"]
        exact = valeurs["exact"]
        erreur = abs(calc - exact)
        print(f"--- {nom} ---")
        print(f"Calculé : {calc:.8f}")
        print(f"Exact   : {exact:.8f}")
        print(f"Erreur  : {erreur:.2e}\n")

    # --- Test 2 : Complexité ---
    print("\n" + "="*50)
    print("TEST DE COMPLEXITÉ : O(N^3)")
    print("="*50)
    
    tester_complexite(tailles_systeme=range(2, 8), P_val=40, ordre_val=53, afficher_graphique=True)