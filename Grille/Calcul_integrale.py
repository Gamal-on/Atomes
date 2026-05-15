import numpy as np 
import matplotlib.pyplot as plt 

#définition des fonctions élémentaires
def p(mu):
    return ((3/2)*mu - (1/2)*(mu**3))
 
def f_k(mu,k):
    resultat = mu 
    for i in range(k):
        resultat=p(resultat)
    return resultat

def s_k(mu,k=3):
    return (1/2)*(1-f_k(mu,k))

#création des valeurs de mu
mu_values = np.linspace(-1,1,500)

plt.figure(figsize=(8, 6))

for k in range(1,6):
    s_values = s_k(mu_values,k)
    plt.plot(mu_values,s_values,label=f'k={k}')
#on affiche la figure     
plt.title('Profils de coupure de Becke $s_k(\mu)$', fontsize=14)
plt.xlabel('$\mu$ (Coordonnée elliptique)', fontsize=12)
plt.ylabel('$s_k(\mu)$ (Poids spatial)', fontsize=12)

plt.xlim(-1, 1)
plt.ylim(0.0, 1.0) 

plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.show()

N = 2 

coords_liste = np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0])
corrds_array = np.array(coords_liste)
coords_rayon = np.linalg.norm(coords_liste)

# matrice des distances
R_ij = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        R_ij[i][j] = np.linalg.norm(coords_liste[i] - coords_liste[j])

#là ou on calcule l'intégrale
#calcul des polynomes p_i


#méthode de gauss-chebychev : on crée les points 
P = 5000
indices = np.arange(1, P + 1)
x_i = np.cos(np.pi * indices / (P + 1))
r_m = 0.5 #valeur an Angstrom
r_i_array = r_m * (1 + x_i) / (1 - x_i)
w_x_cheby = (np.pi / (P + 1)) * np.sin(np.pi * indices / (P + 1))
w_cheby_r = w_x_cheby*(2 * r_m)/(1 - x_i)**2 #calculé avec le jacobien 





#chargement des points de lebdev
ordre_choisi = 53
chemin_fichier = rf"C:\Users\yannc\Documents\Codesprojet\lebedev_{ordre_choisi}.txt"
donnees_brutes = np.loadtxt(chemin_fichier)
lebedev_data = donnees_brutes.tolist() #je convertis lebedev en une liste de liste
lebedev_data_array = np.array(lebedev_data) #c'est une matrice à 3 collones (theta, phi, w) et Nlebedev lignes
theta = lebedev_data_array[:,0]
phi = lebedev_data_array[:,1]

#on crée les points en cartésien
N=len(coords_liste)

#on calcule alors des tableaux contenant coordonnées cartésiennes sur la sphère unité 
X_unit = np.sin(theta) * np.cos(phi)
Y_unit = np.sin(theta) * np.sin(phi)
Z_unit = np.cos(theta)

#on transfome r_i_array d'un vecteur en une matrice (n,1) et X unit en une matrice (1,nlebedev) et on effectue le produit matriciel correspondant, écrasé par ravel en un veteur 1D
X_rel = (r_i_array[:, np.newaxis] * X_unit[np.newaxis, :]).ravel()
Y_rel = (r_i_array[:, np.newaxis] * Y_unit[np.newaxis, :]).ravel()
Z_rel = (r_i_array[:, np.newaxis] * Z_unit[np.newaxis, :]).ravel()

#on crée une grille de 3 colonnes et n points
grille_relative = np.column_stack((X_rel, Y_rel, Z_rel))
# grille_relative a la forme (Points, 3) -> on lui donne la forme (1, Points, 3)
# coords_matrice a la forme (N, 3)      -> on lui donne la forme (N, 1, 3)
# Translation vectorisée de la grille sur chaque atome via broadcasting -> forme (N, points, 3)
cart_array = grille_relative[np.newaxis, :, :] + corrds_array[:, np.newaxis, :]

#poids de lebedev
w_lebedev = lebedev_data_array[:,2]


#Poids de Becke
N_points_par_grille = cart_array.shape[1] #nombre total de points
points_aplatis = cart_array.reshape(-1, 3) #on a un tableau de taille (ligne = atomes x points, 3 colonnes x y  z)
M = points_aplatis.shape[0] #on return atomes x points 

# Distance de chaque point 3D à chaque noyau (M points, N atomes)
r_dist = np.linalg.norm(points_aplatis[:, np.newaxis, :] - corrds_array[np.newaxis, :, :], axis=2)
P_vals = np.ones((M, N))
for i in range(N):
    for j in range(N):
        if i != j:
            mu_ij = (r_dist[:, i] - r_dist[:, j]) / R_ij[i, j]
            P_vals[:, i] *= s_k(mu_ij, k=3)
            
            
# Normalisation
w_tous_points = P_vals / np.sum(P_vals, axis=1, keepdims=True) 
w_brut_reshape = w_tous_points.reshape(N, N_points_par_grille, N) #origine du point, index du point, valeur des poids

# On stocke les poids dans une matrice (N_atomes, N_points)
w_i_becke = np.zeros((N, N_points_par_grille))
for i in range(N):
    w_i_becke[i, :] = w_brut_reshape[i, :, i]

def calcul_integrale(fonction_F, cart_array, w_leb, w_rad, w_becke):
    S = 0.0
    w_angulaire_total = np.tile(w_leb, len(w_rad))
    w_radial_total = np.repeat(w_rad, len(w_leb))
    
    for i in range(N):
        x = cart_array[i, :, 0]
        y = cart_array[i, :, 1]
        z = cart_array[i, :, 2]
        
        valeurs_F = fonction_F(x, y, z)
        produit = valeurs_F * w_radial_total * w_angulaire_total * w_becke[i, :]
        S += np.sum(produit)
        
    return S

def f(x,y,z):
    return np.exp(-(x**2+y**2+z**2))

resultat = calcul_integrale(f, cart_array, w_lebedev, w_cheby_r, w_i_becke)
print(f"L'intégrale de la gaussienne donne : {resultat:.6f}")