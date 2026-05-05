#représenter l'espace en 3D
#convertir les points en cartésien
#coords est une liste qui représente les coordonnées du centre de chaque atom (il y en a N)
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path


def generate_grid(number_radial_points, centers_coordinates, ordre_choisi, r_m=0.5):
    """
    Returns a multicenters 3D grid
    Arguments : 
    - number_radial_points : number of radial points
    - centers_coordinates : list of coordinates of the centers (N x 3)
    - ordre_choisi : order of the Lebedev gri
    """
    # Loading lebedev data
    directory_script = Path(__file__).parent
    file_name = f"lebedev_{ordre_choisi}.txt"
    path_complete = directory_script / file_name

    try:
        lebedev_data_array = np.loadtxt(path_complete)
    except FileNotFoundError:
        print(f"Erreur : Le fichier {file_name} est introuvable.")
        return None

    theta = lebedev_data_array[:, 0]
    phi = lebedev_data_array[:, 1]

    # Generating radial points (Gauss-Chebyshev type)
    i_range = np.arange(number_radial_points)
    x_i = np.cos(np.pi * i_range / (number_radial_points + 1))
    r_i = r_m * (1 + x_i) / (1 - x_i)
    
    # Cartesian coordinates on the unit sphere
    X_unit = np.sin(theta) * np.cos(phi)
    Y_unit = np.sin(theta) * np.sin(phi)
    Z_unit = np.cos(theta)

    # Combinaison radial x angular
    X_rel = (r_i[:, np.newaxis] * X_unit[np.newaxis, :]).ravel()
    Y_rel = (r_i[:, np.newaxis] * Y_unit[np.newaxis, :]).ravel()
    Z_rel = (r_i[:, np.newaxis] * Z_unit[np.newaxis, :]).ravel()
    
    grille_relative = np.column_stack((X_rel, Y_rel, Z_rel))

    # Translation of the grid on each center
    coords_matrice = np.array(centers_coordinates)
    
    # Final resulat
    cart_array = grille_relative[np.newaxis, :, :] + coords_matrice[:, np.newaxis, :]
    
    return cart_array

def plot_grid(cart_array):
    """
    Affiche la grille 3D générée.
    """
    if cart_array is None:
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    n_atomes = cart_array.shape[0]
    couleurs = ['b', 'r', 'g', 'c', 'm', 'y']
        
    for i in range(n_atomes):
        ax.scatter(
            cart_array[i, :, 0], 
            cart_array[i, :, 1], 
            cart_array[i, :, 2], 
            c=couleurs[i % len(couleurs)], 
            marker='.', s=2, alpha=0.3, 
            label=f'Atome {i}'
        )
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title("Grille de points multicentrique")
    ax.legend() 
    plt.show()


def main() : 
    generate_grid(500, [[0, 0, 0], [1, 1, 1]], 15)

if __name__ == "__main__":
    main()



    
