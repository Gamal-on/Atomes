#représenter l'espace en 3D
#convertir les points en cartésien
#coords est une liste qui représente les coordonnées du centre de chaque atom (il y en a N)
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path


def generate_grid(number_radial_points, centers_coordinates, ordre_choisi, r_m=0.5, units='angstroem'):
    """
    Returns a multicenters 3D grid, in Angstrom, of shape (N_atoms, N_points, 3).
    Arguments : 
    - number_radial_points : number of radial points
    - centers_coordinates : list of coordinates of the centers (N x 3)
    - ordre_choisi : order of the Lebedev grid
    """

    BOHR_TO_ANGSTROM = 0.5291772108

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
    i_range = np.arange(1, number_radial_points+1)
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

    if units == 'bohr':
        coords_matrice = coords_matrice * BOHR_TO_ANGSTROM
    
    # Final resulat
    cart_array = grille_relative[np.newaxis, :, :] + coords_matrice[:, np.newaxis, :]
    
    return cart_array


def export_grid_for_multiwfn(cart_array, output_path):
    """
    Converts the grid tensor (N_atoms, N_points, 3) to Multiwfn format.
    Converts Angstrom -> Bohr.
    
    Arguments:
    - cart_array   : output of generate_grid(), shape (N_atoms, N_points, 3)
    - output_path  : path of the output text file
    """
    ANGSTROM_TO_BOHR = 1.8897259886

    # Flatten + conversion
    all_points = cart_array.reshape(-1, 3) * ANGSTROM_TO_BOHR
    
    n_total = len(all_points)
    
    with open(output_path, 'w') as f:
        f.write(f"{n_total}\n")
        for x, y, z in all_points:
            f.write(f"{x:20.10f}  {y:20.10f}  {z:20.10f}\n")
    
    print(f"{n_total} points written to {output_path}")

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
    grid = generate_grid(500, [[0, 0, 0], [500, 500, 0]], 15)
    print(grid)

if __name__ == "__main__":
    main()



    
