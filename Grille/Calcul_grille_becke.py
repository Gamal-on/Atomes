#représenter l'espace en 3D
#convertir les points en cartésien
#coords est une liste qui représente les coordonnées du centre de chaque atom (il y en a N)
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path


def generate_grid(radial_lenght, centers_coordinates, ordre_choisi, r_m=0.5, units='bohr'):
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
    i_range = np.arange(1, radial_lenght+1)
    x_i = np.cos(np.pi * i_range / (radial_lenght + 1))
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


def export_grid_for_multiwfn(cart_array, output_path, centers_bohr, r_cutoff=20.0):
    """
    r_cutoff : distance max en Bohr à partir de l'atome le plus proche
    """
    ANGSTROM_TO_BOHR = 1.8897259886

    all_points = cart_array.reshape(-1, 3) * ANGSTROM_TO_BOHR
    centers = np.array(centers_bohr)  # (N_atoms, 3) en Bohr

    # Distance minimale de chaque point à n'importe quel atome
    dists = np.linalg.norm(
        all_points[:, np.newaxis, :] - centers[np.newaxis, :, :],
        axis=2
    )  # shape (N_points, N_atoms)
    min_dist = dists.min(axis=1)

    mask = min_dist <= r_cutoff
    filtered = all_points[mask]

    print(f"Points avant cutoff : {len(all_points)}, après : {len(filtered)}")

    with open(output_path, 'w') as f:
        f.write(f"{len(filtered)}\n")
        for x, y, z in filtered:
            f.write(f"{x:20.10f}  {y:20.10f}  {z:20.10f}\n")

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
    centers = np.array([[0.00000000,  0.00000000,  0.27127625],  
                        [0.00000000, 0.00000000,  2.57183145], 
                        [0.00000000,  2.19124748, -1.13016901], 
                        [0.00000000, 3.82399079e+00, -1.62836533e-01], 
                        [0.00000000,  2.24084147e+00, -3.02713496e+00],
                        [0.00000000, -2.19124748e+00, -1.13016901e+00],
                        [0.00000000, -2.24084147e+00,-3.02713496e+00], 
                        [0.00000000, -3.82399079e+00, -1.62836533e-01]])
    grid = generate_grid(10, centers, ordre_choisi=29, r_m=0.5, units='bohr')
    export_grid_for_multiwfn(grid, "urea_grid_multiwfn.txt", centers)

if __name__ == "__main__":
    main()