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

    phi = lebedev_data_array[:, 0]
    theta = lebedev_data_array[:, 1]

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


def export_grid_for_multiwfn(cart_array, output_path, centers_bohr, r_cutoff=50.0):
    """
    Export the grid points to a file in a format compatible with Multiwfn, in bohr
    Arguments:
    - cart_array : array of shape (N_atoms, N_points, 3) in Angstrom
    - output_path : path to the output file
    - centers_bohr : list of coordinates of the centers in Bohr
    - r_cutoff : distance max en Bohr à partir de l'atome le plus proche
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

def plot_grid_tenser(cart_array):
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



def plot_multicenter_grid(cart_array, centers_bohr):
    ANGSTROM_TO_BOHR = 1.8897259886
    n_atoms = cart_array.shape[0]
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # On choisit une palette de couleurs (tab10, jet, viridis, etc.)
    colors = plt.cm.get_cmap('tab10', n_atoms)
    
    # Conversion globale pour le calcul des limites (plus simple)
    grid_bohr = cart_array * ANGSTROM_TO_BOHR
    centers = np.array(centers_bohr)

    # Boucle sur chaque atome
    for i in range(n_atoms):
        atom_points = grid_bohr[i]  # Forme (N_points, 3)
        atom_color = colors(i)
        
        # Affichage de la grille de l'atome i
        # Note : stride [::10] recommandé si c'est trop lent
        ax.scatter(atom_points[:, 0], atom_points[:, 1], atom_points[:, 2], 
                   s=1, alpha=0.2, color=atom_color, label=f'Grille Atome {i}')
        
        # Affichage du centre correspondant
        ax.scatter(centers[i, 0], centers[i, 1], centers[i, 2], 
                   s=150, color=atom_color, marker='o', edgecolors='black', linewidth=2)

    # --- Gestion de l'échelle et des limites ---
    # On calcule les limites sur l'ensemble de la grille
    mins = grid_bohr.min(axis=(0, 1))
    maxs = grid_bohr.max(axis=(0, 1))
    mid_points = (maxs + mins) * 0.5
    max_range = (maxs - mins).max() / 2.0

    ax.set_xlim(mid_points[0] - max_range, mid_points[0] + max_range)
    ax.set_ylim(mid_points[1] - max_range, mid_points[1] + max_range)
    ax.set_zlim(mid_points[2] - max_range, mid_points[2] + max_range)

    ax.set_xlabel('X (Bohr)')
    ax.set_ylabel('Y (Bohr)')
    ax.set_zlabel('Z (Bohr)')
    ax.set_title(f'Grille multicentrique colorée par atome ({n_atoms} centres)')
    
    # On n'affiche la légende que si on n'a pas trop d'atomes
    if n_atoms <= 10:
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
    grid = generate_grid(10, centers, ordre_choisi=53, r_m=0.5, units='bohr')
    export_grid_for_multiwfn(grid, "urea_grid_multiwfn.txt", centers)
    plot_multicenter_grid(grid, centers)

if __name__ == "__main__":
    main()