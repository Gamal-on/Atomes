import numpy as np
import argparse
from data_read import get_all_mo

def molecular_orbitals(position, mo_coefficients, gaussian_exponents):
    """
    Compute the molecule Molecular Orbitals 
    from a basis of Gaussian functions, for a given position.
    
    Parameters:
    - position: 1D array of coordinates (Bohr radius)
    - mo_coefficients: 2D array (M, 36)
    - gaussian_exponents: 1D array (36,)
    
    Returns:
    - 2D array of shape (M, len(position))
    """
    position = np.asarray(position, dtype=float)

    # Gaussiennes : exp(-alpha * x^2)
    gaussians = np.exp(-gaussian_exponents[:, np.newaxis] * position[np.newaxis, :]**2)

    molecular_orbitals = mo_coefficients @ gaussians
    return molecular_orbitals

def molecular_orbitals_3d(x_range, y_range, z_range,
                          mo_coefficients, gaussian_exponents, centers):
    """
    - x, y, z : 
    - mo_coefficients
    - gaussian_exp
    - centers : array (36, 3) contenant les positions (x,y,z) des gaussiennes
    """

    x, y, z = np.meshgrid(x_range, y_range, z_range, indexing='ij')

    # Broadcast pour chaque centre
    dx = x[np.newaxis, :, :, :] - centers[:, 0][:, np.newaxis, np.newaxis, np.newaxis]
    dy = y[np.newaxis, :, :, :] - centers[:, 1][:, np.newaxis, np.newaxis, np.newaxis]
    dz = z[np.newaxis, :, :, :] - centers[:, 2][:, np.newaxis, np.newaxis, np.newaxis]

    r2 = dx**2 + dy**2 + dz**2

    gaussians = np.exp(-gaussian_exponents[:, np.newaxis, np.newaxis, np.newaxis] * r2)

    molecular_orbitals = np.tensordot(mo_coefficients, gaussians, axes=(1, 0))

    return molecular_orbitals


def double_factorial(n):
    """Calcul the double factorail (n!!) """
    if n <= 0:
        return 1
    return np.prod(np.arange(n, 0, -2))


def calculate_molecular_orbitals(centers, exponents, powers, coefficients, x, y, z):
    """Calcule la valeur des orbitales moléculaires sur une grille 3D.

    Paramètres:
    -----------
    - centers : np.array (36, 3), coordinates of the atom for each gaussian
    - exponents : np.array de forme (36)
    - powers : np.array de forme (36, 3), the numbers [n,l,m] for each gaussian
    - coefficients : np.array de forme (m, 36) c_ij coefficients for each gaussian
    - x, y, z : np.array de forme (N) ranges to iitialise the grid
    Returns:
    ---------
    np.array de forme (m, Nx, Ny, Nz)
    """

    # Initilising the grid
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    grid_shape = X.shape
    n_gaussians = 36
    n_orbitals = coefficients.shape[0]

    # Initialisation of the output array
    mo_grids = np.zeros((n_orbitals, *grid_shape))

    # Calcul of the normalisation factor for each gaussian
    # N = np.ones(n_gaussians)
    N = np.zeros(n_gaussians)
    for i in range(n_gaussians):
        l, m, n = powers[i]
        alpha = exponents[i]

        num = (2 * alpha / np.pi) ** 0.75 * (4 * alpha) ** (
            (l + m + n) / 2.0)
        den = np.sqrt(
            double_factorial(2 * l - 1)
            * double_factorial(2 * m - 1)
            * double_factorial(2 * n - 1))
        N[i] = num /den

    # Contribution of each gaussian
    for i in range(n_gaussians):
        xc, yc, zc = centers[i]
        alpha = exponents[i]
        l, m, n = powers[i]

        # Distance to the center of the atom
        dx = X - xc
        dy = Y - yc
        dz = Z - zc
        r2 = dx**2 + dy**2 + dz**2

        # Value of the primitivve gaussian
        gaussian_val = N[i] * (dx**l) * (dy**m) * (dz**n) * np.exp(-alpha * r2)

        for j in range(n_orbitals):
            mo_grids[j] += coefficients[j, i] * gaussian_val

    return mo_grids

def main() :
    pass


if __name__ == "__main__"  :
    main()