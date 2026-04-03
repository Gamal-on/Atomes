import numpy as np

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



def main() :
    pass


if __name__ == "__main__"  :
    main()