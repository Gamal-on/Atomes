import numpy as np


def molecular_orbitals(position, mo_coefficients, gaussian_exponents):
    """
    Compute the molecule’s Molecular Orbitals, 
    from a basis of Gaussian functions, for a given position.
    
    Parameters:
    - x, y, z: nd arrays of atomic distances (unit: Bohr radius)
    - mo_coefficients: 2D array (size M, 36), coefficients of the basis functions of the orbitals
    - gaussian_exponents: 1D array (size 36), exponents of the Gaussian basis functions
    
    Returns:
    - 2D array containing the values of the MOs for each position value.
    """
    # Test coherent 
    # if mo_coefficients.shape[1] != len(gaussian_exponents):
      #  raise ValueError(
        #    f"Incohérence : {mo_coefficients.shape[1]} coeffs mais {len(gaussian_exponents)} exposants !"
       # )
    
    # On calcule toutes les gaussiennes pour toutes les positions d'un coup
    gaussians = np.exp(-gaussian_exponents[:, np.newaxis] * (position)[np.newaxis, :])
    
    # On applique les coefficients par multiplication matricielle
    molecular_orbitals = mo_coefficients @ gaussians
    
    return molecular_orbitals

def molecular_orbitals_3d(x_range, y_range, z_range, mo_coefficients, gaussian_exponents):
    """
   Return the molecular orbitals on a 3D grid.
    
    Parameters:
    - x_range, y_range, z_range : 1D array of the coordinates, in Bohr radius
    - mo_coefficients: 2D array (size M, 36), coefficients of the basis functions of the orbitals
    - gaussian_exponents: 1D array (size 36), exponents of the Gaussian basis functions
    
    Returns:
    - 2D array containing the values of the MOs for each position value.
    """
    # réation de la grille 3D
    x, y, z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
    r2 = x**2 + y**2 + z**2  # Forme (Nx, Ny, Nz)
    
    # Calcul des 36 gaussiennes sur toute la grille d'un coup
    gaussians = np.exp(-gaussian_exponents[:, np.newaxis, np.newaxis, np.newaxis] * r2[np.newaxis, :, :, :])
    
    # On somme
    molecular_orbitals = np.tensordot(mo_coefficients, gaussians, axes=(1, 0))
    
    return molecular_orbitals



def main() :
    pass


if __name__ == "__main__"  :
    main()