import numpy as np


def molecular_orbitals(position, mo_coefficients, gaussian_exponents):
    """
    Compute the molecule’s Molecular Orbitals, 
    from a basis of Gaussian functions, for a given position.
    
    Parameters:
    - position: 1D array (size N) of atomic distances (unit: Bohr radius)
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
    gaussians = np.exp(-gaussian_exponents[:, np.newaxis] * (position**2)[np.newaxis, :])
    
    # On applique les coefficients par multiplication matricielle
    molecular_orbitals = mo_coefficients @ gaussians
    
    return molecular_orbitals



def main() :
    pass


if __name__ == "__main__"  :
    main()