import numpy as np


def molecular_orbitals(position, mo_coefficients, gaussian_exponents):
    """
    Returns the OM-value for each molecular orbital.
    
    Parameters:
    - position : tableau 1D (taille N) de distances atomiques  (unité : Bohr radius)
    - mo_coefficients : tableau 2D (taille M, 36), coefficients des fonctions de base des orbitales
    - gaussian_exponents : tableau 1D (taille 36), exposants des fonctions gaussiennes de base
    Returns:
    - 2D-array containing the values of the OM for each position-value.
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