import numpy as np


def molecular_orbitals(position, mo_coefficients, gaussian_exponents):
    """
    Calcule les orbitales moléculaires comme une somme de gaussiennes.
    
    Parameters:
    - position : tableau 1D (taille N)
    - mo_coefficients : tableau 2D (taille M, 36)
    - gaussian_exponents : tableau 1D (taille 36)
    
    Returns:
    - Un tableau 2D de taille (M, N) contenant les valeurs des orbitales.
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