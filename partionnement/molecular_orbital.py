import numpy as np


def molecular_orbitals(position, molecular_orbitals_coefficients, gaussian_exponents):
    """
    Return each electronic orbital of the molecule, as the sum of ..
    r : 1D-nparray of the position
    """
    l = len(gaussian_exponents)
    nb_mo = molecular_orbitals_coefficients.shape[0]
    # Il faut faire un test molecular_coefficients.shape[1] == l

    for i in range(0, nb_mo) :
        for i in range (0,l) :
            for r in position : 
    
    return molecular_orbitals


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
    if mo_coefficients.shape[1] != len(gaussian_exponents):
        raise ValueError(
            f"Incohérence : {mo_coefficients.shape[1]} coeffs mais {len(gaussian_exponents)} exposants !"
        )
    
    # On calcule toutes les gaussiennes pour toutes les positions d'un coup
    gaussians = np.exp(-gaussian_exponents[:, np.newaxis] * (position**2)[np.newaxis, :])
    
    # On applique les coefficients par multiplication matricielle
    molecular_orbitals = mo_coefficients @ gaussians
    
    return molecular_orbitals


def main() :
    position = np.arange(0, 2e-10, 1e-12)
    gaussian = np.array([0, 1e-1])
    mo_coefficients = np.array([1, 1])
    molecular_orbitals = molecular_orbitals(position, mo_coefficients, gaussian)
    print(molecular_orbitals)

if __name__ == "__main__"  :
    main       