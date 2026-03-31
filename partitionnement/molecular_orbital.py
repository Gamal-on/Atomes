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
    position = np.arange(1e-10, 2e-10, 1e-11)
    gaussian = np.array([5.48e3, 8.25e2, 1.88e2])
    mo_coefficients = np.array([[8.25e-1, 1.51, 2.45], [1.74e-1, 3.21e-1, 5.18e-1]])
    result = molecular_orbitals(position, mo_coefficients, gaussian)
    print(result, len(position), len(result), position[0])
    r1 = position[0]
    term1 = (8.25e-1)*np.exp(-(5.48e3)*(r1)**2)
    term2 = (1.51)*np.exp(-(8.25e2)*(r1)**2)
    term3 = (2.45)*np.exp(-(1.88e2)*(r1)**2)
    print(term1+term2+term3)


if __name__ == "__main__"  :
    main()