from molecular_orbital import molecular_orbitals
import numpy as np

def test_molecular_orbitals():
    position = np.arange(1e-10, 2e-10, 1e-11)
    gaussian = np.array([5.48e3, 8.25e2, 1.88e2])
    mo_coefficients = np.array([[8.25e-1, 1.51, 2.45], [1.74e-1, 3.21e-1, 5.18e-1]])
    r1 = position[0]

    # Experimental result
    result_exp = molecular_orbitals(position, mo_coefficients, gaussian)

    # Theorical result
    term1 = (8.25e-1)*np.exp(-(5.48e3)*(r1)**2)
    term2 = (1.51)*np.exp(-(8.25e2)*(r1)**2)
    term3 = (2.45)*np.exp(-(1.88e2)*(r1)**2)
    result_th = term1 + term2 + term3

    assert result_th == result_exp[0,0]