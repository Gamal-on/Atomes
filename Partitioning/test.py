import numpy as np

def sum_weights(atom_dict : dict):
    assert len({len(a['coord']) for a in atom_dict.values()}) == 1, "Tailles de 'coord' incohérentes entre atomes"

    weights_sum = sum(atome['weight'] for atome in atom_dict.values())
    assert np.allclose(weights_sum, 1.0), f"Écart max à 1 : {np.max(np.abs(weights_sum - 1)):.2e}"


