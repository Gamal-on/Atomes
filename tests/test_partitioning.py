import numpy as np
import pytest

def test_coord_lengths_consistent(atom_dict):
    """Toutes les arrays 'coord' doivent avoir la même taille (même grille)."""
    tailles = {nom: len(atome["coords"]) for nom, atome in atom_dict.items()}
    assert len(set(tailles.values())) == 1, f"Tailles de 'coord' incohérentes : {tailles}"


def test_hirshfeld_weights_sum_to_one(atom_dict):
    """En chaque point de la grille où la densité totale n'est pas négligeable,
    la somme des poids de Hirshfeld doit valoir 1."""
    weights_sum = sum(atome["hirshfeld_weight"] for atome in atom_dict.values())
    total_density = sum(atome["rho"] for atome in atom_dict.values())

    mask = total_density > 1e-10  # seuil à ajuster selon l'échelle de tes densités
    ecart = np.max(np.abs(weights_sum[mask] - 1))
    assert np.allclose(weights_sum[mask], 1.0), (
        f"Écart max à 1 (points significatifs) : {ecart:.2e}"
    )