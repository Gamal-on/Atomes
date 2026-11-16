"""Shared state of the project.

Only main.py writes into the state: the other modules receive plain arrays/dicts
and return their results, they never import or modify the state themselves.

Units: atomic units, as in Multiwfn (lengths in bohr, densities in e/bohr^3).
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Grid:
    """Multicenter grid on which every density of the project is evaluated.

    points        : (N_points, 3) cartesian coordinates, in bohr.
    n_radial      : number of radial shells per atom (Gauss-Chebyshev).
    lebedev_order : order of the Lebedev angular grid.
    r_m           : radial mapping parameter.
    weights       : (N_points,) integration weights (radial x angular x Becke),
                    needed to integrate quantities over the grid.
    """

    points: np.ndarray
    n_radial: int | None = None
    lebedev_order: int | None = None
    r_m: float | None = None
    weights: np.ndarray | None = None

    @property
    def n_points(self) -> int:
        return len(self.points)


@dataclass
class ProjectState:
    """Everything computed along the pipeline.

    grid                : the common grid (see Grid).
    centers             : (N_atoms, 3) nuclear positions in bohr, ordered by atom number.
    atom_dict           : per-atom data, keyed by atom number (1, 2, ...):
                          {'element': str, 'Z': float, 'center': (3,) array,
                           'rho': (N_points,) isolated atomic density,
                           'hirshfeld_weight': (N_points,) array, ...}
    molecular_density   : (N_points,) electron density of the molecule.
    promolecule_density : (N_points,) sum of the isolated atomic densities.
    """

    grid: Grid | None = None
    centers: np.ndarray | None = None
    atom_dict: dict[int, dict] = field(default_factory=dict)
    molecular_density: np.ndarray | None = None
    promolecule_density: np.ndarray | None = None
