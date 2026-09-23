from dataclasses import dataclass, field
import numpy as np

# On fait que lire le state, jamais autre chose !!!

@dataclass
class ProjectState:
    centers: np.ndarray = None
    atom_dict: dict = field(default_factory=dict) #{'element, 'coord', 'rho', 'weight}
    promolecule_density: np.ndarray = None
    grid_points: np.ndarray = None