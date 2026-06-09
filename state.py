from dataclasses import dataclass, field
import numpy as np

# On fait que lire le state, jamais autre chose !!!

@dataclass
class ProjectState:
    centers: np.ndarray = None
