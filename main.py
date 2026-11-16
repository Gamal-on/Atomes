from pathlib import Path

from Extract_Parameters.extract_promolecule_density import load_density_promolecule
from state import ProjectState
from Extract_Parameters.extract_centers_atoms import extract_centers_atoms
from Grid.becke_grid import generate_grid, plot_multicenter_grid
from Extract_Parameters.extract_eletronic_densities import load_electronic_density
from Partitioning.calcul_weights import hirshfeld_weight
from Partitioning.plot_weights import plot_weight
from Partitioning.calcul_densities import calcul_final_densities

state = ProjectState()

## --- 1. Parsing the information from Multiwfn files ---

# Molecular (centers, number of atom ...) from urea.wfn

state.centers, state.atom_dict = extract_centers_atoms("Multiwfn_Data/urea.txt")

# Parsing isolated atomic densities: one file per atom in Multiwfn_Data/Isolated_Densities.
# The atom each file belongs to is read from its header (e.g. "N  3"), so the order does not matter.

for density_file in sorted(Path("Multiwfn_Data/Isolated_Densities").glob("*.txt")):
    state.atom_dict = load_electronic_density(density_file, state.atom_dict)

# Parsing promolecule density 

state.promolecule_density = load_density_promolecule("Multiwfn_Data/urea_promolecule_density.txt")


## --- 2. Generating the grid ---

# Do not change the parameters, unless you are willing to generate new density files with Multiwfn

state.grid = generate_grid(19, state.centers, 53)

# Plotting the grid : to visualize the grid, enter the following line of code : 

# plot_multicenter_grid(state.grid, state.centers)



## --- 3. Hirshfeld weights ---

#Calculating weights

state.atom_dict = hirshfeld_weight(state.atom_dict, state.promolecule_density)
#state.atom_dict = calcul_final_densities(state.atom_dict, )

# Plotting Hirshfeld weights : for the  atom number x, enter 

#grid_flat = state.grid.reshape(-1, 3)
#plot_weight(grid_flat,state.atom_dict)