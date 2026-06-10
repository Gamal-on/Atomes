from Extract_Parameters.extract_promolecule_density import load_density_promolecule
from state import ProjectState
from Extract_Parameters.extract_urea_data import extract_urea_data
from Grid.becke_grid import generate_grid, plot_multicenter_grid
from Extract_Parameters.extract_atomic_densities import load_density
from Partitioning.calcul_weights import hirshfeld_weight
from Partitioning.plot_weights import plot_weight

state = ProjectState()

## --- 1. Parsing the information from Multiwfn files ---

# Molecular (centers, number of atom ...) from urea.wfn

state.centers, state.atom_dict = extract_urea_data("Multiwfn_Data/urea.txt")

# Parsing atomic densities info from Multiwfn_Data/Isolated_Densities

state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_3.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/carbon_1.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_4.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_5.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/oxygen_2.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_8.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_7.txt", state.atom_dict)

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

# Plotting Hirshfeld weights : for the  atom number x, enter 

grid_flat = state.grid.reshape(-1, 3)
plot_weight(grid_flat,state.atom_dict)

