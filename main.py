from Extract_Parameters.extract_promolecule_density import load_density_promolecule
from state import ProjectState
from Extract_Parameters.extract_urea_data import extract_urea_data
from Grid.becke_grid import generate_grid
from Extract_Parameters.extract_atomic_densities import load_density
from Partitioning.calcul_weights import hirshfeld_weight
from Partitioning.plot_weights import plot_weight

state = ProjectState()

# Parsing the molecular info from urea.wfn

state.centers, state.atom_dict = extract_urea_data("Multiwfn_Data/urea.txt")

# Grid : do not change the parameters, unless you are willing to generate new density files with Multiwfn

state.grid = generate_grid(19, state.centers, 53)


# Parsing atomic densities info from Multiwfn_Data/Isolated_Densities

state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_3.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/carbon_1.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_4.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_5.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/oxygen_2.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_8.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/hydrogen_7.txt", state.atom_dict)

# Parsing promolecule density from Multiwfn_Data/Promolecule_Density/urea_promolecule_density.txt

state.promolecule_density = load_density_promolecule("Multiwfn_Data/urea_promolecule_density.txt")

# Calculating Hirshfeld weights

state.atom_dict = hirshfeld_weight(state.atom_dict, state.promolecule_density)

# print(state.atom_dict[1]['hirshfeld_weight'])
# print(state.grid)

# Exporting Hirshfedl weights

# Plotting Hirshfeld weights : for the  atom number x, enter 

print(state.atom_dict)

# grid_flat = state.grid.reshape(-1, 3)
# plot_weight(grid_flat, state.atom_dict[1]['hirshfeld_weight'])