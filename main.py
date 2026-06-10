from state import ProjectState
from Extract_Parameters.extract_urea_data import extract_urea_data
from Grid.becke_grid import generate_grid, plot_multicenter_grid, export_grid_for_multiwfn
from Extract_Parameters.extract_atomic_densities import load_density

state = ProjectState()
state.centers = extract_urea_data("Multiwfn_Data/urea.wfn")

# Parsing the molecular info from urea.wfn

state.centers, state.atom_dict = extract_urea_data("Multiwfn_Data/urea.wfn")

# Parsing desnities info from Multiwfn_Data/Isolated_Densities

state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)
state.atom_dict = load_density("Multiwfn_Data/Isolated_Densities/azote_6.txt", state.atom_dict)