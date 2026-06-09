from state import ProjectState
from Extract_Parameters.extract_urea_data import extract_urea_data
from Grid.becke_grid import generate_grid, plot_multicenter_grid, export_grid_for_multiwfn

state = ProjectState()
state.centers = extract_urea_data("Multiwfn_Data/urea.wfn")


