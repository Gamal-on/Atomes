import numpy as np

from plot_weights import plot_weight
from calcul_weights import hirshfeld_weight


# Reading data from .txt files

data_promolecule = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/urea_promolecule_density.txt", skiprows=1)
data_carbon_1 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/carbon_1.txt", skiprows=1)
data_oxygen_2 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/oxygen_2.txt", skiprows=1)
data_azote_3 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/azote_3.txt", skiprows=1)
data_hydrogen_4 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/hydrogen_4.txt", skiprows=1)
data_hydrogen_5 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/hydrogen_5.txt", skiprows=1)
data_azote_6 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/azote_6.txt", skiprows=1)
data_hydrogen_7 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/hydrogen_7.txt", skiprows=1)
data_hydrogen_8 = np.loadtxt("../Multiwfn_data/Weights_Urea_Grid/hydrogen_8.txt", skiprows=1)

# Extracting density values from the data arrays

density_promolecule = data_promolecule[:, 3]
density_carbon_1 = data_carbon_1[:, 3]  
density_oxygen_2 = data_oxygen_2[:, 3]
density_azote_3 = data_azote_3[:, 3]
density_hydrogen_4 = data_hydrogen_4[:, 3]
density_hydrogen_5 = data_hydrogen_5[:, 3]
density_azote_6 = data_azote_6[:, 3]
density_hydrogen_7 = data_hydrogen_7[:, 3]
density_hydrogen_8 = data_hydrogen_8[:, 3]    

# Extracting grid coordinates from the data arrays

grid = data_promolecule[:, :3]

# Calculating Hirshfeld weights for each atom

weight_carbon_1 = hirshfeld_weight(density_carbon_1, density_promolecule)
weight_oxygen_2 = hirshfeld_weight(density_oxygen_2, density_promolecule)
weight_azote_3 = hirshfeld_weight(density_azote_3, density_promolecule)
weight_hydrogen_4 = hirshfeld_weight(density_hydrogen_4, density_promolecule)
weight_hydrogen_5 = hirshfeld_weight(density_hydrogen_5, density_promolecule)
weight_azote_6 = hirshfeld_weight(density_azote_6, density_promolecule)
weight_hydrogen_7 = hirshfeld_weight(density_hydrogen_7, density_promolecule)
weight_hydrogen_8 = hirshfeld_weight(density_hydrogen_8, density_promolecule)



def main() :
    plot_weight(grid, weight_oxygen_2)

if __name__ == "__main__":
    main()