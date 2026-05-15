import numpy as np
from plotly import data


data_promolecule = np.loadtxt("Weights_Urea_Grid/urea_promolecule_density.txt", skiprows=1)
data_carbon_1 = np.loadtxt("Weights_Urea_Grid/carbon_1.txt", skiprows=1)
data_oxygen_2 = np.loadtxt("Weights_Urea_Grid/oxygen_2.txt", skiprows=1)
data_azote_3 = np.loadtxt("Weights_Urea_Grid/azote_3.txt", skiprows=1)
data_hydrogen_4 = np.loadtxt("Weights_Urea_Grid/hydrogen_4.txt", skiprows=1)
data_hydrogen_5 = np.loadtxt("Weights_Urea_Grid/hydrogen_5.txt", skiprows=1)
data_azote_6 = np.loadtxt("Weights_Urea_Grid/azote_6.txt", skiprows=1)
data_hydrogen_7 = np.loadtxt("Weights_Urea_Grid/hydrogen_7.txt", skiprows=1)
data_hydrogen_8 = np.loadtxt("Weights_Urea_Grid/hydrogen_8.txt", skiprows=1)

density_promolecule = data_promolecule[:, 3]
density_carbon_1 = data_carbon_1[:, 3]  
density_oxygen_2 = data_oxygen_2[:, 3]
density_azote_3 = data_azote_3[:, 3]
density_hydrogen_4 = data_hydrogen_4[:, 3]
density_hydrogen_5 = data_hydrogen_5[:, 3]
density_azote_6 = data_azote_6[:, 3]
density_hydrogen_7 = data_hydrogen_7[:, 3]
density_hydrogen_8 = data_hydrogen_8[:, 3]    

grid = data_promolecule[:, :3]

if __name__ == "__main__": 
    print(grid)