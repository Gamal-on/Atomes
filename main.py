from Extract_Parameters.extract_urea_data import extract_urea_data
from Grid.becke_grid import generate_grid, plot_multicenter_grid, export_grid_for_multiwfn



def main():

    # Charge data
    
    centers, atome_dict = extract_urea_data("./Multiwfn_data/urea.wfn")
    
    # Grid, default paramaters
    radial_lenght = 20
    ordre_choisi = 53
    r_m = 0.5
    
    grid = generate_grid(radial_lenght, centers, ordre_choisi, r_m=0.5, units='bohr')
    
    # Si les 
    # export_grid(grid, output_path, centers_bohr, r_cutoff=50.0), si jamais, on a déjà le fichier
# plot_multicenter_grid (optional)

# grid = generate_grid(20, centers, ordre_choisi=53, r_m=0.5, units='bohr')
# export_grid_for_multiwfn(grid, "urea_grid_multiwfn.txt", centers)
# plot_multicenter_grid(grid, centers)

# Partitionnement 

# les fichiers de poids sont dans Data/
# pour plot le poids d'un atome: main_partitioning

if __name__ == "__main__":
    main()