from extract_atomic_densities import load_density
from extract_urea_data import extract_urea_data
from extract_promolecule_density import load_density_promolecule

def main():
    centers, atom_dict = extract_urea_data("./../Multiwfn_Data/urea.txt")
    atom_dict = load_density("./../Multiwfn_Data/Isolated_Densities/azote_6.txt", atom_dict)
    promolecule_density = load_density_promolecule("./../Multiwfn_Data/urea_promolecule_density.txt")

    print(centers, atom_dict, promolecule_density)

if __name__ == "__main__":
    main()