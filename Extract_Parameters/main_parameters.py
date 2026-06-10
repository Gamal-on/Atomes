from extract_atomic_densities import load_density
from extract_urea_data import extract_urea_data

def main():
    centers, atom_dict = extract_urea_data("./../Multiwfn_Data/urea.txt")
    atom_dict = load_density("./../Multiwfn_Data/Isolated_Densities/azote_6.txt", atom_dict)
    print(centers, atom_dict)

if __name__ == "__main__":
    main()