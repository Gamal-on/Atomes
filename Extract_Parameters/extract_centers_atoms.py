import re

import numpy as np

def extract_centers_atoms(file_path):
    # List of lists to store the rows of centers
    centers_table = []
    # Dictionary to group atoms and their corresponding numbers
    atom_dict = {}

    # Regular expression to match only the atom coordinate lines
    pattern = r"^\s*([A-Za-z]+)\s+(\d+)\s+\(CENTRE\s+(\d+)\)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+CHARGE\s*=\s*([\d.]+)"

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            match = re.match(pattern, line)
            if match:
                element = match.group(1)
                atom_num = int(match.group(2))
                x = float(match.group(4))
                y = float(match.group(5))
                z = float(match.group(6))
                charge = float(match.group(7))

                # Build the centers table as a list of lists (table of tables)
                centers_table.append([x, y, z])

                # Build the atom dictionary {atom_num: {'element', 'Z', 'center'}}
                atom_dict[atom_num] = {
                    'element': element,
                    'Z': charge,
                    'center': np.array([x, y, z]),
                }

    return np.array(centers_table), atom_dict


