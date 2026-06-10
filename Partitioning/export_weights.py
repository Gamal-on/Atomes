import numpy as np
from pathlib import Path

def export_hirshfeld_weights(grid_coordinates, weights, filename, folder_name="Weights_Urea_Grid"):
    """Export the Hirshfeld weights to a text file."""

    data_to_save = np.column_stack((grid_coordinates, weights))

    output_dir = Path(folder_name)
    file_path = output_dir / filename
    
    np.savetxt(file_path, data_to_save)

