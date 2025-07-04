#-----------------------------------------
# Generate the event matrix.
#-----------------------------------------
import numpy as np
import spglib
from ase import Atoms
from ase.io import vasp

def find_point_group_operation(atoms, tolerance=1e-5):
    """
    Extract rotational symmetry operations of a crystal structure.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atomic structure representation containing lattice and atomic positions
    tolerance : float, optional
        Symmetry search tolerance (in Angstroms). Controls precision for
        symmetry detection. Default: 1e-5
        
    Returns
    -------
    numpy.ndarray
        3x3 rotation matrices representing crystallographic symmetry operations.
        Shape: (n_operations, 3, 3). Rotations are defined in lattice coordinates.
    
    Notes
    -----
    - Uses spglib library for symmetry analysis
    - Rotation matrices operate on fractional coordinates
    - Combined with translations (not returned) these form full symmetry operations
    """
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec=tolerance)
    symmetry_ops = spglib.get_symmetry(cell, symprec=tolerance)
    rotations = symmetry_ops['rotations']
    return rotations

def find_near_center_site(atoms):
    """
    Search site near to center.
    
    Parameters:
    -----------
        lattice(atoms): The lattice you want to find the center.
    
    Returns:
    --------
        nearest_centre_site_index: The index of the nearest site to the center
    """
    center_site = Atoms('H', scaled_positions=[(0.5, 0.5, 0.5)], cell=atoms.cell)
    atoms = atoms + center_site
    dist = [atoms.get_distance(i, -1, mic=True) 
        for i in range(sum(atoms.numbers))]
    nearest_centre_site_index = np.argsort(dist)[1]
    return nearest_centre_site_index

def main():
    from parse import parse_structure
    atoms = parse_structure('./Fe.vasp')
    rots = find_point_group_operation(atoms,1E-5)
    print(rots)

if __name__ == "__main__":
    main()