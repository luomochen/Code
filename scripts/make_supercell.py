#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------
# A script for creating supercells with enhanced functionality.
#
# Features:
# 1. Outputs supercell matrix in YAML format
# 2. Smart supercell creation: 
#    - Uses existing supercell.yaml if available
#    - Generates optimal cell using ASE when YAML doesn't exist
# 3. Automatic atomic sorting by element order
#
#--------------------------------------------------------------
import os
import copy
import yaml
import numpy as np
from ase import Atoms
from ase.io import vasp
from ase.build import make_supercell, find_optimal_cell_shape

def sort_atoms(atoms: Atoms, order_list: list) -> Atoms:
    """Sort atoms in ASE Atoms object according to specified element order.
    
    Args:
        atoms (Atoms): ASE Atoms object to sort
        order_list (list): Desired element ordering (e.g., ['Si', 'O', 'H'])
        
    Returns:
        Atoms: New sorted Atoms object
    """
    element_num = len(order_list)
    # Create element-specific sub-lattice
    sub_lattice_list = []
    for i in range(element_num):
        sub_lattice = copy.deepcopy(atoms)
        index_list = []
        for ind, iatom in enumerate(atoms):
            if iatom.symbol != order_list[i]: 
                index_list.append(ind)
        del sub_lattice[index_list]
        sub_lattice_list.append(sub_lattice)
    sorted_atoms = Atoms(pbc=True, cell=atoms.cell)
    # Combine sorted components
    for element in sub_lattice_list:
        sorted_atoms = sorted_atoms + element
    return sorted_atoms

def create_supercell(primitive_cell: Atoms, multiple: int) -> Atoms:
    """Generate optimal supercell using ASE's cell shape optimization.
    
    Args:
        primitive_cell (Atoms): Input primitive cell
        multiple (int): Target number of atoms in supercell
        
    Returns:
        Atoms: Generated supercell
    """
    P = find_optimal_cell_shape(primitive_cell.cell, multiple, 'sc')
    return make_supercell(primitive_cell, P), P

def save_supercell_matrix(matrix: np.ndarray, filename: str = 'supercell.yaml'):
    """Save supercell matrix to YAML file with metadata.
    
    Args:
        matrix (np.ndarray): Supercell transformation matrix
        filename (str): Output filename (default: supercell.yaml)
    """
    data = {
        'supercell_matrix': matrix.tolist(),
        'transposed_matrix': matrix.T.tolist(),
        'note': 'Matrix format: ASE uses P matrix, Phonopy/VESTA use P.T'
    }
    
    with open(filename, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

def load_supercell_matrix(filename: str = 'supercell.yaml') -> np.ndarray:
    """Load supercell matrix from YAML file.
    
    Args:
        filename (str): Input filename (default: supercell.yaml)
        
    Returns:
        np.ndarray: Supercell transformation matrix
    """
    with open(filename, 'r') as f:
        data = yaml.safe_load(f)
    return np.array(data['supercell_matrix'])

def process_supercell(input_file: str, multiple: int, element_order: list):
    """Main processing pipeline for supercell creation.
    
    Args:
        input_file (str): Input VASP structure file (POSCAR/CONTCAR)
        multiple (int): Target size multiplier for supercell
        element_order (list): Element sorting order
    """
    primitive = vasp.read_vasp(input_file)
    yaml_file = 'supercell.yaml'
    
    # Supercell matrix handling
    if os.path.exists(yaml_file):
        print(f"Loading existing supercell matrix from {yaml_file}")
        P = load_supercell_matrix(yaml_file)
    else:
        print("Generating new optimal supercell matrix")
        supercell, P = create_supercell(primitive, multiple)
        save_supercell_matrix(P, yaml_file)
    
    # Create final supercell
    final_supercell = make_supercell(primitive, P)
    sorted_supercell = sort_atoms(final_supercell, element_order)
    vasp.write_vasp('supercell.vasp', sorted_supercell, direct=True)

def main():
    process_supercell(
        input_file='coesite_primitive_cell.vasp',
        multiple=6,
        element_order=['Si', 'O']
    )

if __name__ == "__main__":
    main()