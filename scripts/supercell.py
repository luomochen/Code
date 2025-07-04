#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------
# A script for creating supercells with enhanced functionality.
#--------------------------------------------------------------
import yaml
import argparse
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell, find_optimal_cell_shape

def sort_atoms(atoms: Atoms, order_list: list) -> Atoms:
    """Sort atoms in ASE Atoms object according to specified element order."""
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    cell = atoms.cell
    pbc = atoms.pbc

    sorted_indices = []
    for element in order_list:
        sorted_indices.extend([i for i, s in enumerate(symbols) if s == element])

    sorted_atoms = Atoms(
        [symbols[i] for i in sorted_indices],
        positions=[positions[i] for i in sorted_indices],
        cell=cell,
        pbc=pbc
    )
    return sorted_atoms

def save_supercell_data(pstruc_file, multiple, element_order, expand_type, P):
    """Save supercell data to YAML file with metadata."""
    data = {
        "primitive_structure_file": pstruc_file,
        "expand_times": multiple,
        "element_order": element_order,
        "supercell_shape": expand_type,
        "supercell_matrix": P.tolist(),
        "transposed_matrix": P.T.tolist(),
        "note": "Matrix format: ASE uses P matrix, Phonopy/VESTA use P.T"
    }

    with open("supercell.yaml", 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False)

def load_supercell_data(filename="supercell.yaml"):
    """Load supercell data from YAML file."""
    with open(filename, 'r') as f:
        data = yaml.safe_load(f)
        pstruc_file = data["primitive_structure_file"]
        multiple = data["expand_times"]
        element_order = data["element_order"]
        expand_type = data["supercell_shape"]
        P = np.array(data["supercell_matrix"])
    return pstruc_file, multiple, element_order, expand_type, P

def process_supercell(pstruc_file, multiple, element_order, expand_type, parameter_read):
    """Main processing pipeline for supercell creation."""
    if parameter_read:
        try:
            pstruc_file, multiple, element_order, expand_type, P = load_supercell_data("supercell.yaml")
            primitive = read(pstruc_file)
        except FileNotFoundError:
            print("Cannot find the supercell.yaml file.")
            return
    else:
        print("Generating new optimal supercell matrix")
        primitive = read(pstruc_file)

        if expand_type not in ["sc", "fcc"]:
            raise ValueError("expand_type must be 'sc' or 'fcc'.")

        P = find_optimal_cell_shape(primitive.cell, multiple, expand_type)

        # If no element_order is provided, use the default chemical order
        if element_order is None:
            element_order = list(dict.fromkeys(primitive.get_chemical_symbols()))
        
        save_supercell_data(pstruc_file, multiple, element_order, expand_type, P)

    final_supercell = make_supercell(primitive, P, expand_type)
    sorted_supercell = sort_atoms(final_supercell, element_order)
    write('supercell.vasp', sorted_supercell, format='vasp', direct=True)

def main():
    parser = argparse.ArgumentParser(description="Run making supercell using YAML input.")
    parser.add_argument("-p", "--primitive", 
                        type=str, default="primitive_cell.vasp", 
                        help="Path to primitive cell you need to expand.")
    parser.add_argument("-m", "--multiple",
                        type=int, default=1,
                        help="The times of cell expand.")
    parser.add_argument("-o", "--order", 
                        type=str, nargs='+', default=None,
                        help="Designate the order of elements in the supercell.")
    parser.add_argument("-e", "--expand_type",
                        default="sc", type=str,
                        help="The way to expand supercell, sc or fcc.")
    parser.add_argument("-r", "--read",
                        action='store_true',
                        help="Read cell expansion parameters from YAML file.")
    args = parser.parse_args()

    process_supercell(
        pstruc_file=args.primitive,
        multiple=args.multiple,
        element_order=args.order,
        expand_type=args.expand_type,
        parameter_read=args.read
    )

if __name__ == "__main__":
    main()