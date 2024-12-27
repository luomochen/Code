#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for searching suitable supercell and 
# making supercell.
#---------------------------------------------
import copy
from ase import Atoms
from ase.io import vasp
from ase.build import make_supercell, find_optimal_cell_shape

def sort_atoms(atoms, order_list):
    """Sort atom in Atoms with signified element order.

    Args:
        atoms (Atoms): The atoms object which you want to sort.
        order_list (list[string]): The list contain element name with a signified order.

    Returns:
        Atoms: Sorted atoms object.
    """
    element_num = len(order_list)
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
    for element in sub_lattice_list:
        sorted_atoms = sorted_atoms + element
    return sorted_atoms

def creat_supercell(filename, multiple, element_list, supercell_type=0):
    """
    """
    # Determine the supercell type.
    if supercell_type == 0:
        sctype = 'sc'
    else:
        sctype = 'fcc'
    primitivecell = vasp.read_vasp(filename)
    # Determine whether start search suitable supercell or not.
    if isinstance(multiple, int):
        P = find_optimal_cell_shape(primitivecell.cell, multiple, sctype)
    else:
        P = multiple
    supercell = make_supercell(primitivecell, P)
    supercell_sort = sort_atoms(supercell, element_list)
    lattice_parameter = supercell_sort.cell.cellpar()
    vasp.write_vasp('supercell.vasp', supercell_sort, direct=True)
    with open('supercell.dat', 'w') as file:
        file.write(f'ase: \n{P}\n')
        file.write(f'phonopy: \n{P.T}\n')
        file.write(f'lattice parameter: {lattice_parameter}')
    return supercell_sort.cell.cellpar()

def main():
    element_list = ['Si', 'O']
    lattice_parameter = creat_supercell('PbO2.vasp', 9, element_list, 1)
    print(f'lattice parameters: {lattice_parameter}')

if __name__ == "__main__":
    main()