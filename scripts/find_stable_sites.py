#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used here for find the stable configuration of
# interstitial hydrogen in a supercell.
#-------------------------------------------------------------
import os
import re
import copy
import numpy as np
from ase import Atom
from ase import Atoms
from ase.io import vasp
from ase.build import make_supercell

from supercell import sort_atoms

def clean(filepath):
    """Clean old file.
    
    Args:
        filepath(string): The work path.
        
    Returns:
        None.
    """
    # The folder is searched for the existence of the previously generated 
    # POSCAR file, and if it exists, it is deleted to avoid repeated generation.
    files = os.listdir(filepath)
    clean_mark = False
    for file in files:
        if re.search(r"POSCAR", file):
            os.remove(filepath+"/"+file)
            clean_mark = True
    if clean_mark:
        print("Old POSCAR files is already deleted.")

def get_unit(matrix):
    """Convert matrix row vector to unit vector.

    Args:
        matrix (ndarray): The matrix which you want to convert.

    Returns:
        ndarray: Converted matrix.
    """
    row_norm = np.linalg.norm(matrix, axis=1)
    unit_matrix = matrix.T / row_norm
    return unit_matrix.T

def split_sublattice(atoms, order_list):
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
    return sub_lattice_list

def mesh_generation(element_list, primitive_cell, Transform_matrix, scale_matrix, mesh_dist, bond_length):
    # $(a \ b \ c)(\vec{i} \ \vec{j} \ \vec{k}).t = (a \ b \ c)T^{-1}T(\vec{i} \ \vec{j} \ \vec{k}).t$ 
    # Here the T is the transformation matrix. The transfomation should keep the atom vector norm. 
    # Here we only consider a space which can be rotated and generate the whole space.
    element_list.append('H')
    pm_cell = vasp.read_vasp(primitive_cell)
    sub_lattice_list = split_sublattice(pm_cell, element_list)
    latti_vec_matrix = pm_cell.cell
    search_space = scale_matrix@latti_vec_matrix
    directon_matrix = get_unit(latti_vec_matrix)
    norms = np.linalg.norm(search_space, axis=1)
    mesh_size = []
    for norm in norms:
        mesh_size.append(int((norm//mesh_dist)+1))
    print(directon_matrix)
    print(mesh_size)
    mesh_site_list = [mesh_dist*i*directon_matrix[0, :]+
                      mesh_dist*j*directon_matrix[1, :]+
                      mesh_dist*k*directon_matrix[2, :]
                      for i in range(mesh_size[0])
                      for j in range(mesh_size[1])
                      for k in range(mesh_size[2])]
    mesh_excluded_list = []
    H_mesh = Atoms(cell=latti_vec_matrix)
    i = 0
    for site in mesh_site_list:
        H = Atom('H', position=site)
        MH = sub_lattice_list[0] + H
        OH = sub_lattice_list[1] + H
        dist_MH = np.array(
            [MH.get_distance(i,-1,mic=True) 
            for i in range(len(MH)-1)]).min()
        dist_OH = np.array(
            [OH.get_distance(i,-1,mic=True) 
            for i in range(len(OH)-1)]).min()
        if dist_OH > bond_length[0] and dist_MH > bond_length[1]:
            i = i + 1
            supercell = make_supercell(pm_cell, Transform_matrix, wrap=True)
            supercell = sort_atoms(supercell + H, element_list)
            vasp.write_vasp("POSCAR_" + str(i), 
                            supercell,
                            direct=True,
                            wrap=True)
            mesh_excluded_list.append(site)
            H_mesh.extend(H)
    pm_cell_Hmesh = pm_cell + H_mesh
    vasp.write_vasp("Hmesh.vasp", 
                    pm_cell_Hmesh, 
                    direct=True, 
                    wrap=True)
    print(f"{len(H_mesh)} POSCARs is generated")

def main():
    clean('./')
    Transform_matrix = np.array([[0, 2, 3],
                                 [2, 0, 3],
                                 [2, 2, 0]])
    mesh_dist = 0.5
    bond_length = [0.8, 1.0]
    scale_matrix = np.diag([0.5, 0.5, 0.5])
    element_list = ['Si', 'O']
    mesh_generation(element_list, 'stishovite.vasp', Transform_matrix, scale_matrix, mesh_dist, bond_length)

if __name__ == "__main__":
    main()
