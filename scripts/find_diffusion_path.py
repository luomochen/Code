#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to find potential diffusion path.
#-------------------------------------------------------------
import os
import re
import numpy as np
from ase import Atoms
from ase import spacegroup
from ase.io import vasp
from ase.build import make_supercell

from supercell import sort_atoms

def backtocell(Direct):
    """Move the direct coordinate into the cell.

    Args:
        Direct (ndarray): Direct coordinate.

    Returns:
        ndarray: Moved direct coordinate.
    """
    for i in range(3):
        if Direct[i] < 0:
            Direct[i] = Direct[i] + 1
        if Direct[i] > 1:
            Direct[i] = Direct[i] - 1
    return Direct

def find_equivalent_site(site, Transform_matrix, filename):
    """Read the primitive cell and Find the equivalent
    sites of designated site. You need to provide the vasp
    POSCAR file of primitive cell.
    
    Args:
        site(ndarray): Direct coordinate.
        Transform_matrix(ndarray): Cell expand matrix.
        filename(string): The vasp PODCAR file of primitive cell.
        
    Returns:
        ndarray: a set of equivalent sites in primitive cell of 
        input site.
        atoms: an atoms class of primitive cell.
    """
    primitive_cell = vasp.read_vasp(filename)
    sg = spacegroup.get_spacegroup(primitive_cell)
    transformed_site = backtocell(site@Transform_matrix)
    sites, kinds = sg.equivalent_sites(transformed_site)
    return sites, primitive_cell

def cell_expand(sites, transform_matrix, primitive_cell, perfect_element):
    """Creat supercell of primitive cell and marked all H equivalent sites
    supercell.
    
    Args:
        sites(ndarray): the equivalent sites of designated site, which comes
        from find_equivalent_site function.
        transform_matrix(ndarray): Cell expand matrix.
        primitive_cell(atoms): atoms class of primitive cell, which comes from
        find_equivalent_site function.
        perfect_element(list[str]): perfect crystal contained elements using for
        sort expanded supercell atoms. 
        
    Returns:
        atoms: a supercell expanded by transform matrix only contain all equivalent
        sites of defect marked by H atom.
        atoms: a supercell expanded by transform matrix of primitive cell.
    """
    sublattice_element = "H" + str(sites.shape[0])
    sublattice_atoms = Atoms(sublattice_element, scaled_positions=sites, 
                             cell=primitive_cell.cell, pbc=True)
    # Generate supercell.
    sublattice_atoms_supercell = make_supercell(sublattice_atoms, transform_matrix)
    supercell = make_supercell(primitive_cell, transform_matrix)
    supercell = sort_atoms(supercell, perfect_element)
    return sublattice_atoms_supercell, supercell

def find_near_center_site(lattice):
    """Search site near to center.
    
    Args:
        lattice(atoms): The lattice you want to find the center.
    
    Returns:
        integer: The index of the nearest site to the center
    """
    center_site = Atoms('H', scaled_positions=[(0.5, 0.5, 0.5)], cell=lattice.cell)
    lattice = lattice + center_site
    dist = [lattice.get_distance(i, -1, mic=True) 
        for i in range(sum(lattice.numbers))]
    #min_site_center_dist = np.sort(dist)[1]
    min_site_index = np.argsort(dist)[1]
    return min_site_index

def clean_generate(filepath):
    """Clean old file and generate the new needed floder.
    
    Args:
        filepath(string): The work path.
        
    Returns:
        None.
    """
    # Determines if the specified folder exists in the current path, 
    # and creates it if not.
    folder = os.path.exists(filepath)
    if not folder:
        os.mkdir(filepath)
        print(filepath+' does not exist, creat now.')
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

def potential_path(radius, accuracy, supercell, sublattice_atoms_supercell, filepath, pair_distance=False):
    """Searching potential diffusion path under designating rule. 1. In specific cutting
    distance. 2. exclude equal distance path.
    
    Args:
        radius(float): specific cutting distance.
        accuracy(float): the accuracy which you specific to compare two distances are whether
        equivalent.
        supercell(atoms): an ase atoms which contains perfect structre.
        sublattice_atoms_supercell(atoms): an ase atoms which only contains all equivalent sites.
        filepath(str): designate file path for save generate file. 
        pair_distance(float): if pair_distance = float, it mainly is figure out the pair site distance
        which you want to exclude. 
        
    Returns:
        None.
    """
    # Find the site nearest to the center.
    near_center_site = find_near_center_site(sublattice_atoms_supercell)
    # Calculate the distance, and select the sites in the cutoff raidus.
    dist = [sublattice_atoms_supercell.get_distance(i, near_center_site, mic=True) 
            for i in range(sum(sublattice_atoms_supercell.numbers))]
    diffusion_path = []
    for ind, idist in enumerate(dist):
        if idist < radius:
            diffusion_path.append([ind, idist])
    diffusion_path.sort(key=lambda x:x[1])
    # Output 1.
    with open (filepath+"/"+filepath+".dat", 'w') as output:
        output.write("Vesta count: Start from 1!!!\n")
        for path in diffusion_path:
            output.write(f"{path[0]+1}, {path[1]}\n")
    # Second time selection, here we will delete the pair site to decrease the 
    # calculation cost.
    if isinstance(pair_distance, float):
        pair_site_list = []
        #diffusion_path_copy = copy.deepcopy(diffusion_path)
        for path in diffusion_path:
            for index, ipath in enumerate(diffusion_path):
                distance = sublattice_atoms_supercell.get_distance(path[0], ipath[0], mic=True)
                if path[1] < ipath[1] and abs(distance - pair_distance) < accuracy:
                        pair_site_list.append(index)
        diffusion_path = [ipath for index, ipath in enumerate(diffusion_path) 
                          if index not in pair_site_list]
        # Output 2.
        with open (filepath+"/"+filepath+".dat", 'a') as output:
            output.write(f"Accuracy: {accuracy}, Pair_distance: {pair_distance}\n")
            output.write("Selective start and final state: \n")
            for path in diffusion_path:
                output.write(f"{path[0]+1}, {path[1]}\n")
    # Third time selection, the rule is excluding the equal length sites.
    for index_1, ipath_1 in enumerate(diffusion_path):
        for index_2, ipath_2 in enumerate(diffusion_path):
            if abs(ipath_1[1] - ipath_2[1]) < accuracy and  index_1 != index_2:
                diffusion_path.pop(index_2)
    # Output 3.
    i = 0
    with open (filepath+"/"+filepath+".dat", 'a') as output:
        output.write(f"Accuracy: {accuracy}, Cuttoff Radius: {radius}\n")
        output.write("Selective start and final state: \n")
    # Generate the POSCARs for further calculation.
    for path in diffusion_path:
        with open (filepath+"/"+filepath+".dat", 'a') as output:
            output.write(f"state{i}: {path[0]+1}, {path[1]}\n")
        state = supercell + sublattice_atoms_supercell[path[0]]
        vasp.write_vasp(filepath+"/POSCAR_"+str(i), state, direct=True)
        i = i + 1
    # Generate the POSCAR file contained all equivalent sites.
    vasp.write_vasp(filepath+"/"+filepath+".vasp", sublattice_atoms_supercell, direct=True)

def delete_duplicate_atoms(sublattice_atoms_supercell, pair_length, accuracy):
    duplicate_site_list = []
    for i in range(sum(sublattice_atoms_supercell.numbers)):
        if i not in duplicate_site_list:
            for j in range(sum(sublattice_atoms_supercell.numbers)):
                dist = sublattice_atoms_supercell.get_distance(i, j, mic=True)
                if abs(dist - pair_length) < accuracy:
                    duplicate_site_list.append(j)
    del sublattice_atoms_supercell[duplicate_site_list]
    return sublattice_atoms_supercell

def main():
    accuracy = 1E-5
    limit_radius = 4
    pair_distance = 0.05
    #cell_1 = vasp.read_vasp('coesite.vasp')
    #cell_2 = vasp.read_vasp('P21c_0GPa.vasp')
    #lattice_vector_matrix_1 = cell_1.cell
    #lattice_vector_matrix_2 = cell_2.cell
    #B = np.linalg.inv(lattice_vector_matrix_2@np.linalg.inv(lattice_vector_matrix_1))
    #transform_matrix_1 = np.array([[1, 0, 0],
    #                               [0, 1, 0],
    #                               [1, 0, 2]])@B
    transform_matrix_2 = np.array([[0, 2, 3],
                                                         [2, 0, 3],
                                                         [2, 2, 0]])    
    defect_states = 'Hi+'
    perfect_element = ['Si', 'O']
    pair_length = [0.0068405]
    defect_sites = np.array([0.9002059032290829, 0.1002244572524508, 0.1366027217647086])
    clean_generate(defect_states)
    sites, primitive_cell = find_equivalent_site(defect_sites, transform_matrix_2, 'stishovite.vasp')
    sublattice_atoms_supercell, supercell = cell_expand(sites, transform_matrix_2, primitive_cell, perfect_element)
    sublattice_atoms_supercell = delete_duplicate_atoms(sublattice_atoms_supercell, pair_length, accuracy)
    potential_path(limit_radius, accuracy, supercell, sublattice_atoms_supercell, defect_states, 0.64169)
        
if __name__ == "__main__":
    main()
