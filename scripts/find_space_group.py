#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to identify space groups.
#-------------------------------------------------------------
import os
import re
from ase import spacegroup
from ase.io import vasp

precs = {"Hyperfine": 0.00001, 
    "Superfine":0.0001, 
    "Ultrafine":0.001, 
    "Fine":0.01, 
    "Medium":0.1, 
    "Coarse":0.2, 
    "Ultracoarse":0.5}

def find_vasp_type_file():
    prob_vasp_type_files = []
    files = os.listdir()
    for file in files:
        file_exist = re.search(r"POSCAR", file) \
            or re.search(r"CONTCAR", file) \
            or re.search(r".vasp", file)
        if file_exist:
            prob_vasp_type_files.append(file)
    return prob_vasp_type_files

def find_space_group(acc, prob_vasp_type_files):
    if acc in precs.keys():
        prec = precs[acc]
    else:
        print("The specified precise level is not in the list,")
        print("the default level is used.")
        prec = precs['Hyperfine']
    if prob_vasp_type_files is []:
        print("VASP structure file not find!")
    else:    
        for file in prob_vasp_type_files:
            try:
                atoms = vasp.read_vasp(file)
                sg = spacegroup.get_spacegroup(atoms=atoms, symprec=prec)
                print(f"{file}: {sg.symbol}, No. {sg.no}")
            except IndexError:
                print(f"{file} format error!")

def constraint_symmetry(atoms, ):
    aligned_atoms = spacegroup.align(atoms)
    pass

def main():
    for key, value in precs.items():
        print(f"{key}: {value}")
    acc = input("Please choose an accuracy level of space group identification: ")
    prob_vasp_type_files = find_vasp_type_file()
    find_space_group(acc, prob_vasp_type_files)

if __name__ == "__main__":
    main()