#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to unpack the structure file and obtain 
# the structure file of each element for differential charge 
# density calculation.
#-------------------------------------------------------------
import re
import copy
from ase.io import vasp

print("----------------------------------------------------------------------")
file = input("Please specify the VASP structure file that needs to be split: ")
print("----------------------------------------------------------------------")

try:
    atoms = vasp.read_vasp(file)
    order_list = re.findall(r"(\D+)", str(atoms.symbols))

    element_num = len(order_list)
    for i in range(element_num):
        sub_lattice = copy.deepcopy(atoms)
        index_list = []
        for ind, iatom in enumerate(atoms):
            if iatom.symbol != order_list[i]: 
                index_list.append(ind)
        del sub_lattice[index_list]
        vasp.write_vasp(order_list[i]+".vasp", sub_lattice, direct=True)
        print(f"{order_list[i]}.vasp is generated!")
except FileNotFoundError:
    print("VASP structure file not find!")