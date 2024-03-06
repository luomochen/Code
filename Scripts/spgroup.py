#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to identify space groups
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
for key, value in precs.items():
    print(f"{key}: {value}")
print("----------------------------------------------------------------------")
acc = input("Please enter the accuracy level of space group identification: ")
print("----------------------------------------------------------------------")
if acc in precs.keys():
    prec = precs[acc]
else:
    print("The specified precise level is not in the list, the default level is used.")
    prec = precs['Hyperfine']
files = os.listdir()
opjud = []

for file in files:
    if re.search(r"POSCAR", file):
        atoms = vasp.read_vasp(file)
        sg = spacegroup.get_spacegroup(atoms=atoms, symprec=prec)
        print(f"{file}: {sg.symbol}, No. {sg.no}")
        opjud.append(True)
    elif re.search(r"CONTCAR", file):
        atoms = vasp.read_vasp(file)
        sg = spacegroup.get_spacegroup(atoms=atoms, symprec=prec)
        print(f"{file}: {sg.symbol}, No. {sg.no}")
        opjud.append(True)
    elif re.search(r".vasp", file):
        atoms = vasp.read_vasp(file)
        sg = spacegroup.get_spacegroup(atoms=atoms, symprec=prec)
        print(f"{file}: {sg.symbol}, No. {sg.no}")
        opjud.append(True)
    else:
        opjud.append(False)

if True not in opjud:    
    print("VASP structure file not find!")