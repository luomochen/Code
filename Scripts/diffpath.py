#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to obtain the diffusion path and generate
# structure file.
#-------------------------------------------------------------
import os
import re
from ase import Atoms
from ase.io import vasp

# When generating a new structured file with a diffusion path, 
# you need to first place a perfectly structured folder outside 
# the folder. Otherwise the structure file will only has diffusion
# path.
if os.path.exists("../pft.vasp"):
    diffpath = vasp.read_vasp("../pft.vasp")
else:
    print("Perfect crystal structure file not found!")
    atoms = vasp.read_vasp("ini/CONTCAR")
    diffpath = Atoms(pbc=True, cell=atoms.cell)

files = os.listdir()
for file in files:
    if re.match(r"0\d", file):
        if os.path.exists(file+"/CONTCAR"):
            atoms = vasp.read_vasp(file+"/CONTCAR")
        else:
            atoms = vasp.read_vasp(file+"/POSCAR")
        for ind, iatom in enumerate(atoms):
            if iatom.symbol == 'H':
                diffpath.extend(iatom)
vasp.write_vasp("diffpath.vasp", diffpath, direct=True)