#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to ensure that the lattice constants of 
# the structure file do not change.
#-------------------------------------------------------------
import os
import re
from ase.io import vasp

# You need to first place a perfectly structured folder outside 
# the work folder.
if os.path.exists("../pft.vasp"):
    pft = vasp.read_vasp("../pft.vasp")
    files = os.listdir()
    for file in files:
        if re.match(r"0\d", file) and os.path.exists(file+"/POSCAR"):
            atoms = vasp.read_vasp(file+"/POSCAR")
            atoms.set_cell(pft.cell)
            vasp.write_vasp(file+"/POSCAR", atoms, direct=True)
else:
    print("Perfect crystal structure file not found!")