#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------------
# This script has two function: Firstly it can be used to constraint
# the atoms and assign the atom which you want to release, Secondly it 
# can expand the primitive cell and keep its original element order.
# All function is based on the ase package.
#--------------------------------------------------------------------
import re
from ase.io import vasp
from ase.constraints import FixAtoms

#prompt = "Please input the vasp structure file which you want to constraint: "
#filepath = input(prompt)
filepath = "POSCAR"
atoms = vasp.read_vasp(filepath)
a = atoms.get_chemical_formula()
constraint_atoms = FixAtoms(
    indices=[atom.index for atom in atoms if atom.symbol != 'H'])
atoms.set_constraint(constraint_atoms)
#atoms = atoms*[2, 2, 2]

b = re.findall(r"\d+", a)
c = re.findall(r"[^0-9]+", a)
d = []
for i in range(len(b)):
    d = [c[i]]*int(b[i]) + d
print(a)
""" tags = atoms.get_chemical_symbols()
print(tags)
deco = sorted(enumerate(tags))
indices = [i for tag, i in deco]
print(indices) """
vasp.write_vasp(file="POSCAR_new", atoms=atoms, direct=True)