#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to mark all equivalent atoms and regenerate 
# POSCAR in the order of equivalent atoms.
#-------------------------------------------------------------
from ase.io import vasp
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

poscar = Poscar.from_file('POSCAR')
sga = SpacegroupAnalyzer(poscar.structure)
indices = sga.get_symmetrized_structure().equivalent_indices
print(indices)