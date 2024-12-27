#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-----------------------------------------------
# Automatically setting KPOINTS and POTCAR
#-----------------------------------------------
from pymatgen.io.vasp import Poscar, Kpoints, Potcar

def main():
    structure = Poscar.from_file("POSCAR").structure
    symbol = Poscar.from_file("POSCAR").site_symbols
    kpoint_density = [40, 40, 40]  # k点密度, 可以根据需要调整
    kpoints = Kpoints.automatic_density_by_lengths(structure, kpoint_density, force_gamma=True)
    kpoints.write_file("KPOINTS")
    potcar = Potcar(symbol, "PBE")
    potcar.write_file("POTCAR")

if __name__ == "__main__":
    main()