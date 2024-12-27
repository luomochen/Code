#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# This script is used to obtain the diffusion path and generate
# structure file.
#-------------------------------------------------------------

def parse_args():
    # 使用argparse获取参数
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--diffusionpath', type=str, choices=['r', 'w'], default="r",
                   help="Get diffusion in perfect structure file.")
    p.add_argument('-p', '--perfect', type=str, default="../pft.vasp",
                   help="The perfect Vasp structure file for diffusion path generated.")
    p.add_argument('-i', '--interstitial', nargs='+', type=str, default='H',
                   help="The intersitial atom in the host structure.")
    p.add_argument('-v', '--view', type=str, default='./diffpath.vasp',
                   help="View structure file.")
    return p

def get_perfect_structure(pft_path):
    import os
    from ase import Atoms
    from ase.io import vasp
    try:
        if os.path.exists(pft_path):
            pft = vasp.read_vasp(pft_path)
        else:
            print("Perfect crystal structure file not found!")
            print("Initial state lattice constants will be used to generate a box.")
            atoms = vasp.read_vasp("00/POSCAR")
            pft = Atoms(pbc=True, cell=atoms.cell)
        return pft
    except FileNotFoundError:
            print("Perfect crystal structure file not found!")

def get_image_number():
    import re
    import os
    files = os.listdir()
    image_number = 0
    for file in files:
        if re.match(r"0\d", file):
            image_number = image_number + 1
    return image_number

def get_diffusion_path(readmode, image_number, interstitial_atom, pft):
    import os
    from ase.io import vasp
    try:
        if readmode == 'r':
            for i in range(image_number):
                # Use ase to read the OUTCAR.
                if os.path.exists("0" + str(i) + "/OUTCAR"):
                    image = vasp.read_vasp_out("0" + str(i) + "/OUTCAR")
                elif os.path.exists("0" + str(i) + "/CONTCAR"):
                    image = vasp.read_vasp("0" + str(i) + "/CONTCAR")
                else:
                    image = vasp.read_vasp("0" + str(i) + "/POSCAR")
                for ind, iatom in enumerate(image):
                    if iatom.symbol == interstitial_atom:
                        pft.extend(iatom)
            vasp.write_vasp("./diffpath.vasp", pft, direct=True)
            return 0
        else:
            if os.path.exists("./diffpath.vasp"):
                diffpath = vasp.read_vasp("./diffpath.vasp")
                return diffpath
    except FileNotFoundError:
        print("Image file not found!")

def write_image(image_number, diffpath, pft, interstitial_atom):
    import os
    import copy
    from ase import Atoms
    from ase.io import vasp
    atoms = Atoms(pbc=True, cell=pft.cell)
    for ind, iatom in enumerate(diffpath):
        if iatom.symbol == interstitial_atom:
            atoms.extend(iatom)
    try:
        for i in range(image_number):
        # Use ase to read the OUTCAR.
            if i != 0 and i != image_number-1:
                image = vasp.read_vasp("0" + str(i) + "/POSCAR")
                image_copy = copy.deepcopy(image)
                for ind, iatom in enumerate(image_copy):
                    if iatom.symbol == interstitial_atom:
                        del image_copy[ind]
                image_copy.extend(atoms[i])
                os.popen("cp 0"+str(i)+"/POSCAR 0"+str(i)+"/POSCAR.old")
                vasp.write_vasp("0" + str(i) + "/POSCAR", image_copy, direct=True)
                print(f"New image {i} is generated!")
    except FileNotFoundError:
        print("Image file not found!")

def view_structure_file(atoms):
    from ase.visualize import view
    view(atoms)
 
def main():
    p = parse_args()
    args = p.parse_args()
    try:
        pft = get_perfect_structure(args.perfect)
        image_number = get_image_number()
        path = get_diffusion_path(args.diffusionpath, image_number, args.interstitial, pft)
        if path == 0:
            print("diffpath.vasp is generated.")
        else:
            write_image(image_number, path, pft, args.interstitial)
    except:
        p.print_help()
           
if __name__ == "__main__":
    main()