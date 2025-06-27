#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------
# This scripts is used to plot the free energy difference
# and generate a new flod to store the saddle point for
# frequency calculation. It can be also used to merge the
# path. 
#--------------------------------------------------------
import os
import re
import matplotlib.pyplot as plt

def get_directories(args):
    """If don't designate the folders of images, the script will automatically 
    found the propably images in the work folder.
    """
    if not args:
        dirs = [d for d in os.listdir('.') if os.path.isdir(d) and re.match(r'^\d\d$', d)]
        return sorted(dirs, key=int)
    return args

def parse_outcar():
    """Found the force and energy in ionic steps.
    """
    with open("OUTCAR", "r") as f:
        lines = f.readlines()

    energy_lines = [line for line in lines if 'energy  w' in line]
    force_lines = [line for line in lines if 'FORCES: m' in line]

    return force_lines, energy_lines

def extract_data(forces, energies):
    """Extract the data from the raw strings.
    """
    data = []
    e0 = float(energies[0].split()[6]) if energies else 0.0
    for j, (f_line, e_line) in enumerate(zip(forces, energies)):
        f = f_line.strip().split()
        e = e_line.strip().split()
        f_value = float(f[4])
        e_value = float(e[6])
        delta_e = e_value - e0
        data.append((j, f_value, e_value, delta_e))
    return data

def plot_data(data, idx):
    steps = [d[0] for d in data]
    forces = [d[1] for d in data]
    energies = [d[2] for d in data]
    #delta_energies = [d[3] for d in data]

    fig, ax1 = plt.subplots(figsize=(12, 8), dpi=150)
    ax1.set_xlabel('Ionic step', fontsize=20)
    ax1.set_ylabel('Force (eV/$\mathrm{\AA}$)', color='tab:blue', fontsize=20)
    ax1.plot(steps, forces, 'o-', label='Force', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=18)

    ax2 = ax1.twinx()
    ax2.set_ylabel('Energy (eV)', color='tab:red', fontsize=20)
    ax2.plot(steps, energies, 'o-', label='Energy', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=18)

    ax2.set_title('Force and Energy vs Step', fontsize=20)
    ax2.grid(True)
    fig.tight_layout()
    fig.savefig(f"vaspout{idx}.png")

    plt.close()

def process_directory(directory, idx):
    os.chdir(directory)
    print(f"In {directory}...")

    if not os.path.exists("OUTCAR"):
        raise FileNotFoundError("No OUTCAR file found.")

    forces, energies = parse_outcar()
    print(f"{len(forces)}")
    data = extract_data(forces, energies)
    plot_data(data, idx)

    os.makedirs("../vaspgr", exist_ok=True)
    os.rename(f"vaspout{idx}.png", f"../vaspgr/vaspout{idx}.png")
    os.chdir("..")

def main():
    import sys
    args = sys.argv[1:]
    directories = get_directories(args)

    os.makedirs("vaspgr", exist_ok=True)

    for i in range(1, len(directories) - 1):
        process_directory(directories[i], i)

if __name__ == "__main__":
    main()