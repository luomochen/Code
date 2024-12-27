#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#-------------------------------------------------------------
# The script is used to process data after defect calculation.
# It will generate a scatter graph which shows the energy
# difference between each stable sites and a csv file contain
# all processed data.
#-------------------------------------------------------------
import os
import re
import numpy as np
import pandas as pd
from ase.io import vasp
from matplotlib import pyplot as plt

def get_energy_sites():
    # List all the files in the directory.
    # Use the regular expersion to match the filename
    # and sum the file number.
    floders = os.listdir()
    floder_list = []
    for floder in floders:
        if re.match(r"\d", floder):
            floder_list.append(int(floder))
    floder_list.sort()
    coor_en_list = []
    # Use the ase get the total energy and coordinates
    # of the whole sites.
    for i in floder_list:
        coor_en = []
        filepath = str(i) + "/OUTCAR"
        atoms = vasp.read_vasp_out(filepath)
        coor_en.append(i)
        coor_en.append(atoms.get_total_energy())
        coor_en.append(atoms.get_scaled_positions(wrap=True)[-1])
        coor_en_list.append(coor_en)
    # Use energy to sort the list.
    # Calculate the energy difference.
    sorted_list = sorted(coor_en_list, key=lambda x:x[1])
    en_diff_list = [[outcome[0], outcome[1]-sorted_list[0][1], outcome[2]] 
                for outcome in sorted_list]
    return en_diff_list

def data_process(en_diff_list):
    # Plot scatter plots.
    x = [outcome[1] for outcome in en_diff_list]
    y = [outcome[0] for outcome in en_diff_list]
    fig, ax = plt.subplots(figsize=(8, 15), dpi=150)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_xlim(-0.05, 1)
    #ax.set_xticks(np.arange(-0.05, 0.505, 0.05)) 
    ax.set_xlabel(r"$\mathrm{\Delta}$H (eV)", fontsize=20)
    ax.grid(axis="x")
    ax.scatter(x, y, edgecolors='red', facecolors='white', s=150)
    # Lable on each points.
    for i, energy in enumerate(x):
        ax.text(energy, y[i], y[i],
                ha='center', va='center', fontsize=10, color='black')
    fig.savefig("outcome.png")
    # Use csv format to save the data.
    name = ["No", "energy difference", "direct"]
    data = pd.DataFrame(columns=name, data=en_diff_list)
    data.to_csv("outcome.csv", index=False)

def main():
    en_diff_list = get_energy_sites()
    data_process(en_diff_list)

if __name__ == "__main__":
    main()