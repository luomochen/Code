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
import pandas as pd
from ase.io import vasp
from matplotlib import pyplot as plt

# List all the files in the directory.
# Use the regular expersion to match the filename
# and sum the file number.
files = os.listdir()
file_number = 0
for file in files:
    if re.match(r"\d", file):
        file_number = file_number + 1
coor_en_list = []

# Use the ase get the total energy and coordinates
# of the whole sites.
for i in range(file_number):
    coor_en = []
    filepath = str(i+1) + "/OUTCAR"
    ZrOH = vasp.read_vasp_out(filepath)
    coor_en.append(i+1)
    coor_en.append(ZrOH.get_total_energy())
    coor_en.append(ZrOH.get_scaled_positions(wrap=True)[-1])
    coor_en_list.append(coor_en)
# Use energy to sort the list.
# Calculate the energy difference.
sorted_list = sorted(coor_en_list, key=lambda x:x[1])
en_diff_list = [[outcome[0], outcome[1]-sorted_list[0][1], outcome[2]] 
               for outcome in sorted_list]

# Plot scatter plots.
x = [outcome[0] for outcome in en_diff_list]
y = [outcome[1] for outcome in en_diff_list]
plt.xlabel("Search sites")
plt.ylabel("Energy difference (eV)")
plt.scatter(x, y)
for i, label in enumerate(x):
    plt.annotate(label, (x[i], y[i]))
plt.savefig("outcome.png")

# Use csv format to save the data.
name = ["No", "energy difference", "direct"]
data = pd.DataFrame(columns=name, data=en_diff_list)
data.to_csv("outcome.csv", index=False)