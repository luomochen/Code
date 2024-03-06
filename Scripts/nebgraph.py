#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------
# This scripts is used to plot the free energy difference
# and generate a new flod to store the saddle point for
# frequency calculation.
#--------------------------------------------------------
import os
import re
import numpy as np
import pandas as pd
from ase.io import vasp
from scipy import interpolate
from matplotlib import pyplot as plt

def image_distance(position1, position2):
    """ Calculate the distance between each image.
        The distance was defined as $\sum_{i=1}^n(x^\prime_n-x_n)$.
    Args:
        position1 (list[list[]]): high dimension coordiante of image 1.
        position2 (_type_): high dimension coordiante of image 2.

    Returns:
        float: distance.
    """
    differences = np.linalg.norm(position1-position2)
    return differences

print("--------------------------------------------------------------------------------")
paths = list(input(
    "If multiple diffusion paths need to be merged, specify their numbers: ").split(" "))
print("--------------------------------------------------------------------------------")

# List all the files in the directory.
# Use the regular expersion to match the filename
# and find the image created in the path.
free_energy_list = []
dist_list = []
if paths == [""]:
    files = os.listdir()
    image_number = 0
    for file in files:
        if re.match(r"0\d", file):
            image_number = image_number + 1

    for i in range(image_number):
        # Use ase to read the OUTCAR.
        image = vasp.read_vasp_out("0" + str(i) + "/OUTCAR")
        free_energy = image.get_total_energy()
        carts = image.get_positions()
        if i == 0:
            dist = 0
        else:
            dist = image_distance(carts, cartsp)
        cartsp = carts
        free_energy_list.append(free_energy)
        dist_list.append(dist)
else:
    for path in paths:
        path_free_energy_list = []
        path_dist_list = []
        files = os.listdir("path"+path)
        image_number = 0
        for file in files:
            if re.match(r"0\d", file):
                image_number = image_number + 1

        for i in range(image_number):
            # Use ase to read the OUTCAR.
            image = vasp.read_vasp_out("path"+path+ "/0" + str(i) + "/OUTCAR")
            free_energy = image.get_total_energy()
            carts = image.get_positions()
            if i == 0:
                dist = 0
            else:
                dist = image_distance(carts, cartsp)
            cartsp = carts
            path_free_energy_list.append(free_energy)
            path_dist_list.append(dist)
        
        if re.search(r".1", path) == None:
            path_free_energy_list.pop(0)
            path_dist_list.pop(0)
        
        free_energy_list.extend(path_free_energy_list)
        dist_list.extend(path_dist_list)
# Compare the initial state and final state 
# to find the lower one as the energy baseline.
if free_energy_list[0] > free_energy_list[-1]:
    free_energy_diff_list = [free_energy - free_energy_list[-1] 
                             for free_energy in free_energy_list]
else:
    free_energy_diff_list = [free_energy - free_energy_list[0] 
                             for free_energy in free_energy_list]
# When calculating reaction coordinates, 
# you need to calculate the distance 
# between the previous image and the next image.
# and then add them together.
for i in range(len(dist_list)):
    if i == 0:
        dist_list[i] = dist_list[i]
    else:
        dist_list[i] = dist_list[i] + dist_list[i-1]   
data = pd.DataFrame({"Reaction coordinate": dist_list, 
                     "Energy difference": free_energy_diff_list})
data.to_csv("neboutcome.csv", index=False)
print(data)

# Interpolate to get smooth curve.
model = interpolate.interp1d(
    dist_list, free_energy_diff_list, kind="quadratic")
xs = np.linspace(0, dist_list[-1], 500)
ys = model(xs)

# Plot the free energy - reaction coordianate diagram.
plt.xlabel("Reaction coordiante ($\AA$)")
plt.ylabel("Energy difference (eV)")
plt.grid(axis="y")
plt.scatter(dist_list, free_energy_diff_list, c="red")
plt.plot(xs, ys, c="black")
plt.show()
plt.savefig("neboutcome.png")