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
import numpy as np
import pandas as pd
from ase.io import vasp
from scipy import interpolate
from matplotlib import pyplot as plt

def image_distance(directs1, directs2, latt_vec_matrix):
    """ Calculate the distance between each image using pbc.
        The distance was defined as $\sum_{i=1}^n(x^\prime_n-x_n)$.
    Args:
        position1 (list[list[]]): high dimension coordiante of image 1.
        position2 (_type_): high dimension coordiante of image 2.

    Returns:
        float: distance.
    """
    delta_dirs = directs1 - directs2
    size = delta_dirs.shape
    for i in range(size[0]):
        for j in range(size[1]):
            if delta_dirs[i, j] <= -0.5:
                delta_dirs[i, j] += 1.0
            elif delta_dirs[i, j] > 0.5:
                delta_dirs[i, j] -= 1.0
    differences = np.linalg.norm(delta_dirs@latt_vec_matrix)
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
        latt_vec_matrix = image.cell
        free_energy = image.get_total_energy()
        directs = image.get_scaled_positions()
        if i == 0:
            dist = 0
            latt_vec_matrix_p = latt_vec_matrix
        else:
            dist = image_distance(directs, directsp, latt_vec_matrix_p)
            if np.linalg.norm(latt_vec_matrix_p-latt_vec_matrix) > 1E-6:
                print(f"Image{i} lattice vector matrix is changed!!!")
        directsp = directs
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
            latt_vec_matrix = image.cell
            free_energy = image.get_total_energy()
            directs = image.get_scaled_positions()
            if i == 0:
                dist = 0
                latt_vec_matrix_p = latt_vec_matrix
            else:
                dist = image_distance(directs, directsp, latt_vec_matrix_p)
                if np.linalg.norm(latt_vec_matrix_p-latt_vec_matrix) > 1E-6:
                    print(f"Image{i} lattice vector matrix is changed!!!")
            directsp = directs
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
plt.savefig("neboutcome.png")
plt.show()