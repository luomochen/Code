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
from ase.io import vasp
from scipy import interpolate
from matplotlib import pyplot as plt

filepath = "/root/code/Scripts/Test/path2.1/"
os.chdir(filepath)

def image_distance(position1, position2):
    """ Calculate the distance between each image.
        The distance was defined as $\sum_{i=1}^n(x^\prime_n-x_n)$.
    Args:
        position1 (list[list[]]): high dimension coordiante of image 1.
        position2 (list[list[]]): high dimension coordiante of image 2.

    Returns:
        float: distance.
    """
    differences = np.linalg.norm(position1-position2)
    return differences

files = os.listdir()
image_number = 0
for file in files:
    if re.match(r"0\d", file):
        image_number = image_number + 1
     
for i in range(image_number):
    image = vasp.read_vasp_out("0" + str(i) + "/OUTCAR")
    Cart = image.get_positions()
    if i == 0:
        Cart0 = Cart
    reaction_coordiante = image_distance(Cart, Cart0)
    print(reaction_coordiante)