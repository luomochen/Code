#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------------------------
# 该脚本用于随晶格常数变化自动选取合适的N_1,N_2和N_3并生成KPOINTS。
#---------------------------------------------------------------

from subprocess import getstatusoutput

filename = "KPOINTS"
print("Tips: If you want to use single Gamma point in calculation, please input 'G' or 'g'!!!")
print("Otherwise you can use the empirical equation Ka = 30-40 to set KPOINTS, which 'a' is the lattice constant.")
print("Here we only use gamma center.")
Ka = input("Please input the Ka = ")
rm = getstatusoutput("rm KPOINTS")
if rm[0] == 0:
    print("Old KPOINTS is removed! New KPOINTS will be generated...")
if Ka == "G" or Ka == "g":
    with open(filename, 'a') as file_object:
        output = "Gamma_only\n 0\nGamma\n1 1 1\n0 0 0"
        file_object.write(output)
else:
    lattice_vec = []
    N_list = []
    for i in range(3):
        for j in range(3):    
            get_vec = "sed -n '" + str(i+3) + "p' POSCAR | gawk '{print $" + str(j+1) + "}'"
            vec = getstatusoutput(get_vec)
            if vec[0] == 0:
                lattice_vec.append(float(vec[1].rstrip()))
            else:
                lattice_vec.append("NaN")
    try:
        for i in range(3, 10, 3):
            N = int(Ka) // (lattice_vec[i-3]**2 + lattice_vec[i-2]**2 + lattice_vec[i-1]**2)**0.5
            N_list.append(int(N))
    except ValueError:
        print("Some error happen, the lattice vector is not read")
    with open(filename, 'a') as file_object:
        output = "KPOINTS\n 0\nGamma\n" + str(N_list[0]) + " " + str(N_list[1]) + " " + str(N_list[2]) + "\n0 0 0"
        file_object.write(output)