#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------
#该脚本用于根据POSCAR添加FM和AFM磁序。
#------------------------------------------------------------

import numpy as np
from subprocess import getstatusoutput

# 获取晶格矢量矩阵。
lattice_vec = []
for i in range(3):
    for j in range(3):    
        get_vec = "sed -n '" + str(i+3) + "p' POSCAR | gawk '{print $" + str(j+1) + "}'"
        vec = getstatusoutput(get_vec)
        if vec[0] == 0:
            lattice_vec.append(float(vec[1].rstrip()))
        else:
            print("The format of POSCAR has some problem! Can not read lattice vector!")
lattice_vec_matric = np.array([lattice_vec[0:3], lattice_vec[3:6], lattice_vec[6:9]])
# 求晶格矢量矩阵的逆矩阵。
lattice_vec_matric_inv = np.linalg.inv(lattice_vec_matric)
# 计算内积。
def inner_product(direct1, direct2):
    # 计算度规矩阵(Metric Matrix)，Inner_product = (a b c)*(i j k)T*(i j k)*(a b c)T，MM = (i j k)T*(i j k)
    Metric_Matrix = lattice_vec_matric.dot(lattice_vec_matric.T)
    product = direct1.dot(Metric_Matrix).dot(direct2.T)
    return product

# 获取原子的分数坐标，并按照不同元素分别存放在不同的列表中。
# 确定各个种类原子数和总原子数。
get_U_number = getstatusoutput("sed -n '7p' POSCAR | gawk '{print $1}'")
get_O_number = getstatusoutput("sed -n '7p' POSCAR | gawk '{print $2}'")
if get_U_number[0] == 0 and get_O_number[0] == 0:
    U_number = int(get_U_number[1].rstrip())
    O_number = int(get_O_number[1].rstrip())
    atom_number = U_number + O_number
else:
    print("The format of POSCAR has some problem! Can not read atom number!")
# 读取每种原子分数坐标并储存在列表中。
U_coordi_list = []
O_coordi_list = []
for i in range(1, atom_number+1):
    U_coordinate=[]
    O_coordinate=[]
    line_number = i + 8
    if i <= U_number:
        for j in range(1, 4):
            get_coordinate = getstatusoutput("sed -n " + str(line_number) + "p POSCAR | gawk '{print $"+ str(j) +"}'")
            if get_coordinate[0] == 0:
                U_coordinate.append(float(get_coordinate[1].rstrip()))
            else:
                print("The format of POSCAR has some problem! Can not read U direct coordinate!")
        U_coordi_list.append(U_coordinate)
    else:
        for j in range(1, 4):
            get_coordinate = getstatusoutput("sed -n " + str(line_number) + "p POSCAR | gawk '{print $"+ str(j) +"}'")
            if get_coordinate[0] == 0:
                O_coordinate.append(float(get_coordinate[1].rstrip()))
            else:
                print("The format of POSCAR has some problem! Can not read O direct coordinate!")  
        O_coordi_list.append(O_coordinate)

# 确定沿着轴a同一平面上的原子。
""" U_coordi_list = sorted(U_coordi_list, key=(lambda x:[x[2], x[1]]))
print(U_coordi_list, len(U_coordi_list)) """
accuracy = 1E-4

projection_list = []
a = np.array([0, 0, 1])

for i in range(len(U_coordi_list)):
    U_coordinate = np.array(U_coordi_list[i])
    d = inner_product(a, U_coordinate)
    projection_list.append([i, d])
proj_select_a = []
for i in range(len(projection_list)):
    projection_select_plane = []
    for U_coordinate in projection_list:
        if abs(projection_list[i][1] - U_coordinate[1]) < accuracy:
            projection_select_plane.append(U_coordinate)
    if projection_select_plane not in proj_select_a:
        proj_select_a.append(projection_select_plane)
# 确定沿着轴b同一平面的原子。
projection_list = []
b = np.array([0, 1, 0])

for i in range(len(U_coordi_list)):
    U_coordinate = np.array(U_coordi_list[i])
    d = inner_product(a, U_coordinate)
    projection_list.append([i, d])
proj_select_b = []
for i in range(len(projection_list)):
    projection_select_plane = []
    for U_coordinate in projection_list:
        if abs(projection_list[i][1] - U_coordinate[1]) < accuracy:
            projection_select_plane.append(U_coordinate)
    if projection_select_plane not in proj_select_b:
        proj_select_b.append(projection_select_plane)
# 确定沿着轴c同一平面的原子。
projection_list = []
b = np.array([1, 0, 0])

for i in range(len(U_coordi_list)):
    U_coordinate = np.array(U_coordi_list[i])
    d = inner_product(a, U_coordinate)
    projection_list.append([i, d])
proj_select_c = []
for i in range(len(projection_list)):
    projection_select_plane = []
    for U_coordinate in projection_list:
        if abs(projection_list[i][1] - U_coordinate[1]) < accuracy:
            projection_select_plane.append(U_coordinate)
    if projection_select_plane not in proj_select_c:
        proj_select_c.append(projection_select_plane)
print(proj_select_a, "\n", proj_select_b, "\n", proj_select_c)