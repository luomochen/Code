#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------------------------------
# 该脚本用于自动提取数据得到U值，对VASP生成文件的命名方式有特定的格式要求。
# 该脚本用于线性响应U值获取。
#--------------------------------------------------------------------------------------

import subprocess
import numpy as np
from scipy.optimize import curve_fit

# 从OUTCAR中读取电子占据数。
charge_nscf_list = []
charge_scf_list = []
str_V_list = ['0.00', '0.02', '-0.02', '0.04', '-0.04', '0.06', '-0.06']
rm = subprocess.getstatusoutput("rm Outcome.txt")
if rm[0] == 0:
    print("Already clean old Outcome file.")
for v in str_V_list:
    get_DAV_nscf = "tail -n 2 OSZICAR.V=" + v + ".ICHARG=11 | grep DAV | gawk '{print $2}'"
    get_DAV_scf = "tail -n 2 OSZICAR.V=" + v + " | grep DAV | gawk '{print $2}'"
    get_charge_nscf = "grep -A 4 'total charge' OUTCAR.V=" + v + ".ICHARG=11 | tail -n 1 | gawk '{print $5}'"
    get_charge_scf = "grep -A 4 'total charge' OUTCAR.V=" + v + " | tail -n 1 | gawk '{print $5}'"
    # 获取OSZICAR的电子步数并判断是否收敛。如果执行时候出现问题输出的数据均为NaN。
    DAV_nscf = subprocess.getstatusoutput(get_DAV_nscf)
    DAV_scf = subprocess.getstatusoutput(get_DAV_scf)
    
    if DAV_nscf[0] == 0 and int(DAV_nscf[1]) < 150:
        charge_nscf = subprocess.getstatusoutput(get_charge_nscf)
        if charge_nscf[0] == 0:
            charge_nscf_list.append(float(charge_nscf[1]))
        else:
            charge_nscf_list.append('NaN')
            print("Can't get V =" + v + "nscf data!")
    else:
        charge_nscf_list.append('NaN')
        print("The nscf calculation in V =" + v + "is unconvergence!")
    
    if DAV_scf[0] == 0 and int(DAV_scf[1]) < 150:
        charge_scf = subprocess.getstatusoutput(get_charge_scf)
        if charge_scf[0] == 0:
            charge_scf_list.append(float(charge_scf[1]))
        else:
            charge_scf_list.append('NaN')
            print("Can't get V =" + v + "scf data!")
    else:
        charge_scf_list.append('NaN')
        print("The scf calculation in V =" + v + "is unconvergence!")

# 将三个列表合并为一个dataframe，按照微扰大小进行排序。
filename = "Outcome.txt"
V_list = [float(v) for v in str_V_list]
data = np.array(list(zip(V_list, charge_nscf_list, charge_scf_list)))

# 拟合用函数。
def func(x, a):
    # 读取V=0时的电子占据数。
    b = data[0, 2]
    return  a * x + b
# 计算R^2函数。
def r_squared(slope, xdata, ydata):
    residuals = ydata - func(xdata, slope)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

popt, pcov = curve_fit(func, data[:, 0], data[:, 1])
slope_nscf = popt[0]
r_squared_nscf = r_squared(slope_nscf, data[:, 0], data[:, 1])
popt, pcov = curve_fit(func, data[:, 0], data[:, 2])
slope_scf = popt[0]
r_squared_scf = r_squared(slope_scf, data[:, 0], data[:, 2])
# 计算U值，如果拟合度较差将不返回U值。
if r_squared_scf > 0.90 and r_squared_nscf > 0.90:
    U = 1 / slope_scf - 1 / slope_nscf
    U = round(U,2)
else:
    U = "NaN"

data = data[data[:, 0].argsort()]
number = np.size(data, 0)
with open(filename, 'a') as file_object:
    output = "V(eV), N(nscf), N(scf) \n"
    file_object.write(output)
    for i in range(number):
        for j in range(3):
            output = str(round(data[i, j], 2)) + ", "
            file_object.write(output)
        file_object.write("\n")
    output = "U: "+ str(U) + "\n" + "slope_nscf: " + str(round(slope_nscf, 2)) + "\n" + "R^2_nscf: " + str(round(r_squared_nscf, 2)) \
    + "\n" + "slope_scf: " + str(round(slope_scf, 2)) + "\n" + "R^2_scf: " + str(round(r_squared_scf, 2))
    file_object.write(output)
print(output)
