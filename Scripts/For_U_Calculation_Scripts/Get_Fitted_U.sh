#!/bin/bash
#--------------------------------------------------------------------------------------
# 该脚本用于自动提取数据得到U值，对VASP生成文件有特定的格式要求。
# 该脚本用于d电子线性响应U值获取。
# 该脚本是同时使用python和Shell的混合脚本。
# 请阅读注释并根据实际需求或选择命令。
#--------------------------------------------------------------------------------------

rm Outcome.txt
# 获取所有自洽和非自洽计算得到的f电子数并生成数据文件。注意，如果不收敛将不读取f电子数。
V_list="0.00 0.01 -0.01 0.02 -0.02 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.1 -0.1"
echo "V(eV),N(nscf),N(scf)" >> Outcome.txt
for v in $V_list
do
    DAV_nscf=$(tail -n 2 OSZICAR.V="$v".ICHARG=11 | grep DAV | gawk '{print $2}')
    DAV_scf=$(tail -n 2 OSZICAR.V="$v" | grep DAV | gawk '{print $2}')
    if [ "$DAV_nscf" -le 150 ]
    then
        # 用于ase扩胞且排序后d电子获取U值：
        #charge_nscf=$(grep -A 100 "total charge" OUTCAR.V="$v".ICHARG=11 | grep "96" | tail -n 1 | gawk '{print $4}')
        # 用于ase扩胞且排序后f电子获取U值：
        #charge_nscf=$(grep -A 100 "total charge" OUTCAR.V="$v".ICHARG=11 | grep "96" | tail -n 1 | gawk '{print $5}')
        # 用于d电子：
        #charge_nscf=$(grep -A 4 "total charge" OUTCAR.V="$v".ICHARG=11 | tail -n 1 | gawk '{print $4}')
        # 用于f电子：
        charge_nscf=$(grep -A 4 "total charge" OUTCAR.V="$v".ICHARG=11 | tail -n 1 | gawk '{print $5}')
    else
        charge_nscf="NaN"
    fi
    if [ "$DAV_scf" -le 150 ] 
    then
        # 用于ase扩胞且排序后d电子获取U值：
        #charge_scf=$(grep -A 100 "total charge" OUTCAR.V="$v" | grep "96" | tail -n 1 | gawk '{print $4}')
        # 用于ase扩胞且排序后f电子获取U值：
        #charge_scf=$(grep -A 100 "total charge" OUTCAR.V="$v" | grep "96" | tail -n 1 | gawk '{print $5}')
        # 用于d电子：
        #charge_scf=$(grep -A 4 "total charge" OUTCAR.V="$v" | tail -n 1 | gawk '{print $4}')
        # 用于f电子：
        charge_scf=$(grep -A 4 "total charge" OUTCAR.V="$v" | tail -n 1 | gawk '{print $5}')
    else
        charge_scf="NaN"
    fi
    echo "$v,$charge_nscf,$charge_scf" >> Outcome.txt
done
# 自动生成用于线性拟合的python脚本
cat >> Get_Fitted_U.py <<!
# 该脚本用于线性拟合电子占据数和施加的微扰获取U值。
# 相较于之前的线性拟合脚本，该脚本拟合的时候能够固定截距。
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# 读取数据并按照微扰大小进行排序。
filename = 'Outcome.txt'
data = pd.read_csv(filename)
data.sort_values(by='V(eV)').to_csv(filename, index=False)
# 将DataFrame格式的数据单列抽取出来转换为一维数组。
nscf = np.array(data['N(nscf)'].tolist())
scf = np.array(data['N(scf)'].tolist())
V = np.array(data['V(eV)'].tolist())
# 分别进行线性拟合

# 拟合用函数。
def func(x, a):
    # 读取V=0时的电子占据数。
    row_V_eq_0 = data[ data['V(eV)'] == 0 ]
    N_V_eq_0 = row_V_eq_0['N(scf)'].tolist()
    b = N_V_eq_0[0]
    return  a * x + b
# 计算R^2函数。
def r_squared(slope, xdata, ydata):
    residuals = ydata - func(xdata, slope)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

popt, pcov = curve_fit(func, V, nscf)
slope_nscf = popt[0]
r_squared_nscf = r_squared(slope_nscf, V, nscf)
popt, pcov = curve_fit(func, V, scf)
slope_scf = popt[0]
r_squared_scf = r_squared(slope_scf, V, scf)
# 计算U值，如果拟合度较差将不返回U值。
if r_squared_scf > 0.90 and r_squared_nscf > 0.90:
    U = 1 / slope_scf - 1 / slope_nscf
    U = round(U,2)
else:
    U = "NaN"
print(U, round(slope_nscf, 2), round(r_squared_nscf, 2), round(slope_scf, 2), round(r_squared_scf, 2))

# 格式化输出，所有值只保留两位小数。
output = "U: "+ str(U) + "\n" + "slope_nscf: " + str(round(slope_nscf, 2)) + "\n" + "R^2_nscf: " + str(round(r_squared_nscf, 2)) \
+ "\n" + "slope_scf: " + str(round(slope_scf, 2)) + "\n" + "R^2_scf: " + str(round(r_squared_scf, 2))
with open(filename, 'a') as file_object:
    file_object.write(output)
!
python Get_Fitted_U.py
rm Get_Fitted_U.py