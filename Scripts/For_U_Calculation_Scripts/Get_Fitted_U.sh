#!/bin/bash
#--------------------------------------------------------------------------------------
# 该脚本用于自动提取数据得到U值，对VASP生成文件有特定的格式要求。
# 该脚本是同时使用python和Shell的混合脚本。
#--------------------------------------------------------------------------------------

# 获取所有自洽和非自洽计算得到的f电子数并生成数据文件。注意，如果不收敛将不读取f电子数。
V_list="0.001 -0.001 0.002 -0.002 0.003 -0.003 0.004 -0.004 0.006 -0.006 0.008 -0.008 \
0.01 -0.01 0.02 -0.02 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.1 -0.1 0.2 -0.2 0.4 -0.4 0.6 \
-0.6 0.8 -0.8"
echo "V(eV),N(nscf),N(scf)" >> Outcome.txt
for v in $V_list
do
    DAV_nscf=$(tail -n 2 OSZICAR.V="$v".ICHARG=11 | grep DAV | gawk '{print $2}')
    DAV_scf=$(tail -n 2 OSZICAR.V="$v" | grep DAV | gawk '{print $2}')
    if [ "$DAV_nscf" -le 150 ]
    then
        charge_nscf=$(grep -A 4 "total charge" OUTCAR.V="$v".ICHARG=11 | tail -n 1 | gawk '{print $5}')
    else
        charge_nscf="NaN"
    fi
    if [ "$DAV_scf" -le 150 ] 
    then
        charge_scf=$(grep -A 4 "total charge" OUTCAR.V="$v" | tail -n 1 | gawk '{print $5}')
    else
        charge_scf="NaN"
    fi
    echo "$v,$charge_nscf,$charge_scf" >> Outcome.txt
done
# 自动生成用于线性拟合的python脚本
cat >> Get_Fitted_U.py <<!
import pandas as pd
import numpy as np

# 读取数据并按照微扰大小进行排序。
filename = 'Outcome.txt'
data = pd.read_csv(filename)
data.sort_values(by='V(eV)').to_csv(filename, index=False)
# 将DataFrame格式的数据单列抽取出来转换为一维数组。
nscf = np.array(data['N(nscf)'].tolist())
scf = np.array(data['N(scf)'].tolist())
V = np.array(data['V(eV)'].tolist())
# 分别进行线性拟合。
coeff_1 = np.polyfit(V, nscf, 1)
coeff_2 = np.polyfit(V, scf, 1)
# 分别计算r_squared。
correlation_matrix_1 = np.corrcoef(V, nscf)
correlation_V_nscf = correlation_matrix_1[0,1]
r_squared_1 = correlation_V_nscf**2
correlation_matrix_2 = np.corrcoef(V, scf)
correlation_V_nscf = correlation_matrix_2[0,1]
r_squared_2 = correlation_V_nscf**2
# 计算U值，如果拟合度较差将不返回U值。
if r_squared_1 > 0.90 and r_squared_2 > 0.90:
    U = 1 / coeff_2[0] - 1 / coeff_1[0]
    U = round(U,2)
else:
    U = "NaN"
# 格式化输出，所有值只保留两位小数。
output = "U: "+ str(U) + "\n" + "slope_nscf: " + str(round(coeff_1[0], 2)) + "\n" + "R^2_nscf: " + str(round(r_squared_1, 2)) \
+ "\n" + "slope_scf: " + str(round(coeff_2[0], 2)) + "\n" + "R^2_scf: " + str(round(r_squared_2, 2))
with open(filename, 'a') as file_object:
    file_object.write(output)
print(U,coeff_1,r_squared_1,coeff_2,r_squared_2)
!
python Get_Fitted_U.py
rm Get_Fitted_U.py