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