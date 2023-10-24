# 该脚本用于线性拟合电子占据数和施加的微扰获取U值。
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