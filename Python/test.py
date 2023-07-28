import linecache as lc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D    #这里只是用Axes3D函数，所以只导入了Axes3D

filename = 'Python\m_data.txt'
x_lattice_vetor = [5.1894846975942732, 0.0000000000000000, -0.0206394262052893]
y_lattice_vetor = [0.0000000000000000, 5.2437694475432117, 0.0000000000000000]
z_lattice_vetor = [-0.8803837383861647, 0.0000000000000000, 5.3044071932761900]
xcoordinates = []
ycoordinates = []
zcoordinates = []
with open(filename) as file_object:
    for line in file_object:
        xcoordinate = float(line.split()[3]) * x_lattice_vetor[0] 
        + float(line.split()[4]) * y_lattice_vetor[0] 
        + float(line.split()[5]) * z_lattice_vetor[0]
        xcoordinates.append(xcoordinate)
        ycoordinate = float(line.split()[3]) * x_lattice_vetor[1] 
        + float(line.split()[4]) * y_lattice_vetor[1] 
        + float(line.split()[5]) * z_lattice_vetor[1]
        ycoordinates.append(ycoordinate)
        zcoordinate = float(line.split()[3]) * x_lattice_vetor[2] 
        + float(line.split()[4]) * y_lattice_vetor[2] 
        + float(line.split()[5]) * z_lattice_vetor[2] 
        zcoordinates.append(zcoordinate)

file_name = 'Python\CONTCAR.txt'
with open(filename) as file_object:
    latt_ba_atoms_nums = lc.getline(filename, 9)
    Zr_atoms_num = int(latt_ba_atoms_nums.split()[0])
    O_atoms_num = int(line.split()[1])
    line_Zr = Zr_atoms_num + 9
    for value in range(9, line_Zr):
        Zr_fcoordinate = lc.getline(filename, value)
        Zr_xcoordinate = float(Zr_fcoordinate.split()[0]) * x_lattice_vetor[0] 
        + float(Zr_fcoordinate.split()[1]) * y_lattice_vetor[0] 
        + float(Zr_fcoordinate.split()[2]) * z_lattice_vetor[0]
        Zr_ycoordinate = float(Zr_fcoordinate.split()[0]) * x_lattice_vetor[1] 
        + float(Zr_fcoordinate.split()[1]) * y_lattice_vetor[1] 
        + float(Zr_fcoordinate.split()[2]) * z_lattice_vetor[1]
        Zr_zcoordinate = float(Zr_fcoordinate.split()[0]) * x_lattice_vetor[2] 
        + float(Zr_fcoordinate.split()[1]) * y_lattice_vetor[2] 
        + float(Zr_fcoordinate.split()[2]) * z_lattice_vetor[2] 
        print(coordinate)

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(xcoordinates, ycoordinates, zcoordinates, c='r')
plt.show()
