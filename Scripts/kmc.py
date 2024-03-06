#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#----------------------------------------------------------------
# 该脚本用于进行氢原子间隙扩散的Mote Carlo模拟
#----------------------------------------------------------------

import math
import numpy as np
import matplotlib.pyplot as plt

lattice_vector = np.array([[4.194, 0, 0],
                           [0, 4.194, 0],
                           [0, 0, 4.194]])
c = np.array([0.25, 0.25, 0.25])
i1 = [0.5, 0., 0.]
i2 = [-0.5, 0., 0.]
i3 = [0., 0.5, 0.]
i4 = [0.5, -0.5, 0.]
i5 = [0., 0., 0.5]
i6 = [0., 0., -0.5]

def jump_rate(T):
    prefactor = 1.22E12
    Energy_barrier = 0.5*1.60217662E-19
    k_B = 1.380649E-23
    rate = prefactor*math.exp(-Energy_barrier/k_B/T)
    return rate
def KMC(Time):
    t = 0
    step = []
    while t < Time:
        rate = jump_rate(500)
        k_sum = np.arange(rate, rate*6, rate)
        r_1 = np.random.random()
        r_2 = np.random.random()
        t = t - math.log(r_1) / k_sum[5]
        sampling = r_2 * k_sum[5]
        event = np.argwhere(k_sum >= sampling)[0, 0]
        if event == 0:
            step.append(i1)
        if event == 1:
            step.append(i2)
        if event == 2:
            step.append(i3)
        if event == 3:
            step.append(i4)
        if event == 4:
            step.append(i5)
        if event == 5:
            step.append(i6)
    return step

i = 0
displacesqare_sum = 0
while i < 100:    
    step = KMC(1E-2)
    displace_sum = np.zeros(3)
    for displace in step:
        displace_sum = displace_sum + np.array(displace)
    dispalce_Car = displace_sum.dot(lattice_vector)
    displacesqare = dispalce_Car.dot(dispalce_Car.T) * 1E-20
    displacesqare_sum = displacesqare_sum + displacesqare
    i = i + 1
D = displacesqare_sum / 100 / 2 / 1E-2
print(D)
print(len(step))
""" steps = KMC(1E-1)
x = np.zeros(len(steps))
y = np.zeros(len(steps))
z = np.zeros(len(steps))
for i in range(len(steps)):
    c = c + np.array(steps[i])
    x[i] = c[0]
    y[i] = c[1]
    z[i] = c[2]
x = x[0:10000]
y = y[0:10000]
z = z[0:10000]
end_x = x[9999]
end_y = y[9999]
end_z = z[9999]

fig = plt.figure(figsize=(8, 8), dpi=100)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='purple', alpha=0.5, s=0.5)
ax.plot(x, y, z, c='blue', alpha=0.75, lw=0.5)
ax.plot(0, 0, 0, c='red', marker='D')
ax.plot(end_x, end_y, end_z, c='black', marker='o')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.title('3D Random Walk')
plt.show()   """