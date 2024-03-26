#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#--------------------------------------------------------------------------------------
# 该脚本用于EOS拟合。
#--------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from subprocess import getstatusoutput
from scipy.optimize import curve_fit

Scale = np.arange(0.95, 1.05, 0.01)

LattVec = []
for i in range(3):
    for j in range(3):    
        GetVec = "sed -n '" + str(i+3) + "p' 1.0/POSCAR | gawk '{print $" \
                + str(j+1) + "}'"
        Vec = getstatusoutput(GetVec)
        if Vec[0] == 0:
            LattVec.append(float(Vec[1].rstrip()))
LattMatrix = np.array([LattVec[0:3], LattVec[3:6], LattVec[6:9]])
V = LattMatrix[0, :]@np.cross(LattMatrix[1, :], LattMatrix[2, :])
Volume = V*(Scale**3)
print(Volume)
#Volume = Volume**(-2/3)

Energy = []
for i in range(11):
    GetE = "grep enthalpy " + str(round(Scale[i], 2)) \
            + "/OUTCAR | tail -n 1 | gawk '{print $5}'"
    E = getstatusoutput(GetE)
    if E[0] == 0:
        Energy.append(float(E[1]))
print(Energy)

def BM_3order_EOS(x, a, b, c, d):
    return a * (x**-2) + b * (x**(-4/3)) + c * (x**(-2/3)) + d
def M_EOS(x, a, b, c, d):
    return a + b*c*((1/d/(d-1))*((x/b)**(1-b))+(1/b)*(x/b)-1/(d-1))
""" def r_squared(coeff, x, y):
    residuals = y - func(x, coeff[0], coeff[1], coeff[2], coeff[3])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared """

popt, pcov = curve_fit(BM_3order_EOS, Volume, Energy)
""" Rsquared = r_squared(popt, Volume, Energy) """

EnergyFit = BM_3order_EOS(Volume, popt[0], popt[1], popt[2], popt[3])
print(EnergyFit)
xmin = (-2 * popt[1] + (4 * (popt[1]**2) - 12 * popt[0] * popt[2])**0.5) / (6 * popt[0])
Vreal = xmin**(-3/2)
print(Vreal)
plt.scatter(Volume, Energy, c="blue")
plt.plot(Volume, EnergyFit, c="red")
plt.show()

popt, pcov = curve_fit(M_EOS, Volume, Energy)
""" Rsquared = r_squared(popt, Volume, Energy) """