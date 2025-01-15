#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#----------------------------------------------------
# Process the kmc outcome.
# Function:
# 1. Calculate the diffusion coefficient in a crystal
# orientation.
# 2. Calculate the total diffusion coefficient.
#----------------------------------------------------
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from ase.io import vasp
from scipy.optimize import curve_fit
from scipy.constants import physical_constants

import input
def convert_ndarray(data, column_name):
    numeric = pd.to_numeric(data[column_name])
    ndarray = np.array(numeric.tolist())
    return ndarray

def r_squared(slope, intercept, xdata, ydata):
    residuals = ydata - linear_func(xdata, slope, intercept)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def linear_func(x, a, b):
    return  a * x + b

def linear_fit(xdata, ydata):
    """Use Arrhenius equation to fit the data.
    """
    kb = physical_constants['Boltzmann constant in eV/K'][0]
    popt, pcov = curve_fit(linear_func, xdata, ydata)
    barrier = -popt[0]*kb*1000*np.log(10)
    D0 = 10**popt[1]
    return popt[0], popt[1], barrier, D0

def diffusion_coefficient_in_crystal_orientation(T, crystal_orientation, vasp_filename):
    atom = vasp.read_vasp(vasp_filename)
    displacement_list = ['d_x', 'd_y', 'd_z']
    data = pd.read_csv(str(T)+".csv")
    t = convert_ndarray(data, 't')
    displacement_matrix = []
    for displacement_direction in displacement_list:
        displacement = convert_ndarray(data, displacement_direction)
        displacement_matrix.append(displacement)
    displacement_matrix = np.array(displacement_matrix)
    # Project on the designated crystal orientation.
    # transform to unit vector.
    crystal_vector = crystal_orientation@atom.cell
    unit_crystal_vector = crystal_vector/np.linalg.norm(crystal_vector)
    projection = unit_crystal_vector@displacement_matrix
    # Calculate the diffusion coefficient.
    projection_square_sum = np.sum(np.square(projection)*10**-20)
    mean_diffusion_coeff = projection_square_sum / (2*np.sum(t))
    return mean_diffusion_coeff

def diffusion_coefficient(T):
    displacement_list = ['d_x', 'd_y', 'd_z']
    data = pd.read_csv(str(T)+".csv")
    t = convert_ndarray(data, 't')
    displacement_matrix = []
    for displacement_direction in displacement_list:
        displacement = convert_ndarray(data, displacement_direction)
        displacement_matrix.append(displacement)
    displacement_matrix = np.array(displacement_matrix)
    # Calculate the diffusion coefficient in each direction.
    displacement_square_sum = np.sum(np.square(displacement_matrix)*10**-20, axis=1)
    mean_diffusion_coeff = displacement_square_sum / (2*np.sum(t))
    mean_diffusion_coeff = np.append(mean_diffusion_coeff, np.sum(mean_diffusion_coeff)/3)
    return mean_diffusion_coeff

def processed_data_output(slope, intercept, Rsquared, barrier, D0, diffusion_coeffs, reciprocal_T_list) -> None:
    # Plot
    x_data = np.linspace(min(reciprocal_T_list)-0.1, max(reciprocal_T_list)+0.1, 1000)
    y_data = linear_func(x_data, slope, intercept)
    fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
    ax.set_xlabel(r"1000/T ($\mathrm{K^{-1}}$)", fontsize=20)
    ax.set_ylabel("log10(D)", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.scatter(reciprocal_T_list, diffusion_coeffs, edgecolors='red', 
               facecolors='white', s=100, linewidths=2, zorder=2)
    ax.plot(x_data, y_data, c="black", linewidth=3, zorder=1)
    ax.set_title(f"Barrier_fit = {round(barrier, 3)} eV", fontsize=20)
    func_text = f"y = {slope:.2f}x + {intercept:.2f}\n$R^2$ = {Rsquared:.3f}"
    ax.text(np.median(x_data)+0.1, np.median(y_data)+0.1, func_text, fontsize=25, color='red')
    ax.grid(True)
    fig.savefig("linear_fit.png")
    # Output the fitted data.
    with open ("process_kmc.dat", 'w') as output:
        output.write(f"Linear_fit: {slope}*x+{intercept}\n")
        output.write(f"R_squared: {Rsquared}\n")
        output.write(f"Diffusion_barrier: {barrier} eV\n")
        output.write(f"Prefactor: {D0}")
    plt.show()

def main():
    _, _, _, _, _, T_list, _, _, _= input.read_input()
    diffusion_coeffs = []
    for T in T_list:
        coeff = diffusion_coefficient(T)
        diffusion_coeffs.append(np.log10(coeff))
    reciprocal_T_list = 1000 / np.array(T_list)
    diffusion_coeffs = np.array(diffusion_coeffs)
    slope, intercept, barrier, D0 = linear_fit(reciprocal_T_list, diffusion_coeffs[:,3])
    Rsquared = r_squared(slope, intercept, reciprocal_T_list, diffusion_coeffs[:,3])
    processed_data_output(slope, intercept, Rsquared, barrier, D0,
                               diffusion_coeffs[:,3], reciprocal_T_list)
    #crystal_orientation = [1, 0, 0]
    #diffusion_coeffs = []
    #for T in T_list:
    #    coeff = diffusion_coefficient_in_crystal_orientation(T, crystal_orientation, 'sites.vasp')
    #    diffusion_coeffs.append(np.log10(coeff))
    #reciprocal_T_list = 1000 / np.array(T_list)
    #diffusion_coeffs = np.array(diffusion_coeffs)
    #slope, intercept, barrier, D0 = linear_fit(reciprocal_T_list, diffusion_coeffs)
    #Rsquared = r_squared(slope, intercept, reciprocal_T_list, diffusion_coeffs)
    #processed_data_output(slope, intercept, Rsquared, barrier, D0,
    #                           diffusion_coeffs, reciprocal_T_list)
    
if __name__ == "__main__":
    main()