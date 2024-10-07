#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for ploting the fat band structure.
#---------------------------------------------
import math
import numpy as np
from doped import analysis
from pymatgen.io.vasp import Vasprun
from matplotlib import pyplot as plt

def charge_correction(filepath, stastic_dielec):
    defect_entry = analysis.defect_entry_from_paths(defect_path=filepath, 
                                                    bulk_path="pft", 
                                                    dielectric=stastic_dielec)
    print(f"Charge: {defect_entry.charge_state} at site: {defect_entry.defect_supercell_site.frac_coords}")
    correction, error = defect_entry.get_kumagai_correction(return_correction_error=True)
    print(f"error range: {error}")
    correction = defect_entry.corrections['kumagai_charge_correction']
    return correction, error

def defect_data_porcess(stastic_dielec, band_structure, enthalpy_list):
    correction_Hip, error_Hip = charge_correction('Hi+', stastic_dielec)
    intercept_Hip = enthalpy_list[0] - enthalpy_list[3] - enthalpy_list[4] + band_structure[0] + correction_Hip
    correction_Him, error_Him = charge_correction('Hi-', stastic_dielec)
    intercept_Him = enthalpy_list[1] - enthalpy_list[3] - enthalpy_list[4] - band_structure[0] + correction_Him
    intercept_Hin = enthalpy_list[2] - enthalpy_list[3] - enthalpy_list[4]
    Formation_Enthalpy_Hi = intercept_Hin
    pm_level = (intercept_Him - intercept_Hip) / 2
    pn_level = Formation_Enthalpy_Hi - intercept_Hip
    nm_level = -Formation_Enthalpy_Hi + intercept_Him
    pm_Formation_enthalpy = pm_level + intercept_Hip
    U = 2 * (pm_Formation_enthalpy - Formation_Enthalpy_Hi)
    processed_data = [intercept_Hip, intercept_Him, intercept_Hin, pm_level, pn_level, nm_level, U]
    output = f"Formation Enthalpy of Hi+: E_f + {intercept_Hip}\
        \nFormation Enthalpy of Hi-: -E_f + {intercept_Him}\
        \nFormation Enthalpy of Hi: {Formation_Enthalpy_Hi}\
        \n+/- CTL Fermi level: {pm_level}\
        \n+/0 CTL Fermi level: {pn_level}\
        \n0/- CTL Fermi level: {nm_level}\
        \nCTL Formation enthalpy: {pm_Formation_enthalpy}\
        \nU: {U}"
    print(output)
    with open('defect.dat', 'w') as file:
        file.write(output)
    return processed_data

def plot_formation_enthalpy_fermi_level_diagram(defect_data, band_structure):
    level_span = max(defect_data[3:6])+0.5
    enthalpy_span = np.linspace(defect_data[0]-0.5, defect_data[1]+0.5, 120)
    E_fermi = np.arange(0, level_span+0.5, 0.1)
    E_formation_Hip = E_fermi + defect_data[0]
    E_formation_Him = -E_fermi + defect_data[1]
    E_formation_Hin = 0 * E_fermi + defect_data[2]
    E_formation = [E_formation_Hip, E_formation_Him, E_formation_Hin]
    gap = np.ones(enthalpy_span.size)*band_structure[2]
    # Plot
    colour = ['blue', 'red', 'black']
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    ax.set_xlabel("Fermi level (eV)", fontsize=22)
    ax.set_ylabel("Formation Enthalpy (eV)", fontsize=22)
    ax.set_xlim((0, level_span))
    ax.set_ylim((defect_data[0]-0.5, defect_data[1]+0.5))
    ax.set_xticks(np.arange(0, math.floor(level_span)+1, 1))
    ax.set_xticklabels(np.arange(0, math.floor(level_span)+1, 1), 
                       fontsize=18)
    ax.set_yticks(np.arange(math.ceil(defect_data[0]-0.5), 
                            math.floor(defect_data[1]+0.5)+1, 2))
    ax.set_yticklabels(np.arange(math.ceil(defect_data[0]-0.5), 
                                 math.floor(defect_data[1]+0.5)+1, 2), 
                       fontsize=18)
    for i in range(len(colour)):
        ax.plot(E_fermi, E_formation[i], c=colour[i], linewidth=4)
    ax.plot(gap, enthalpy_span, linewidth=2, color="black", linestyle="--")
    ax.legend(labels=[r'$\mathrm{H}_\mathrm{i}^\mathrm{{\cdot}}$', 
                    r'$\mathrm{H}_\mathrm{i}^\mathrm{{\prime}}$', 
                    r'$\mathrm{H}_\mathrm{i}$'], 
            loc=1,
            framealpha=1,
            fontsize=22)
    fig.savefig('Formation_enthalpy-Fermi_level')

def main():
    VBM = 8.8992
    CBM = 15.1630
    BAND_GAP = 6.2637
    band_structure = [VBM, CBM, BAND_GAP]
    enthalpy_Hip = -637.76984459
    enthalpy_Him = -611.10019699
    enthalpy_Hin = -623.56955505
    enthalpy_pft = -625.45286166
    chemical_potential = -40.49991932 / 16
    enthalpy_list = [enthalpy_Hip, enthalpy_Him, enthalpy_Hin, enthalpy_pft, chemical_potential]
    high_freq_dielec = np.array([[3.302206, 0.000000, 0.000000],
                                [0.000000, 3.302206, 0.000000],
                                [0.000000, 0.000000, 3.486863]])
    ionic_dielec = np.array([[6.552192, 0.000000, 0.000000],
                            [0.000000, 6.552192, 0.000000],
                            [0.000000, 0.000000, 5.038697]])
    stastic_dielec = high_freq_dielec + ionic_dielec
    defect_data = defect_data_porcess(stastic_dielec, band_structure, enthalpy_list)
    plot_formation_enthalpy_fermi_level_diagram(defect_data, band_structure)

if __name__ == "__main__":
    main()