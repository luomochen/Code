#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#---------------------------------------------
# A script for ploting the fat band structure.
#---------------------------------------------
import re
import math
import numpy as np
from subprocess import getstatusoutput
from doped import analysis
from pymatgen.io.vasp import Vasprun
from matplotlib import pyplot as plt

def get_free_energy(OUTCAR_file_path):
    """Get free energy from the defect syetem.
    """
    input = "gawk '/TOTEN/' " + OUTCAR_file_path
    get_free_energy = getstatusoutput(input)
    if get_free_energy[0] == 0:
        free_energy_list = re.findall(r"=(.+?)eV", get_free_energy[1])
        try:
            free_energy_list = [float(free_energy.strip()) 
                                for free_energy in free_energy_list]
        except ValueError:
            free_energy_list = ['NaN']
        return free_energy_list
    else:
        return ['NaN']

def charge_correction(defect_file_path, stastic_dielec, chemical_potential, defect_charge_state):
    """Get energy correction using eFNV method.
    """
    defect_entry = analysis.defect_entry_from_paths(defect_path=defect_file_path,
                                                    charge_state=defect_charge_state, 
                                                    bulk_path="pft", 
                                                    dielectric=stastic_dielec)
    print(f"Charge: {defect_entry.charge_state} at site: {defect_entry.defect_supercell_site.frac_coords}")
    correction, error = defect_entry.get_kumagai_correction(return_correction_error=True)
    print(f"error range: {error}")
    correction = defect_entry.corrections['kumagai_charge_correction']
    defect_free_energy = get_free_energy(defect_file_path)
    pft_free_energy = get_free_energy("pft")
    intercept = defect_free_energy - pft_free_energy + chemical_potential + correction
    return intercept, [correction, error]

def transition_level(charge_state_1, charge_state_2, intercept_1, intercept_2):
    trans_level = (intercept_1 - intercept_2) / (charge_state_2 - charge_state_1)
    trans_energy = trans_level * charge_state_1 + intercept_1
    return trans_level, trans_energy

def piecewise_function(defect_founction_list):
    transition_data_list = []
    for defect in pairs:
        trans_level, trans_energy = transition_level(pair[0][0], pair[0][1], pair[1][0], pair[1][1])
        transition_data_list.append([trans_level, trans_energy])
    

def plot_formation_enthalpy_fermi_level_diagram(defect_data, band_structure):
    level_span = band_structure[2]+0.5
    enthalpy_span = np.linspace(defect_data[0]-0.5, defect_data[1]+0.5, 120)
    E_fermi = np.arange(0, band_structure[2]+0.5, 0.1)
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