#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#----------------------------------------------------------------------
# The script is used to calculate the reaction frequency.
# Before starting to work, it needs the OUTCAR of VASP after finite
# displacement or DFPT calculation. It will automatically exclude the
# imaginary frequency of saddle points and use formula of the TST theory:
# miu_1*mu_2*...miu_n-1*miu_n/miu'_1*miu'_2...miu_n-1 to calculate the 
# reaction frequency.
# New function: Calculate the zero point energy.
#----------------------------------------------------------------------
import os
import re
import numpy as np
import pandas as pd
from subprocess import getstatusoutput
from scipy.constants import physical_constants

def get_saddle_point_number():
    files = os.listdir()
    saddle_number = 0
    stable_number = 0
    for file in files:
        if re.match(r"saddle_\d", file):
            saddle_number = saddle_number + 1
    return saddle_number

def get_frequency(filename):
    """Get frequency from VASP OUTCAR
    """
    input = "gawk '/THz/' " + filename
    get_frequency = getstatusoutput(input)
    if get_frequency[0] == 0:
        f_list = re.findall(r"f  =(.+?)THz", get_frequency[1])
        try:
            f_list = np.array([float(frequency.strip()) for frequency in f_list])
        except ValueError as e:
            raise e
        fi_list = re.findall(r"f/i=(.+?)THz", get_frequency[1])
        try:
            fi_list = np.array([float(frequency.strip()) for frequency in fi_list])
        except ValueError as e:
            raise e
        return f_list, fi_list
    else:
        raise ValueError

def zero_energy(f_list):
    """Calculate the zero energy.
    """
    f_list = f_list * (10**12) * physical_constants["Planck constant in eV s"][0] * 0.5
    total_zero_energy = np.sum(f_list)
    return total_zero_energy

def process_freq_calculation() -> None:
    stable_site_path = input("Please designate the path of stable site: ")
    if not stable_site_path.strip(): 
        stable_site_path = 'stable'
    point_number = get_saddle_point_number()
    jump_probability_list = []
    zero_energy_list = []
    img_freq_list = []
    real_freq_list = []
    f_list_stable, fi_list_stable  = get_frequency(stable_site_path+'/OUTCAR')
    img_freq_list.append(fi_list_stable)
    real_freq_list.append(f_list_stable)
    stable_freq_number= len(f_list_stable)
    zero_energy_stable = zero_energy(f_list_stable)
    for saddle in range(point_number):
        f_list_saddle, fi_list_saddle = get_frequency('saddle_'+str(saddle+1)+'/OUTCAR')
        img_freq_list.append(fi_list_saddle)
        real_freq_list.append(f_list_saddle)
        zero_energy_saddle = zero_energy(f_list_saddle)
        saddle_freq_number = len(f_list_saddle)
        jump_probability = 1
        for i in range(saddle_freq_number):
            jump_probability = jump_probability * f_list_stable[i] / f_list_saddle[i]
        number_diff = stable_freq_number - saddle_freq_number
        if number_diff > 0:
            for j in range(abs(number_diff)):
                jump_probability = jump_probability * f_list_stable[-(j+1)]
        elif number_diff < 0:
            for j in range(abs(number_diff)):
                jump_probability = jump_probability * f_list_saddle[-(j+1)]
        zero_energy_difference = zero_energy_saddle - zero_energy_stable
        zero_energy_list.append(zero_energy_difference)
        jump_probability_list.append(round(jump_probability,5))
    reaction_data = pd.DataFrame({"Jump_probability (THz)": jump_probability_list, 
                                  "Zero_energy_difference (eV)": zero_energy_list})
    reaction_data.to_csv('jump_freq.csv', index=False)
    img_freq_data = pd.DataFrame(img_freq_list).transpose()
    img_freq_data.to_csv('img_freq.csv', index=False)
    real_freq_data = pd.DataFrame(real_freq_list).transpose()
    real_freq_data.to_csv('real_freq.csv', index=False)
    print('Jump&img&real_freq.csv is generated!')

def main() -> None:
    process_freq_calculation()

if __name__ == "__main__":
    main()