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

def get_calculated_point():
    files = os.listdir()
    saddle_file_name_list = []
    stable_file_name_list = []
    for file in files:
        if re.match(r"^saddle(_\d+(\.\d+)?)?$", file):
            saddle_file_name_list.append(file)
        elif re.match(r"^stable(_\d+(\.\d+)?)?$", file):
            stable_file_name_list.append(file)
    return stable_file_name_list, saddle_file_name_list

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

def process_freq_calculation(stable_file_name_list, saddle_file_name_list) -> None:
    #img_freq_data = pd.DataFrame()
    #real_freq_data = pd.DataFrame()
    for stable_point in stable_file_name_list:
        jump_probability_list = []
        zero_energy_list = []
        f_list_stable, fi_list_stable  = get_frequency(stable_point+'/OUTCAR')
        #img_freq_data[stable_point] = fi_list_stable
        #real_freq_data[stable_point] = f_list_stable
        stable_freq_number= len(f_list_stable)
        zero_energy_stable = zero_energy(f_list_stable)
        for saddle_point in saddle_file_name_list:
            f_list_saddle, fi_list_saddle = get_frequency(saddle_point+'/OUTCAR')
            #if stable_point == stable_file_name_list[0]:
                #img_freq_data[saddle_point] = fi_list_saddle
                #real_freq_data[saddle_point] = f_list_saddle
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
        reaction_data = pd.DataFrame({"Saddle_point_name": saddle_file_name_list,
                                      "Jump_probability (THz)": jump_probability_list, 
                                      "Zero_energy_difference (eV)": zero_energy_list})
        reaction_data.to_csv(stable_point+'_jump_freq.csv', index=False)
    #img_freq_data.to_csv('img_freq.csv', index=False)
    #real_freq_data.to_csv('real_freq.csv', index=False)

def main() -> None:
    stable_file_name_list, saddle_file_name_list = get_calculated_point()
    process_freq_calculation(stable_file_name_list, saddle_file_name_list)

if __name__ == "__main__":
    main()