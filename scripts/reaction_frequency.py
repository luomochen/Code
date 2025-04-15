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
# New function2: Calculate the correction coefficient of Quantum harmonic 
# oscillator.
# New function3: Calculate the  crossover temperature, Tc, Below this 
# crossover temperature, tunneling can make a signiﬁcant contribution 
# to the total rate, while above this temperature tunneling can be
# neglected.
#----------------------------------------------------------------------
import os
import re
import yaml
import numpy as np
from scipy.constants import physical_constants

def get_calculated_point():
    """
    Identify folders matching 'stable', 'stable_#', 'saddle', or 'saddle_#'.
    These are considered calculation directories.
    """
    files = os.listdir()
    saddle_list = [f for f in files if re.match(r"^saddle(_\d+(\.\d+)?)?$", f)]
    stable_list = [f for f in files if re.match(r"^stable(_\d+(\.\d+)?)?$", f)]
    return stable_list, saddle_list

def get_frequency(filename, max_lines=100000):
    """
    Extract real and imaginary frequencies from OUTCAR (in THz).
    Reads only the last `max_lines` lines for performance.
    """
    f_list, fi_list = [], []
    try:
        with open(os.path.join(filename, "OUTCAR"), 'r', encoding='utf-8', errors='ignore') as f:
            f.seek(0, os.SEEK_END)
            file_size = f.tell()
            block_size = 1024
            blocks, lines_found = [], 0

            while file_size > 0 and lines_found < max_lines:
                file_size = max(0, file_size - block_size)
                f.seek(file_size)
                block = f.read(block_size)
                blocks.append(block)
                lines_found = sum(block.count('\n') for block in blocks)

            tail = ''.join(reversed(blocks)).splitlines()[-max_lines:]
    except Exception:
        with open(os.path.join(filename, "OUTCAR"), 'r', encoding='utf-8', errors='ignore') as f:
            tail = f.readlines()

    for line in tail:
        if "THz" in line and "f" in line:
            match_real = re.search(r"f\s+=\s*([-+]?[0-9]*\.?[0-9]+)", line)
            match_img = re.search(r"f/i=\s*([-+]?[0-9]*\.?[0-9]+)", line)
            if match_real:
                f_list.append(float(match_real.group(1)))
            if match_img:
                fi_list.append(float(match_img.group(1)))
    return np.array(f_list), np.array(fi_list)

def zero_point_energy(freq_THz):
    """
    Compute total zero-point energy (ZPE) from vibrational frequencies (THz).
    ZPE = 0.5 * h * nu
    """
    freq_Hz = freq_THz * 1e12
    zpe = 0.5 * physical_constants["Planck constant in eV s"][0] * freq_Hz
    return zpe.sum()

def quantum_correction(freq_THz, T):
    """
    Compute quantum statistical correction factor at temperature T (K).
    sinh(h$/mu/$2kT) / (h$/mu$/2kT)
    """
    freq_Hz = freq_THz * 1e12
    h = physical_constants["Planck constant"][0]
    kB = physical_constants["Boltzmann constant"][0]
    x = 0.5 * h * freq_Hz / (kB * T)
    return np.sinh(x) / x

def jump_freq(stable_f, saddle_f):
    """
    Calculate raw TST reaction rate constant k0 using frequency ratios.
    """
    # Here can not to use np.prod function, because it number is too huge,
    # which will bring overflow. 
    log_k0 = np.sum(np.log(stable_f)) - np.sum(np.log(saddle_f))
    return np.exp(log_k0)

def jump_freq_zpe_correction(stable_f, saddle_f, T):
    """
    Calculate quantum correction factor to k0 at temperature T.
    """
    stable_corr = quantum_correction(stable_f, T)
    saddle_corr = quantum_correction(saddle_f, T)  
    log_corr = np.sum(np.log(stable_corr)) - np.sum(np.log(saddle_corr))
    return np.exp(log_corr)

def crossover_tem():
    pass

def main():
    stable_dirs, saddle_dirs = get_calculated_point()
    if not stable_dirs or not saddle_dirs:
        raise RuntimeError("Folders named 'stable' or 'saddle' not found.")

    temperatures = [200, 300, 400, 500, 600, 800, 1000, 1400, 1800, 2200]
    output_data = []

    for stable in stable_dirs:
        stable_f, _ = get_frequency(stable)
        zpe_stable = zero_point_energy(stable_f)

        for saddle in saddle_dirs:
            saddle_f, fi_saddle = get_frequency(saddle)
            zpe_saddle = zero_point_energy(saddle_f)

            zpe_corr = zpe_saddle - zpe_stable
            k0 = jump_freq(stable_f, saddle_f)
            correction_factors = {}

            for T in temperatures:
                corr = jump_freq_zpe_correction(stable_f, saddle_f, T)
                correction_factors[T] = float(corr)

            entry = {
                "stable": stable,
                "saddle": saddle,
                "zpe_correction (eV)": float(zpe_corr),
                "k0 (TST rate constant)": float(k0),
                "quantum_corrections": correction_factors
            }
            output_data.append(entry)

    # Write to YAML file
    with open("reaction_rates.yaml", "w", encoding='utf-8') as f:
        yaml.dump(output_data, f, allow_unicode=True, sort_keys=False)

    print("Results saved to reaction_rates.yaml")

if __name__ == "__main__":
    main()