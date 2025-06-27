#!/usr/bin/env python
# -*- coding: utf-8 -*-
#--------------------------------------------------
# This script is used for calculating the reaction
# rate from processed rate yaml file.
#--------------------------------------------------
import os
import re
import yaml
import argparse
import numpy as np
from scipy.constants import physical_constants

# Physical constants.
h = physical_constants["Planck constant"][0]
h_in_eV_s = physical_constants["Planck constant in eV s"][0]
kB = physical_constants["Boltzmann constant"][0]

def load_input_from_yaml(filepath):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    # Read the input parameters.
    temperatures = data.get("temperatures", [300])
    output_file = data.get("output_file", "reaction_rates.yaml")

    return temperatures, output_file

def imag_freq_check(saddle_fi):
    """Avoid incorrect saddle points with more than one large imaginary frequency.
    """
    second_large_fi = np.partition(np.unique(saddle_fi), -2)[-2]
    if second_large_fi > 0.6:
        raise ValueError(
            f"The second largest imaginary frequency is {second_large_fi} THz, "
            "which is larger than 0.6 THz. Terminating the script."
        )

def zero_point_energy(freq_THz):
    """
    Compute total zero-point energy (ZPE) from vibrational frequencies (THz).
    ZPE = 0.5 * h * nu
    """
    freq_Hz = freq_THz * 1e12
    zpe = 0.5 * h_in_eV_s * freq_Hz
    return zpe.sum()

def quantum_correction(freq_THz, T):
    """
    Compute quantum statistical correction factor at temperature T (K).
    sinh(h$/mu/$2kT) / (h$/mu$/2kT)
    """
    freq_Hz = freq_THz * 1e12
    x = 0.5 * h * freq_Hz / (kB * T)
    return np.sinh(x) / x

def jump_freq(stable_f, saddle_f):
    """
    Calculate raw TST reaction rate constant k0 using frequency ratios.
    """
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

def get_calculated_point():
    """Find the frequecny file in the folder, such as stabel_f.yaml, stable_1_f.yaml,
    saddle_f.yaml, saddle_fi_1.yaml
    """
    files = os.listdir()
    stable_pattern = re.compile(r"^(stable.*?)_f\.ya?ml$")
    saddle_pattern = re.compile(r"^(saddle.*?)_f\.ya?ml$")

    stable_list = []
    saddle_list = []

    for f in files:
        if match := stable_pattern.match(f):
            stable_list.append(match.group(1))
        elif match := saddle_pattern.match(f):
            saddle_list.append(match.group(1))

    return stable_list, saddle_list

def load_freq(name):
    """Load processed frequency file.
    """
    with open(f"{name}_f.yaml", 'r', encoding="utf-8") as f:
        f_data = np.array(yaml.safe_load(f))
    with open(f"{name}_fi.yaml", 'r', encoding="utf-8") as f:
        fi_data = np.array(yaml.safe_load(f))
    return f_data, fi_data

def main():
    # The user can disignate the input file name.
    parser = argparse.ArgumentParser(description="Run tunneling correction using YAML input.")
    parser.add_argument("--input", type=str, default="input.yaml", help="Path to input YAML file")
    args = parser.parse_args()

    temperatures, output_file = load_input_from_yaml(args.input)
    stable_dirs, saddle_dirs = get_calculated_point()
    print(stable_dirs, saddle_dirs)
    if not stable_dirs or not saddle_dirs:
        raise RuntimeError("No *_f.yaml files found.")

    output_data = []

    for stable in stable_dirs:
        stable_f, stable_fi = load_freq(stable)
        zpe_stable = zero_point_energy(stable_f)
        imag_freq_check(stable_fi)
        for saddle in saddle_dirs:
            saddle_f, saddle_fi = load_freq(saddle)
            imag_freq_check(saddle_fi)
            zpe_saddle = zero_point_energy(saddle_f)
            zpe_corr = zpe_saddle - zpe_stable

            k0 = jump_freq(stable_f, saddle_f)
            correction_factors = {}
            for T in temperatures:
                corr = jump_freq_zpe_correction(stable_f, saddle_f, T)
                correction_factors[T] = float(corr)
            print(correction_factors)
            entry = {
                "stable": stable,
                "saddle": saddle,
                "zpe_correction (eV)": float(zpe_corr),
                "k0 (Classical TST rate constant (THz))": float(k0),
                "quantum_corrections": correction_factors
            }
            output_data.append(entry)

    with open(output_file, "w", encoding='utf-8') as f:
        yaml.dump(output_data, f, allow_unicode=True, sort_keys=False)
    print("Results saved to reaction_rates.yaml")

if __name__ == "__main__":
    main()