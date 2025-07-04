#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------
# This script is used for calculating the quantum 
# tunneling correction factors of reaction rates.
# J. T. Fermann and S. Auerbach, J. Chem. Phys. 112, 6787 (2000).
#----------------------------------------------------------------
import numpy as np
import argparse
import yaml
import warnings
from scipy.constants import physical_constants
from scipy.integrate import quad

# Physical constants
h = physical_constants["Planck constant"][0]
kB = physical_constants["Boltzmann constant"][0]
e = physical_constants["elementary charge"][0]
# The difference of the result between the article and our script 
# may come from the rounding error.
#h = float(f"{h:.2g}")
#kB = float(f"{kB:.2g}")
#e = float(f"{e:.2g}")

def load_input_from_yaml(filepath):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    # Read the input parameters.
    deltaE0 = data["deltaE0"]
    i_omega = np.float64(data["i_omega"]*10**12)
    down_limit = data.get("down_limit", -np.inf)
    temperatures = data.get("temperatures", [300])
    output_file = data.get("output_file", "qtst_tunning.yaml")

    return deltaE0, i_omega, down_limit, temperatures, output_file

def tunneling_correction(T, deltaE0_eV, i_omega, down_limit):
    of_warning = False
    beta = 1 / (T * kB)
    i_energy = h * i_omega
    deltaE0 = deltaE0_eV * e

    theta_up_limit = np.pi * deltaE0 / i_energy

    def integrand(theta):
        return 0.5 * np.exp(beta * i_energy * theta / np.pi) * (np.cosh(theta)**-2)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result, abserr = quad(integrand, down_limit, theta_up_limit, epsabs=1e-10, epsrel=1e-10)
        for warn in w:
            if issubclass(warn.category, RuntimeWarning) and "overflow encountered in cosh" in str(warn.message):
                of_warning = True
                
    surf_term = np.exp(beta * deltaE0) / (1 + np.exp(2 * np.pi * deltaE0 / i_energy))
    gamma = result + surf_term

    return gamma, abserr, of_warning

def crossover_temperature(deltaE0_eV, i_omega):
    i_energy = h * i_omega
    deltaE0 = deltaE0_eV * e
    numerator = i_energy * deltaE0
    denominator = kB * (2 * np.pi * deltaE0 - i_energy * np.log(2))
    return numerator / denominator

def save_to_yaml(file_name, deltaE0, i_omega, down_limit, temperatures, corrections, errors, Tc):
    data = {
        "input_parameters": {
            "deltaE0 (eV)": str(deltaE0),
            "i_omega (THz)": str(i_omega/(10**12)),
            "down_limit": str(down_limit)
        },
        "results": {
            "tunneling_corrections": {str(T): str(val) for T, val in zip(temperatures, corrections)},
            "integration_errors": {str(T): str(val) for T, val in zip(temperatures, errors)},
            "crossover_temperature_K": str(Tc)
        }
    }

    with open(file_name, "w") as f:
        yaml.dump(data, f, sort_keys=False)

def main():
    print("-----------------Tips-----------------------")
    print("------------I hope you know-----------------")
    print("----------what you are doing----------------")
    print("--The tunnling correction is only need for--")
    print("--------Hydrogen and its Isotopes-----------")
    
    parser = argparse.ArgumentParser(description="Run tunneling correction using YAML input.")
    parser.add_argument("--input", type=str, default="qtst_input.yaml", help="Path to input YAML file")
    args = parser.parse_args()

    deltaE0, i_omega, down_limit, temperatures, output_file = load_input_from_yaml(args.input)

    corrections = []
    errors = []
    overflows = []
    for T in temperatures:
        gamma, err, of_warning = tunneling_correction(T, deltaE0, i_omega, down_limit)
        corrections.append(gamma)
        errors.append(err)
        overflows.append(of_warning)
        
    Tc = crossover_temperature(deltaE0, i_omega)
    save_to_yaml(output_file, deltaE0, i_omega, down_limit, temperatures, corrections, errors, Tc)

    print(f"\n++++Results written to {output_file}++++")
    
    if True in overflows:
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("------------overflow warning----------------")
        print("------Please increase the down_limit--------")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

if __name__ == "__main__":
    main()