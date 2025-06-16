#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------
# The script is used to calculate the formation enthalpy of defects.
#----------------------------------------------------------------------
import yaml
import numpy as np
from doped import analysis

def charge_correction(defect_path, dielectric, host_path):
    """Get the kumagai charge correction.
    """
    defect_entry = analysis.defect_entry_from_paths(defect_path=defect_path,
                                                    bulk_path=host_path,
                                                    dielectric=dielectric)
    print(f"Charge: {defect_entry.charge_state} at site: {defect_entry.defect_supercell_site.frac_coords}")
    correction, error = defect_entry.get_kumagai_correction(return_correction_error=True)
    print(f"Error range: {error}")
    correction = defect_entry.corrections['kumagai_charge_correction']
    return correction, error

def calculate_expression(defect, VBM, host_energy, host_path, dielectric):
    """E_f = E_{total}(X^q) - E_perfect - \sum_i n_i\mu_i + q E_fermi + E_coor
    In this formula, the E_f is depent variable and E_fermi is independent variable.
    """
    q = defect['charge']
    if q != 0:
        correction, error = charge_correction(defect['path'], dielectric, host_path)
        correction = float(correction)
        error = float(error)
    else:
        correction = 0
        error = 0
    energy = defect['energy']
    slope = q
    intercept = energy - host_energy + q * VBM - defect['chemical_potential'] + correction
    intercept = float(intercept)
    return {
        'name': defect['name'],
        'slope': slope,
        'intercept': intercept,
        'kumagai_charge_correction': correction,
        'error': error
    }

def safe_load_yaml(defect_data='defect_data.yaml'):
    """Load defect properties form defect_data.yaml.
    """
    with open(defect_data, 'r') as f:
        data = yaml.safe_load(f)

    VBM = data['VBM']
    high_frequency_dielectric = np.array(data['high_frequency_dielectric'])
    ionic_dielectric = np.array(data['ionic_dielectric'])
    dielectric = high_frequency_dielectric + ionic_dielectric

    host_energy = data['host_energy']
    host_path = data['host_path']
    defects = data['defects']

    return VBM, dielectric, host_energy, host_path, defects

def main():
    VBM, dielectric, host_energy, host_path, defects = safe_load_yaml()
    results = []
    for defect in defects:
        result = calculate_expression(defect, VBM, host_energy, host_path, dielectric)
        results.append(result)
        print(f"{result['name']}: ΔH(E_f) = {result['slope']}*E_f + {result['intercept']:.3f}")
    with open('defect_outcome.yaml', 'w') as f:
        yaml.dump(results, f, sort_keys=False)

if __name__ == '__main__':
    main()