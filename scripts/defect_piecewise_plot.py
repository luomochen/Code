#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------
# The script is calculate the pairwise function of formation enthalpy
# vs fermi level.
#----------------------------------------------------------------------
import yaml
import numpy as np
import matplotlib.pyplot as plt

def load_thermal_data(thermal_data='defect_outcome.yaml'):
    """Load defect thermal properties from yaml file.
    """
    with open(thermal_data, "r") as f:
        data = yaml.safe_load(f)
        
    defects = []
    for item in data:
        slope = item["slope"]
        intercept = item["intercept"]
        name = item["name"]
        defects.append({
            "name": name,
            "q": slope,
            "E0": intercept
        })
    return defects

def load_host_data(host_data="defect_data.yaml"):
    """Load host band properties from yaml file.
    """
    with open(host_data, "r") as f:
        data = yaml.safe_load(f)
    
    Gap = data['Gap']
    return Gap

def find_intersection(line1, line2):
    """Calculate the intersection of two lines.
    """
    q1, E01 = line1['q'], line1['E0']
    q2, E02 = line2['q'], line2['E0']
    if q1 == q2:
        return None
    return (E02 - E01) / (q1 - q2)

def get_piecewise_segments(defects, fermi_max=1):
    """Calculate the piece functions.
    """
    # Get all the intersections.
    crossings = set()
    for i in range(len(defects)):
        for j in range(i + 1, len(defects)):
            ef = find_intersection(defects[i], defects[j])
            if ef is not None and 0 <= ef <= fermi_max:
                crossings.add(ef)
    ef_points = sorted(list(crossings) + [0, fermi_max])

    segments = []
    # Divide sorted ef_points into intervals.
    for i in range(len(ef_points) - 1):
        ef_left = ef_points[i]
        ef_right = ef_points[i + 1]
        ef_mid = (ef_left + ef_right) / 2 # Using mid points to represent the intervals.
        # Find the lowest line.
        min_energy = float("inf")
        best_line = None
        for line in defects:
            energy = line["E0"] + line["q"] * ef_mid
            if energy < min_energy:
                min_energy = energy
                best_line = line

        segments.append({
            "Ef_range": [round(ef_left, 6), round(ef_right, 6)],
            "name": best_line["name"],
            "q": best_line["q"],
            "E0": round(best_line["E0"], 6)
        })
    return segments

def plot_segments(segments, fermi_max=1.0, filename="piecewise_plot.png"):
    """Plot.
    """
    ef_values = np.linspace(0, fermi_max, 1000)
    energies = []

    for ef in ef_values:
        for seg in segments:
            ef0, ef1 = seg["Ef_range"]
            if ef0 <= ef <= ef1:
                E = seg["E0"] + seg["q"] * ef
                energies.append(E)
                break
    
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.plot(ef_values, energies, color='red', linewidth=3)
    for seg in segments:
        ef0, ef1 = seg["Ef_range"]
        ef_mid = (ef0 + ef1) / 2
        e_mid = seg["E0"] + seg["q"] * ef_mid
        ax.text(ef_mid, 
                e_mid - 0.2, 
                f"${seg['name']}$", 
                fontsize=15, 
                ha='center', 
                va='bottom', 
                rotation=15)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.set_xlabel("Fermi level (eV)", fontsize=20)
    ax.set_ylabel("$\mathrm{\Delta H}$ (eV)", fontsize=20)
    ax.set_title("Defect formation enthalpy (minimum envelope)", fontsize=20)
    ax.grid(True)
    ax.legend()
    
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    
    plt.show()

def save_segments_to_yaml(segments, out_path="piecewise_result.yaml"):
    with open(out_path, "w") as f:
        yaml.safe_dump(segments, f, sort_keys=False)

def main():
    defects = load_thermal_data()
    Gap = load_host_data()
    segments = get_piecewise_segments(defects, Gap)
    save_segments_to_yaml(segments, "defect_piecewise_result.yaml")
    plot_segments(segments, Gap, "defect_piecewise_plot.png")

    print("defect_piecewise_result.yaml is saved")
    print("defect_piecewise_plot.png is saved")

if __name__ == "__main__":
    main()