#-----------------------------------------
# Parse the input file.
# 1. The basic information of simulation.
# 2. The outcome of frequency analysis.
# 3. The structure file in vasp type.
#-----------------------------------------
import yaml
import numpy as np
from ase.io import vasp

def parse_input(yaml_path='input.yaml'):
    """
    Read simulation input parameters from a YAML file and return core parameters.

    Parameters
    ----------
    yaml_path : str, optional
        Path to the input YAML file (default 'input.yaml').

    Returns
    -------
    pw : any
        PWKMC simulation parameter.
    events_mode : any
        Event list generation mode.
    control_step : int
        Control step for the simulation.
    repeat_run : int
        Number of repeats for the simulation.
    nsteps : int
        Total number of simulation steps.
    T_list : list of float
        List of temperatures at which to simulate.
    barriers : any
        Diffusion barriers data.
    select_path : any
        Distance or path selection parameter.

    Notes
    -----
    Jump frequencies (k0_list) are obtained separately from rate data and
    should not be read from this input file.
    """
    # Load YAML
    with open(yaml_path, "r") as file:
        input_data = yaml.safe_load(file)

    # Extract simulation parameters
    pw = input_data["simulation_parameters"]["PWKMC"]
    events_mode = input_data["simulation_parameters"]["event_list_generation_mode"]
    control_step = int(float(input_data["simulation_parameters"]["control_step"]))
    repeat_run = int(float(input_data["simulation_parameters"]["repeat_run"]))
    nsteps = int(float(input_data["simulation_parameters"]["nsteps"]))
    quantum_correction = input_data["simulation_parameters"]["quantum_correction"]
    T_list = input_data["simulation_parameters"]["temperatures"]
    barriers = input_data["simulation_parameters"]["diffusion_barriers"]
    select_path = input_data["simulation_parameters"]["distances"]

    return pw, events_mode, control_step, repeat_run, nsteps, \
        quantum_correction, T_list, barriers, select_path

def parse_rate_data(yaml_path='reaction_rates.yaml', desired_temps=None):
    """
    Parse rate data from a YAML file and filter quantum corrections by desired temperatures.

    Parameters
    ----------
    yaml_path : str, optional
        Path to the YAML file containing rate data (default 'reaction_rates.yaml').
    desired_temps : list of float or int
        Temperatures to filter from the YAML data. Values will be rounded to integers.

    Returns
    -------
    qc_matrix : list of list of float
        Quantum correction coefficients matrix with shape (len(used_temps), num_paths).
        Each row corresponds to a temperature in used_temps and each column to a path.
    k0_list : list of float
        List of TST rate constants for each path, in the same order as in the YAML.
    used_temps : list of int
        List of temperatures actually found and used from the YAML, in the order of desired_temps.
    """
    # Load the YAML file containing rate data
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)

    # Extract list of TST rate constants (k0) for each path
    k0_list = [entry['k0 (TST rate constant)'] for entry in data]

    if desired_temps is not None:
        # Extract all available temperatures from the first entry and sort them
        available = sorted(float(t) for t in data[0]['quantum_corrections'].keys())

        # Filter out temperatures not present in the YAML, preserving order
        used_temps = [T for T in desired_temps if T in available]
        if len(used_temps) < len(desired_temps):
            missing = set(desired_temps) - set(used_temps)
            print(f"Warning: The following temperatures were not found and will be ignored: {sorted(missing)}")

        # Build the quantum correction matrix for each used temperature
        qc_matrix = []
        for T in used_temps:
            # For each temperature, collect corrections for all paths
            row = [entry['quantum_corrections'][T] for entry in data]
            qc_matrix.append(row)

        return qc_matrix, k0_list, used_temps
    else:
        return k0_list

def parse_structure(vasp_file_path='sites.vasp'):
    """
    Parse structure file containing all the sites.
    
    Parameters
    ----------
    vasp_file_path : str, optional
        Path to the vasp type file containing all the considered sites.
        
    Returns
    -------
    atoms: A ase 
    
    Notes
    -----
    Only support vasp type structure files.
    """
    atoms = vasp.read_vasp(vasp_file_path)
    return atoms





# --------------------------------------------------------------------------------
def main():
    import rate
    # Read core simulation parameters
    pw, events_mode, control_step, repeat_run, nsteps, quantum_correction, \
        T_list, barriers, select_path = parse_input('input.yaml')
    if quantum_correction == True:
        # Parse rate data and align jump frequencies with temperatures
        qc_matrix, k0_list, used_temps = parse_rate_data('reaction_rates.yaml', T_list)
        total_rate, rates = rate.rate(T_list, barriers, k0_list)
        qc_rates = rate.quantum_correction(rates, qc_matrix)
        print(f"{qc_rates}")

if __name__ == '__main__':
    main()