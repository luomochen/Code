#-----------------------------------------
# Calculate the rate.
#-----------------------------------------
import numpy as np
from scipy import constants as C

def rate(T, barriers, k0):
    """
    Calculate reaction rates and total rate based on DFT results and Arrhenius equation.

    Parameters
    ----------
        T (array_like): Temperature(s) in Kelvin. Can be a scalar or array-like.
        barriers (array_like): Energy barriers (eV) for each reaction pathway.
        k0 (array_like): Pre-exponential factor(s) for each pathway. Can be a scalar or array-like.

    Returns
    -------
        tuple:
            total_rate (float): Total reaction rate summed over all barriers.
            rates (numpy.ndarray): Array of individual rates. 

    Notes:
    ------
        Uses the Arrhenius formula: rate = k0 * exp(-barrier / (k_B * T))
    """
    # Convert inputs to numpy arrays
    T = np.asarray(T)
    barriers = np.asarray(barriers)
    k0 = np.asarray(k0) * 10**12
    # Calculate Boltzmann constant and beta values
    # Betas is seen as a row vector
    k_b = C.physical_constants['Boltzmann constant in eV/K'][0]
    Betas = 1.0/(k_b * T)
    # Outer Product. A new two-dimensional matrix
    # like np.array([Beta1barrier1, Beta1barrier2, ...],[Beta2barrier1, Beta2barrier2, ...])
    # This matrix will be element-wise multiplied with a row vector 
    # after taking the natural exponent
    rates = np.exp(-np.outer(Betas, barriers)) * k0
    total_rate = np.sum(rates)
    return total_rate, rates

def quantum_correction(rates, qc_matrix):
    """
    Apply quantum correction by computing the Hadamard product of rates and qc_matrix.
    
    Parameters
    ----------
        rates (numpy.ndarray or list): Array/list of individual rates.
        qc_matrix (numpy.ndarray or list): Quantum correction matrix (same shape as rates).
        
    Returns
    -------
        qc_rates (numpy.ndarray): Rates after quantum correction (Hadamard product).
        
    Notes:
    ------
        ValueError: If input matrices have different shapes.
    """
    # Convert inputs to numpy arrays
    rates = np.asarray(rates)
    qc_matrix = np.asarray(qc_matrix)
    
    # Check shape consistency
    if rates.shape != qc_matrix.shape:
        raise ValueError(f"Shape mismatch: rates {rates.shape} vs qc_matrix {qc_matrix.shape}")
    
    # Compute Hadamard product
    qc_rates = rates * qc_matrix  # Element-wise multiplication
    
    return qc_rates

def main():
    T = np.array([100, 200])
    barriers = [1, 2]
    jumpfreqs = [1E12, 2E12]
    print(rate(T, barriers, jumpfreqs))

if __name__ == "__main__":
    main()