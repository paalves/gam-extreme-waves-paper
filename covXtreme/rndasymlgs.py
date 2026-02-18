import numpy as np
from itertools import combinations
import math
import rndsymlgs.py

def rndasymlgs(alp, theta, nD, n):
    """
    Generate random values from nD dimension asymmetric logistic function with 
    Frechet margins.
    
    Algorithm 2.1 Alec Stephenson, "Simulating Multivariate Extreme Value
    Distributions of Logistic Type", Extremes 6, 49-59, 2003.

    Args:
        alp (list or np.ndarray): Dependency parameter(s).
        theta (list of lists/arrays): Weighting parameters. Structure mirrors the set B.
        nD (int): Number of dimensions.
        n (int): Number of realisations.

    Returns:
        F (np.ndarray): [nD x n] Random values on Frechet scale.
    """
    
    # --- Find all sets B ---
    # In MATLAB code: B matrix rows contain the elements of the subset.
    # We can replicate this structure or handle it more dynamically in Python.
    # Let's replicate the structure to stay close to the algorithm.
    
    # Calculate total number of non-empty subsets: 2^nD - 1
    nB = 2**nD - 1 
    
    # B will be a list of lists (or arrays) containing the indices of each subset
    B = []
    
    # Generate all subsets of size 1 to nD
    # itertools.combinations returns tuples
    for i in range(1, nD + 1):
        combs = list(combinations(range(1, nD + 1), i)) # 1-based indices to match MATLAB logic
        for c in combs:
            B.append(np.array(c))
            
    # --- Generate Symmetric Values (Z) ---
    Z = [None] * nB
    
    for b in range(nB):
        current_subset = B[b]
        dim_subset = len(current_subset)
        
        # Get specific alpha for this subset
        # Note: alp input in MATLAB seems to align with the index b
        # Assuming alp is a vector of length nB
        current_alp = alp[b]
        
        if dim_subset == 1:
            # 1D is unit Frechet: F = -1 / log(U)
            # MATLAB: 1./(-(log(rand(1,n))))
            Z[b] = 1.0 / (-np.log(np.random.rand(1, n)))
        else:
            # Call symmetric generator
            # Note: rndsymlgs returns [dim_subset x n]
            Z[b] = rndsymlgs(current_alp, dim_subset, n)

    # --- Find Max over Sets for Asymmetric Values ---
    F = np.full((nD, n), np.nan)
    
    # theta is expected to be a list where theta[b] contains weights for subset B[b]
    
    for iD in range(1, nD + 1): # For each dimension 1..nD
        # Find all subsets B containing dimension iD
        # I contains indices of subsets in B that include iD
        I = [b_idx for b_idx, subset in enumerate(B) if iD in subset]
        
        if not I:
            continue
            
        # Initialize F[iD-1, :] with the first relevant subset
        b_idx_first = I[0]
        subset_first = B[b_idx_first]
        
        # Find index of iD within the subset (tI)
        # Python 0-based index
        tI_first = np.where(subset_first == iD)[0][0]
        
        # Term: theta * Z
        # theta[b_idx] should correspond to the subset
        term = theta[b_idx_first][tI_first] * Z[b_idx_first][tI_first, :]
        F[iD-1, :] = term
        
        # Loop over remaining subsets containing iD and take max
        for i in range(1, len(I)):
            b_idx = I[i]
            subset = B[b_idx]
            tI = np.where(subset == iD)[0][0]
            
            term = theta[b_idx][tI] * Z[b_idx][tI, :]
            F[iD-1, :] = np.maximum(F[iD-1, :], term)
            
    return F