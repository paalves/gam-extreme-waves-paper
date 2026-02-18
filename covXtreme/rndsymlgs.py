import numpy as np
from scipy.special import gamma

def rndsymlgs(alp, nD, n):
    """
    Generate random values from nD dimension symmetric logistic function with
    Frechet margins.
    
    Algorithm 1.1 Alec Stephenson, "Simulating Multivariate Extreme Value
    Distributions of Logistic Type", Extremes 6, 49-59, 2003.

    Args:
        alp (float or np.ndarray): 1x1 or nx1 vector of the dependency parameter (0 < alp <= 1).
        nD (int): Number of dimensions.
        n (int): Number of realisations.

    Returns:
        F (np.ndarray): [nD x n] Random values on Frechet margins.
    """
    
    # Handle alpha input (scalar vs vector)
    alp = np.array(alp).flatten()
    
    # Check if alpha is unique scalar or vector
    if np.unique(alp).size > 1:
        # If alphas vary per realization, ensure input length matches n
        if alp.size != n:
             # In MATLAB code: n=numel(alp) effectively overrides n input
             pass 
        nalp = n
        # Replicate alpha to shape [1, n] for broadcasting
        alp_broad = alp.reshape(1, n)
    else:
        nalp = 1
        # Use scalar alpha
        alp_scalar = alp[0]
        alp_broad = alp_scalar # For broadcasting later if needed

    # --- Step 1: Get T from unit exponentials ---
    # W ~ Exp(1)
    W = np.random.exponential(1.0, size=(nD, n))
    # T = W / sum(W) -> Dirichlet related
    T = W / np.sum(W, axis=0)

    # --- Find mixture probabilities from recurrence relationship ---
    # p will be [nD x nD] or [nD x nD x n] if alpha varies
    
    if nalp == 1:
        p = np.zeros((nD, nD))
        p[0, 0] = 1.0 # p(1,1) in 1-based indexing
        
        for iD in range(1, nD): # iD is 0-based index of row (so this is row 2 to nD)
            row = iD + 1 # 1-based row number
            
            # First column (col 0)
            # gamma(iD - alp) / (gamma(iD) * gamma(1 - alp))
            # Note: iD in formula is the row index 2..nD.
            # Using 'row' variable for math consistency with paper/MATLAB
            term1 = gamma(row - alp_scalar) / (gamma(row) * gamma(1 - alp_scalar))
            p[iD, 0] = term1
            
            # Middle columns
            for jD in range(1, row - 1): # jD is 0-based col index (cols 2 to row-1)
                col = jD + 1 # 1-based col number
                
                # Recurrence formula
                # p(iD-1, jD) -> p[iD-1, jD]
                # p(iD-1, jD-1) -> p[iD-1, jD-1]
                
                # MATLAB: t1 = (row-1 - alp*col) * p(row-1, col) + alp*(col-1) * p(row-1, col-1)
                t1 = (row - 1 - alp_scalar * col) * p[iD-1, jD] + \
                     alp_scalar * (col - 1) * p[iD-1, jD-1]
                
                p[iD, jD] = t1 / (row - 1)
            
            # Diagonal (col = row)
            # p(iD, iD) = alp^(iD-1)
            p[iD, iD] = alp_scalar**(row - 1)
            
        # Cumulative probability P (take last row)
        # MATLAB: P = cumsum(shiftdim(p(end,:,:), 1))
        probs = p[nD-1, :]
        P = np.cumsum(probs)
        # P is [nD]
        
    else:
        # Vectorized version for variable alpha
        p = np.zeros((nD, nD, nalp))
        p[0, 0, :] = 1.0
        
        for iD in range(1, nD):
            row = iD + 1
            
            # First column
            num = gamma(row - alp)
            den = gamma(row) * gamma(1 - alp)
            p[iD, 0, :] = num / den
            
            # Middle columns
            for jD in range(1, row - 1):
                col = jD + 1
                t1 = (row - 1 - alp * col) * p[iD-1, jD, :] + \
                     alp * (col - 1) * p[iD-1, jD-1, :]
                p[iD, jD, :] = t1 / (row - 1)
                
            # Diagonal
            p[iD, iD, :] = alp**(row - 1)
            
        # P: [nD x n]
        probs = p[nD-1, :, :]
        P = np.cumsum(probs, axis=0) 

    # --- Step 2: Find k ---
    U = np.random.rand(1, n)
    
    if nalp == 1:
        # P is [nD]. Reshape to [nD, 1] to broadcast against U [1, n]
        # logic: sum(U > P) counts how many CDF steps U passed.
        # P needs to be column vector for broadcasting?
        # MATLAB: sum(U > P, 1). P is vector.
        # In MATLAB U>P creates matrix.
        
        # Python: P[:, None] is [nD, 1]. U is [1, n].
        # P[:, None] < U checks if Probability threshold is less than random draw?
        # MATLAB code: sum(U > P, 1). U is 1xn, P is nDx1 (due to shiftdim).
        # This implies checking U against thresholds.
        
        P_col = P.reshape(-1, 1)
        # Compare: U > P. 
        # If P is CDF [p1, p1+p2, ... 1], finding interval.
        k = np.sum(U > P_col, axis=0) + 1
    else:
        # P is [nD x n]. U is [1 x n].
        k = np.sum(U > P, axis=0) + 1

    # --- Step 3: Get Z from gamma(k, 1) ---
    # MATLAB: gamrnd(k, 1, [1, n]) -> Shape=k, Scale=1
    # Numpy: gamma(shape, scale)
    Z = np.random.gamma(shape=k, scale=1.0, size=(1, n))

    # --- Step 4: Find F on Frechet margins ---
    # F = 1 ./ (Z .* (T .^ alp))
    
    # Broadcast alp if scalar
    if nalp == 1:
        alp_use = alp_scalar
    else:
        alp_use = alp
        
    denom = Z * (T ** alp_use)
    F = 1.0 / denom
    
    return F