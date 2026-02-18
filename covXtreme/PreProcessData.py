import numpy as np

def PreProcessData(RspOrg, CvrOrg, AscOrg, VslHd):
    """
    Pre-process case study data: remove "0" responses due to heading measurement error 
    giving 90deg.

    Args:
        RspOrg (np.ndarray): Original main response vector.
        CvrOrg (np.ndarray): Original covariate matrix.
        AscOrg (np.ndarray): Original associated response matrix.
        VslHd (np.ndarray): Vessel heading data.

    Returns:
        tuple: (Rsp, Cvr, Asc)
            - Rsp: Filtered main response.
            - Cvr: Filtered covariate.
            - Asc: Filtered associated response.
    """
    # Ensure inputs are numpy arrays
    RspOrg = np.array(RspOrg)
    CvrOrg = np.array(CvrOrg)
    AscOrg = np.array(AscOrg)
    VslHd = np.array(VslHd)

    # Ensure VslHd is shaped correctly for broadcasting against matrices (N x 1)
    if VslHd.ndim == 1:
        VslHd = VslHd.reshape(-1, 1)

    # Create logical mask for rows to keep
    # MATLAB: KpInd = ~any(VslHd==90 | AscOrg== 0, 2);
    # Logic: Mark row if VslHd is 90 OR if any value in AscOrg row is 0.
    # axis=1 performs the 'any' check across columns (rows in the output)
    
    error_condition = (VslHd == 90) | (AscOrg == 0)
    
    # np.any(..., axis=1) reduces rows. We use ~ to negate (keep valid rows).
    if error_condition.ndim > 1:
        KpInd = ~np.any(error_condition, axis=1)
    else:
        KpInd = ~error_condition.flatten()

    # Filter data
    # Note: Variable assignments follow the input logic (Asc gets AscOrg, etc.)
    # despite the potentially conflicting inline comments in the original MATLAB source.
    Asc = AscOrg[KpInd]
    Rsp = RspOrg[KpInd]
    Cvr = CvrOrg[KpInd]

    return Rsp, Cvr, Asc