import numpy as np
import matplotlib.pyplot as plt
import os

class PeakPickedData:
    """
    Class to hold the output data structure similar to the MATLAB struct.
    """
    def __init__(self):
        self.Y = None       # [nPk x (1+nAsc)] Peak response + associated vars
        self.X = None       # [nPk x nCvr] Covariate data at peaks
        self.RspLbl = None  # Labels
        self.CvrLbl = None  # Labels
        self.IsPrd = None   # Periodic flags

def PeakPicking(Rsp, Cvr, Asc, IsPrd, NEP, RspLbl, CvrLbl):
    """
    Perform Peak Over Threshold (POT) analysis to pick storm peaks.

    Functionality:
    1. Filters NaN values from inputs.
    2. Determines a threshold based on the Non-Exceedance Probability (NEP).
    3. Identifies "storm periods" where the response continuously exceeds the threshold.
    4. Extracts the single maximum peak (and corresponding covariates/associated vars) 
       from each storm period.
    5. Generates diagnostic plots (Marginal and Joint).

    Args:
        Rsp (np.ndarray): [n x 1] Vector of main response.
        Cvr (np.ndarray): [n x nCvr] Matrix of covariates (e.g., direction, season).
        Asc (np.ndarray): [n x nAsc] Matrix of other associated variables.
        IsPrd (np.ndarray): [nCvr] Boolean/Int array, flag if covariate is periodic.
        NEP (float): Quantile level (0-1) used to find threshold.
        RspLbl (list): List of strings, labels for Response (first) and Associated variables.
        CvrLbl (list): List of strings, labels for Covariates.

    Returns:
        Dat (PeakPickedData): Object containing:
            - Dat.Y: [nPk x (1+nAsc)] Matrix with Peak Response in col 0, Associated vars in cols 1+.
            - Dat.X: [nPk x nCvr] Matrix of Covariates at the peak locations.
            - Metadata (labels, flags).
    """
    
    # --- 1. Input Validation and Setup ---
    # Ensure inputs are numpy arrays
    Rsp = np.array(Rsp).flatten()
    Cvr = np.array(Cvr)
    Asc = np.array(Asc)
    IsPrd = np.array(IsPrd).flatten()
    
    n = Rsp.size
    
    # Reshape Cvr/Asc if they are 1D arrays but intended as matrices
    if Cvr.ndim == 1: Cvr = Cvr.reshape(-1, 1)
    if Asc.ndim == 1: Asc = Asc.reshape(-1, 1)

    nCvr = Cvr.shape[1]
    nAsc = Asc.shape[1]

    # Validate shapes
    if Cvr.shape[0] != n: raise ValueError("Cvr rows must match Rsp length")
    if Asc.shape[0] != n: raise ValueError("Asc rows must match Rsp length")
    if IsPrd.size != nCvr: raise ValueError("IsPrd must match number of covariates")
    if not (0 <= NEP <= 1): raise ValueError("NEP must be between 0 and 1")
    if len(RspLbl) != nAsc + 1: raise ValueError("RspLbl length mismatch")
    if len(CvrLbl) != nCvr: raise ValueError("CvrLbl length mismatch")

    # --- 2. Remove NaNs ---
    # MATLAB: I=isnan(Rsp) | any(isnan(Cvr),2) | any(isnan(Asc),2);
    mask_nan = np.isnan(Rsp) | np.any(np.isnan(Cvr), axis=1) | np.any(np.isnan(Asc), axis=1)
    
    # Keep only valid data
    Rsp = Rsp[~mask_nan]
    Cvr = Cvr[~mask_nan, :]
    Asc = Asc[~mask_nan, :]
    
    # Update n after filtering
    n = Rsp.size

    # --- 3. Exceedance Threshold ---
    # Find threshold corresponding to NEP
    Thr = np.quantile(Rsp, NEP)
    
    # Index for exceedances
    IExc = Rsp > Thr
    
    # Observation numbers (indices) for exceedance data
    # np.flatnonzero returns indices where condition is true
    ObsExc = np.flatnonzero(IExc) 

    if len(ObsExc) == 0:
        print("No exceedances found.")
        Dat = PeakPickedData()
        Dat.Y = np.empty((0, 1 + nAsc))
        Dat.X = np.empty((0, nCvr))
        return Dat

    # --- 4. Cluster Storms (De-clustering) ---
    # MATLAB logic uses diff > 1 to find gaps.
    # In Python, we split the array of indices wherever the jump is > 1.
    
    # Find where the difference between consecutive indices is > 1
    gaps = np.where(np.diff(ObsExc) > 1)[0]
    
    # Split ObsExc into subarrays based on gaps. 
    # np.split requires indices to split *at*, so we shift gaps by +1
    storm_clusters = np.split(ObsExc, gaps + 1)
    
    nExc = len(storm_clusters) # Number of independent storms

    # --- 5. Find Peaks within Storms ---
    # We need the global index of the maximum Rsp within each storm cluster
    maxIndOrg = np.zeros(nExc, dtype=int)
    
    for i, cluster_indices in enumerate(storm_clusters):
        # Extract Rsp values for this specific storm
        storm_rsp = Rsp[cluster_indices]
        
        # Find index of max value within this storm (local index)
        local_max_idx = np.argmax(storm_rsp)
        
        # Map back to global index
        maxIndOrg[i] = cluster_indices[local_max_idx]

    # --- 6. Populate Output ---
    Dat = PeakPickedData()
    Dat.Y = np.full((nExc, 1 + nAsc), np.nan)
    
    # Col 0: Maxima response
    Dat.Y[:, 0] = Rsp[maxIndOrg]
    
    # Col 1+: Associated variables at the location of maxima
    if nAsc > 0:
        Dat.Y[:, 1:] = Asc[maxIndOrg, :]
        
    # Covariates at the location of maxima
    Dat.X = Cvr[maxIndOrg, :]
    
    Dat.RspLbl = RspLbl
    Dat.CvrLbl = CvrLbl
    Dat.IsPrd = IsPrd

    nDmn = Dat.Y.shape[1]

    # --- 7. Plotting ---
    if not os.path.exists('Figures'):
        os.makedirs('Figures')

    # Figure 1: Marginal Plot
    # Grid of subplots: Rows = Response/Associated, Cols = Covariates
    fig1 = plt.figure(1, figsize=(10, 8)) # Adjust size as needed
    plt.clf()
    
    c = 0
    # Loop over dimensions (Response + Associated)
    for i in range(nDmn):
        # Loop over Covariates
        for iC in range(nCvr):
            c += 1
            ax = plt.subplot(nDmn, nCvr, c)
            
            # Plot background (all data) in grey
            if i == 0:
                # Main response vs Covariate
                plt.plot(Cvr[:, iC], Rsp, '.', color=[0.7, 0.7, 0.7])
            else:
                # Associated var vs Covariate
                # Note: Asc index is i-1 because Y includes Rsp at 0
                plt.plot(Cvr[:, iC], Asc[:, i-1], '.', color=[0.7, 0.7, 0.7])
            
            # Plot Peaks in black
            plt.plot(Dat.X[:, iC], Dat.Y[:, i], 'k.')
            
            plt.xlabel(CvrLbl[iC])
            plt.ylabel(RspLbl[i])
            plt.grid(True)
            plt.autoscale(enable=True, axis='both', tight=True)
            
            # Periodic axis handling
            if IsPrd[iC]:
                plt.xlim([0, 360])
                plt.xticks(np.arange(0, 361, 45))

    plt.tight_layout()
    plt.savefig('C:/these_docs/mon_papier/covXtreme/figs/Stg1_Data_Margins.png')

    # Figure 2: Joint Plot
    # Scatter of Main Response vs Associated Variables (at peaks)
    if nDmn > 1:
        fig2 = plt.figure(2, figsize=(10, 4))
        plt.clf()
        
        # Number of associated variables
        nAsc_actual = nAsc 
        
        for iA in range(nAsc_actual):
            plt.subplot(1, nAsc_actual, iA + 1)
            
            # Background: All Rsp vs All Asc
            plt.plot(Rsp, Asc[:, iA], '.', color=[0.7, 0.7, 0.7])
            
            # Peaks: Peak Rsp vs Peak Asc
            # Dat.Y[:, 0] is Peak Rsp
            # Dat.Y[:, iA+1] is Peak Asc
            plt.plot(Dat.Y[:, 0], Dat.Y[:, iA+1], 'k.')
            
            plt.xlabel(RspLbl[0])
            plt.ylabel(RspLbl[iA+1])
            plt.grid(True)
            
        plt.tight_layout()
        plt.savefig('C:/these_docs/mon_papier/covXtreme/figs/Stg1_Data_Joint.png')

    return Dat