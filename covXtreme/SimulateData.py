import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal
import os

class SimulatedData:
    """Class to hold simulated data structure."""
    def __init__(self):
        self.Y = None
        self.X = None
        self.IsPrd = True
        self.RspLbl = ['Response', 'Associated']
        self.CvrLbl = ['Direction']

def SimulateData(Mdl):
    """
    Simulate model data from a directional model.

    Args:
        Mdl: Populated model object containing:
             - nBin, nDat, nDmn
             - DrcEdg (np.ndarray): [nBin+1 x 1] Direction edges
             - Rat (np.ndarray): [nBin x 1] Rates
             - MM (list of objects): Marginal models for each dimension with Scl, Shp, Thr
             - Jnt (object): Joint model parameters (Mth, Alp, Theta, Rho)

    Returns:
        Dat: Populated SimulatedData object.
    """
    
    # --- 1. Check Input Parameter Consistency ---
    # Python uses 0-based indexing, dimensions checked via shape[0]
    
    # Check dimensions of DrcEdg and Rat
    # Note: DrcEdg usually has nBin+1 edges or Mdl.DrcEdg might be defined differently.
    # MATLAB code implies size(Mdl.DrcEdg,1) == Mdl.nBin which suggests DrcEdg might be bin centers or starts?
    # Or strict equality check. Let's stick to the logic provided.
    if Mdl.nBin != Mdl.DrcEdg.shape[0] or Mdl.nBin != Mdl.Rat.shape[0]:
        raise ValueError('Dimension of DrcEdg and/or Rate inconsistent with nBins')

    # Check Marginal Models (MM)
    if Mdl.nDmn > 1:
        if Mdl.nBin != Mdl.MM[1].Scl.shape[0]: # Index 1 for 2nd dimension
            raise ValueError('Dimension of Scl inconsistent with nBins in dimension 2')
        if np.array(Mdl.MM[1].Shp).size > 1:
            raise ValueError('GP shape assumed constant, Mdl.MM[1].Shp should be a scalar')
    else:
        if Mdl.nBin != Mdl.MM[0].Scl.shape[0]:
            raise ValueError('Dimension of Scl inconsistent with nBins in dimension 1')
        if np.array(Mdl.MM[0].Shp).size > 1:
            raise ValueError('GP shape assumed constant, Mdl.MM[0].Shp should be a scalar')

    if hasattr(Mdl.Jnt, 'Alp') and np.size(Mdl.Jnt.Alp) > 1:
        if Mdl.Jnt.Alp.shape[0] != Mdl.nBin:
            raise ValueError('Alpha should be [1 x 1] or [Mdl.nBin x 1]')

    # --- 2. Simulate Covariate (Directions) ---
    # Rate CDF
    RatSum = np.sum(Mdl.Rat)
    if RatSum == 0: RatSum = 1.0 # Avoid div by zero
    RatCdf = np.cumsum(Mdl.Rat) / RatSum
    
    # Simulate Uniform values [0,1]
    rand_vals = np.random.rand(Mdl.nDat)
    
    # Allocate to bins 0 to nBin-1 based on RatCdf
    # np.searchsorted finds indices where elements should be inserted to maintain order.
    # This is equivalent to discretize/histc
    A = np.searchsorted(RatCdf, rand_vals)
    # Ensure indices are within bounds (searchsorted returns len(RatCdf) if val > max)
    A = np.clip(A, 0, Mdl.nBin - 1)

    # --- 3. Simulate Joint Tail Behavior ---
    U = np.zeros((Mdl.nDat, Mdl.nDmn))
    
    if Mdl.nDmn > 1:
        method = Mdl.Jnt.Mth
        
        if method == 'ASL': # Asymmetric Logistic
            # 
            if np.max(Mdl.Jnt.Alp) > 1 or np.min(Mdl.Jnt.Alp) < 0:
                raise ValueError('Alpha must be in [0,1] for ASL')
            if np.any(np.array(Mdl.Jnt.Theta) > 1) or np.any(np.array(Mdl.Jnt.Theta) < 0):
                raise ValueError('Theta must be in [0,1] for ASL')
            
            # NOTE: rndasymlgs is a custom function not provided. 
            # Placeholder call assuming it exists or needs implementation.
            # Alp = [NaN, NaN, Mdl.Jnt.Alp]
            # Th = {Theta1, Theta2, [1-Theta1, 1-Theta2]}
            # F = rndasymlgs(Alp, Th, Mdl.nDmn, Mdl.nDat).T
            # U = np.exp(-(1.0 / F))
            raise NotImplementedError("ASL method requires external function 'rndasymlgs'")

        elif method == 'LGS': # Logistic
            # 
            if np.max(Mdl.Jnt.Alp) > 1 or np.min(Mdl.Jnt.Alp) < 0:
                raise ValueError('Alpha must be in [0,1] for LGS')
            
            # Use Alpha specific to the bin A
            # NOTE: rndsymlgs is a custom function.
            # F = rndsymlgs(Mdl.Jnt.Alp[A], Mdl.nDmn, Mdl.nDat).T
            # U = np.exp(-(1.0 / F))
            raise NotImplementedError("LGS method requires external function 'rndsymlgs'")

        elif method == 'MVN': # Multivariate Normal
            # 
            rho_vals = np.array(Mdl.Jnt.Rho).flatten()
            
            if np.max(rho_vals) > 1 or np.min(rho_vals) < -1:
                raise ValueError('Rho must be in [-1, 1]')
            if np.any(rho_vals < 0) or np.any(rho_vals >= 1):
                 # Keeping MATLAB warning logic regarding negative dependence setup
                 pass 

            nRho = len(rho_vals)
            Z = np.full((Mdl.nDat, Mdl.nDmn), np.nan)
            
            if nRho == 1:
                rho = rho_vals[0]
                cov = np.array([[1, rho], [rho, 1]])
                Z = multivariate_normal.rvs(mean=np.zeros(Mdl.nDmn), cov=cov, size=Mdl.nDat)
            else:
                # Loop through unique bins/Rhos if Rho varies by bin
                # Assumes Mdl.Jnt.Rho matches bins if nRho > 1? 
                # MATLAB code: loops iRho=1:NoRho. Assumes A corresponds to indices of Rho?
                # Actually, MATLAB logic: `tabulate(A)` implies A maps to Rho indices?
                # If Rho is [nBin x 1], then A (bin index) maps directly to Rho.
                
                # Iterate over unique bins present in A
                unique_bins = np.unique(A)
                for iBin in unique_bins:
                    # Map bin index to Rho index (assuming 1-to-1 or Rho is [nBin])
                    if iBin < nRho:
                        rho = rho_vals[iBin]
                        cov = np.array([[1, rho], [rho, 1]])
                        mask = (A == iBin)
                        count = np.sum(mask)
                        Z[mask, :] = multivariate_normal.rvs(mean=np.zeros(Mdl.nDmn), cov=cov, size=count)
            
            U = norm.cdf(Z) # Transform Normal to Uniform

        else:
            raise ValueError('Joint tail method not recognised (ASL, LGS, MVN)')
            
    else: # 1D
        U[:, 0] = np.random.rand(Mdl.nDat)

    # --- 4. Transform data to Generalised Pareto margins ---
    Dat = SimulatedData()
    Dat.Y = np.full((Mdl.nDat, Mdl.nDmn), np.nan)
    Dat.X = np.full((Mdl.nDat, 1), np.nan)

    # Bin Edges logic: Add final bin for wrapping (circular 0-360)
    # MATLAB: Edg=[Mdl.DrcEdg(end)-360; Mdl.DrcEdg]; 
    # Assumes DrcEdg is sorted. If DrcEdg(0) is 0, DrcEdg(end) is 360?
    # Python equivalent:
    Edg = np.concatenate(([Mdl.DrcEdg[-1] - 360], Mdl.DrcEdg))
    BinSze = np.diff(Edg)

    # Simulate covariate (directions)
    # A is 0-based index. In MATLAB A was 1-based index into DrcEdg?
    # MATLAB: BinSze(A) + Edg(A). 
    # Since we prepended one edge to Edg, the "original" bins 1..N now correspond to indices 1..N in new Edg array.
    # If A is 0..N-1, we likely need to shift by 1 to access the "main" bins in the modified array, 
    # or just use Mdl.DrcEdg directly if A corresponds to that.
    
    # Let's check logic:
    # MATLAB: A is from 1 to nBin.
    # Edg has nBin+1 elements. BinSze has nBin elements.
    # BinSze(A) picks width. Edg(A) picks start.
    # If A=1, picks Edg(1) (which is DrcEdg(end)-360). This implies bin 1 wraps around?
    # Let's assume standard logic: A maps to indices of BinSze.
    
    # In Python, A is 0..nBin-1.
    Dat.X = (np.random.rand(Mdl.nDat, 1).flatten() * BinSze[A] + Edg[A]) % 360.0

    # Transform to GP Margins
    for i in range(Mdl.nDmn):
        # Retrieve parameters for this dimension
        shp = Mdl.MM[i].Shp
        # Parameters vary by bin A
        scl = Mdl.MM[i].Scl[A].flatten()
        thr = Mdl.MM[i].Thr[A].flatten()
        
        # Inverse GP: u = GP(y) -> y = GP_inv(u)
        # y = thr + (scl/shp) * ((1-u)^(-shp) - 1)
        # Using the static method from previous MarginalModel class if available, or inline:
        
        u_vec = U[:, i]
        
        # Handle scalar shape
        if np.isscalar(shp) or np.size(shp) == 1:
            xi = float(shp)
            if abs(xi) < 1e-5: # Gumbel limit
                Dat.Y[:, i] = thr - scl * np.log(1 - u_vec)
            else:
                Dat.Y[:, i] = thr + (scl / xi) * ((1 - u_vec)**(-xi) - 1)
        
    Dat.IsPrd = True
    
    # --- 5. Plotting ---
    if not os.path.exists('Figures'):
        os.makedirs('Figures')

    # Marginal Plot 
    plt.figure(1)
    plt.clf()
    for i in range(Mdl.nDmn):
        plt.subplot(Mdl.nDmn, 1, i + 1)
        plt.plot(Dat.X, Dat.Y[:, i], 'k.')
        plt.xlabel(Dat.CvrLbl[0])
        plt.xlim([0, 360])
        plt.xticks(np.arange(0, 361, 45))
        plt.ylabel(Dat.RspLbl[i])
    plt.tight_layout()
    plt.savefig('C:/these_docs/mon_papier/covXtreme/figs/Stg1_Data_Simulated_Margins.png')

    # Joint Plot 
    if Mdl.nDmn > 1:
        plt.figure(2)
        plt.clf()
        plt.plot(Dat.Y[:, 0], Dat.Y[:, 1], 'k.')
        plt.xlabel(Dat.RspLbl[0])
        plt.ylabel(Dat.RspLbl[1])
        plt.title('Original Margins')
        plt.grid(True)
        plt.savefig('C:/these_docs/mon_papier/covXtreme/figs/Stg1_Data_Simulated_Joint.png')

    return Dat