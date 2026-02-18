import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt

# Import the MarginalModel class (assuming it's in MarginalModel.py)
from MarginalModel import MarginalModel
# Import the CovariateBinning class (needed to reconstruct Bn object)
from CovariateBinning import CovariateBinning

def main():
    # Clear console and close plots
    plt.close('all')

    print("Stage 3: Fit Marginal Generalised Pareto")
    print("Fit Piecewise constant extreme value model for margin / dimension 'iDmn'")

    # %% Inputs
    # Load data from Stage 1
    data_path = os.path.join('C:/these_docs/mon_papier/covXtreme/data', 'Data.mat')
    if not os.path.exists(data_path):
        print(f"Error: {data_path} not found. Please run Stage 1 first.")
        return
        
    try:
        mat_dat = sio.loadmat(data_path, squeeze_me=True, struct_as_record=False)
        Dat = mat_dat['Dat']
        
        # Verify Dat integrity (e.g. check if X became 1D)
        # The MarginalModel class constructor expects Dat.X to be [nDat x nCvr]
        if hasattr(Dat, 'X') and Dat.X.ndim == 1:
            Dat.X = Dat.X.reshape(-1, 1)
        if hasattr(Dat, 'Y') and Dat.Y.ndim == 1:
            Dat.Y = Dat.Y.reshape(-1, 1)
            
    except Exception as e:
        print(f"Error loading Data: {e}")
        return

    # Load bins from Stage 2
    bin_path = os.path.join('C:/these_docs/mon_papier/covXtreme/data', 'Bin.mat')
    if not os.path.exists(bin_path):
        print(f"Error: {bin_path} not found. Please run Stage 2 first.")
        return

    try:
        mat_bn = sio.loadmat(bin_path, squeeze_me=True, struct_as_record=False)
        Bn_struct = mat_bn['Bn']
        
        # Reconstruct the CovariateBinning object from the loaded struct/dictionary
        # Since we can't save the Python class instance directly to .mat easily,
        # we often need to manually re-instantiate or just pass the struct if the
        # MarginalModel accepts it. 
        # However, MarginalModel checks `if ~isa(Bn,'CovariateBinning')`.
        # So we must recreate a CovariateBinning-like object or use the actual class.
        
        # Option A: Re-instantiate CovariateBinning (requires original inputs X, Edges, etc.)
        # Ideally, we should pickle the python objects instead of using .mat for intermediate python steps.
        # But sticking to the workflow: we will attach the loaded attributes to a dummy CovariateBinning.
        
        class BnWrapper:
            pass
        
        # A workaround to make it look like a CovariateBinning object if we don't want to re-run the binning logic
        Bn = BnWrapper()
        Bn.__class__ = CovariateBinning # Force class type
        Bn.nBin = Bn_struct.nBin
        Bn.nCvr = Bn_struct.nCvr
        Bn.A = Bn_struct.A
        Bn.Edg = Bn_struct.Edg
        if isinstance(Bn.Edg, np.ndarray) and Bn.Edg.dtype == object:
             # Fix object array of arrays if necessary
             Bn.Edg = [Bn.Edg[i] for i in range(len(Bn.Edg))]
        Bn.BinLbl = Bn_struct.BinLbl
        Bn.CvrLbl = Bn_struct.CvrLbl
        Bn.nDat = len(Bn.A)
        
        # Important: Ensure Bn.A is 0-based for Python
        # If it came from MATLAB (1-based), subtract 1. 
        # If it came from Python Stage 2 (0-based), keep it.
        # Check min value to guess.
        if np.min(Bn.A) == 1:
            Bn.A = Bn.A - 1
            
    except Exception as e:
        print(f"Error loading Bin: {e}")
        return

    # %% Dimension to fit
    iDmn = 0  # 0-based index for python (corresponds to MATLAB iDmn=1)
    
    NEP = [0.7, 0.9] # GP non exceedence probability range
    nB = 100         # number bootstrap resamples
    Yrs = 34         # number of years of data
    RtrPrd = [10, 100] # vector of return Periods

    # %% Cross Validation defaults
    CV = {}
    CV['CVMth'] = 0    # 0: Only Cross Validate smoothness for original dataset (fast)
    CV['nCV'] = 10     # number cross-validation groups
    CV['nSmth'] = 10   # number smoothnesses tried in CV
    CV['SmthLB'] = -4  # lower bound (log10) for smoothness range
    CV['SmthUB'] = 4   # upper bound (log10) for smoothness range

    # %% Transformation to Standard-margins 
    MarginType = 'Laplace' # 'Laplace' or 'Gumbel'

    # %% Fit model
    # Initialize and fit the model
    # Note: iDmn is passed as integer. The class should handle the column selection.
    print(f"Fitting Marginal Model for Dimension {iDmn}...")
    MM = MarginalModel(Dat, iDmn, NEP, Bn, nB, Yrs, RtrPrd, CV, MarginType)

    # %% Save result
    # We use pickle here because saving a complex custom class instance (MM) to .mat is difficult/impossible.
    # If downstream MATLAB compatibility is strictly required, you must manually export MM attributes to a dict.
    import pickle
    
        
    output_filename = f'C:/these_docs/mon_papier/covXtreme/data/MM{iDmn+1}.pkl' # Using +1 to match MATLAB naming convention (MM1, MM2)
    
    with open(output_filename, 'wb') as f:
        pickle.dump(MM, f)
        
    print(f"Marginal Model saved to {output_filename}")

if __name__ == "__main__":
    main()