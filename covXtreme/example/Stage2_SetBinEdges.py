import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt

# Import the CovariateBinning class
# Ensure this matches your file structure. 
# If getting ModuleNotFoundError, see the sys.path fix from the previous step.
from CovariateBinning import CovariateBinning

def main():
    # Clear console and close plots
    plt.close('all')

    print("Stage 2: Bin Data")
    print("Choose covariate bins")

    # %% Inputs 
    # Load 'Output/Data.mat' containing 'Dat'
    # Fixed absolute path as per your code
    input_path = os.path.join('C:/these_docs/mon_papier/covXtreme/data', 'Data.mat')
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found. Please run Stage 1 first.")
        return

    try:
        # squeeze_me=True is convenient but dangerous for matrices that must stay 2D
        # struct_as_record=False allows attribute access (Dat.X instead of Dat['X'])
        mat_contents = sio.loadmat(input_path, squeeze_me=True, struct_as_record=False)
        Dat = mat_contents['Dat']
        
        # Access fields
        X = Dat.X
        Y = Dat.Y
        IsPrd = Dat.IsPrd
        RspLbl = Dat.RspLbl
        CvrLbl = Dat.CvrLbl
        
        # --- CRITICAL FIX: Ensure 2D arrays ---
        # CovariateBinning expects [Rows, Cols]. squeeze_me turned [N, 1] into [N,]
        
        # Fix X (Covariates)
        if X.ndim == 1:
            X = X.reshape(-1, 1)
            
        # Fix Y (Response)
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)
            
        # Ensure IsPrd is iterable/array if it became a scalar boolean
        if np.ndim(IsPrd) == 0:
            IsPrd = np.array([IsPrd])
            
        # Ensure Labels are lists/arrays (squeeze_me turns single cell to string)
        if isinstance(RspLbl, str): RspLbl = [RspLbl]
        if isinstance(CvrLbl, str): CvrLbl = [CvrLbl]

    except Exception as e:
        print(f"Error loading Data: {e}")
        # It's helpful to print the shape to debug
        # print(f"X shape: {Dat.X.shape}") 
        return

    # %% Choose bins for each covariate dimension
    # Python: List of numpy arrays.
    BinEdg = [np.array([0, 20, 60, 225, 270, 315])]

    # %% Allocate Data to Bins and Plot
    print(" allocating data to bins...")
    
    # Now X is guaranteed to be (N, 1), so x.shape will return 2 values
    Bn = CovariateBinning(X, BinEdg, IsPrd, Y, RspLbl, CvrLbl)

    # %% Save bin Edges
    output_path = os.path.join('C:/these_docs/mon_papier/covXtreme/data', 'Bin.mat')
    
    # --- Helper function to robustly get attributes ---
    # This handles the mismatch between Python lowercase attributes (self.a) 
    # and the expected Capitalized keys for the .mat file (A).
    def get_attr(obj, names):
        for name in names:
            if hasattr(obj, name):
                return getattr(obj, name)
        return None

    # Construct the dictionary using robust lookups
    bn_dict = {
        'A':      get_attr(Bn, ['A', 'a', 'allocation']),
        'nBin':   get_attr(Bn, ['nBin', 'n_bin', 'n_bins']),
        'Edg':    get_attr(Bn, ['Edg', 'edg', 'edges']),
        'BinLbl': get_attr(Bn, ['BinLbl', 'bin_lbl', 'bin_labels']),
        'nCvr':   get_attr(Bn, ['nCvr', 'n_cvr', 'n_covariates']),
        'CvrLbl': get_attr(Bn, ['CvrLbl', 'cvr_lbl', 'cvr_labels']),
        'IsPrd':  get_attr(Bn, ['IsPrd', 'is_prd', 'is_periodic'])
    }
    
    # Debugging: Check if any key failed to load
    for key, val in bn_dict.items():
        if val is None:
            print(f"Warning: Could not find attribute for '{key}' in CovariateBinning object.")

    sio.savemat(output_path, {'Bn': bn_dict})
    print(f"Binning data saved to {output_path}")

if __name__ == "__main__":
    main()