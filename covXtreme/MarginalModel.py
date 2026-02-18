import numpy as np
import scipy.stats as stats
from scipy.optimize import fmin
from scipy.special import gamma, gammaln, gammainc, gammaincinv
import os
import matplotlib.pyplot as plt

class MarginalModel:
    """
    Fit Piecewise constant EVA (Extreme Value Analysis) Model.
    
    Functionality:
    - Estimates piecewise constant threshold using local bin quantile given non-exceedance probability (NEP).
    - Annual rate calculated using local Poisson MLE for each bin.
    - Generalised Pareto distribution fitted above the threshold:
        - Shape assumed constant across all bins.
        - Scale varies across covariate bins (constant within bins).
        - Variation in scale parameter across bins controlled by smoothness penalty estimated using cross validation.
    
    Attributes:
        X (np.ndarray): nDat x nCvr, covariate direction data.
        Y (np.ndarray): nDat x 1, response data.
        Yrs (float): Number of years of data.
        RspLbl (str): Response label.
        RspSavLbl (str): Response label for saving plots.
        CvrLbl (list): Covariate label vector.
        nBoot (int): Number of bootstrap resamples.
        RtrPrd (np.ndarray): Return Period in years.
        Bn (object): Bin class defining binning properties (CovariateBinning).
        Scl (np.ndarray): Fitted Generalised Pareto Scale (nBin x nBoot).
        Shp (np.ndarray): Fitted Generalised Pareto Shape (nBoot x 1).
        Omg (np.ndarray): Gamma shape parameter (nBin x nBoot).
        Kpp (np.ndarray): Gamma scale parameter (nBin x nBoot).
        GmmLct (np.ndarray): Gamma location parameter (nBin x 1).
        NEP (np.ndarray): Non Exceedence Probability Quantile threshold (nBoot x 1).
        Thr (np.ndarray): Extreme value threshold (nBin x nBoot).
        Rat (np.ndarray): Count no. observations in each bin (nBin x nBoot).
        BSInd (np.ndarray): Bootstrap index (nDat x nBoot).
        nCvr (int): Number of covariates.
        nDat (int): Number of data observations.
        RVSml (dict): Return value simulation results.
        MarginType (str): 'Laplace' or 'Gumbel'.
        FigureFolder (str): Folder to save figures.
        BnMax (np.ndarray): Bin max (used in check for upper endpoint).
        IsPrd (np.ndarray): Flag for if it is a periodic covariate.
        CVMth (int): CV method (0=Original only, 1=Every bootstrap).
        nCV (int): No. cross-validation groups.
        nSmth (int): No. smoothness values tried in CV.
        SmthLB (float): Lower bound (log10) for smoothness range.
        SmthUB (float): Upper bound (log10) for smoothness range.
        SmthSet (np.ndarray): Set of candidate smoothness parameters.
        OptSmth (np.ndarray): Optimal smoothness from SmthSet.
        CVLackOfFit (np.ndarray): Lack of Fit of candidate smoothness parameters.
    """

    def __init__(self, Dat=None, iDmn=None, NEP=None, Bn=None, nB=1, Yrs=1, RtrPrd=None, CV=None, MarginType='Laplace'):
        """
        Constructor for the MarginalModel class.

        Args:
            Dat (object): Data object containing:
                - X: nDat x nCvr covariate values.
                - Y: nDat x nDmn response data.
                - IsPrd: nCvr x 1 vector (boolean/int) indicating periodic covariates.
                - CvrLbl: List of strings for covariate labels.
                - RspLbl: List of strings for response labels.
            iDmn (int): Dimension of Y to use as response (0-based index recommended for Python, 
                        but logic below adapts to 1-based if strictly following MATLAB or 0-based).
            NEP (float or np.ndarray): Non Exceedence Probability Quantile threshold.
            Bn (object): CovariateBinning class instance.
            nB (int, optional): Number of bootstrap resamples. Defaults to 1.
            Yrs (float, optional): Number of years of data. Defaults to 1.
            RtrPrd (list/array, optional): Return period in years. Defaults to 100.
            CV (dict, optional): Cross validation control parameters.
            MarginType (str, optional): 'Gumbel' or 'Laplace'. Defaults to 'Laplace'.

        Returns:
            MarginalModel: Initialized object.
        """
        # Default property initialization
        self.Yrs = 1
        self.nBoot = 100
        self.RtrPrd = np.array([100])
        self.MarginType = 'Laplace'
        self.FigureFolder = 'C:/these_docs/mon_papier/covXtreme/figs'
        
        # CV Defaults
        self.CVMth = 0
        self.nCV = 10
        self.nSmth = 10
        self.SmthLB = -4
        self.SmthUB = 4
        
        # Empty Constructor Check
        if Dat is None:
            return

        # Check inputs
        self.nDat, self.nCvr = Dat.X.shape
        
        # Validate IsPrd
        if not (np.issubdtype(np.array(Dat.IsPrd).dtype, np.number) or np.issubdtype(np.array(Dat.IsPrd).dtype, np.bool_)):
             raise ValueError("IsPrd must be numeric or logical")
        self.IsPrd = np.array(Dat.IsPrd).flatten()

        for iC in range(self.nCvr):
            if self.IsPrd[iC]: # periodic check
                if np.any(Dat.X[:, iC] > 360) or np.any(Dat.X[:, iC] < 0):
                    raise ValueError(f"Covariate {iC} is periodic but values are outside [0, 360]")
        
        self.X = Dat.X

        # Dimension check and Y assignment
        # Note: Assuming iDmn is passed as 0-based index for Python. 
        # If passed as 1-based from MATLAB context, subtract 1: iDmn - 1
        idx = int(iDmn) 
        if not isinstance(idx, int) or idx < 0:
             raise ValueError("iDmn must be a positive integer index")
             
        self.Y = Dat.Y[:, idx]
        self.RspLbl = Dat.RspLbl[idx]
        self.RspSavLbl = Dat.RspLbl[idx] # Logic for removing special chars can be added here
        self.CvrLbl = Dat.CvrLbl

        # Validate NEP
        NEP_arr = np.array(NEP)
        if np.any(NEP_arr > 1) or np.any(NEP_arr < 0):
            raise ValueError("NEP must be between 0 and 1")

        # Validate Binning Object
        # In Python, we check if the object has specific attributes expected of CovariateBinning
        if not hasattr(Bn, 'nBin'): 
             raise TypeError("Input Bn should be of class: CovariateBinning")
        self.Bn = Bn

        # Optional Inputs processing
        if nB is not None:
            if not isinstance(nB, (int, float)) or nB < 0:
                raise ValueError("nB must be a non-negative scalar")
            self.nBoot = int(nB)
            if self.nBoot == 0:
                self.nBoot = 1
        
        # NEP processing for bootstraps
        # range(NEP)/2 + min(NEP) ; sort(rand...)*range + min
        nep_min = np.min(NEP_arr)
        nep_range = np.ptp(NEP_arr) # Peak to peak (max - min)
        
        first_nep = nep_range / 2.0 + nep_min
        rand_neps = np.sort(np.random.rand(self.nBoot - 1)) * nep_range + nep_min
        self.NEP = np.concatenate(([first_nep], rand_neps))

        if Yrs is not None:
            if Yrs <= 0:
                raise ValueError("Yrs must be positive")
            self.Yrs = Yrs
        
        if RtrPrd is not None:
            self.RtrPrd = np.sort(np.array(RtrPrd))
        else:
            self.RtrPrd = np.array([100])

        if MarginType is not None:
            if MarginType not in ['Gumbel', 'Laplace']:
                raise ValueError("MarginType must be 'Gumbel' or 'Laplace'")
            self.MarginType = MarginType

        # Cross Validation Input Processing
        if CV is not None:
            if 'CVMth' in CV: self.CVMth = CV['CVMth']
            if 'nCV' in CV: self.nCV = CV['nCV']
            if 'nSmth' in CV: self.nSmth = CV['nSmth']
            if 'SmthLB' in CV: self.SmthLB = CV['SmthLB']
            if 'SmthUB' in CV: self.SmthUB = CV['SmthUB']

        # Preallocate Output Arrays
        self.Scl = np.full((self.Bn.nBin, self.nBoot), np.nan)
        self.Shp = np.full((self.nBoot, 1), np.nan)
        self.Omg = np.full((self.Bn.nBin, self.nBoot), np.nan)
        self.Kpp = np.full((self.Bn.nBin, self.nBoot), np.nan)
        self.GmmLct = np.full((self.Bn.nBin, 1), np.nan)
        self.Thr = np.full((self.Bn.nBin, self.nBoot), np.nan)
        self.Rat = np.full((self.Bn.nBin, self.nBoot), np.nan)
        self.OptSmth = np.full((self.nBoot, 1), np.nan)
        self.CVLackOfFit = np.full((self.nSmth, self.nBoot), np.nan)
        self.BSInd = np.full((self.nDat, self.nBoot), np.nan)
        
        # logspace in numpy includes the endpoint by default, similar to MATLAB
        self.SmthSet = np.logspace(self.SmthLB, self.SmthUB, self.nSmth)

        # Fit model
        self.Fit()

        # Return Value
        RVMCRls = 10000
        # Passing 0-based indices for bins if nBin > 0
        bin_indices = np.arange(self.Bn.nBin)
        self.RVSml = self.sample_RV_MC(RVMCRls, bin_indices, self.Bn.nBin)

        # Plots
        if not os.path.exists('Figures'):
            os.makedirs('Figures')
        self.Plot()

    @property
    def nRtr(self):
        """Dependent property: Number of return periods."""
        return self.RtrPrd.size

    def Fit(self):
        """
        Bootstrap fit the threshold, rate, and generalised Pareto model.
        Generalised Pareto estimates smoothness of scale using k-fold cross validation.
        """
        # Reset random seed to get same bootstrap sample for all marginals (Reproducibility)
        np.random.seed(1)
        
        # Create Bootstrap Indices [nDat x nBoot]
        # Column 0 is original data (indices 0 to N-1)
        # Columns 1 to nBoot-1 are resampled with replacement
        original_indices = np.arange(self.nDat).reshape(-1, 1)
        if self.nBoot > 1:
            resampled_indices = np.random.randint(0, self.nDat, size=(self.nDat, self.nBoot - 1))
            self.BSInd = np.hstack((original_indices, resampled_indices))
        else:
            self.BSInd = original_indices

        # Reset RNG to shuffle for subsequent operations
        np.random.seed(None) 

        # Calculate Max per Bin (for upper endpoint check) using the original bin allocation
        # Assuming self.Bn.A is 0-based indices. If 1-based (from MATLAB), subtract 1.
        # Here we assume the CovariateBinning class provides a numpy array 'A' of bin indices.
        bin_indices = self.Bn.A.astype(int)
        
        # Determine if indices are 0-based or 1-based
        if np.min(bin_indices) == 1:
            bin_indices -= 1
            
        self.BnMax = np.full(self.Bn.nBin, np.nan)
        # Calculate max Y in each bin
        for i in range(self.Bn.nBin):
            mask = (bin_indices == i)
            if np.any(mask):
                self.BnMax[i] = np.max(self.Y[mask])

        # Gamma Location Parameter (GmmLct) - Simple quantile or min
        # MATLAB: min(x) - range(x)*0.001
        self.GmmLct = np.full(self.Bn.nBin, 0.0)
        for i in range(self.Bn.nBin):
            mask = (bin_indices == i)
            if np.any(mask):
                vals = self.Y[mask]
                self.GmmLct[i] = np.min(vals) - np.ptp(vals) * 0.001

        # Bootstrap Loop
        for iBt in range(self.nBoot):
            if self.nBoot > 1:
                print(f'Fitting for bootstrap sample {iBt+1} of {self.nBoot}')
            else:
                print('Fitting sample')

            # Get Bootstrap Sample
            tY, A = self.GetBootstrapSample(self.BSInd[:, iBt])
            
            # Ensure A is 0-based for indexing
            if np.min(A) == 1: 
                A = A - 1
            A = A.astype(int)

            # --- Fit Gamma Distribution ---
            for iBn in range(self.Bn.nBin):
                I = (A == iBn)
                if np.any(I):
                    # Fit gamma. MATLAB gamfit returns [shape, scale]. 
                    # Scipy gamma.fit returns (shape, loc, scale).
                    # We fix loc=0 relative to the shifted data (tY - GmmLct).
                    data_shifted = tY[I] - self.GmmLct[iBn]
                    
                    # We force floc=0 to match MATLAB's 2-parameter fit on shifted data
                    # alpha (shape), _, beta (scale)
                    alpha, _, beta = stats.gamma.fit(data_shifted, floc=0)
                    
                    self.Omg[iBn, iBt] = alpha
                    # Orthogonal parameterisation k = alpha * beta? 
                    # MATLAB code: Kpp = p(1).*p(2)
                    self.Kpp[iBn, iBt] = alpha * beta

            # --- Threshold Calculation ---
            # Calculate threshold based on NEP (Non-Exceedance Probability)
            # Uses the inverse Gamma CDF (gaminv)
            self.Thr[:, iBt] = self.gaminv(self.NEP[iBt], 
                                           self.Omg[:, iBt], 
                                           self.Kpp[:, iBt], 
                                           self.GmmLct)

            # --- Get Exceedances ---
            # Threshold for every data point based on its bin
            thr_per_point = self.Thr[A, iBt]
            IExc = (tY > thr_per_point)
            AExc = A[IExc]

            # Rsd: Residuals (Y - u) above threshold
            Rsd = tY[IExc] - self.Thr[AExc, iBt]
            
            # ObsMax used for upper endpoint constraint (MaxObs - Threshold)
            # Note: AExc are indices of the bins where exceedances occurred
            ObsMax = self.BnMax[AExc] - self.Thr[AExc, iBt]

            # --- Rate of occurrence in each bin ---
            # Count occurrences in A (bootstrap sample)
            counts = np.bincount(A, minlength=self.Bn.nBin)
            self.Rat[:, iBt] = counts / self.Yrs

            # --- Generalised Pareto Fit ---
            self.GPCrossValidation(Rsd, ObsMax, AExc, iBt)

        return self

    def GetBootstrapSample(self, indices):
        """
        Extract Y and BinAllocation (A) for a given bootstrap index set.
        """
        # Convert floats to integers for indexing if necessary
        idx = indices.astype(int)
        Y_sample = self.Y[idx]
        A_sample = self.Bn.A[idx]
        return Y_sample, A_sample

    def GPCrossValidation(self, Rsd, ObsMax, AExc, iBt):
        """
        Fit piecewise constant GP to threshold exceedances.
        Includes Cross Validation logic for smoothness parameter.
        """
        # --- Constant Starting Solution (Method of Moments) ---
        if len(Rsd) == 0:
            return # Handle case with no exceedances

        xbar = np.mean(Rsd)
        s2 = np.var(Rsd, ddof=1) # ddof=1 for sample variance
        
        # MoM estimates
        xi0 = 0.0
        # Formula from code: 0.5 * xbar * (xbar^2 / s2 + 1) * 2 
        # (The *2 seems to cancel the 0.5, usually sigma = 0.5 * xbar * (xbar^2/s2 + 1))
        sigma0 = 0.5 * xbar * ((xbar**2 / s2) + 1) * 2
        
        # Initial parameters: [xi, log(sigma_bin1), log(sigma_bin2), ...]
        p0 = np.concatenate(([xi0], np.log(np.full(self.Bn.nBin, sigma0))))

        nExc = len(Rsd)
        
        # --- Cross Validation Logic ---
        # Do CV if: (CV method is 'Original' AND this is the first bootstrap) OR (CV method is 'Always')
        do_cv = (self.CVMth == 0 and iBt > 0) or (self.nSmth == 1)
        
        # Note: Logic in MATLAB: if (obj.CVMth==0 && iBt>1) || obj.nSmth==1 -> CV OFF
        # So: Do CV if (CVMth != 0 OR iBt == 0) AND nSmth > 1
        
        # Logic translated:
        if (self.CVMth == 0 and iBt > 0) or self.nSmth == 1:
            # Cross validation OFF
            if self.nSmth == 1:
                self.OptSmth[iBt] = self.SmthSet[0]
            else:
                # Use first bootstrap's optimal smoothness
                self.OptSmth[iBt] = self.OptSmth[0]
        else:
            print('Starting Cross Validation to estimate GP Scale smoothness:')
            LackOfFit = np.full((self.nCV, self.nSmth), np.nan)
            
            # Split data into nCV groups
            # MATLAB: randi(obj.nCV, nExc, 1) -> Integers 1 to nCV
            ICV = np.random.randint(0, self.nCV, size=nExc)

            for iCV in range(self.nCV):
                print('.', end='')
                # Fit set
                mask_fit = (ICV != iCV)
                ObsMaxFit = ObsMax[mask_fit]
                RsdFit = Rsd[mask_fit]
                AExcFit = AExc[mask_fit]
                
                # Prediction set
                mask_prd = (ICV == iCV)
                ObsMaxPrd = ObsMax[mask_prd]
                RsdPrd = Rsd[mask_prd]
                AExcPrd = AExc[mask_prd]

                for iLam in range(self.nSmth):
                    curr_smth = self.SmthSet[iLam]
                    
                    # Define Objective Function
                    def obj_func(p):
                        return self.GPLikeNonStationary(p, RsdFit, AExcFit, ObsMaxFit, curr_smth)

                    # Check initial NLL
                    if np.isinf(obj_func(p0)):
                        # If initial guess is bad, this might crash, but following MATLAB logic:
                         raise ValueError('Bad starting value in CV')

                    # Optimize
                    # MATLAB fminsearch -> scipy.optimize.fmin (Nelder-Mead)
                    PrmHat = fmin(obj_func, p0, disp=False)
                    
                    # Calculate Lack of Fit on Prediction set (unpenalised NLL)
                    LackOfFit[iCV, iLam] = self.GPLikeNonStationary(PrmHat, RsdPrd, AExcPrd, ObsMaxPrd, 0)
            
            print('') # Newline
            
            # Sum Lack of Fit across folds
            self.CVLackOfFit[:, iBt] = np.sum(LackOfFit, axis=0)
            
            # Find index of min Lack of Fit
            tI = np.argmin(self.CVLackOfFit[:, iBt])
            self.OptSmth[iBt] = self.SmthSet[tI]

        # --- Final Fit Model using optimal smoothness ---
        final_smth = self.OptSmth[iBt]
        
        def final_obj_func(p):
            return self.GPLikeNonStationary(p, Rsd, AExc, ObsMax, final_smth)

        if np.isinf(final_obj_func(p0)):
             raise ValueError('Bad starting value in Final Fit')

        PrmHat = fmin(final_obj_func, p0, disp=False)
        
        # Store parameters
        self.Shp[iBt] = PrmHat[0] # Xi
        self.Scl[:, iBt] = np.exp(PrmHat[1:]) # Sigmas (converted back from log)
        
        return self

    @staticmethod
    def GPLikeNonStationary(PARAMS, Dat, BinAlc, ObsMax, SigPen=0):
        """
        Calculate Penalised Negative Log Likelihood for Piecewise Constant Marginal Model.
        
        Args:
            PARAMS: [(nBins+1)] [xi, log(sig_1), ..., log(sig_nBins)]
            Dat: Threshold exceedances (Y - u)
            BinAlc: Bin allocation indices
            ObsMax: Max observed in each bin (for upper endpoint constraint)
            SigPen: Smoothness penalty
        """
        # Unpack parameters
        Xi = PARAMS[0]
        log_Sig = PARAMS[1:]
        Sig = np.exp(log_Sig)
        nBins = len(Sig)
        
        NLOGL = 0
        
        # Calculate likelihood for each bin
        # We can vectorize this instead of looping for performance in Python, 
        # but looping matches MATLAB structure and is safe for bin-wise constraints.
        
        # Identify unique bins present in the data to avoid looping all if sparse
        unique_bins = np.unique(BinAlc)
        
        for iBin in unique_bins:
            I = (BinAlc == iBin)
            # subset data
            dat_i = Dat[I]
            sig_i = Sig[iBin]
            obs_max_i = ObsMax[I]
            
            # GP Negative Log Likelihood
            # Z = x / sigma
            Z = dat_i / sig_i
            
            # Constraints:
            # 1. Xi < -0.5 (Variance infinite/unstable) -> Infinite NLL
            # 2. Upper endpoint violation: if Xi < 0, max(x) must be < -sigma/Xi
            
            if Xi < -0.5:
                NLOGL = np.inf
                break
            
            if Xi < 0 and np.any(obs_max_i > (-sig_i / Xi)):
                NLOGL = np.inf
                break
            
            # Calculate term (1 + Xi * Z)
            term = 1 + Xi * Z
            
            if np.any(term <= 0):
                NLOGL = np.inf
                break
                
            # Log Likelihood formula:
            # -n*log(sigma) - (1 + 1/xi) * sum(log(1 + xi*z))
            # However, we want Negative Log Likelihood:
            # n*log(sigma) + (1 + 1/xi) * sum(log(1 + xi*z))
            
            # Handle limit xi -> 0 (Gumbel)
            if abs(Xi) < 1e-5:
                # Gumbel: sum(log(sigma) + Z)
                # NLL = n*log(sigma) + sum(Z)
                nll_i = len(dat_i) * np.log(sig_i) + np.sum(Z)
            else:
                nll_i = len(dat_i) * np.log(sig_i) + (1 + 1/Xi) * np.sum(np.log(term))
                
            if np.isinf(nll_i):
                NLOGL = np.inf
                break
                
            NLOGL += nll_i
            
        # Add Smoothness Penalty
        # SigPen * sum((Sig - mean(Sig))^2)
        if not np.isinf(NLOGL):
            penalty = SigPen * np.sum((Sig - np.mean(Sig))**2)
            NLOGL += penalty
            
        return NLOGL


    def Margins(self, iBt=None):
        """
        Transform response Y to Uniform margins (GP -> Unif -> Standard margins).
        
        Args:
            iBt (int or list, optional): Index on bootstrap resample. Defaults to all bootstraps.

        Returns:
            YMrg (np.ndarray): Response data (Y) on Standard margins (Laplace/Gumbel).
            YUnif (np.ndarray): Response data (Y) on Uniform margins.
        """
        if iBt is None:
            # Default to original data (index 0) + all bootstraps
            # MATLAB: iBt=1:obj.nBoot -> indices 0 to nBoot-1
            iBt = np.arange(self.nBoot)
        
        # Ensure iBt is iterable
        if np.isscalar(iBt):
            iBt = np.array([iBt])
            
        YUnif = np.full((len(self.Y), len(iBt)), np.nan)
        YMrg = np.full((len(self.Y), len(iBt)), np.nan)
        
        for i, jBt in enumerate(iBt):
            # Get indices for this bootstrap sample
            # self.BSInd is [nDat x nBoot]
            bs_indices = self.BSInd[:, jBt].astype(int)
            
            tY = self.Y[bs_indices]
            # Assumes Bn.A is 0-based bin indices corresponding to original data
            tA = self.Bn.A[bs_indices].astype(int)
            
            # Check Upper Endpoint
            # If shape (xi) < 0, check upper bound: u - sigma/xi
            shp_val = self.Shp[jBt]
            if shp_val < 0:
                # Vectorized check
                # Thr[tA, jBt] gets threshold for each point's bin
                # Scl[tA, jBt] gets scale for each point's bin
                upper_bound_check = tY - self.Thr[tA, jBt] + self.Scl[tA, jBt] / shp_val
                if np.max(upper_bound_check) > 0:
                     raise ValueError('Upper end point invalid in Margins transformation')

            # Calculate CDF (Uniform Margins)
            # We pass jBt as a scalar to CDF, but it expects specific format usually.
            # Here we follow the logic: CDF(X, A, I)
            YUnif[:, i] = self.CDF(tY, tA, jBt)
            
            # Transform from Uniform to Standard Margins
            YMrg[:, i] = self.INV_Standard(YUnif[:, i])
            
            # Handle numerical edge case where YUnif == 1 (Log of 0 issue)
            # Use Survivor function for better precision at the tail
            J = (YUnif[:, i] == 1)
            if np.any(J):
                Q = self.SURVIVOR(tY[J], tA[J], jBt)
                YMrg[J, i] = self.INV_Standard_survivor(Q)

        if np.any(np.isinf(YMrg)):
             raise ValueError(f'Bad transformation to {self.MarginType}')

        return YMrg, YUnif

    def sample_MC(self, nRls, A=None, nA=None):
        """
        Simulate Monte Carlo draws from marginal model.
        
        Args:
            nRls (int): Number of realizations.
            A (np.ndarray, optional): Index of super bins/bins to simulate.
            nA (int, optional): Number of super bins.

        Returns:
            dict: Dictionary containing 'A', 'I', 'Unf', 'Org'.
        """
        # Default A (Omni situation)
        if A is None:
            # Assuming bin indices 0 to nBin-1. 
            # If 'super bins' aggregate multiple original bins, A maps original bins to super bins.
            # Default: map all bins to super-bin 0? 
            # MATLAB: A=ones(obj.Bn.nBin,1) -> All bins map to super-bin 1.
            # Python: All map to 0.
            A = np.zeros(self.Bn.nBin, dtype=int)
            
        if nA is None:
            nA = np.max(A) + 1 # +1 for 0-based count

        uA = np.unique(A[A >= 0]) # Unique super bins
        
        Sml = {}
        
        # Choose bootstraps (random integer 0 to nBoot-1)
        # MATLAB: randi(obj.nBoot, 1, nRls)
        Sml['I'] = np.random.randint(0, self.nBoot, size=(1, nRls))
        
        # Repmat to fit nA rows
        Sml['I'] = np.tile(Sml['I'], (nA, 1))
        # Flatten for vectorized INV call later
        Sml_I_flat = Sml['I'].flatten(order='F') # MATLAB uses column-major

        # Choose bins from Poisson rate model
        Sml['A'] = np.full((nA, nRls), np.nan)
        Sml['nRls'] = nRls
        
        if self.Bn.nBin > 1:
            for iA in uA:
                # Find original bins that map to super bin iA
                tA = np.where(A == iA)[0] 
                
                if len(tA) == 1:
                    Sml['A'][iA, :] = tA[0]
                else:
                    # Simulate covariate with right rate
                    # Rate of all observations in these bins across selected bootstraps
                    # tRat: [len(tA) x nRls]
                    # We need to pick specific bootstraps for specific columns
                    # Sml['I'][iA, :] is the bootstrap index for each realization
                    
                    # Create index grid
                    # rows: bins in tA, cols: realizations
                    bs_indices = Sml['I'][iA, :]
                    tRat = self.Rat[tA][:, bs_indices]
                    
                    # Cumulative Sum for weighted random choice
                    RatSum = np.sum(tRat, axis=0)
                    # Avoid division by zero
                    RatSum[RatSum == 0] = 1 
                    RatCdf = np.cumsum(tRat, axis=0) / RatSum
                    
                    # Random draw
                    rand_vals = np.random.rand(1, nRls)
                    
                    # Determine which bin is selected (sum(rand > cdf))
                    # This counts how many CDF thresholds the random value passed
                    # tJ is 0-based index into tA
                    tJ = np.sum(rand_vals > RatCdf, axis=0)
                    
                    # Handle edge case where tJ might exceed index due to float precision (unlikely with >)
                    tJ = np.minimum(tJ, len(tA) - 1)
                    
                    Sml['A'][iA, :] = tA[tJ]
        else:
            Sml['A'] = np.zeros((1, nRls))

        # Simulate data on Uniform margins
        Sml['Unf'] = np.random.rand(nA, nRls)
        
        # Probability integral transform to get data onto original margins
        # INV(P, I, A)
        # Reshape to vector for processing
        P_vec = Sml['Unf'].flatten(order='F')
        A_vec = Sml['A'].flatten(order='F').astype(int)
        I_vec = Sml_I_flat.astype(int)
        
        org_vec = self.INV(P_vec, I_vec, A_vec)
        Sml['Org'] = org_vec.reshape((nA, nRls), order='F')
        
        return Sml

    def sample_RV_MC(self, nRls, A=None, nA=None):
        """
        Simulate Monte Carlo draws from Return Value distribution.
        
        Args:
            nRls (int): Number of realizations.
            A (np.ndarray): Super bin mapping.
            nA (int): Number of super bins.
            
        Returns:
            dict: Simulation results.
        """
        if A is None:
            A = np.zeros(self.Bn.nBin, dtype=int)
        if nA is None:
            nA = np.max(A) + 1

        uA = np.unique(A[A >= 0])
        Sml = {'nRls': nRls}
        
        # Choose bootstraps
        # [1 x nRls]
        bs_row = np.random.randint(0, self.nBoot, size=nRls)
        
        # Preallocate
        Sml['A'] = np.full((nA, nRls, self.nRtr), np.nan)
        Sml['Unf'] = np.full((nA, nRls, self.nRtr), np.nan)
        Sml['Org'] = np.full((nA, nRls, self.nRtr), np.nan)
        
        for iRtr in range(self.nRtr):
            # Adjust for return period (inv of exp(-L*T*(1-C)))
            # rho: annual rate [nBin x nRls] based on chosen bootstraps
            rho = self.Rat[:, bs_row]
            
            # LT: Poisson Rate = rho * ReturnPeriod
            LT = rho * self.RtrPrd[iRtr]
            
            # UX: Random Uniform [nBin x nRls]
            UX = np.random.rand(self.Bn.nBin, nRls)
            
            # U: Transformed Probability
            # P_max = exp(-LT(1-U_orig)) -> U_orig = 1 + log(P_max)/LT ??
            # Code: U = 1 + log(UX) ./ LT
            # Note: UX here represents the cumulative probability of the maximum over T years?
            # Actually, this looks like the inverse of the Poisson-Process Maxima CDF.
            with np.errstate(divide='ignore', invalid='ignore'):
                U = 1 + np.log(UX) / LT
            
            # U < 0 implies non-occurrence (below lowest threshold effectively)
            U[U < 0] = np.nan
            
            # Transform to Original Scale for ALL bins
            # We need to broadcast bs_row (1xNRls) to (nBin x NRls)
            bs_grid = np.tile(bs_row, (self.Bn.nBin, 1))
            
            # Flatten for INV
            U_flat = U.flatten(order='F')
            bs_flat = bs_grid.flatten(order='F')
            # Bin indices 0..nBin-1 repeated
            bin_flat = np.tile(np.arange(self.Bn.nBin), nRls) # Order needs checking
            # Standard flatten is row-major in numpy (C-style), MATLAB is F-style.
            # Construct bin_flat to match F-style flattening of U
            bin_flat = np.tile(np.arange(self.Bn.nBin)[:, None], (1, nRls)).flatten(order='F')

            # Get Original Values
            tOrgAllBin_flat = self.INV(U_flat, bs_flat, bin_flat)
            tOrgAllBin = tOrgAllBin_flat.reshape((self.Bn.nBin, nRls), order='F')
            
            # Max over super-bins
            for iA in uA:
                tA = np.where(A == iA)[0]
                
                # Extract subset for this super-bin
                subset_org = tOrgAllBin[tA, :]
                
                # Find max along axis 0 (bins within super-bin)
                # nanmax handles NaNs (non-occurrences)
                # If all are NaN, result is NaN
                max_vals = np.nanmax(subset_org, axis=0)
                
                # Find argmax (which bin produced the max)
                # We fill NaNs with -inf to find argmax safely, then restore
                subset_filled = subset_org.copy()
                subset_filled[np.isnan(subset_filled)] = -np.inf
                argmax_local = np.argmax(subset_filled, axis=0)
                
                Sml['Org'][iA, :, iRtr] = max_vals
                Sml['A'][iA, :, iRtr] = tA[argmax_local]
                
                # Get corresponding U value
                # joint index logic:
                # row indices: tA[argmax_local]
                # col indices: 0..nRls-1
                Sml['Unf'][iA, :, iRtr] = U[tA[argmax_local], np.arange(nRls)]

        # Expand I to match output shape
        Sml['I'] = np.tile(bs_row, (nA, 1, self.nRtr))
        
        return Sml

    def sample_MC_CondX(self, nRls, A=None, nA=None, RspCond=None):
        """
        Simulate Monte Carlo draws conditional on chosen input XCond (via RspCond).
        
        Args:
            nRls (int): Number of realizations.
            A (np.ndarray): Super bin mapping.
            nA (int): Number of super bins.
            RspCond (np.ndarray): [nA x nBoot x nRtr] Response values to condition on.
        """
        if A is None:
            A = np.zeros(self.Bn.nBin, dtype=int)
        if nA is None:
            nA = np.max(A) + 1
        
        uA = np.unique(A[A >= 0])
        Sml = {'nRls': nRls}
        
        # Decide bootstrap sample
        Sml['I'] = np.random.randint(0, self.nBoot, size=nRls)
        
        if RspCond is None:
             raise ValueError("RspCond is required for sample_MC_CondX")

        tnRtr = RspCond.shape[2]
        
        Sml['Unf'] = np.full((nA, nRls, tnRtr), np.nan)
        Sml['A'] = np.full((nA, nRls, tnRtr), np.nan)
        
        # Assign response from RspCond based on selected bootstraps
        # RspCond is [nA x nBoot x nRtr]
        # We want [nA x nRls x nRtr]
        # We need to pick columns from axis 1 of RspCond corresponding to Sml['I']
        Sml['Org'] = RspCond[:, Sml['I'], :]
        
        for iA in uA:
            tA = np.where(A == iA)[0]
            
            if len(tA) > 0:
                # [numel(tA) x nBoot x tnRtr]
                tP = np.full((len(tA), self.nBoot, tnRtr), np.nan)
                tf = np.full((len(tA), self.nBoot, tnRtr), np.nan)
                
                # Calculate PDF and CDF for the conditional response values
                for iRtr in range(tnRtr):
                    # RspCond value for this super-bin, all bootstraps, this return period
                    val = RspCond[iA, :, iRtr] # [nBoot]
                    
                    # We need to evaluate CDF/PDF at these values for EACH sub-bin in tA
                    # This implies broadcasting val across tA?
                    # MATLAB: CDF(val, tA) -> implies val matches structure or broadcasts
                    # Here val is [nBoot], tA is [nSubBins]
                    # We loop subbins for clarity or use broadcasting inside CDF/PDF
                    
                    # Let's loop subbins to be safe given complex broadcasting
                    for idx_sub, sub_bin in enumerate(tA):
                         # Pass val [nBoot], sub_bin (scalar), indices 0..nBoot (implicit)
                         tP[idx_sub, :, iRtr] = self.CDF(val, sub_bin, np.arange(self.nBoot))
                         tf[idx_sub, :, iRtr] = self.PDF(val, sub_bin, np.arange(self.nBoot))

                # Relative weights across each bin (PDF density)
                # Sum over axis 0 (sub-bins)
                SumTf = np.sum(tf, axis=0) # [nBoot x tnRtr]
                SumTf[SumTf == 0] = 1
                W = np.cumsum(tf, axis=0) / SumTf
                
                # Sample bin based on weights W
                # W is [nSubBins x nBoot x tnRtr]
                # We need to sample for each Realization (nRls)
                # Each realization uses a specific Bootstrap (Sml['I'])
                
                # Extract W for the chosen bootstraps: [nSubBins x nRls x tnRtr]
                W_selected = W[:, Sml['I'], :]
                
                # Random draw [1 x nRls x tnRtr]
                # Or [1 x nRls x 1] broadcast? 
                # MATLAB: rand(1, nRls) > W(:, Sml.I, :)
                # Dimensions: (1, nRls, 1) vs (nSubBins, nRls, tnRtr)
                rand_vals = np.random.rand(1, nRls, 1) # Broadcast across tnRtr?
                # Actually MATLAB rand(1, nRls) creates one random number per realization
                # and compares it against the weights for all return periods?
                # "rand(1,nRls)>W(:,Sml.I,:)" -> W is 3D.
                # numpy broadcasting (1, nRls, 1) vs (nSub, nRls, tnRtr) works.
                
                tJ = np.sum(rand_vals > W_selected, axis=0) # Index of selected bin
                tJ = np.minimum(tJ, len(tA) - 1)
                
                # Assign
                Sml['A'][iA, :, :] = tA[tJ]
                
                # Assign Unif
                # tP is [nSubBins x nBoot x tnRtr]
                # tP_selected is [nSubBins x nRls x tnRtr]
                tP_selected = tP[:, Sml['I'], :]
                
                # We want values from tP_selected at indices tJ
                # advanced indexing:
                # dim0: tJ [nRls x tnRtr]
                # dim1: 0..nRls [nRls x tnRtr]
                # dim2: 0..tnRtr [nRls x tnRtr]
                m_grid, r_grid = np.meshgrid(np.arange(nRls), np.arange(tnRtr), indexing='ij')
                Sml['Unf'][iA, :, :] = tP_selected[tJ, m_grid, r_grid]

        Sml['I'] = np.tile(Sml['I'][:, None, None], (nA, 1, tnRtr))
        
        return Sml

        
    def CDF(self, X, A=None, I=None):
        """
        CDF function for marginal model (empirical below threshold - GP above).
        
        Args:
            X (np.ndarray): Locations to compute CDF at.
            A (int or np.ndarray, optional): Bin indices.
            I (int or np.ndarray, optional): Bootstrap indices.
            
        Returns:
            P (np.ndarray): Probability.
        """
        if A is not None and I is None:
            # Default to all bootstraps if A is specified but I is not
            I = np.arange(self.nBoot)
        
        if A is None:
            # CDF for everything (all bins, all bootstraps)
            # Transpose Shape and NEP to match broadcasting (1 x nBoot)
            return self.gamgpcdf(X, self.Shp.T, self.Scl, self.Thr, 
                                 self.Omg, self.Kpp, self.GmmLct, self.NEP.T)
        else:
            # Specified bins and bootstraps
            A = np.array(A, dtype=int)
            I = np.array(I, dtype=int)
            
            if A.size == I.size and A.size > 1:
                # 1-to-1 mapping: Compute for specific bin A[k] and bootstrap I[k]
                # In Python, fancy indexing Scl[A, I] retrieves the specific elements
                # Note: GmmLct is [nBin x 1], so we index by A only
                
                # Reshape for broadcasting if necessary, though 1-to-1 usually returns vector
                return self.gamgpcdf(X, self.Shp[I].T, self.Scl[A, I], self.Thr[A, I], 
                                     self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A].T, self.NEP[I].T)
            else:
                # Compute for subset defined by bins A and bootstraps I (Matrix output likely)
                # Reshape I to row vector, A to column vector for broadcasting if needed, 
                # or rely on numpy broadcasting rules for Scl[A, :][:, I]
                
                # Scl[A][:, I] fetches submatrix
                # Shp[I] fetches subset
                
                # Use np.ix_ to construct open meshes for correct slicing
                sub_scl = self.Scl[np.ix_(A.flatten(), I.flatten())]
                sub_thr = self.Thr[np.ix_(A.flatten(), I.flatten())]
                sub_omg = self.Omg[np.ix_(A.flatten(), I.flatten())]
                sub_kpp = self.Kpp[np.ix_(A.flatten(), I.flatten())]
                sub_gmmlct = self.GmmLct[A] # [nA x 1]
                
                return self.gamgpcdf(X, self.Shp[I].T, sub_scl, sub_thr, 
                                     sub_omg, sub_kpp, sub_gmmlct, self.NEP[I].T)

    def SURVIVOR(self, X, A=None, I=None):
        """Survivor function for marginal model."""
        if A is not None and I is None:
            I = np.arange(self.nBoot)
            
        if A is None:
            return self.gamgpsurvivor(X, self.Shp.T, self.Scl, self.Thr, 
                                      self.Omg, self.Kpp, self.GmmLct, self.NEP.T)
        else:
            A = np.array(A, dtype=int)
            I = np.array(I, dtype=int)
            
            if A.size == I.size and A.size > 1:
                return self.gamgpsurvivor(X, self.Shp[I].T, self.Scl[A, I], self.Thr[A, I], 
                                          self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A].T, self.NEP[I].T)
            else:
                sub_scl = self.Scl[np.ix_(A.flatten(), I.flatten())]
                sub_thr = self.Thr[np.ix_(A.flatten(), I.flatten())]
                sub_omg = self.Omg[np.ix_(A.flatten(), I.flatten())]
                sub_kpp = self.Kpp[np.ix_(A.flatten(), I.flatten())]
                sub_gmmlct = self.GmmLct[A]
                
                return self.gamgpsurvivor(X, self.Shp[I].T, sub_scl, sub_thr, 
                                          sub_omg, sub_kpp, sub_gmmlct, self.NEP[I].T)

    def PDF(self, X, A=None, I=None):
        """PDF function for marginal model."""
        if A is not None and I is None:
            I = np.arange(self.nBoot)
            
        if A is None:
            return self.gamgppdf(X, self.Shp.T, self.Scl, self.Thr, 
                                 self.Omg, self.Kpp, self.GmmLct, self.NEP.T)
        else:
            A = np.array(A, dtype=int)
            I = np.array(I, dtype=int)
            
            if A.size == I.size and A.size > 1:
                return self.gamgppdf(X, self.Shp[I].T, self.Scl[A, I], self.Thr[A, I], 
                                     self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A].T, self.NEP[I].T)
            else:
                sub_scl = self.Scl[np.ix_(A.flatten(), I.flatten())]
                sub_thr = self.Thr[np.ix_(A.flatten(), I.flatten())]
                sub_omg = self.Omg[np.ix_(A.flatten(), I.flatten())]
                sub_kpp = self.Kpp[np.ix_(A.flatten(), I.flatten())]
                sub_gmmlct = self.GmmLct[A]
                return self.gamgppdf(X, self.Shp[I].T, sub_scl, sub_thr, 
                                     sub_omg, sub_kpp, sub_gmmlct, self.NEP[I].T)

    def LogPDF(self, X, A=None, I=None):
        """Log PDF function for marginal model."""
        if A is not None and I is None:
            I = np.arange(self.nBoot)
            
        if A is None:
            return self.loggamgppdf(X, self.Shp.T, self.Scl, self.Thr, 
                                    self.Omg, self.Kpp, self.GmmLct, self.NEP.T)
        else:
            A = np.array(A, dtype=int)
            I = np.array(I, dtype=int)
            
            if A.size == I.size and A.size > 1:
                return self.loggamgppdf(X, self.Shp[I].T, self.Scl[A, I], self.Thr[A, I], 
                                        self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A].T, self.NEP[I].T)
            else:
                sub_scl = self.Scl[np.ix_(A.flatten(), I.flatten())]
                sub_thr = self.Thr[np.ix_(A.flatten(), I.flatten())]
                sub_omg = self.Omg[np.ix_(A.flatten(), I.flatten())]
                sub_kpp = self.Kpp[np.ix_(A.flatten(), I.flatten())]
                sub_gmmlct = self.GmmLct[A]
                return self.loggamgppdf(X, self.Shp[I].T, sub_scl, sub_thr, 
                                        sub_omg, sub_kpp, sub_gmmlct, self.NEP[I].T)

    def INV(self, P, I, A=None):
        """
        Inverse CDF for marginal model.
        
        Args:
            P: Probability.
            I: Index of bootstraps to use.
            A: Index of bins to use.
        """
        if A is None:
            # Implicit assumption: if A is missing, we might be doing scalar or matrix case based on P/I
            # But Python requires explicit logic. Assuming 3 cases based on I size.
            pass
            
        I = np.array(I, dtype=int)
        
        # Case 1: Scalar I (Single bin/bootstrap context)
        if I.size == 1:
            A = np.array(A, dtype=int)
            # Need to handle if A is scalar or vector
            # Accessing via [A, I] works for both if A is array
            return self.gamgpinv(P, self.Shp[I], self.Scl[A, I], self.Thr[A, I], 
                                 self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A], self.NEP[I])

        # Check dimensions of P to determine if Matrix or Vector case
        # MATLAB: prod(p(2:end)) > 1 -> Matrix
        p_shape = np.array(P).shape
        is_matrix = len(p_shape) > 1 and np.prod(p_shape[1:]) > 1
        
        if is_matrix:
            # Case 3: Matrix P, Matrix I/A context (All bootstraps/bins)
            # I acts as columns?
            return self.gamgpinv(P, self.Shp[I].T, self.Scl[:, I], self.Thr[:, I], 
                                 self.Omg[:, I], self.Kpp[:, I], self.GmmLct, self.NEP[I].T)
        else:
            # Case 2: Vector P/I (Across sampled bins and bootstraps)
            A = np.array(A, dtype=int)
            if self.Bn.nBin == 1:
                return self.gamgpinv(P, self.Shp[I], self.Scl[I].T, self.Thr[I].T, 
                                     self.Omg[I].T, self.Kpp[I].T, self.GmmLct, self.NEP[I])
            else:
                # 1-to-1 mapping A[k] with I[k]
                return self.gamgpinv(P, self.Shp[I], self.Scl[A, I], self.Thr[A, I], 
                                     self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A], self.NEP[I])

    def INV_survivor(self, Q, I, A=None):
        """Inverse CDF using survival probability."""
        I = np.array(I, dtype=int)
        
        # Case 1: Scalar
        if I.size == 1:
            A = np.array(A, dtype=int)
            return self.gamgpinvsurvivor(Q, self.Shp[I], self.Scl[A, I], self.Thr[A, I], 
                                         self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A], self.NEP[I])

        p_shape = np.array(Q).shape
        is_matrix = len(p_shape) > 1 and np.prod(p_shape[1:]) > 1
        
        if is_matrix:
            # Case 3
            return self.gamgpinvsurvivor(Q, self.Shp[I].T, self.Scl[:, I], self.Thr[:, I], 
                                         self.Omg[:, I], self.Kpp[:, I], self.GmmLct, self.NEP[I].T)
        else:
            # Case 2
            A = np.array(A, dtype=int)
            if self.Bn.nBin == 1:
                return self.gamgpinvsurvivor(Q, self.Shp[I], self.Scl[I].T, self.Thr[I].T, 
                                             self.Omg[I].T, self.Kpp[I].T, self.GmmLct, self.NEP[I])
            else:
                return self.gamgpinvsurvivor(Q, self.Shp[I], self.Scl[A, I], self.Thr[A, I], 
                                             self.Omg[A, I], self.Kpp[A, I], self.GmmLct[A], self.NEP[I])

    def INV_Standard(self, P):
        """Transform from uniform to standard margins using inverse CDF."""
        if self.MarginType == 'Gumbel':
            X = -np.log(-np.log(P))
        elif self.MarginType == 'Laplace':
            # sign(0.5-P) * log(2*min(1-P, P))
            # Note: P=0.5 yields log(1)=0, sign=0 -> 0. Correct.
            X = np.sign(0.5 - P) * np.log(2 * np.minimum(1 - P, P))
        else:
            raise ValueError('Margin Type not recognised')
        return X

    def INV_Standard_survivor(self, Q):
        """Transform to standard margins using survival probability."""
        if self.MarginType == 'Gumbel':
            # -log(-log1p(-Q))
            X = -np.log(-np.log1p(-Q))
        elif self.MarginType == 'Laplace':
            X = np.sign(Q - 0.5) * np.log(2 * np.minimum(Q, 1 - Q))
        else:
            raise ValueError('Margin Type not recognised')
        return X

    def CDF_Standard(self, X):
        """Transform from standard to uniform margins using CDF."""
        if self.MarginType == 'Gumbel':
            P = np.exp(-np.exp(-X))
        elif self.MarginType == 'Laplace':
            P = (X > 0).astype(float) - 0.5 * np.sign(X) * np.exp(-np.abs(X))
        else:
            raise ValueError('Margin Type not recognised')
        return P

    def Survivor_Standard(self, X):
        """Transform from standard to uniform margins using survival probability."""
        if self.MarginType == 'Gumbel':
            # -expm1(-exp(-X)) -> -(exp(-exp(-X)) - 1) -> 1 - P
            Q = -np.expm1(-np.exp(-X))
        elif self.MarginType == 'Laplace':
            Q = (X <= 0).astype(float) + 0.5 * np.sign(X) * np.exp(-np.abs(X))
        else:
            raise ValueError('Margin Type not recognised')
        return Q

    def PDF_Standard(self, X):
        """Get probability density on standard margins."""
        if self.MarginType == 'Gumbel':
            F = np.exp(-(X + np.exp(-X)))
        elif self.MarginType == 'Laplace':
            F = 0.5 * np.exp(-np.abs(X))
        else:
            raise ValueError('Margin Type not recognised')
        
        F[np.isnan(F)] = 0
        return F

    def LogPDF_Standard(self, X):
        """Get log probability density on standard margins."""
        if self.MarginType == 'Gumbel':
            F = -(X + np.exp(-X))
        elif self.MarginType == 'Laplace':
            F = np.log(0.5) - np.abs(X)
        else:
            raise ValueError('Margin Type not recognised')
            
        F[np.isnan(F)] = 0
        return F


        # -------------------------------------------------------------------------
    # Static Distribution Methods
    # -------------------------------------------------------------------------

# -------------------------------------------------------------------------
    # Static Distribution Methods (Fixed Broadcasting)
    # -------------------------------------------------------------------------

    @staticmethod
    def gpinv(P, Xi, Sgm, Thr):
        """Inverse of Generalized Pareto (GP) CDF."""
        # Broadcast all inputs to a common shape
        # This handles cases where P is scalar but parameters are vectors (per bin)
        P, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(P), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        
        X = np.full(P.shape, np.nan)
        
        # Gumbel case (Xi close to 0)
        I_gumbel = np.abs(Xi) < 1e-5
        
        # Standard GP case
        mask_gp = ~I_gumbel
        if np.any(mask_gp):
            with np.errstate(divide='ignore', invalid='ignore'):
                t1 = ((1 - P[mask_gp])**(-Xi[mask_gp]) - 1) / Xi[mask_gp]
            X[mask_gp] = Sgm[mask_gp] * t1 + Thr[mask_gp]

        # Apply mask for Gumbel
        if np.any(I_gumbel):
            t1 = -np.log(1 - P[I_gumbel])
            X[I_gumbel] = Sgm[I_gumbel] * t1 + Thr[I_gumbel]

        return X

    @staticmethod
    def gpinvsurvivor(Q, Xi, Sgm, Thr):
        """Inverse of GP using survivor probability Q."""
        Q, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(Q), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        X = np.full(Q.shape, np.nan)

        I_gumbel = np.abs(Xi) < 1e-5
        
        mask_gp = ~I_gumbel
        if np.any(mask_gp):
            with np.errstate(divide='ignore', invalid='ignore'):
                t1 = (Q[mask_gp]**(-Xi[mask_gp]) - 1) / Xi[mask_gp]
            X[mask_gp] = Sgm[mask_gp] * t1 + Thr[mask_gp]

        if np.any(I_gumbel):
            t1 = -np.log(Q[I_gumbel])
            X[I_gumbel] = Sgm[I_gumbel] * t1 + Thr[I_gumbel]

        return X

    @staticmethod
    def gppdf(X, Xi, Sgm, Thr):
        """PDF of Generalized Pareto."""
        X, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        F = np.zeros(X.shape)
        
        Z = (X - Thr) / Sgm
        
        # Standard Case
        I_gp = np.abs(Xi) >= 1e-5
        
        # We need a combined mask for valid domain
        mask_gp = I_gp & (1 + Xi * Z > 0)
        
        if np.any(mask_gp):
            t1_valid = 1 + Xi[mask_gp] * Z[mask_gp]
            F[mask_gp] = (1.0 / Sgm[mask_gp]) * (t1_valid ** (-1.0/Xi[mask_gp] - 1.0))

        # Gumbel case
        I_gum = np.abs(Xi) < 1e-5
        if np.any(I_gum):
            F[I_gum] = (1.0 / Sgm[I_gum]) * np.exp(-Z[I_gum])

        # Zero below threshold (already 0 init, but ensures logic)
        F[X < Thr] = 0
        return F

    @staticmethod
    def loggppdf(X, Xi, Sgm, Thr):
        """Log PDF of Generalized Pareto."""
        X, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        logF = np.full(X.shape, -np.inf)
        
        Z = (X - Thr) / Sgm
        
        I_gp = np.abs(Xi) >= 1e-5
        mask_gp = I_gp & (1 + Xi * Z > 0) & (X >= Thr)
        if np.any(mask_gp):
            t1 = 1 + Xi[mask_gp] * Z[mask_gp]
            logF[mask_gp] = -np.log(Sgm[mask_gp]) - (1.0/Xi[mask_gp] + 1.0) * np.log(t1)

        I_gum = np.abs(Xi) < 1e-5
        mask_gum = I_gum & (X >= Thr)
        if np.any(mask_gum):
            logF[mask_gum] = -np.log(Sgm[mask_gum]) - Z[mask_gum]

        return logF

    @staticmethod
    def gpcdf(X, Xi, Sgm, Thr):
        """CDF of Generalized Pareto."""
        X, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        P = np.zeros(X.shape)
        
        Z = (X - Thr) / Sgm
        
        mask_gp = (np.abs(Xi) >= 1e-5) & (X > Thr)
        # Check upper bound 1+xi*z > 0
        mask_gp &= (1 + Xi * Z > 0)
        
        if np.any(mask_gp):
            t1 = 1 + Xi[mask_gp] * Z[mask_gp]
            P[mask_gp] = 1 - t1**(-1.0/Xi[mask_gp])
            
        # Upper bound handling: if 1+xi*z <= 0 and xi<0, P=1
        mask_ub = (np.abs(Xi) >= 1e-5) & (X > Thr) & (1 + Xi * Z <= 0) & (Xi < 0)
        P[mask_ub] = 1.0

        mask_gum = (np.abs(Xi) < 1e-5) & (X > Thr)
        if np.any(mask_gum):
            P[mask_gum] = 1 - np.exp(-Z[mask_gum])

        return P

    @staticmethod
    def gpsurvivor(X, Xi, Sgm, Thr):
        """Survivor function of Generalized Pareto."""
        X, Xi, Sgm, Thr = np.broadcast_arrays(np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr))
        Q = np.ones(X.shape) # Default 1 (below threshold)
        
        Z = (X - Thr) / Sgm
        
        mask_exc = (X > Thr)
        
        mask_gp = mask_exc & (np.abs(Xi) >= 1e-5)
        # Check upper bound
        mask_gp_valid = mask_gp & (1 + Xi * Z > 0)
        
        if np.any(mask_gp_valid):
            t1 = 1 + Xi[mask_gp_valid] * Z[mask_gp_valid]
            Q[mask_gp_valid] = t1**(-1.0/Xi[mask_gp_valid])
            
        # If beyond upper bound, survivor is 0
        mask_gp_invalid = mask_gp & (1 + Xi * Z <= 0)
        Q[mask_gp_invalid] = 0

        mask_gum = mask_exc & (np.abs(Xi) < 1e-5)
        if np.any(mask_gum):
            Q[mask_gum] = np.exp(-Z[mask_gum])

        return Q

    @staticmethod
    def gampdf(X, Alp, Bet, GmmLct):
        """Gamma PDF using orthogonal parameterisation."""
        X, Alp, Bet, GmmLct = np.broadcast_arrays(np.asarray(X), np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct))
        Z = X - GmmLct
        Z[Z < 0] = 0
        
        scale = Bet / Alp
        
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = (scale)**(-Alp) / gamma(Alp)
            term2 = Z**(Alp - 1)
            term3 = np.exp(-Z / scale)
            F = term1 * term2 * term3
            
        F[np.isnan(F)] = 0
        return F

    @staticmethod
    def loggampdf(X, Alp, Bet, GmmLct):
        """Log Gamma PDF."""
        X, Alp, Bet, GmmLct = np.broadcast_arrays(np.asarray(X), np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct))
        Z = X - GmmLct
        Z[Z < 0] = 0
        
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = -Alp * (np.log(Bet) - np.log(Alp))
            term2 = -gammaln(Alp)
            term3 = (Alp - 1) * np.log(Z)
            term4 = -Alp * Z / Bet
            logF = term1 + term2 + term3 + term4
            
        logF[Z <= 0] = -np.inf
        return logF

    @staticmethod
    def gamcdf(X, Alp, Bet, GmmLct):
        """Gamma CDF."""
        X, Alp, Bet, GmmLct = np.broadcast_arrays(np.asarray(X), np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct))
        Z = X - GmmLct
        Z[Z < 0] = 0
        
        Z_scaled = Z * (Alp / Bet)
        return gammainc(Alp, Z_scaled)

    @staticmethod
    def gamsurvivor(X, Alp, Bet, GmmLct):
        """Gamma Survivor function."""
        # Broadcasting handled in gamcdf
        return 1.0 - MarginalModel.gamcdf(X, Alp, Bet, GmmLct)

    @staticmethod
    def gaminv(P, Alp, Bet, GmmLct):
        """Inverse Gamma CDF (Fixed for Broadcasting)."""
        # Ensure all inputs are broadcast to the same shape 
        # This fixes the TypeError when P is scalar and parameters are vectors
        P, Alp, Bet, GmmLct = np.broadcast_arrays(np.asarray(P), np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct))
        
        q = np.full(P.shape, np.nan)
        
        # Lower tail
        # Now P and Alp have same shape, so boolean indexing I works for both
        I = ~np.isnan(P) & (P >= 1e-2)
        if np.any(I):
            q[I] = gammaincinv(Alp[I], P[I])
            
        # Upper tail logic
        I_low = ~np.isnan(P) & (P < 1e-2)
        if np.any(I_low):
            q[I_low] = gammaincinv(Alp[I_low], P[I_low])

        X = q * (Bet / Alp) + GmmLct
        return X

    @staticmethod
    def gaminvsurvivor(Q, Alp, Bet, GmmLct):
        """Inverse Gamma using Survival Probability Q."""
        from scipy.special import gammainccinv
        
        Q, Alp, Bet, GmmLct = np.broadcast_arrays(np.asarray(Q), np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct))
        
        q = np.full(Q.shape, np.nan)
        
        I = ~np.isnan(Q)
        if np.any(I):
            # gammainccinv is Inverse of the upper incomplete gamma function (regularized)
            q[I] = gammainccinv(Alp[I], Q[I])
            
        X = q * (Bet / Alp) + GmmLct
        return X

    @staticmethod
    def gamgpcdf(X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """CDF of Gamma-GP Mixture."""
        # Broadcast all
        X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        IBlw = X < Thr
        P = np.full(X.shape, np.nan)
        
        if np.any(IBlw):
            gam_X = MarginalModel.gamcdf(X[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            gam_Thr = MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            P[IBlw] = gam_X * Tau[IBlw] / gam_Thr
            
        if np.any(~IBlw):
            gp_X = MarginalModel.gpcdf(X[~IBlw], Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            P[~IBlw] = gp_X * (1 - Tau[~IBlw]) + Tau[~IBlw]
            
        return P

    @staticmethod
    def gamgpsurvivor(X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """Survivor function of Gamma-GP Mixture."""
        X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        IBlw = X < Thr
        Q = np.full(X.shape, np.nan)
        
        if np.any(IBlw):
            gam_surv_X = 1 - MarginalModel.gamcdf(X[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            gam_cdf_Thr = MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            Q[IBlw] = gam_surv_X * Tau[IBlw] / gam_cdf_Thr
            
        if np.any(~IBlw):
            gp_surv_X = MarginalModel.gpsurvivor(X[~IBlw], Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            Q[~IBlw] = gp_surv_X * (1 - Tau[~IBlw])
            
        return Q

    @staticmethod
    def gamgppdf(X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """PDF of Gamma-GP Mixture."""
        X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        IBlw = (X <= Thr)
        f = np.full(X.shape, np.nan)
        
        if np.any(IBlw):
            f1 = MarginalModel.gampdf(X[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            cdf_thr = MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            f[IBlw] = f1 * Tau[IBlw] / cdf_thr
            
        if np.any(~IBlw):
            f2 = MarginalModel.gppdf(X[~IBlw], Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            f[~IBlw] = f2 * (1 - Tau[~IBlw])
            
        return f

    @staticmethod
    def loggamgppdf(X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """Log PDF of Gamma-GP Mixture."""
        X, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(X), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        IBlw = (X <= Thr)
        f = np.full(X.shape, np.nan)
        
        if np.any(IBlw):
            log_f1 = MarginalModel.loggampdf(X[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            log_cdf = np.log(MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw]))
            f[IBlw] = log_f1 + np.log(Tau[IBlw]) - log_cdf
            
        if np.any(~IBlw):
            log_f2 = MarginalModel.loggppdf(X[~IBlw], Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            f[~IBlw] = log_f2 + np.log(1 - Tau[~IBlw])
            
        return f

    @staticmethod
    def gamgpinv(P, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """Inverse CDF of Gamma-GP Mixture."""
        # Broadcast all inputs
        P, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(P), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        X = np.full(P.shape, np.nan)
        
        IBlw = P < Tau
        
        if np.any(IBlw):
            gam_cdf_Thr = MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            PGam = P[IBlw] * gam_cdf_Thr / Tau[IBlw]
            PGam = np.clip(PGam, 0, 1)
            X[IBlw] = MarginalModel.gaminv(PGam, Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            
        if np.any(~IBlw):
            PGP = (P[~IBlw] - Tau[~IBlw]) / (1 - Tau[~IBlw])
            PGP = np.clip(PGP, 0, 1)
            X[~IBlw] = MarginalModel.gpinv(PGP, Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            
        return X

    @staticmethod
    def gamgpinvsurvivor(Q, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau):
        """Inverse CDF using Survival Probability."""
        Q, Xi, Sgm, Thr, Alp, Bet, GmmLct, Tau = np.broadcast_arrays(
            np.asarray(Q), np.asarray(Xi), np.asarray(Sgm), np.asarray(Thr), 
            np.asarray(Alp), np.asarray(Bet), np.asarray(GmmLct), np.asarray(Tau)
        )
        
        P = 1 - Q
        X = np.full(Q.shape, np.nan)
        
        IBlw = P < Tau
        
        if np.any(IBlw):
            gam_cdf_Thr = MarginalModel.gamcdf(Thr[IBlw], Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            QGam = Q[IBlw] * gam_cdf_Thr / Tau[IBlw] + (Tau[IBlw] - gam_cdf_Thr) / Tau[IBlw]
            QGam = np.clip(QGam, 0, 1)
            X[IBlw] = MarginalModel.gaminvsurvivor(QGam, Alp[IBlw], Bet[IBlw], GmmLct[IBlw])
            
        if np.any(~IBlw):
            QGP = Q[~IBlw] / (1 - Tau[~IBlw])
            X[~IBlw] = MarginalModel.gpinvsurvivor(QGP, Xi[~IBlw], Sgm[~IBlw], Thr[~IBlw])
            
        return X
    # -------------------------------------------------------------------------
    # Plotting Methods
    # -------------------------------------------------------------------------

    def Plot(self):
        """Plot results of PPC marginal fit."""
        # Check/Create directory handled in init
        
        # Figure 1: Margin Check
        plt.figure(1)
        plt.clf()
        self.PlotMarginCheck()
        
        # Figure 2: Non Stationary Parameters (Scale, Gamma Shape, Gamma Scale)
        plt.figure(2)
        plt.clf()
        self.PlotNonStatParam()
        
        # Figure 3: GP Shape
        plt.figure(3)
        plt.clf()
        self.PlotGPShape()
        
        # Figure 4: Smoothness
        if not np.all(np.isnan(self.CVLackOfFit)) and self.Bn.nBin > 1:
            plt.figure(4)
            plt.clf()
            self.PlotSmoothness()
        
        # Figure 5 & 6: QQ Plots
        self.PlotQQ()
        
        # Figure 7: Threshold Stability
        if self.nBoot > 1 and len(np.unique(self.NEP)) > 1:
            plt.figure(7)
            plt.clf()
            self.PlotThresholdStability()
            
        # Figure 8: Return Value CDF
        plt.figure(8)
        plt.clf()
        self.PlotRV()

def PlotMarginCheck(self):
        """
        Plot of data on standard and uniform margin as check of the goodness of fit.
        """
        # Calculate median threshold over bootstraps
        # self.Thr is [nBin x nBoot]
        mThr = np.nanmedian(self.Thr, axis=1)
        
        # Determine threshold for every data point based on its bin
        # self.Bn.A contains bin indices for the data
        thr_per_point = mThr[self.Bn.A.astype(int)]
        
        # Identify exceedances
        IExc = (self.Y > thr_per_point).flatten()
        
        # Get margins for the first bootstrap (original data)
        YMrg, YUnif = self.Margins([0]) 
        YMrg = YMrg.flatten()
        YUnif = YUnif.flatten()
        
        # Determine number of subplots
        nSubPlt = 1
        if not np.all(np.isnan(YUnif)):
            nSubPlt = 3
            
        for iC in range(self.nCvr):
            # --- Subplot 1: Raw Data ---
            plt.subplot(self.nCvr, nSubPlt, 1 + nSubPlt*iC)
            plt.plot(self.X[IExc, iC], self.Y[IExc], 'k.')
            plt.plot(self.X[~IExc, iC], self.Y[~IExc], '.', color=[0.7, 0.7, 0.7])
            
            if iC == 0:
                plt.title(f'{self.RspLbl}: Raw data')
            
            # --- Calls to CovariateBinning Class ---
            # MATLAB: PlotBinEdge(obj.Bn, iC)
            self.Bn.PlotBinEdge(iC)
            
            # MATLAB: PlotParameter(obj.Bn, obj.GmmLct, iC, 'color', 'r', 'linewidth', 2)
            self.Bn.PlotParameter(self.GmmLct, iC, color='r', linewidth=2)
            
            # MATLAB: PlotParameter(obj.Bn, obj.Thr, iC, 'color', 'b', 'linewidth', 2)
            self.Bn.PlotParameter(self.Thr, iC, color='b', linewidth=2)
            # ---------------------------------------

            plt.xlabel(self.CvrLbl[iC])
            plt.ylabel(self.RspLbl)
            plt.axis('tight')
            plt.grid(True)
            
            if nSubPlt > 1:
                # --- Subplot 2: Uniform Margins ---
                plt.subplot(self.nCvr, nSubPlt, 2 + nSubPlt*iC)
                plt.plot(self.X[:, iC], YUnif, 'k.')
                if iC == 0:
                    plt.title(f'On uniform margins\nNEP = {self.NEP[0]:.2f}')
                
                plt.xlabel(self.CvrLbl[iC])
                plt.ylabel(f'{self.RspLbl}-uniform scale')
                
                # MATLAB: PlotBinEdge(obj.Bn, iC)
                self.Bn.PlotBinEdge(iC)
                
                plt.ylim([0, 1])
                plt.xlim([0, 360])
                plt.xticks(np.arange(0, 361, 45))
                plt.grid(True)
                
                # --- Subplot 3: Standard Margins ---
                plt.subplot(self.nCvr, nSubPlt, 3 + nSubPlt*iC)
                plt.plot(self.X[:, iC], YMrg, 'k.')
                if iC == 0:
                    plt.title(f'On {self.MarginType} margins\nNEP = {self.NEP[0]:.2f}')
                
                plt.xlabel(self.CvrLbl[iC])
                plt.ylabel(f'{self.RspLbl}-{self.MarginType} scale')
                
                # MATLAB: PlotBinEdge(obj.Bn, iC)
                self.Bn.PlotBinEdge(iC)

                plt.xlim([0, 360])
                plt.xticks(np.arange(0, 361, 45))
                plt.grid(True)

        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_1_DataTransform.png'))

def PlotNonStatParam(self):
        """
        Piecewise Sigma: fitted GP scale (black)
        """
        for iC in range(self.nCvr):
            # --- 1. GP Scale ---
            plt.subplot(self.nCvr, 3, iC*3 + 1)
            
            if self.Bn.nBin > 1:
                # MATLAB: PlotParameter(obj.Bn, obj.Scl, iC, 'color', 'k', 'linewidth', 2)
                self.Bn.PlotParameter(self.Scl, iC, color='k', linewidth=2)
                
                # MATLAB: PlotBinEdge(obj.Bn, iC)
                self.Bn.PlotBinEdge(iC)
                
                # MATLAB: PlotParameter(obj.Bn, obj.Scl, iC, 'color', 'k', 'linewidth', 2)
                # (Repeated in original code, likely to ensure overlay on top of grid/edges)
                self.Bn.PlotParameter(self.Scl, iC, color='k', linewidth=2)
                
                plt.xlabel(self.CvrLbl[iC])
                plt.ylabel(r'$\nu$') # LaTeX for 'nu'
            else:
                # Stationary histogram
                plt.hist(self.Scl.flatten(), density=True, facecolor='k', edgecolor='none')
                plt.xlabel(r'$\nu$')
                plt.ylabel('Empirical density')
            
            plt.title(f'{self.RspLbl}: GP scale')
            plt.grid(True)
            
            # --- 2. Gamma Shape ---
            plt.subplot(self.nCvr, 3, iC*3 + 2)
            
            if self.Bn.nBin > 1:
                # MATLAB: PlotParameter(obj.Bn, obj.Omg, iC, 'color', 'k', 'linewidth', 2)
                self.Bn.PlotParameter(self.Omg, iC, color='k', linewidth=2)
                
                # MATLAB: PlotBinEdge(obj.Bn, iC)
                self.Bn.PlotBinEdge(iC)
                
                # MATLAB: PlotParameter(obj.Bn, obj.Omg, iC, 'color', 'k', 'linewidth', 2)
                self.Bn.PlotParameter(self.Omg, iC, color='k', linewidth=2)
                
                plt.xlabel(self.CvrLbl[iC])
                plt.ylabel(r'$\omega$')
            else:
                plt.hist(self.Omg.flatten(), density=True, facecolor='k', edgecolor='none')
                plt.xlabel(r'$\omega$')
                plt.ylabel('Empirical density')
                
            plt.title(f'{self.RspLbl}: Gamma shape')
            plt.grid(True)
            
            # --- 3. Gamma Scale ---
            plt.subplot(self.nCvr, 3, iC*3 + 3)
            
            if self.Bn.nBin > 1:
                # MATLAB: PlotParameter(obj.Bn, obj.Kpp, iC, 'color', 'k', 'linewidth', 2)
                self.Bn.PlotParameter(self.Kpp, iC, color='k', linewidth=2)
                
                # MATLAB: PlotBinEdge(obj.Bn, iC)
                self.Bn.PlotBinEdge(iC)
                
                # MATLAB: PlotParameter(obj.Bn, obj.Kpp, iC, 'color', 'k', 'linewidth', 2)
                self.Bn.PlotParameter(self.Kpp, iC, color='k', linewidth=2)
                
                plt.xlabel(self.CvrLbl[iC])
                plt.ylabel(r'$\kappa$')
            else:
                plt.hist(self.Kpp.flatten(), density=True, facecolor='k', edgecolor='none')
                plt.xlabel(r'$\kappa$')
                plt.ylabel('Empirical density')
                
            plt.title(f'{self.RspLbl}: Gamma scale')
            plt.grid(True)
            
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_2_Parameters.png'))
        
def PlotGPShape(self):
        """Plot GP Shape."""
        plt.hist(self.Shp.flatten(), density=True, facecolor='k')
        plt.xlabel('xi')
        plt.ylabel('Empirical density')
        plt.title(f'{self.RspLbl}: GP shape')
        plt.grid(True)
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_3_ParametersShape.png'))

def PlotSmoothness(self):
        """Plot CV Smoothness results."""
        if self.nBoot > 1 and self.CVMth == 1:
            median_LOF = np.nanmedian(self.CVLackOfFit, axis=1)
            lb_LOF = np.nanquantile(self.CVLackOfFit, 0.025, axis=1)
            ub_LOF = np.nanquantile(self.CVLackOfFit, 0.975, axis=1)
            plt.plot(self.SmthSet, median_LOF, 'k-', linewidth=2)
            plt.plot(self.SmthSet, lb_LOF, 'k--', linewidth=2)
            plt.plot(self.SmthSet, ub_LOF, 'k--', linewidth=2)
        else:
            plt.plot(self.SmthSet, self.CVLackOfFit[:, 0], 'k-', linewidth=2)
            
        # Plot optimal vertical line
        med_opt = np.median(self.OptSmth)
        plt.axvline(med_opt, color='r', linestyle='--')
        
        plt.ylabel('Lack of fit')
        plt.xlabel('lambda')
        plt.xscale('log')
        plt.title(f'{self.RspLbl}: GP cross-validation lack of fit')
        plt.grid(True)
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_4_CV.png'))

def PlotQQ(self):
        """QQ plot of goodness of fit."""
        nQ = 100
        Q = np.linspace(np.min(self.Y), np.max(self.Y)*1.2, nQ)
        
        # Calculate fitted CDF for these Q values across bootstraps
        # Reshape Q for broadcasting: [nQ x 1] vs [1 x nBoot] parameters
        Q_col = Q[:, None]
        
        # This will return [nQ x nBin x nBoot] or similar depending on implementation
        # Assuming gamgpcdf broadcasts correctly
        C = self.gamgpcdf(Q_col, self.Shp.T, self.Scl, self.Thr, 
                          self.Omg, self.Kpp, self.GmmLct, self.NEP.T)
        # C is likely [nQ x nBin x nBoot]
        
        # Figure 5: Sector Goodness of Fit
        if self.Bn.nBin > 1:
            plt.figure(5)
            plt.clf()
            nPlt1 = int(np.ceil(np.sqrt(self.Bn.nBin)))
            nPlt2 = int(np.ceil(self.Bn.nBin / nPlt1))
            
            for iB in range(self.Bn.nBin):
                plt.subplot(nPlt2, nPlt1, iB + 1)
                I = (self.Bn.A == iB).flatten()
                
                if np.any(I):
                    nI = np.sum(I)
                    P = np.arange(nI) / nI
                    plt.plot(np.sort(self.Y[I]), np.log10(1 - P), 'r.')
                    
                    if self.nBoot > 1:
                        # Quantiles over bootstraps (axis 2 of C)
                        # C slice: [nQ x nBin x nBoot] -> [nQ x nBoot]
                        C_bin = C[:, iB, :]
                        qC = np.quantile(C_bin, [0.025, 0.5, 0.975], axis=1).T
                        plt.plot(Q, np.log10(1 - qC[:, 1]), 'k-') # Median
                        plt.plot(Q, np.log10(1 - qC[:, 0]), 'k--')
                        plt.plot(Q, np.log10(1 - qC[:, 2]), 'k--')
                    else:
                        plt.plot(Q, np.log10(1 - C[:, iB, 0]), 'k-')
                        
                    plt.xlim([np.min(self.Y[I]), np.max(self.Y[I])])
                    plt.ylabel('log(1-p)')
                    plt.grid(True)
            plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_5_SectorGoodnessOfFit.png'))

        # Figure 6: Overall Goodness of Fit
        plt.figure(6)
        plt.clf()
        P_all = np.arange(self.nDat) / self.nDat
        plt.plot(np.sort(self.Y.flatten()), np.log10(1 - P_all), 'r.')
        
        # Calculate Omni CDF (weighted sum of bin CDFs)
        # Rat is [nBin x nBoot]
        w = self.Rat / np.sum(self.Rat, axis=0) # Weights
        
        # C is [nQ x nBin x nBoot]
        # Weighted sum over bins (axis 1)
        # Reshape w to [1 x nBin x nBoot]
        w_reshaped = w[None, :, :]
        COmni = np.sum(C * w_reshaped, axis=1) # [nQ x nBoot]
        
        COmni = np.clip(COmni, 0, 1)
        
        if self.nBoot > 1:
            qCOmni = np.quantile(COmni, [0.025, 0.5, 0.975], axis=1).T
            plt.plot(Q, np.log10(1 - qCOmni[:, 1]), 'k-')
            plt.plot(Q, np.log10(1 - qCOmni[:, 0]), 'k--')
            plt.plot(Q, np.log10(1 - qCOmni[:, 2]), 'k--')
        else:
            plt.plot(Q, np.log10(1 - COmni[:, 0]), 'k-')
            
        plt.ylabel('log(1-p)')
        plt.title(f'{self.RspLbl}: Overall goodness of fit')
        plt.grid(True)
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_6_OverallGoodnessOfFit.png'))

def PlotThresholdStability(self):
        """Plot threshold stability."""
        # Simple plot of Shp vs NEP
        plt.plot(self.NEP, self.Shp, 'k.', markersize=10)
        plt.xlabel('NEP')
        plt.ylabel('GP shape xi')
        plt.grid(True)
        plt.title(f'{self.RspLbl}: GP shape stability by threshold')
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_7_ThresholdStability.png'))

def PlotRV(self):
        """Plot Return Value CDFs."""
        # Using HSV colormap
        colors = plt.cm.hsv(np.linspace(0, 1, self.Bn.nBin))
        
        for iRtr in range(self.nRtr):
            plt.subplot(1, self.nRtr, iRtr + 1)
            
            if self.Bn.nBin > 1:
                for iBin in range(self.Bn.nBin):
                    # Sort simulations
                    tX = np.sort(self.RVSml['Org'][iBin, :, iRtr])
                    # Remove NaNs
                    tX = tX[~np.isnan(tX)]
                    plt.plot(tX, np.linspace(0, 1, len(tX)), color=colors[iBin], linewidth=2)
            
            # Omni (Max over bins)
            # RVSml['Org'] is [nBin x nRls x nRtr]
            # Max over axis 0
            XOmni = np.nanmax(self.RVSml['Org'][:, :, iRtr], axis=0)
            XOmni = np.sort(XOmni)
            plt.plot(XOmni, np.linspace(0, 1, len(XOmni)), 'k', linewidth=2)
            
            # Lines for 0.5 and 1/e
            plt.axhline(0.5, color='k', linestyle='--')
            plt.axhline(np.exp(-1), color='k', linestyle='--')
            
            plt.ylabel('Cumulative probability')
            plt.xlabel(self.RspLbl)
            plt.title(f'Maximum over {self.RtrPrd[iRtr]} years')
            plt.grid(True)
            
        plt.savefig(os.path.join(self.FigureFolder, f'Stg3_{self.RspSavLbl}_8_ReturnValueCDF.png'))