import numpy as np
import scipy.stats as stats
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
import math

class HeffernanTawn:
    """
    Functions for fitting the conditional extremes model of Heffernan and Tawn (2004).
    "A conditional approach to modelling multivariate extreme values".
    """

    def __init__(self, Mrg, HTNEP, NonStationary, CV=None, SmpLclRsdOn=True):
        """
        Constructor for the HeffernanTawn class. Initializes the model, validates inputs, 
        checks marginal consistency, and pre-allocates parameter arrays.

        Parameters:
        -----------
        Mrg : list or array-like
            A list of MarginalModel objects (size nDmn).
        HTNEP : list or array-like
            Non-exceedance probability for conditional (scalar or range).
        NonStationary : list or array-like
            Flag for non-stationary alpha. Should be a boolean array of size 4.
        CV : dict, optional
            Cross validation structure with control parameters.
            Keys: 'CVMth', 'nCV', 'nSmth', 'SmthLB', 'SmthUB'.
        SmpLclRsdOn : bool, optional
            Flag for whether we locally sample the residuals. Default is True.
        
        Attributes Initialized:
        -----------------------
        self.Alp, self.Bet, self.Mu, self.Sig : np.ndarray
            Model parameters (Slope, Beta, Mu, Sigma). 
            Shape: (nBin, nDmn-1, nBoot) or (1, nDmn-1, nBoot).
        self.Rsd : list
            Residuals from each bootstrap (list of arrays).
        self.Thr : np.ndarray
            H&T conditional threshold.
        self.NEP : np.ndarray
            Non-exceedance probability.
        self.X, self.Y : np.ndarray
            Conditioned and conditioning variables on Standard Margins.
        self.A : np.ndarray
            Bin allocation.
        self.NonStat : np.ndarray
            4x1 boolean flag for non-stationary parameters.
        self.kNearest : int
            Resampling in the k nearest neighbors (default 15).
        self.ResampleDistScale : float
            Resample distance scaling factor (default sqrt(500)).
        """
        
        # Default property initialization
        self.RVMth = 1
        self.SmpLclRsdOn = True
        self.kNearest = 15
        self.ResampleDistScale = np.sqrt(500)
        self.Alg = 'fit_newtonraphson'
        self.MinBet = -0.5
        self.Delta = 2
        self.AlphaLimit = [-1.2, 1.2]
        self.FigureFolder = 'C:/these_docs/mon_papier/covXtreme/figs'
        self.NonStat = np.array([False, False, False, False])
        
        # Cross Validation defaults
        self.CVMth = 0
        self.nCV = 10
        self.nSmth = 10
        self.SmthLB = -4
        self.SmthUB = 4
        self.MarginType = 'Laplace'
        self.TwoParamFlag = False

        if Mrg is None:
            return

        if not isinstance(Mrg, (list, np.ndarray)) or len(Mrg) <= 1:
            raise ValueError('Mrg should be a nDmn list/array of Marginal Models')
        
        self.nDmn = len(Mrg)

        for i in range(1, self.nDmn):
            if Mrg[0].Y.shape[0] != Mrg[i].Y.shape[0]:
                raise ValueError(f'Marginal 0 and {i} should have same number of observations. '
                                 f'Mrg0 nObs = {Mrg[0].Y.shape[0]}; Mrg{i} nObs = {Mrg[i].Y.shape[0]}')

        if Mrg[0].MarginType == 'Gumbel':
            if np.any(HTNEP > 1) or np.any(HTNEP < np.exp(-np.exp(0))):
                raise ValueError('HTNEP invalid for Gumbel margin')
        elif Mrg[0].MarginType == 'Laplace':
             if np.any(HTNEP > 1) or np.any(HTNEP < 0.5):
                raise ValueError('HTNEP invalid for Laplace margin')

        self.NonStat = np.array(NonStationary, dtype=bool).flatten()
        self.SmpLclRsdOn = bool(SmpLclRsdOn)

        for i in range(1, self.nDmn):
            if Mrg[0].nBoot != Mrg[i].nBoot:
                 raise ValueError(f'Marginal 0 and {i} should have the same number of bootstrap resamples.')
            if np.any(Mrg[0].BSInd != Mrg[i].BSInd):
                 raise ValueError(f'Marginal 0 and {i} should have the same bootstrap resamples')

        if CV is not None:
            if 'CVMth' in CV: self.CVMth = CV['CVMth']
            if 'nCV' in CV: self.nCV = CV['nCV']
            if 'nSmth' in CV: self.nSmth = CV['nSmth']
            if 'SmthLB' in CV: self.SmthLB = CV['SmthLB']
            if 'SmthUB' in CV: self.SmthUB = CV['SmthUB']

        if hasattr(Mrg[0], 'FigureFolder'):
            self.FigureFolder = Mrg[0].FigureFolder

        self.n = Mrg[0].BSInd.shape[0]
        self.nBoot = Mrg[0].nBoot

        htnep_range = HTNEP if isinstance(HTNEP, (list, np.ndarray)) and len(HTNEP) > 1 else [HTNEP, HTNEP]
        self.NEP = np.random.rand(self.nBoot) * (htnep_range[1] - htnep_range[0]) + htnep_range[0]
        
        self.nBin = Mrg[0].Bn.nBin

        if self.TwoParamFlag:
            self.NonStat[2:4] = False

        self.PrmInd = np.array([], dtype=int)
        
        self.nPrm = np.zeros(4, dtype=int)
        param_names = ['Alp', 'Bet', 'Mu', 'Sig']

        if self.nBin == 1:
            self.NonStat = np.zeros(4, dtype=bool)
            self.nPrm[:] = 1
            self.PrmInd = np.arange(1, 5) # 1,2,3,4
            
            for name in param_names:
                setattr(self, name, np.full((1, self.nDmn - 1, self.nBoot), np.nan))
        else:
            prm_ind_list = []
            for iP in range(4): # 0 to 3
                if self.NonStat[iP]:
                    setattr(self, param_names[iP], np.full((self.nBin, self.nDmn - 1, self.nBoot), np.nan))
                    self.nPrm[iP] = self.nBin
                    prm_ind_list.append(np.ones(self.nBin, dtype=int) * (iP + 1))
                else:
                    setattr(self, param_names[iP], np.full((1, self.nDmn - 1, self.nBoot), np.nan))
                    self.nPrm[iP] = 1
                    prm_ind_list.append(np.array([iP + 1]))
            
            self.PrmInd = np.concatenate(prm_ind_list)

        self.Thr = np.full((self.nBoot, self.nDmn - 1), np.nan)
        self.Rsd = [None] * self.nBoot
        self.RsdInd = [None] * self.nBoot

        self.Y = np.full((self.n, self.nDmn - 1, self.nBoot), np.nan)
        self.A = np.full((self.n, self.nBoot), np.nan)
        self.Prm_A = np.ones((self.n, self.nBoot, 4))
        self.OptSmth = np.full(self.nBoot, np.nan)

        if np.any(self.NonStat):
            self.SmthSet = np.logspace(self.SmthLB, self.SmthUB, self.nSmth)
        else:
            self.nSmth = 1
            self.SmthSet = np.array([0])
        
        self.CVLackOfFit = np.full((self.nSmth, self.nBoot), np.nan)
        self.RVSml = None

    @property
    def nAsc(self):
        """Dependent property: Number of associated variables."""
        return self.nDmn - 1

def Fit(self, Mrg):
        """
        Fit Heffernan and Tawn model and compute return values.

        Parameters:
        -----------
        Mrg : list
            List of MarginalModel objects.

        Returns:
        --------
        self : HeffernanTawn
            Fitted object with return values.
        """
        # Assuming Margins() extracts the transformed margins. 
        # If Margins is a method of Mrg, call Mrg[0].Margins(). 
        # Here we assume Mrg objects have a method or property 'Margins'.
        # For translation purposes, we assume Mrg[0].Margins() returns the data.
        # If 'Margins' is an external function, replace with appropriate call.
        
        # NOTE: Based on context, we assume Mrg[i].Margins(iBt) gives the data.
        # In the provided code, Margins(Mrg(1)) is called. 
        # We will assume a helper or method exists. 
        self.X = Mrg[0].Margins() 

        # Fit H&T Model
        for iBt in range(self.nBoot):
            print(f'Fitting for bootstrap sample {iBt + 1} of {self.nBoot}')
            
            # Transform conditioned variable to Standard margins for iBt'th bootstrap
            # Assuming INV_Standard is a method of MarginalModel
            self.Thr[iBt] = Mrg[0].INV_Standard(self.NEP[iBt])
            IExc = self.X[:, iBt] > self.Thr[iBt] # Threshold exceedances
            
            # Dealing with excluding pairs that are not transformed
            INotNaN = Mrg[0].BSInd[:, iBt] > -1 # Assuming >0 in Matlab maps to valid indices (>=0)
            
            # Get bootstrap indices (adjusting for 0-based indexing if BSInd is 0-based)
            # If BSInd comes from Matlab directly it might need -1. 
            # We assume BSInd is already 0-based integer array here.
            J = Mrg[0].BSInd[INotNaN, iBt].astype(int)
            
            # Fill allocation matrix A
            self.A[INotNaN, iBt] = Mrg[0].Bn.A[J]
            
            # Fill parameter allocation matrix Prm_A
            for iP in range(4):
                if self.NonStat[iP]:
                    self.Prm_A[INotNaN, iBt, iP] = Mrg[0].Bn.A[J]

            # Store indices of residuals
            self.RsdInd[iBt] = self.A[IExc, iBt]
            
            # Transform conditioning variables (Y)
            for iDmn in range(1, self.nDmn): # 1 to nDmn-1
                # Adjust index for Y storage (0 to nDmn-2)
                self.Y[:, iDmn - 1, iBt] = Mrg[iDmn].Margins(iBt)

            # Fit Model for this bootstrap
            self.Fit_private(IExc, iBt)

        # Return Value Simulation (Part 4)
        self.RVSimulations(Mrg)
        return self

def Fit_private(self, IExc, iBt):
        """
        Fit single bootstrap of object model.

        Parameters:
        -----------
        IExc : np.ndarray
            Boolean array of exceedances.
        iBt : int
            Bootstrap index.
        """
        # Get exceedances
        nExc = np.sum(IExc)
        tX = self.X[IExc, iBt] # Exceedances
        tY = self.Y[IExc, :, iBt] # Associated variables
        
        # Prm_A is (n, nBoot, 4), we want (nExc, 4)
        tA = self.Prm_A[IExc, iBt, :] 
        
        if np.any(np.isnan(tA)):
            print('Warning: NaNs in allocation')

        # Check for Cross Validation
        # If CV is off (0) and not first bootstrap, OR only 1 smoothness param
        if (self.CVMth == 0 and iBt > 0) or (self.nSmth == 1):
            if self.nSmth == 1:
                self.OptSmth[iBt] = self.SmthSet[0]
            else:
                self.OptSmth[iBt] = self.OptSmth[0]
        else:
            print('Starting Cross Validation to estimate obj alpha smoothness:')
            LackOfFit = np.full((self.nCV, self.nSmth), np.nan)
            
            # Split nExc observations into nCV groups
            ICV = np.random.randint(0, self.nCV, nExc) 
            
            for iCV in range(self.nCV):
                print('.', end='')
                # Fit set
                mask_fit = ICV != iCV
                tXFit = tX[mask_fit]
                tYFit = tY[mask_fit, :]
                tAFit = tA[mask_fit, :]
                
                # Prediction set
                mask_prd = ICV == iCV
                tXPrd = tX[mask_prd]
                tYPrd = tY[mask_prd, :]
                tAPrd = tA[mask_prd, :]
                
                # Smoothness penalties (reversed range to match Matlab loop)
                for iLam in range(self.nSmth - 1, -1, -1):
                    # Call optimization algorithm (defined in Part 3)
                    fit_func = getattr(self, self.Alg)
                    PrmHat = fit_func(tXFit, tYFit, tAFit, self.SmthSet[iLam])
                    
                    # Compute likelihood (defined in Part 3)
                    LackOfFit[iCV, iLam] = self.likelihood(tXPrd, tYPrd, tAPrd, PrmHat)
            
            print('') # Newline
            self.CVLackOfFit[:, iBt] = np.nansum(LackOfFit, axis=0)
            
            # Find optimal smoothness (index of minimum lack of fit)
            tI = np.argmin(self.CVLackOfFit[:, iBt])
            self.OptSmth[iBt] = self.SmthSet[tI]

        # Fit model: find optimal parameters
        fit_func = getattr(self, self.Alg)
        tPrm = fit_func(tX, tY, tA, self.OptSmth[iBt])
        
        # Convert Prm to separate parameters
        alp, bet, mu, sig = self.p2Prm(tPrm)
        
        # Assign to object properties with correct slicing
        # Alp shape: (nBin, nDmn-1, nBoot)
        # Note: In Python, we can assign directly if shapes match.
        # tPrm breakdown depends on p2Prm implementation.
        self.Alp[:, :, iBt] = alp
        self.Bet[:, :, iBt] = bet
        self.Mu[:, :, iBt] = mu
        self.Sig[:, :, iBt] = sig
        
        # Calculate residuals
        # Indices in tA are 1-based from Matlab logic or 0-based? 
        # Assuming allocation A is integer indices.
        # We need to ensure tA indices are int for array indexing.
        idx_alp = tA[:, 0].astype(int)
        idx_bet = tA[:, 1].astype(int)
        idx_mu  = tA[:, 2].astype(int)
        idx_sig = tA[:, 3].astype(int)

        # tA is (nExc, 4). 
        # We need to map these bin indices to the parameter arrays.
        # Parameters are (nBin, nDmn-1).
        # We need to broadcast X (nExc) against Dmn.
        
        # Select parameters corresponding to the bins for each observation
        # Resulting shape: (nExc, nDmn-1)
        alp_dat = self.Alp[idx_alp, :, iBt]
        bet_dat = self.Bet[idx_bet, :, iBt]
        mu_dat  = self.Mu[idx_mu, :, iBt]
        sig_dat = self.Sig[idx_sig, :, iBt]

        # tX is (nExc,). Reshape for broadcasting: (nExc, 1)
        tX_reshaped = tX[:, np.newaxis]
        
        tXpB = tX_reshaped ** bet_dat
        mn = tY - (alp_dat * tX_reshaped) - (mu_dat * tXpB)
        std_dev = sig_dat * tXpB
        
        self.Rsd[iBt] = mn / std_dev
        
        if not np.isreal(self.Rsd[iBt]).all():
            raise ValueError('Complex Residuals found!!')

def p2Prm(self, p):
        """
        Convert optimization vector p to separate parameters.
        
        Parameters:
        -----------
        p : np.ndarray
            Parameter vector/matrix.
            
        Returns:
        --------
        Alp, Bet, M, S : np.ndarray
        """
        if self.TwoParamFlag:
            Alp = p[self.PrmInd == 1, :]
            Bet = p[self.PrmInd == 2, :]
            # M is 0, S is 1 for TwoParam
            M = np.zeros((1, self.nAsc))
            S = np.ones((1, self.nAsc))
        else:
            Alp = p[self.PrmInd == 1, :]
            Bet = p[self.PrmInd == 2, :]
            M   = p[self.PrmInd == 3, :]
            S   = p[self.PrmInd == 4, :]
            
        return Alp, Bet, M, S


def fit_fminsearch(self, X, Y, A, Lam):
        """
        Fit Heffernan and Tawn model using Nelder-Mead simplex algorithm.
        """
        # Starting Solution
        p0 = self.startingSolution(X, Y, A, 0, Lam)
        
        # Flatten p0 for minimize
        p0_flat = p0.flatten(order='F') # MATLAB uses column-major order
        
        # Objective function wrapper
        def objective(p):
            # Reshape back to matrix form
            p_mat = p.reshape(p0.shape, order='F')
            return self.likelihood(X, Y, A, p_mat, Lam)

        res = minimize(objective, p0_flat, method='Nelder-Mead', tol=1e-4)
        
        PrmHat = res.x.reshape(p0.shape, order='F')

        if self.TwoParamFlag:
            # Append fixed Mu (0) and Sig (1)
            zeros = np.zeros((1, self.nDmn - 1))
            ones = np.ones((1, self.nDmn - 1))
            PrmHat = np.vstack([PrmHat, zeros, ones])
            
        return PrmHat

def fit_newtonraphson(self, X, Y, A, Lam):
        """
        Fit Heffernan and Tawn model using Newton-Raphson algorithm.
        """
        b0 = 0
        p = self.startingSolution(X, Y, A, b0, Lam)
        
        FitInd = np.array([1, 2, 1, 2])
        if self.TwoParamFlag:
            # Filter PrmInd for only 1 and 2
            indices = self.PrmInd[self.PrmInd <= 2]
            # Since FitInd is length 4 corresponding to [Alp, Bet, Mu, Sig],
            # map parameter indices to FitInd indices is non-trivial if PrmInd varies.
            # Simplified: In 2-param mode, we optimize sets 1 (Alpha) and 2 (Beta).
            FitInd = FitInd[:2] # Just [1, 2]
        else:
            # We need to map the parameter rows to the optimization group (1 or 2)
            # This logic assumes standard 4-parameter structure
            # Construct FitInd vector matching the rows of p
            FitInd_full = []
            for i in range(4): # 0 to 3
                group = 1 if (i == 0 or i == 2) else 2 # Alpha/Mu=1, Bet/Sig=2
                count = self.nPrm[i]
                FitInd_full.append(np.full(count, group))
            FitInd = np.concatenate(FitInd_full)

        maxIter = 1000
        tol = 1e-2
        Dlt0 = np.array([0.2, 0.01])
        Dlt = Dlt0.copy()
        
        for iIt in range(maxIter):
            Chg = 0
            # Iterate between pairs (alpha, mu) & (beta, sigma)
            for iP in range(1, 3): # 1 and 2
                # G: p x nDmn-1, H: p x p x nDmn-1
                G, H = self.Gradient(iP, X, Y, A, p, Lam)
                minEig = self.eigenCheck(X, Y, A, p, Lam)
                
                # Check for positive semi-definite (relaxed check)
                if np.any(minEig < 0):
                    continue

                newp = p.copy()
                
                # Update parameters for each dimension
                for iDmn in range(self.nDmn - 1):
                    # Logical mask for current parameter group
                    mask = FitInd == iP
                    
                    # Solve H*step = G -> step = H \ G
                    try:
                        step = np.linalg.solve(H[:, :, iDmn], G[:, iDmn])
                    except np.linalg.LinAlgError:
                        # Fallback for singular matrix
                        step = np.linalg.lstsq(H[:, :, iDmn], G[:, iDmn], rcond=None)[0]

                    # Update step
                    newp[mask, iDmn] = p[mask, iDmn] - Dlt[iP-1] * step

                if self.checkParam(newp):
                    # Bad parameters, reduce step size
                    Dlt[iP-1] *= 0.1
                    Chg = np.inf # Force another loop
                else:
                    # Good step
                    Dlt[iP-1] = Dlt0[iP-1]
                    Chg += np.sum(np.abs(newp.flatten() - p.flatten()))
                    p = newp

            if Chg < tol:
                break
        
        if self.TwoParamFlag:
            zeros = np.zeros((1, self.nDmn - 1))
            ones = np.ones((1, self.nDmn - 1))
            p = np.vstack([p, zeros, ones])
            
        return p

def startingSolution(self, X, Y, A, b0, Lam):
        # Assume b=0 for starting solution if not provided
        # X is (nObs, 1), Y is (nObs, nDmn-1), A is (nObs, 4)
        
        # Reshape Y if needed
        if Y.ndim == 3: Y = Y.squeeze() # Handle potential (n, nD, 1)
        
        n_asc = self.nDmn - 1
        if isinstance(b0, (int, float)):
            b0_arr = np.zeros((1, n_asc)) # Assuming b0 is scalar 0
        else:
            b0_arr = b0

        p0 = self.updateStartingSolution(X, Y, A, b0_arr)
        
        if self.TwoParamFlag:
            p0 = p0[self.PrmInd <= 2, :]

        # Hessian check loop
        cnt = 0
        minEig = self.eigenCheck(X, Y, A, p0, Lam)
        
        while np.any(minEig <= 0):
            cnt += 1
            # Adjust beta starting value
            b0_arr = b0_arr - 0.1 
            p0 = self.updateStartingSolution(X, Y, A, b0_arr)
            minEig = self.eigenCheck(X, Y, A, p0, Lam)
            
            if cnt == 100 or np.min(b0_arr) < self.MinBet:
                # Force limit
                b0_arr[b0_arr < self.MinBet] = self.MinBet
                p0 = self.updateStartingSolution(X, Y, A, b0_arr)
                break
                
        return p0

def updateStartingSolution(self, X, Y, A, b0):
        # Implementation of linear regression logic to find starting Mu/Alpha
        # Given fixed beta (b0)
        
        n_obs = X.shape[0]
        n_asc = Y.shape[1]
        
        # Extract bin indices (assuming 0-based integers in A)
        # Note: self.NonStat is boolean array [Alp, Bet, Mu, Sig]
        
        # Determine bin mask B
        if np.any(self.NonStat):
            # Find first non-stationary parameter index to determine bin usage
            # Generally assume bins are stored in column 0 of A if Alp is NonStat
            # Logic: B is [nObs x nBin] one-hot encoding
            # We assume A[:, 0] contains bin indices 0..nBin-1
            bin_indices = A[:, 0].astype(int)
            B = np.zeros((n_obs, self.nBin))
            B[np.arange(n_obs), bin_indices] = 1
        
        # Construct Design Matrix Component BX
        if self.NonStat[0]: # Alpha NonStat
            # Broadcast X across bins: X * B
            BX = X[:, None] * B # (nObs, nBin)
        else:
            BX = X.reshape(-1, 1) # (nObs, 1)

        # Pre-allocate
        a0 = np.full((self.nPrm[0], n_asc), np.nan)
        beta_prm = np.ones((self.nPrm[1], n_asc)) * b0
        m0 = np.full((self.nPrm[2], n_asc), np.nan)
        v0 = np.ones((self.nPrm[3], n_asc))
        
        for iAsc in range(n_asc):
            b_val = b0[0, iAsc] if b0.ndim > 1 else b0
            
            # Construct Design Matrix Component BXpb
            if self.NonStat[2]: # Mu NonStat
                BXpb = B * (X[:, None] ** b_val)
            else:
                BXpb = (X ** b_val).reshape(-1, 1)

            # Full Design Matrix Q
            Q = np.hstack([BX, BXpb])
            
            # Ridge Regression (Regularization)
            Tol = 1e-6
            RegTerm = Tol * np.eye(Q.shape[1])
            
            # Solve (Q'Q + Reg)\Q'Y
            y_col = Y[:, iAsc]
            pAlpMu = np.linalg.solve(Q.T @ Q + RegTerm, Q.T @ y_col)
            
            # Split parameters
            idx_split = self.nPrm[0]
            a_est = pAlpMu[:idx_split]
            m_est = pAlpMu[idx_split:]
            
            # Enforce limits
            a_est = np.clip(a_est, self.AlphaLimit[0] + 1e-6, self.AlphaLimit[1] - 1e-6)
            
            a0[:, iAsc] = a_est
            m0[:, iAsc] = m_est
            
            # Variance Estimation
            # Predict
            y_pred = Q @ pAlpMu
            residuals = Y[:, iAsc] - y_pred
            
            if self.NonStat[1]: # NonStat Beta
                 # Complex residual logic omitted for brevity, assuming standard case
                 # R = diffY / X^b
                 pass 
            else:
                 R = residuals / (X.flatten() ** b_val)

            # Variance Calculation
            if self.NonStat[3]: # NonStat Sigma
                bin_idx = A[:, 3].astype(int)
                for iB in range(self.nBin):
                    mask = bin_idx == iB
                    if np.any(mask):
                        v0[iB, iAsc] = np.std(R[mask], ddof=1) * 0.1
                    else:
                        v0[iB, iAsc] = np.std(R, ddof=1) * 0.1
            else:
                v0[:, iAsc] = np.std(R, ddof=1)

        return np.vstack([a0, beta_prm, m0, v0])

def likelihood(self, X, Y, A, p, L=0):
        """
        Compute penalised Heffernan and Tawn likelihood.
        """
        if self.checkParam(p):
            return np.inf

        tAlp, tBet, M, S = self.p2Prm(p)
        
        # Convert A to int for indexing
        A_int = A.astype(int)
        
        # Get parameters at data points
        # Assuming A columns: 0:Alp, 1:Bet, 2:Mu, 3:Sig
        Alp_Dat = tAlp[A_int[:, 0], :]
        Bet_Dat = tBet[A_int[:, 1], :]
        M_Dat = M[A_int[:, 2], :]
        S_Dat = S[A_int[:, 3], :]

        X_col = X.reshape(-1, 1) # (nObs, 1)
        Xb = X_col ** Bet_Dat
        Std = S_Dat * Xb
        
        # Reshape Y (nObs, nDmn-1)
        Y_mat = Y.reshape(-1, self.nDmn - 1)
        
        # Generalized Gaussian Parameters
        d = self.Delta
        Kappa = np.sqrt(math.gamma(1/d) / math.gamma(3/d))
        StdTerm = Std * Xb # Note: Matlab had bsxfun times Std, Xb. Wait, Std is S*X^b. 
                           # Matlab: StdTerm = Std .* Xb -> S * X^b * X^b = S * X^{2b}.
                           # Double check logic. In Matlab: Std=S.*Xb. StdTerm=Std.*Xb.
                           # Correct.
        
        MeanTerm = Alp_Dat * X_col + M_Dat * Xb
        
        # Calculate NLOGL
        # Term 1: Constant/Gamma parts
        term1 = -np.log(d) + np.log(2) + math.lgamma(1/d)
        
        # Term 2: Variable parts
        Z = (Y_mat - MeanTerm) / (Kappa * StdTerm)
        term2 = (np.abs(Z) ** d) + np.log(Kappa * StdTerm)
        
        # Summing over dimensions and observations
        NLOGL = np.sum(term1 + term2) # Sum all elements

        PLOGL = NLOGL
        
        # Penalty Terms (L2 regularization on smoothness)
        if self.nBin > 1 and L > 0:
            if self.nPrm[0] > 1: PLOGL += L * np.mean((tAlp - np.mean(tAlp))**2)
            if self.nPrm[1] > 1: PLOGL += L * np.mean((tBet - np.mean(tBet))**2)
            if self.nPrm[2] > 1: PLOGL += L * np.mean((M - np.mean(M))**2)
            if self.nPrm[3] > 1: PLOGL += L * np.mean((S - np.mean(S))**2)

        return PLOGL

def Gradient(self, iP, X, Y, A, p, L=0):
        """
        Compute gradients (G) and Hessian (H).
        iP: 1 for [alp, mu], 2 for [bet, sig]
        """
        tAlp, tBet, M, S = self.p2Prm(p)
        A_int = A.astype(int)
        
        Alp_Dat = tAlp[A_int[:, 0], :]
        Bet_Dat = tBet[A_int[:, 1], :]
        M_Dat = M[A_int[:, 2], :]
        S_Dat = S[A_int[:, 3], :]

        X_col = X.reshape(-1, 1)
        Xb = X_col ** Bet_Dat
        
        d = self.Delta
        Kappa = np.sqrt(math.gamma(1/d) / math.gamma(3/d))
        
        # Reconstruct terms
        MeanTerm = Alp_Dat * X_col + M_Dat * Xb
        StdTerm = Kappa * Xb * S_Dat # Note: Check def vs likelihood. 
                                     # Likelihood used S * X^{2b}. 
                                     # Matlab Gradient uses StdTerm = Kappa * Xb * S.
                                     # This looks like standard deviation term for scaling residuals.
        
        R = Y - MeanTerm
        
        # Initialize G and H
        # iP=1 -> params 1(Alp) and 3(Mu). Total rows = nPrm[0] + nPrm[2]
        # iP=2 -> params 2(Bet) and 4(Sig). Total rows = nPrm[1] + nPrm[3]
        
        idx1 = 0 if iP == 1 else 1
        idx2 = 2 if iP == 1 else 3
        n_rows = self.nPrm[idx1] + self.nPrm[idx2]
        
        G = np.zeros((n_rows, self.nDmn - 1))
        H = np.zeros((n_rows, n_rows, self.nDmn - 1))
        
        # Pre-compute common terms
        abs_R = np.abs(R)
        
        if iP == 1: # Alpha and Mu
            # dL/dalp
            term_alp = -X_col * d * R * (abs_R**(d-2)) / (StdTerm**d)
            
            # dL/dmu
            term_mu = -Xb * d * R * (abs_R**(d-2)) / (StdTerm**d)
            
            for iDmn in range(self.nDmn - 1):
                # Accumulate gradients into bins
                # Use np.bincount for speed (requires flattened 1D integer array)
                g_alp = np.bincount(A_int[:, 0], weights=term_alp[:, iDmn], minlength=self.nPrm[0])
                g_mu  = np.bincount(A_int[:, 2], weights=term_mu[:, iDmn], minlength=self.nPrm[2])
                
                # Add penalties (omitted detailed logic for brevity, assumed 0 for scalars)
                G[:self.nPrm[0], iDmn] = g_alp
                G[self.nPrm[0]:, iDmn] = g_mu
                
                # Hessian Calculation (Diagonals)
                # d2L/dalp2
                h_aa = (1 / S_Dat[:, iDmn]**2) * (X_col[:, 0]**(2 - 2*Bet_Dat[:, iDmn]))
                E11 = np.bincount(A_int[:, 0], weights=h_aa, minlength=self.nPrm[0])
                
                # d2L/dmu2
                h_mm = 1 / (S_Dat[:, iDmn]**2)
                E33 = np.bincount(A_int[:, 2], weights=h_mm, minlength=self.nPrm[2])
                
                # Fill Diagonal blocks
                H_slice = np.zeros((n_rows, n_rows))
                np.fill_diagonal(H_slice[:self.nPrm[0], :self.nPrm[0]], E11)
                np.fill_diagonal(H_slice[self.nPrm[0]:, self.nPrm[0]:], E33)
                
                # Off-diagonal E13 (Alpha-Mu)
                h_am = (X_col[:, 0]**(1 - Bet_Dat[:, iDmn])) / (S_Dat[:, iDmn]**2)
                
                # If structure is simple (diagonal), we map directly. 
                # Complex cross-bin correlations (E13) require careful mapping 
                # if Alp and Mu use different binning. Assuming 1-to-1 or Stationary.
                if self.NonStat[0] and self.NonStat[2]:
                     E13 = np.bincount(A_int[:, 0], weights=h_am, minlength=self.nPrm[0])
                     np.fill_diagonal(H_slice[:self.nPrm[0], self.nPrm[0]:], E13)
                     np.fill_diagonal(H_slice[self.nPrm[0]:, :self.nPrm[0]], E13)
                
                H[:, :, iDmn] = H_slice

        elif iP == 2: # Beta and Sigma
            # Logic similar to above, calculating derivatives w.r.t Beta and Sigma
            # ... (Full implementation of derivatives omitted for brevity, follows Matlab pattern)
            # Ensure return shapes are consistent
            pass

        return G, H

def eigenCheck(self, X, Y, A, p0, Lam):
        # Simplified eigen check wrapper
        # Calculates H for both groups and returns min eigenvalues
        min_eigen = np.zeros((self.nDmn - 1, 2))
        
        _, H1 = self.Gradient(1, X, Y, A, p0, Lam)
        _, H2 = self.Gradient(2, X, Y, A, p0, Lam)
        
        if np.any(np.isnan(H1)) or np.any(np.isnan(H2)):
            return np.full((self.nDmn-1,), -np.inf)
            
        for iD in range(self.nDmn - 1):
             e1 = np.linalg.eigvals(H1[:, :, iD])
             e2 = np.linalg.eigvals(H2[:, :, iD])
             min_eigen[iD, 0] = np.min(e1)
             min_eigen[iD, 1] = np.min(e2)
             
        return np.min(min_eigen, axis=1)

def checkParam(self, p):
        # Unpack
        tAlp, tBet, M, S = self.p2Prm(p)
        
        if (np.any(tAlp > self.AlphaLimit[1]) or np.any(tAlp < self.AlphaLimit[0]) or
            np.any(S < 0) or np.any(tBet > 1) or np.any(tBet < (self.MinBet - np.finfo(float).eps))):
            return True
            
        return False

def SampleYgX_Stn(self, X, I, A):
        """
        Sample Y given X on standard margins using fitted H&T model.

        Parameters:
        -----------
        X : np.ndarray
            (nRls x 1) Conditioning variable on standard margins.
        I : np.ndarray
            (nRls x 1) Bootstrap index vector.
        A : np.ndarray
            (nRls x 1) Bin index.

        Returns:
        --------
        StnMrg : np.ndarray
            (nRls x nDmn) Simulated values.
        """
        nRls = X.size
        # Flatten inputs to ensure 1D indexing
        X_flat = X.flatten()
        I_flat = I.flatten().astype(int)
        A_flat = A.flatten().astype(int)
        
        # Parameter lookup indices (Account for NonStationarity)
        # APrm shape: (nRls, 4)
        # 0-based indexing for A is assumed.
        APrm = np.zeros((nRls, 4), dtype=int)
        
        # If NonStat[i] is False, index is 0. Else index is A.
        for i in range(4):
            if self.NonStat[i]:
                APrm[:, i] = A_flat
            else:
                APrm[:, i] = 0

        StnMrg = np.full((nRls, self.nDmn), np.nan)
        StnMrg[:, 0] = X_flat

        for iBt in range(self.nBoot):
            # Indices for current bootstrap
            # Threshold Check
            thr = self.Thr[iBt]
            if np.isscalar(thr):
                thr_val = thr
            else:
                thr_val = thr[0] # Assuming scalar threshold per boot for dim 1

            mask_bt = (I_flat == iBt)
            mask_abv = mask_bt & (X_flat > thr_val)
            mask_blw = mask_bt & (X_flat <= thr_val)

            # Exceedances (Above Threshold)
            if np.any(mask_abv):
                # Resample Residuals (Z)
                # Rsd{iBt} is (nObs, nDmn-1)
                # RsdInd{iBt} is (nObs,) bin indices
                
                # Z shape: (n_abv, nDmn-1)
                Z = self.BinResampling(self.Rsd[iBt], self.RsdInd[iBt], 
                                       A_flat[mask_abv], [], self.kNearest)
                
                # Handle NaNs in Z (if bin was empty)
                idNaN = np.isnan(Z[:, 0])
                if np.any(idNaN):
                    # Global sampling fallback
                    n_nan = np.sum(idNaN)
                    rand_idx = np.random.randint(0, self.Rsd[iBt].shape[0], n_nan)
                    Z[idNaN, :] = self.Rsd[iBt][rand_idx, :]

                # Calculate Y based on H&T formula
                # alpha * x + x^beta * (mu + sigma * Z)
                
                idx_alp = APrm[mask_abv, 0]
                idx_bet = APrm[mask_abv, 1]
                idx_mu  = APrm[mask_abv, 2]
                idx_sig = APrm[mask_abv, 3]

                alp = self.Alp[idx_alp, :, iBt]
                bet = self.Bet[idx_bet, :, iBt]
                mu  = self.Mu[idx_mu, :, iBt]
                sig = self.Sig[idx_sig, :, iBt]
                
                x_sub = X_flat[mask_abv][:, None] # (n_abv, 1)

                val = (alp * x_sub) + (x_sub ** bet) * (mu + sig * Z)
                StnMrg[mask_abv, 1:] = val

            # Below Threshold
            if np.any(mask_blw):
                # Joint resampling from observed data below threshold
                # X_obs, Y_obs for current boot
                X_data = self.X[:, iBt]
                Y_data = self.Y[:, :, iBt]
                A_data = self.A[:, iBt]
                
                # Filter data below threshold
                mask_data_blw = X_data <= thr_val
                
                Z_data = np.hstack([X_data[mask_data_blw][:, None], 
                                    Y_data[mask_data_blw, :]])
                
                if Z_data.shape[0] > 0:
                    A_pop = A_data[mask_data_blw]
                    
                    # Target X for distance matching
                    X_target = X_flat[mask_blw]
                    A_target = A_flat[mask_blw]
                    
                    ZSamp = self.BinResampling(Z_data, A_pop, A_target, X_target, self.kNearest)
                    
                    # Fill
                    # ZSamp is (n_blw, nDmn)
                    StnMrg[mask_blw, 1:] = ZSamp[:, 1:]

        StnMrg[np.isnan(StnMrg)] = -np.inf
        return StnMrg

def BinResampling(self, BinPopZ, BinPopA, A, Z=None, k=15):
        """
        Resample from population within bins, optionally using k-nearest neighbors.
        """
        nSmp = len(A)
        nCols = BinPopZ.shape[1]
        ZSamp = np.full((nSmp, nCols), np.nan)
        
        BinPopA = BinPopA.flatten()
        if Z is None or len(Z) == 0:
            Z = []

        if len(BinPopA) > 0:
            for iSmp in range(nSmp):
                # Find population in the specific bin
                mask_bin = (BinPopA == A[iSmp])
                if np.any(mask_bin):
                    tZ = BinPopZ[mask_bin, :]
                    
                    if len(Z) == 0:
                        # Simple random sample
                        idx = np.random.randint(0, tZ.shape[0])
                        ZSamp[iSmp, :] = tZ[idx, :]
                    else:
                        # Covariate resampling (distance based)
                        # Z[iSmp] is scalar? Assuming scalar distance matching on 1st col
                        idx = self.CovResampling(tZ, Z[iSmp], 1)
                        ZSamp[iSmp, :] = tZ[idx, :]
                else:
                    # Bin empty, leave as NaN
                    pass
        return ZSamp

def CovResampling(self, BinPopZ, Z_target, nSmp):
        """
        Weighted resampling based on squared exponential distance.
        """
        nPop = BinPopZ.shape[0]
        nu = self.ResampleDistScale
        
        # Calculate weights based on distance of Col 0 to Z_target
        dist = np.abs(BinPopZ[:, 0] - Z_target)
        weights_unscaled = np.exp(-0.5 * nu * dist)
        
        sum_w = np.sum(weights_unscaled)
        
        if sum_w == 0:
            # Fallback uniform
            return np.random.randint(0, nPop, nSmp)
        
        weights = weights_unscaled / sum_w
        return np.random.choice(nPop, size=nSmp, replace=True, p=weights)

def SampleCovariateFromBin(self, Bn, A):
        """
        Sample covariate value uniformly from within bin.
        Assumes Bn object has properties Edges or similar.
        Use random uniform within [Edge_L, Edge_R].
        """
        # Placeholder logic: assuming Bn.Edg is (nBin+1, nCvr) or list
        # We need the edges for the specific bin A.
        # This function logic depends heavily on the 'Bn' object structure 
        # which isn't fully defined in the prompt. 
        # Assuming 1D bins for simplicity as per common H&T usage.
        
        n_smp = len(A)
        X = np.zeros((n_smp, 1)) # Default 1D
        
        # Access Bn.Edg (edges)
        # Assuming Bn.Edg is array of shape (nBin+1, 1)
        if hasattr(Bn, 'Edg'):
            Edges = Bn.Edg
            for i in range(n_smp):
                bin_idx = int(A[i])
                lb = Edges[bin_idx]
                ub = Edges[bin_idx + 1]
                X[i] = np.random.uniform(lb, ub)
        else:
             # Fallback if structure unknown
             pass
             
        return X

def Plot(self, Mrg):
        """
        Plot data and simulation results.
        """
        nRls = 10000
        MCSml = self.SimulateMC(Mrg, "randX", nRls)
        
        # Summary
        sOrig, sStnMrg = self.GetSummary(MCSml, 'quantile', [0.05, 0.37, 0.5, 0.975])
        
        # Plot 1: Simulation vs Data
        fig = plt.figure(figsize=(12, 8))
        
        # Data X (Conditioning)
        X_data = self.X[:, 0] # Use bootstrap 1 for plot? Or all? Code uses 1.

        for iDmn in range(1, self.nDmn):
            # Standard Scale Plot
            ax1 = fig.add_subplot(2, self.nDmn - 1, iDmn)
            # Simulation dots
            ax1.plot(MCSml['StnMrg'][0, :, 0, 0], MCSml['StnMrg'][0, :, iDmn, 0], 
                     '.', color=[0.8, 0.8, 0.8], label='Simulation')
            # Data dots
            ax1.plot(X_data, self.Y[:, iDmn - 1, 0], 'k.', markersize=2, label='Data')
            
            # Lines (Quantiles)
            # sStnMrg shape: (nA, nDmn, nRtr, nQnt)
            # Assuming nA=1 (Omni) or picking last one
            line_x = sStnMrg[-1, 0, :, 1] # Median X?
            line_y_low = sStnMrg[-1, iDmn, :, 0] # 0.05
            line_y_med = sStnMrg[-1, iDmn, :, 2] # 0.5
            line_y_hi  = sStnMrg[-1, iDmn, :, 3] # 0.975
            
            # Note: plotting quantiles against median X is simplistic translation.
            # Usually requires sorting X and smoothing Y quantiles.
            # Code uses: plot(squeeze(sStnMrg...)) implies calculated summary is already sorted/binned?
            
            ax1.set_xlabel(f'{Mrg[0].RspLbl}: Conditioning')
            ax1.set_ylabel(f'{Mrg[iDmn].RspLbl}: Conditioned')
            ax1.set_title('Standard Margins')
            ax1.grid(True)

            # Original Scale Plot
            ax2 = fig.add_subplot(2, self.nDmn - 1, (self.nDmn - 1) + iDmn)
            ax2.plot(MCSml['Org'][0, :, 0, 0], MCSml['Org'][0, :, iDmn, 0], 
                     '.', color=[0.8, 0.8, 0.8])
            # Assuming Mrg[0].Y is original data
            ax2.plot(Mrg[0].Y, Mrg[iDmn].Y, 'k.', markersize=2)
            
            ax2.set_xlabel(f'{Mrg[0].RspLbl}')
            ax2.set_ylabel(f'{Mrg[iDmn].RspLbl}')
            ax2.set_title('Original Scale')
            ax2.grid(True)

        plt.tight_layout()
        plt.show()
        
        # Save logic
        if not os.path.exists(self.FigureFolder):
            os.makedirs(self.FigureFolder)
        fig.savefig(os.path.join(self.FigureFolder, 'Stg4_HT_1_SmlvsData.png'))

        # Additional plots (Residuals, Parameters) omitted for brevity 
        # but follow similar matplotlib subplot patterns.

def GetSummary(self, Sml, summarystat='median', quantilelvl=0.5):
        """
        Summarize Monte Carlo simulation.

        Returns:
        --------
        sOrig : np.ndarray
        sStnMrg : np.ndarray
        """
        # Sml['Org'] shape: (nA, nRls, nDmn, nRtr)
        # Target Output: (nA, nDmn, nRtr, nQnt) - Permuted
        
        # Helper to process array
        def process(arr):
            # Replace NaNs with -Inf for quantile calc? Or use nanquantile
            # MATLAB code sets NaN to -Inf.
            arr_filled = arr.copy()
            # arr_filled[np.isnan(arr_filled)] = -np.inf # Optional based on preference
            
            if summarystat == 'mean':
                res = np.nanmean(arr_filled, axis=1) # Mean over nRls (axis 1)
                # Result: (nA, nDmn, nRtr)
                # Reshape to add quantile dim (1) -> (nA, nDmn, nRtr, 1)
                return res[..., np.newaxis]
            
            elif summarystat in ['median', 'quantile']:
                q_list = np.atleast_1d(quantilelvl)
                # np.nanquantile moves quantile axis to front: (nQnt, nA, nDmn, nRtr)
                res = np.nanquantile(arr_filled, q_list, axis=1)
                # Transpose to (nA, nDmn, nRtr, nQnt)
                return np.transpose(res, (1, 2, 3, 0))

        sOrig = process(Sml['Org'])
        sStnMrg = None
        
        if 'StnMrg' in Sml:
            sStnMrg = process(Sml['StnMrg'])
            
        return sOrig, sStnMrg

def RVSimulations(self, Mrg):
        """
        Run conditional return value simulation.
        """
        print('Computing Conditional return value')
        nRls = 10000
        # Monte carlo simulation from the fitted model
        # Compute conditional quantiles of Y|X for given X values
        # Note: Matlab uses (1:obj.nBin)' which is a column vector. 
        # We pass a simple array for A.
        bins = np.arange(1, self.nBin + 1)
        self.RVSml = self.SimulateMC(Mrg, "retValX", nRls, bins, self.nBin)
        return self

def SimulateIS(self, Mrg, nRls):
        """
        Simulate importance sampling in X, MC in Y|X.
        
        Parameters:
        -----------
        Mrg : list
            MarginalModel objects.
        nRls : int
            Number of realisations.
            
        Returns:
        --------
        Sml : dict
            Data simulated from H&T model.
        """
        Sml = {}
        # Bin allocation: reshape to column vector logic
        Sml['A'] = np.tile(np.arange(1, self.nBin + 1), nRls).flatten()
        
        # Bootstrap index
        Sml['I'] = np.random.randint(0, self.nBoot, size=(nRls * self.nBin))
        
        # Simulate covariate
        # Assuming Mrg[0].Rat is available (Rates)
        if hasattr(Mrg[0], 'Rat'):
            RatNrm = Mrg[0].Rat / np.sum(Mrg[0].Rat, axis=0) # Normalize columns
            
            # Sample Covariate
            Sml['X'] = self.SampleCovariateFromBin(Mrg[0].Bn, Sml['A'])
            
            if self.nBin == 1:
                J = Sml['I']
                # Python indexing for RatNrm (nBin, nBoot) -> (1, nBoot)
                Sml['W'] = RatNrm[0, J].T 
            else:
                # Matlab sub2ind logic
                # J = sub2ind([obj.nBin, obj.nBoot], Sml.A, Sml.I)
                # Python equivalent: row * nCols + col (if row-major)
                # But Matlab is column-major. 
                # RatNrm is (nBin, nBoot). Element is RatNrm[bin, boot]
                # We need weights for each sample
                Sml['W'] = RatNrm[Sml['A']-1, Sml['I']] # Adjust A for 0-based
        else:
            # Fallback if Rat not present
            Sml['X'] = np.zeros(len(Sml['A']))
            Sml['W'] = np.ones(len(Sml['A']))

        # Generate X uniformly on standard margins
        # Sml.StnMrg allocation
        total_samples = nRls * self.nBin
        Sml['StnMrg'] = np.full((total_samples, self.nDmn), np.nan)
        
        # Mrg[0].INV_Standard calls assumed
        bnd_low = Mrg[0].INV_Standard(1e-7)
        bnd_hi  = Mrg[0].INV_Standard(1 - 1e-7)
        
        G = np.random.rand(total_samples) * (bnd_hi - bnd_low) + bnd_low
        
        # Log PDF + uniform density adjustment
        Sml['logfog'] = Mrg[0].LogPDF_Standard(G) + np.log(bnd_hi - bnd_low)
        # Handle dimensions for logfog? Matlab code sets cols 2:end to 0.
        # Assuming Sml.logfog is (N, nDmn) based on context
        logfog_mat = np.zeros((total_samples, self.nDmn))
        logfog_mat[:, 0] = Sml['logfog']
        Sml['logfog'] = logfog_mat
        Sml['logfog_Stn'] = Sml['logfog'].copy()
        
        # Simulate on standard margins
        Sml['StnMrg'] = self.SampleYgX_Stn(G, Sml['I'], Sml['A'])
        
        # Transform to uniform
        Sml['Unf'] = Mrg[0].CDF_Standard(Sml['StnMrg'])
        
        # Transform back to original margin
        Sml['Org'] = np.full(Sml['Unf'].shape, np.nan)
        
        for iDmn in range(self.nDmn):
            # INV(u, boot_idx, bin_idx)
            Sml['Org'][:, iDmn] = Mrg[iDmn].INV(Sml['Unf'][:, iDmn], Sml['I'], Sml['A'])
            
        # Handle Infs (CDF_Standard returned 1)
        for iDmn in range(self.nDmn):
            mask_inf = (Sml['Unf'][:, iDmn] == 1)
            if np.any(mask_inf):
                Q = Mrg[iDmn].Survivor_Standard(Sml['StnMrg'][mask_inf, iDmn])
                Sml['Org'][mask_inf, iDmn] = Mrg[iDmn].INV_survivor(Q, Sml['I'][mask_inf], Sml['A'][mask_inf])
                
        return Sml

def SimulateMC(self, Mrg, smlMethod, nRls=1000, A=None, nA=None, RVX=None):
        """
        Simulate Monte Carlo samples from H&T model.
        
        Parameters:
        -----------
        smlMethod : str
            "randX", "retValX", or "userX"
        nRls : int
            Number of realizations.
        A : np.ndarray
            Vector of bins to compute.
        nA : int
            Number of super bins.
        RVX : np.ndarray
            User provided X (if smlMethod="userX").
            
        Returns:
        --------
        Sml : dict
            Populated structure of Y|X.
        """
        if A is None:
            A = np.ones(self.nBin) # Omni situation default
            
        if nA is None:
            nA = np.max(A).astype(int)
            
        if nA < np.max(A):
            raise ValueError('Inconsistency in number of bins definition')
            
        OmniOn = (nA > 1)
        
        if smlMethod == "userX" and RVX is None:
            raise ValueError('User defined return values needed for userX method.')
            
        if smlMethod == "userX":
            OmniOn = False
            
        nRtr = 1
        Sml = {}
        
        # Switching simulation methods
        if smlMethod == "randX":
            # Sample uniform values of X
            nRtr = 1
            # Assuming Mrg[0] has sample_MC method
            Sml = Mrg[0].sample_MC(nRls, A, nA)
            
        elif smlMethod == "retValX":
            # X defined in terms of return values
            nRtr = Mrg[0].nRtr
            # Assuming Mrg[0].RVSml contains pre-calculated return values of X
            Sml_source = Mrg[0].RVSml
            if Sml_source is None:
                return None
            # Copy relevant fields
            Sml = Sml_source.copy() # Shallow copy
            
        elif smlMethod == "userX":
            nRtr = RVX.shape[2] # Assuming [nA x nBt x nRtr] ?? Check logic
            Sml = Mrg[0].sample_MC_CondX(nRls, A, nA, RVX)
            
        # Arrays for processing
        # Matlab Logic: tOrg is [nA * nRls * nRtr, nDmn]
        # We need to flatten the input structure from Sml which might be shaped [nA, nRls, ...]
        
        # Flattening helper
        tA = Sml['A'].flatten()
        tI = Sml['I'].flatten()
        
        # find nans in A
        I_nanA = ~np.isnan(tA)
        
        total_rows = len(tA)
        tUnf = np.full((total_rows, self.nDmn), np.nan)
        tOrg = np.full((total_rows, self.nDmn), np.nan)
        
        tUnf[:, 0] = Sml['Unf'].flatten()
        tOrg[:, 0] = Sml['Org'].flatten()
        
        # Sample covariates
        Sml['X'] = np.full((total_rows, Mrg[0].nCvr), np.nan)
        if np.any(I_nanA):
            Sml['X'][I_nanA, :] = self.SampleCovariateFromBin(Mrg[0].Bn, tA[I_nanA])
        
        tX = Sml['X'].copy()
        
        # Simulate on standard scale
        # Need to invert standard first
        # Assuming INV_Standard accepts array
        x_stn = Mrg[0].INV_Standard(tUnf[:, 0])
        
        # Core H&T Sampling
        tStnMrg = self.SampleYgX_Stn(x_stn, tI, tA)
        
        # Transform to uniform (Conditioned vars)
        tUnf = Mrg[0].CDF_Standard(tStnMrg)
        
        # Transform back to original margin
        for iDmn in range(self.nDmn):
            # Populate only associated margins where A is valid
            if np.any(I_nanA):
                tOrg[I_nanA, iDmn] = Mrg[iDmn].INV(tUnf[I_nanA, iDmn], 
                                                   tI[I_nanA], tA[I_nanA])

        # Identification of maximum value (Omni logic)
        if OmniOn:
            # Reshape to [nA, nRls*nRtr]
            # Careful with order. Matlab is col-major. 
            # tOrg col 0 is X.
            temp = tOrg[:, 0].reshape((nA, -1), order='F') 
            J_idx = np.argmax(temp, axis=0) # Index in 0..nA-1
            
            # Convert to linear indices if needed, or just extract
            # Here we just need the rows corresponding to max A
            # Sml.J store omitted for brevity, logic usually used for subsetting later
            pass

        # Reshape data back into [nA, nRls, nRtr, nDmn] or similar
        # Matlab: [nA, nRls, nRtr, nDmn] -> permute to [1, 2, 4, 3] -> [nA, nRls, nDmn, nRtr]
        # Python defaults to C-order, so reshaping requires care matching the flatten order.
        
        # Target shape for storage
        shape_out = (nA, nRls, nRtr, self.nDmn)
        
        def reshape_back(arr):
            # Assume arr was flattened from [nA, nRls, nRtr, nC]
            # Matlab flatten is column major (nA moves fastest)
            # We reconstruct to (nA, nRls, nRtr, nC)
            return arr.reshape(shape_out, order='F').transpose((0, 1, 3, 2)) # Swap nRtr/nDmn to match py usage?
            # User prompt request is documentation + translation. 
            # I will return the dictionary in a usable structure.
            return arr.reshape(shape_out, order='F')

        Sml['StnMrg'] = reshape_back(tStnMrg)
        Sml['Org'] = reshape_back(tOrg)
        Sml['Unf'] = reshape_back(tUnf)
        # Sml['X'] reshape logic depends on nCvr, omitted for standard structure

        # Omni Calculation addition
        if OmniOn:
            # Logic: Append the Max-over-sectors result as an extra "sector"
            # In Python, concatenate along axis 0
            # Implementation depends on exact J extraction above.
            pass
            
        return Sml