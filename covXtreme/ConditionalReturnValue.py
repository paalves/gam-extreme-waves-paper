import numpy as np

def conditional_return_value(ht, mrg, n_rls=1000, boot_indices=None):
    """
    Computes the conditional return value Y | X using Monte Carlo simulation based on 
    a precomputed Heffernan and Tawn (HT) model and marginal models.

    This function simulates values of dependent variables (Y) given a conditioning 
    variable (X) using the semi-parametric Heffernan and Tawn approach. It handles 
    residual sampling, standard scale transformation, and re-transformation to 
    original scales.

    Parameters
    ----------
    ht : object
        A structured object containing the precomputed Heffernan and Tawn model parameters.
        Expected attributes:
        - n_bin (int): Number of bins (covariates).
        - n_dmn (int): Number of dimensions.
        - rvmth (int): Method for X generation (1=Random T-yr Max, 2=Fixed, 3=Fixed+Boot).
        - n_boot (int): Number of bootstrap samples available.
        - rv_fix (array): Fixed return values for methods 2 & 3.
        - smp_lcl_rsd_on (bool): Flag for local residual sampling.
        - rsd (list/array): Residuals for simulation.
        - rsd_ind (list/array): Indices for residuals.
        - nep (array): Non-exceedance probabilities.
        - alp, bet, mu, sig (arrays): HT model parameters (Alpha, Beta, Mu, Sigma).
        - rv (object): Placeholder for output simulation structure.
        
    mrg : list of objects
        A list of precomputed marginal model objects. 
        - mrg[0]: Represents the conditioning variable X.
        - mrg[1:]: Represent the dependent variables Y.
        Expected attributes/methods:
        - rtr_prd (array): Return periods.
        - rat (array): Annual rate of occurrence.
        - inv(p, idx): Inverse CDF method.
        - cdf(x, cov, idx): CDF method.
        - cdf_standard(x): Method to convert to standard scale.
        - inv_standard(p): Method to convert from standard scale.

    n_rls : int, optional
        Number of realizations in the Monte Carlo sample. Default is 1000.

    boot_indices : array_like, optional
        Indices for the bootstrap samples to use.

    Returns
    -------
    ht : object
        The input `ht` object with the `rv` attribute populated containing:
        - x, y: Simulated values on original scale.
        - x_stn, y_stn: Simulated values on standard scale.
        - n_rls: Number of realizations used.
    """
    
    n_rtr = len(mrg[0].rtr_prd)
    n_asc = ht.n_dmn - 1
    
    trv_x_stn = np.full((ht.n_bin, n_rls, n_rtr), np.nan)
    trv_y_stn = np.full((ht.n_bin, n_asc, n_rls, n_rtr), np.nan)
    trv_y = np.full((ht.n_bin, n_asc, n_rls, n_rtr), np.nan)
    trv_x = np.full((ht.n_bin, n_rls, n_rtr), np.nan)

    if ht.rvmth == 1:
        if boot_indices is None:
            boot_indices = np.random.randint(0, ht.n_boot, n_rls)
    elif ht.rvmth == 2:
        x_fixed = np.transpose(ht.rv_fix, (0, 2, 1))
        trv_x = np.tile(x_fixed, (1, n_rls, 1))
        if boot_indices is None:
            boot_indices = np.random.randint(0, ht.n_boot, n_rls)
    elif ht.rvmth == 3:
        x_fixed = np.transpose(ht.rv_fix, (0, 2, 1))
        trv_x = np.tile(x_fixed, (1, n_rls, 1))

    residuals = np.zeros((n_rls, n_asc, ht.n_bin))

    if not ht.smp_lcl_rsd_on:
        for i_bin in range(ht.n_bin):
            valid_res = ht.rsd[boot_indices] 
            picked = [x[np.random.randint(len(x))] for x in valid_res] 
            residuals[:, :, i_bin] = np.array(picked)
    else:
        for i_bin in range(ht.n_bin):
            for i_rls in range(n_rls):
                boot_idx = boot_indices[i_rls]
                curr_rsd_ind = ht.rsd_ind[boot_idx]
                curr_rsd = ht.rsd[boot_idx]
                
                mask = (curr_rsd_ind == i_bin)
                if np.any(mask):
                    bin_residuals = curr_rsd[mask]
                    rand_idx = np.random.randint(bin_residuals.shape[0])
                    residuals[i_rls, :, i_bin] = bin_residuals[rand_idx]
                else:
                    residuals[i_rls, :, i_bin] = np.nan

    residuals = np.transpose(residuals, (2, 1, 0))

    for i_rtr in range(n_rtr):
        rho = mrg[0].rat[:, boot_indices]
        lt = rho * mrg[0].rtr_prd[i_rtr]
        
        if ht.rvmth == 1:
            ux = np.random.rand(ht.n_bin, n_rls)
            p = 1 + np.log(ux) / lt
            
            nep_threshold = ht.nep[boot_indices]
            p[p < nep_threshold.T] = np.nan
            
            trv_x[:, :, i_rtr] = mrg[0].inv(p, boot_indices)
        elif ht.rvmth in [2, 3]:
            covariates = np.arange(ht.n_bin).reshape(-1, 1)
            p = mrg[0].cdf(trv_x[:, 0, i_rtr][:, None], covariates, boot_indices)

        trv_x_stn[:, :, i_rtr] = mrg[0].inv_standard(p)
        
        tx = trv_x_stn[:, :, i_rtr][:, None, :]
        
        alpha = ht.alp[:, :, boot_indices]
        beta = ht.bet[:, :, boot_indices]
        mu = ht.mu[:, :, boot_indices]
        sigma = ht.sig[:, :, boot_indices]
        z_res = residuals
        
        term1 = alpha * tx
        term2 = (tx ** beta) * (mu + sigma * z_res)
        trv_y_stn[:, :, :, i_rtr] = term1 + term2

        for i_asc in range(n_asc):
            uy = mrg[i_asc + 1].cdf_standard(trv_y_stn[:, i_asc, :, i_rtr])
            trv_y[:, i_asc, :, i_rtr] = mrg[i_asc + 1].inv(uy, boot_indices)

    if ht.n_bin > 1:
        x_omni_stn = np.max(trv_x_stn, axis=0)
        j_indices = np.argmax(trv_x, axis=0)
        x_omni = np.max(trv_x, axis=0)
        
        y_omni = np.full((1, n_asc, n_rls, n_rtr), np.nan)
        y_omni_stn = np.full((1, n_asc, n_rls, n_rtr), np.nan)
        
        for i_rtr in range(n_rtr):
            flat_indices = np.ravel_multi_index(
                (j_indices[:, i_rtr], np.arange(n_rls)), 
                (ht.n_bin, n_rls)
            )
            
            ty_reshaped = np.transpose(trv_y[:, :, :, i_rtr], (0, 2, 1)).reshape(-1, n_asc)
            y_omni[0, :, :, i_rtr] = ty_reshaped[flat_indices, :].T
            
            ty_stn_reshaped = np.transpose(trv_y_stn[:, :, :, i_rtr], (0, 2, 1)).reshape(-1, n_asc)
            y_omni_stn[0, :, :, i_rtr] = ty_stn_reshaped[flat_indices, :].T

        trv_y = np.concatenate((trv_y, y_omni), axis=0)
        trv_y_stn = np.concatenate((trv_y_stn, y_omni_stn), axis=0)
        trv_x = np.concatenate((trv_x, x_omni[None, :, :]), axis=0)
        trv_x_stn = np.concatenate((trv_x_stn, x_omni_stn[None, :, :]), axis=0)

    class SimulationResult:
        pass
    
    sim_res = SimulationResult()
    sim_res.x = trv_x
    sim_res.y = trv_y
    sim_res.x_stn = trv_x_stn
    sim_res.y_stn = trv_y_stn
    sim_res.i = boot_indices
    sim_res.n_rls = n_rls
    
    ht.rv = sim_res
    return ht