
import numpy as np
import matplotlib.pyplot as plt

class CovariateBinning:
    """
    Allocates data to covariate bins and visualizes the binning strategy.

    This class handles the discretization of covariate data into bins, supporting
    both periodic (e.g., angular) and non-periodic variables. It manages bin edges,
    assigns observations to bins, and provides visualization tools for marginal
    and joint distributions.

    Attributes:
        n_cvr (int): Number of covariates.
        n (int): Number of observations.
        n_bin (int): Total number of unique bin combinations across all covariates.
        n_bin_cvr (np.ndarray): Array containing the number of bins for each covariate.
        edg (list): List of arrays defining the bin edges for each covariate.
        a (np.ndarray): Bin allocation index for each observation (linear index).
        cnt (np.ndarray): Count of observations in each bin.
        bin_lbl (list): String labels describing each bin.
        left_edge (np.ndarray): The starting value (left edge) of the bin for each dimension.
        length (np.ndarray): The width/length of the bin for each dimension.
        edg_exd (list): Extended bin edges (padded for periodic variables).
        is_prd (np.ndarray): Boolean array indicating if a covariate is periodic.
        a_prj (np.ndarray): Mapping matrix to project linear bin indices back to individual covariate bin indices.
        rsp_lbl (list): Labels for response variables.
        cvr_lbl (list): Labels for covariate variables.
    """

    def __init__(self, x, edg, is_prd, y=None, rsp_lbl=None, cvr_lbl=None):
        """
        Initializes the CovariateBinning object, performs validation, and executes bin allocation.

        Parameters:
            x (np.ndarray): Input covariate data of shape [n, n_cvr]. Periodic values must be [0, 360).
            edg (list): List of arrays defining bin edges. Periodic bins wrap around 0/360.
            is_prd (list or np.ndarray): Boolean vector indicating periodic covariates.
            y (np.ndarray, optional): Response variables [n, n_dmn] used for plotting.
            rsp_lbl (list, optional): Labels for response variables.
            cvr_lbl (list, optional): Labels for covariate variables.
        """
        self.n, self.n_cvr = x.shape
        self.edg = [np.array(e) for e in edg]
        self.is_prd = np.array(is_prd, dtype=bool)
        self.y = y
        self.rsp_lbl = rsp_lbl if rsp_lbl is not None else [f"Y{i+1}" for i in range(y.shape[1])] if y is not None else []
        self.cvr_lbl = cvr_lbl if cvr_lbl is not None else [f"X{i+1}" for i in range(self.n_cvr)]

        self.n_bin = 0
        self.n_bin_cvr = np.zeros(self.n_cvr, dtype=int)
        self.a = None
        self.cnt = None
        self.bin_lbl = []
        
        self.left_edge = None
        self.length = None
        self.edg_exd = [None] * self.n_cvr
        self.a_prj = None

        self.bin_allocation(x)

        if y is not None:
            self.plot_bins(x, y)

    def bin_allocation(self, x):
        """
        Performs the discretization of covariates into bins.

        This method iterates over each covariate, determines the appropriate bin edges 
        (handling wrapping for periodic variables), and assigns a linear bin index to 
        each observation. It also populates bin metadata like labels, lengths, and counts.

        Parameters:
            x (np.ndarray): The covariate data matrix.
        """
        a_dmn = np.zeros((self.n, self.n_cvr), dtype=int)

        for i_c in range(self.n_cvr):
            self.edg[i_c] = np.sort(self.edg[i_c])

            if self.is_prd[i_c]:
                if np.any(x[:, i_c] > 360) or np.any(x[:, i_c] < 0):
                    raise ValueError(f"Periodic X[:,{i_c}] expected on range [0,360]")
                
                if np.unique(np.mod(self.edg[i_c], 360)).size < self.edg[i_c].size:
                     raise ValueError(f"Covariate {i_c}: Supplied bin edges correspond to the same point on a circle.")

                self.n_bin_cvr[i_c] = self.edg[i_c].size
                
                # Padding edges for periodic handling
                self.edg_exd[i_c] = np.unique(np.concatenate(([0], self.edg[i_c], [360])))
            else:
                self.edg_exd[i_c] = self.edg[i_c]
                self.n_bin_cvr[i_c] = self.edg[i_c].size - 1
                
                if np.any(x[:, i_c] < self.edg[i_c][0]) or np.any(x[:, i_c] > self.edg[i_c][-1]):
                    raise ValueError(f"Non Periodic X[:,{i_c}] outside of range defined by bins")

            # Discretize
            # np.digitize returns 1-based indices for bins defined by bins[i-1] <= x < bins[i]
            # We treat the extended edges as the boundaries.
            inds = np.digitize(x[:, i_c], self.edg_exd[i_c])
            
            # Handle right-edge inclusivity matching MATLAB behavior usually implies 
            # checking bounds carefully. For now, assuming standard left-closed intervals.
            # Fix indices to be 1-based relative to the defined bins.
            
            if self.is_prd[i_c]:
                # If periodic, we need to map the padding bins back to original bins.
                # If edg_exd is [0, 90, 180, 270, 360] (original [90, 180, 270])
                # Bins are: 0-90, 90-180, 180-270, 270-360.
                # If indices go beyond the count, wrap them.
                # Specific logic: wrap the last bin (360) back to 1 if it exceeds count
                pass 
                
            a_dmn[:, i_c] = inds

            if self.is_prd[i_c]:
                # Wrap periodic bins: if it falls in the last padded bin, it might belong to the first logical bin
                # depending on edge definitions.
                # Matching MATLAB logic: ADmn(ADmn==nBinCvr(iC)+1,iC)=1
                mask_wrap = a_dmn[:, i_c] == (self.n_bin_cvr[i_c] + 1)
                a_dmn[mask_wrap, i_c] = 1

        # Combine bins over dimensions to create a unique linear index
        k = np.cumprod(self.n_bin_cvr)
        self.a = a_dmn[:, 0] - 1 # Convert to 0-based for calculation, then will shift back or keep 0-based
        
        for i_c in range(1, self.n_cvr):
            self.a = self.a + (a_dmn[:, i_c] - 1) * k[i_c - 1]
        
        self.a = self.a.astype(int) # This is now a 0-based index
        
        self.n_bin = np.prod(self.n_bin_cvr)
        
        # Unravel indices to map back from linear index to individual covariate indices
        # MATLAB ind2sub equivalent
        linear_indices = np.arange(self.n_bin)
        
        # Note: MATLAB ind2sub and Python unravel_index use different dimension ordering (F vs C)
        # We need to be careful to match the cumprod logic used above.
        # The accumulation above suggests: A = (d1) + (d2)*n1 + (d3)*n1*n2... (Fortran/MATLAB style)
        
        self.a_prj = np.zeros((self.n_bin, self.n_cvr), dtype=int)
        
        temp_ind = linear_indices
        for i_c in range(self.n_cvr):
            self.a_prj[:, i_c] = temp_ind % self.n_bin_cvr[i_c]
            temp_ind = temp_ind // self.n_bin_cvr[i_c]
            
        # Add 1 to make it compatible with 1-based bin counts for logic below if needed, 
        # but we will stick to 0-based indexing for arrays.

        self.bin_lbl = [""] * self.n_bin
        self.left_edge = np.full((self.n_bin, self.n_cvr), np.nan)
        self.length = np.full((self.n_bin, self.n_cvr), np.nan)

        for i_c in range(self.n_cvr):
            if self.is_prd[i_c]:
                bin_st = self.edg[i_c]
                bin_end = np.roll(self.edg[i_c], -1)
                
                if self.edg[i_c][0] > 0:
                    bin_st = np.roll(bin_st, 1)
                    bin_end = np.roll(bin_end, 1)
            else:
                bin_st = self.edg[i_c][:-1]
                bin_end = self.edg[i_c][1:]
            
            for i_b in range(self.n_bin):
                idx = self.a_prj[i_b, i_c]
                
                curr_st = bin_st[idx]
                curr_end = bin_end[idx]
                
                if curr_end < curr_st and self.is_prd[i_c]:
                    self.length[i_b, i_c] = curr_end + 360 - curr_st
                else:
                    self.length[i_b, i_c] = curr_end - curr_st
                
                self.left_edge[i_b, i_c] = curr_st
                
                lbl_part = f"{self.cvr_lbl[i_c][0]}[{curr_st}, {curr_end})"
                
                if i_c == 0:
                    self.bin_lbl[i_b] = lbl_part
                else:
                    self.bin_lbl[i_b] += " x " + lbl_part

        self.cnt = np.bincount(self.a, minlength=self.n_bin)
        
        if np.any(self.cnt < 30):
            print("Warning: Too few exceedances in one or more bins. Consider reducing number of bins.")
            print(self.cnt)

    def plot_bins(self, x, y):
        """
        Generates visualizations of the binning strategy.

        Produces two figures:
        1. Marginal Plots: Scatter plots of Y vs X for each dimension, with vertical lines indicating bin edges.
        2. Joint Plots: Scatter plots of response variables against each other, colored/separated by bin.

        Parameters:
            x (np.ndarray): Covariate values.
            y (np.ndarray): Response values.
        """
        n_dmn = y.shape[1]
        
        # Marginal Plot
        fig1 = plt.figure(figsize=(10, 2 * n_dmn))
        c = 1
        for i in range(n_dmn):
            for i_c in range(self.n_cvr):
                ax = fig1.add_subplot(n_dmn, self.n_cvr, c)
                ax.plot(x[:, i_c], y[:, i], 'k.', markersize=2)
                ax.grid(True)
                
                # Plot edges
                ylim = ax.get_ylim()
                for edge in self.edg[i_c]:
                    ax.plot([edge, edge], ylim, 'r--', linewidth=1)
                
                ax.set_xlabel(self.cvr_lbl[i_c])
                ax.set_ylabel(self.rsp_lbl[i])
                
                if self.is_prd[i_c]:
                    ax.set_xlim([0, 360])
                    ax.set_xticks(np.arange(0, 361, 45))
                
                c += 1
        plt.tight_layout()
        plt.show()

        # Joint Plot
        if n_dmn > 1:
            for i_d in range(1, n_dmn):
                fig = plt.figure(figsize=(12, 8))
                
                rows = 2
                cols = int(np.ceil(self.n_bin / 2)) if self.n_bin > 1 else 1
                
                for i_c in range(self.n_bin):
                    if self.n_bin > 1:
                        ax = fig.add_subplot(rows, cols, i_c + 1)
                    else:
                        ax = fig.add_subplot(1, 1, 1)
                    
                    mask = (self.a == i_c)
                    if np.any(mask):
                        ax.plot(y[mask, 0], y[mask, i_d], 'k.', markersize=2)
                    
                    ax.grid(True)
                    ax.set_xlabel(self.rsp_lbl[0])
                    ax.set_ylabel(self.rsp_lbl[i_d])
                    
                    if self.n_bin > 1:
                        ax.set_title(self.bin_lbl[i_c], fontsize=8)
                
                plt.tight_layout()
                plt.show()

    def plot_bin_edge(self, i_c):
        """
        Adds vertical lines to the current active plot representing bin boundaries for a specific covariate.

        Parameters:
            i_c (int): Index of the covariate (0-based) to plot edges for.
        """
        ax = plt.gca()
        ylim = ax.get_ylim()
        
        # Note: self.edg is stored as 0-based lists in Python init
        # nBinCvr corresponds to edges or edges-1 depending on periodicity
        # Here we just plot all defined edges for that covariate
        for edge in self.edg[i_c]:
            ax.plot([edge, edge], ylim, 'r--', linewidth=1)
        
        ax.set_ylim(ylim)

    def sample_covariate_from_bin(self, a):
        """
        Samples covariate values uniformly from within the specified bins.

        Parameters:
            a (np.ndarray): Array of linear bin indices (0-based) for which to generate samples.

        Returns:
            np.ndarray: Matrix of sampled covariate values [len(a), n_cvr].
        """
        n_samples = len(a)
        x_new = np.random.rand(n_samples, self.n_cvr) * self.length[a, :] + self.left_edge[a, :]
        
        for i_c in range(self.n_cvr):
            if self.is_prd[i_c]:
                x_new[:, i_c] = np.mod(x_new[:, i_c], 360)
                
        return x_new

    def plot_parameter(self, p, i_c, **kwargs):
        """
        Projects parameters onto a specific covariate dimension and plots the trend.

        Calculates the mean and quantiles (2.5%, 97.5%) of the parameters `p` 
        grouped by the bins of the specified covariate `i_c`.

        Parameters:
            p (np.ndarray): Parameter values to be projected [n, n_params] or [n, 1].
            i_c (int): Index of the covariate dimension to project onto.
            **kwargs: Additional keyword arguments passed to the matplotlib plot function.
        """
        # Determine unique bins in the specific covariate dimension
        # a_prj maps linear bin index -> covariate specific bin index
        
        # We need to map the parameter P (which corresponds to bins) to the covariate dimension
        # If P is provided per bin (size n_bin), we group by a_prj[:, i_c]
        
        unique_cvr_bins = np.unique(self.a_prj[:, i_c])
        n_groups = len(unique_cvr_bins)
        
        # Calculate stats
        p_mean = np.zeros((n_groups, p.shape[1] if p.ndim > 1 else 1))
        p_quant = np.zeros((n_groups, 3)) if p.ndim == 1 else np.zeros((n_groups, 3, p.shape[1]))
        
        # Assuming P is aligned with the TOTAL bins (self.n_bin), not observations
        # If P is aligned with observations, logic would differ slightly. 
        # Based on context "Project and plot parameter", usually P matches nBin rows.
        
        for idx, bin_idx in enumerate(unique_cvr_bins):
            # Find all total-bins that correspond to this specific covariate bin
            mask = self.a_prj[:, i_c] == bin_idx
            p_subset = p[mask]
            
            p_mean[idx] = np.mean(p_subset, axis=0)
            if p.ndim > 1:
                p_quant[idx] = np.percentile(p_subset, [2.5, 50, 97.5], axis=0)
            else:
                p_quant[idx] = np.percentile(p_subset, [2.5, 50, 97.5])

        # Create X axis for plotting (approximate 0 to 360 or min/max edge)
        x_plot = np.linspace(np.min(self.edg_exd[i_c]), np.max(self.edg_exd[i_c]), 360)
        
        if self.is_prd[i_c]:
            t_x = np.mod(x_plot, 360)
        else:
            t_x = x_plot
            
        # Discretize plot points to find which stats to plot
        # Using self.edg_exd[i_c] which is 1D array
        t_a = np.digitize(t_x, self.edg_exd[i_c]) - 1 # 0-based
        
        # Handle periodic wrapping for indices
        if self.is_prd[i_c]:
            # n_bin_cvr is logical bins, edg_exd has padding
            # Map the indices to the logical bin range
             t_a[t_a >= self.n_bin_cvr[i_c]] = 0 # Simple wrap logic for visualization

        # Clip to ensure valid indices
        t_a = np.clip(t_a, 0, n_groups - 1)

        ax = plt.gca()
        
        if p.ndim > 1:
            # If P is multidimensional, plot index 1 (median) as solid, others dashed
            # The structure of qP in MATLAB was [0.025, 0.5, 0.975]
            # Here we follow that logic per parameter column if needed, 
            # but MATLAB code implies PPrj is (N, 3) after quantile.
            
            # Simplified plotting for 1D parameter distribution behavior
            ax.plot(x_plot, p_quant[t_a, 1], '-', **kwargs)
            ax.plot(x_plot, p_quant[t_a, 0], '--', **kwargs)
            ax.plot(x_plot, p_quant[t_a, 2], '--', **kwargs)
        else:
            ax.plot(x_plot, p_mean[t_a, 0], '-', **kwargs)
            
        if self.is_prd[i_c]:
            ax.set_xlim([0, 360])
            ax.set_xticks(np.arange(0, 361, 90))