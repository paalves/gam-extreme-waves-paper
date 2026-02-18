import numpy as np

def cdf_interp1(Z, Zi, P=None):
    """
    Performs 'next' neighbor interpolation for a Cumulative Distribution Function (CDF).

    This function computes probabilities for interpolation points based on a sorted 
    dataset Z. It assigns the cumulative probability of the nearest preceding value 
    in Z to each point in Zi. Effectively, this implements a step function where 
    the probability remains constant until the next value in Z is reached.

    Parameters:
    -----------
    Z : array_like
        A sorted input vector of shape (n,) or (n, 1) representing the domain of the CDF 
        (e.g., sample quantiles).
    Zi : array_like
        The points at which to evaluate the CDF. Can be of any shape (p, q).
    P : array_like, optional
        A vector of probabilities corresponding to Z, of shape (n,) or (n, 1). 
        If not provided, P defaults to equally spaced probabilities [1/n, 2/n, ..., 1].

    Returns:
    --------
    Pi : numpy.ndarray
        An array of probabilities with the same shape as Zi, containing values in the range [0, 1].
    """
    Z = np.asarray(Z).flatten()
    Zi = np.asarray(Zi)
    n = Z.size

    if n == 0:
        return np.ones_like(Zi, dtype=float)
    
    if n == 1:
        return (Zi >= Z[0]).astype(float)

    if P is None:
        P = np.arange(1, n + 1) / n
    else:
        P = np.asarray(P).flatten()

    indices = np.searchsorted(Z, Zi, side='right') - 1
    
    Pi = np.zeros_like(Zi, dtype=float)
    
    valid_mask = indices >= 0
    Pi[valid_mask] = P[indices[valid_mask]]

    return Pi