import numpy as np
import pandas as pd

def SsnDgr(DatJln):
    """
    For a given vector of Julian dates (MATLAB datenum format), finds the seasonal degree.
    All years have seasonal degrees in [0, 360) regardless of leap years.

    Args:
        DatJln (np.ndarray or list): [n x 1] Julian dates (MATLAB serial date numbers) 
                                     for storm peak data.

    Returns:
        np.ndarray: [n x 1] Seasonal degree (in [0, 360)) corresponding to Julian dates.
    """
    # Ensure input is a flat numpy array
    DatJln = np.array(DatJln).flatten()

    # Convert MATLAB datenum to Pandas DatetimeIndex
    # MATLAB datenum 719529 is equal to 1970-01-01 (Unix Epoch)
    # We subtract this offset to get days since epoch, then convert to datetime
    timestamps = pd.to_datetime(DatJln - 719529, unit='D')

    # Extract Day of Year (1 to 365 or 366)
    # Note: Pandas .dayofyear is 1-indexed (Jan 1st = 1)
    doy = timestamps.dayofyear.to_numpy()

    # Determine total number of days in the specific year (365 or 366)
    is_leap = timestamps.is_leap_year.to_numpy()
    days_in_year = np.where(is_leap, 366, 365)

    # Calculate Seasonal Degree
    # Formula: 360 * (DayOfYear / DaysInYear)
    SsnDgr = 360.0 * doy / days_in_year

    return SsnDgr