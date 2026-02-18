import numpy as np

class OptionsContours:
    """
    Class of options for estimation of contours using PPC (Predictive Probability Contours).
    
    Four contour estimation methods are supported:
      - 'Exc': Exceedance contour
      - 'HTDns': Heffernan and Tawn Density contour
      - 'Hus': Huseby contour
      - 'HusOld': Huseby Cleaned contour (referred to as HusCln in comments, HusOld in code)
    
    Attributes:
        nGrd (int): Number of grid points for each of x and y when rectangular gridding needed. Default 50.
        nPnt (int): Number of points on the contour. Default 200.
        SmtWdtC (int): Smoothing width for the Huseby contour C function. Default 5.
        BndWdtScl (float): Band width scale for HTDns contour. Default 0.02.
        Mth (list): List of contour Methods to use. Default empty list.
        nSml (int): Number to simulate (number of importance samples). Default 1e6.
    """

    def __init__(self):
        # Initialize with default values using setters to ensure validation
        self.nGrd = 50
        self.nPnt = 200
        self.SmtWdtC = 5
        self.BndWdtScl = 0.02
        self.Mth = []
        self.nSml = 1000000

    @property
    def nSml(self):
        return self._nSml

    @nSml.setter
    def nSml(self, value):
        if not (isinstance(value, (int, float, np.number)) and value > 0 and value % 1 == 0):
             raise ValueError("nSml must be a positive integer scalar.")
        self._nSml = int(value)

    @property
    def nGrd(self):
        return self._nGrd

    @nGrd.setter
    def nGrd(self, value):
        if not (isinstance(value, (int, float, np.number)) and value > 0 and value % 1 == 0):
             raise ValueError("nGrd must be a positive integer scalar.")
        self._nGrd = int(value)

    @property
    def nPnt(self):
        return self._nPnt

    @nPnt.setter
    def nPnt(self, value):
        if not (isinstance(value, (int, float, np.number)) and value > 0 and value % 1 == 0):
             raise ValueError("nPnt must be a positive integer scalar.")
        self._nPnt = int(value)

    @property
    def SmtWdtC(self):
        return self._SmtWdtC

    @SmtWdtC.setter
    def SmtWdtC(self, value):
        if not (isinstance(value, (int, float, np.number)) and value > 0 and value % 1 == 0):
             raise ValueError("SmtWdtC must be a positive integer scalar.")
        self._SmtWdtC = int(value)

    @property
    def BndWdtScl(self):
        return self._BndWdtScl

    @BndWdtScl.setter
    def BndWdtScl(self, value):
        if not (isinstance(value, (int, float, np.number)) and value > 0):
             raise ValueError("BndWdtScl must be a positive scalar.")
        self._BndWdtScl = float(value)

    @property
    def Mth(self):
        return self._Mth

    @Mth.setter
    def Mth(self, value):
        # Handle single string input by wrapping in list
        if isinstance(value, str):
            value = [value]
        
        # Define valid options
        valid_options = {'Exc', 'HTDns', 'Hus', 'HusOld'}
        
        validated_list = []
        if value is not None:
            for item in value:
                # Case-insensitive match check logic could go here, 
                # but following strict MATLAB validatestring behavior for exact matches provided in code
                if item not in valid_options:
                     raise ValueError(f"Method '{item}' is not valid. Allowed: {valid_options}")
                validated_list.append(item)
        
        # Remove duplicates while preserving order (or just use set if order doesn't matter)
        # Using sorted(set()) to mimic unique() behavior
        self._Mth = sorted(list(set(validated_list)), key=validated_list.index)