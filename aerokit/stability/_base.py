import numpy as np
import aerokit.common.numspectral as ns

class LinOperator():

    def __init__(self, n=100) -> None:
        self._diffop = ns.ChebCollocation(n)

    @property
    def dim(self):
        return self._diffop.npts
    
    @property
    def x(self):
        return self._diffop.x
    
    def set_basestate(self, state):
        """set base state as it will be used by derived class"""
        self._basestate = state



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()