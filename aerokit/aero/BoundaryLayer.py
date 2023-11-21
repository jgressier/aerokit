"""
    The ``BoundaryLayer`` module
    =========================
 
    Provides integral models for boundary layer computations
  
    Available functions
    -------------------
 
"""

class base():
    def __init__(self, Re, L=1.):
        self._Re = Re
        self._L  = L

    def delta1(self):
        raise NameError("not implemented function")

    def theta(self):
        raise NameError("not implemented function")

    def Cf(self):
        raise NameError("not implemented function")

class Blasius(base):
    def __init__(self, Re, m, L=1.):
        self._Re = Re
        self._L  = L
        
class FalknerSkan(base):
    def __init__(self, Re, m, L=1.):
        self._Re = Re
        self._L  = L
