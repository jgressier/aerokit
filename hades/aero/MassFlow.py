"""@package MassFlow
  Compressible flow functions
"""

import math
import numpy as np
from . import IterativeSolve

# -- Compressible flow functions  --

def WeightMassFlow(Mach, r=287.1, gamma=1.4):
    return np.sqrt(gamma/r)*Mach*(1.+.5*(gamma-1)*Mach**2)**(-.5*(gamma+1.)/(gamma-1.))

def Sigma_Mach(Mach, gamma=1.4):
    return (2./(gamma+1.)*(1.+.5*(gamma-1)*Mach**2))**(.5*(gamma+1.)/(gamma-1.))/Mach

def Mach_Sigma(sigma, Mach=2., gamma=1.4):
    def sigma_of_mach(m):
        return Sigma_Mach(m, gamma)
    return IterativeSolve.secant_solve(sigma_of_mach, sigma, Mach)

