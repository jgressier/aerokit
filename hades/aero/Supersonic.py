"""@package Supersonic
  supersonic flow functions
"""

import math
import numpy as np
from . import IterativeSolve
from . import Isentropic as Is
from scipy.optimize import fsolve

# -- 2D supersonic invariants --

def PrandtlMeyer_Mach(Mach, gamma=1.4):
    g    = gamma
    gm1  = gamma - 1.
    beta = np.sqrt(Mach**2-1.)
    return (np.sqrt((g+1.)/gm1)*np.arctan(np.sqrt(gm1/(g+1.))*beta) - np.arctan(beta))*180./math.pi

def old_Mach_PrandtlMeyer(omega, gamma=1.4):
    def omega_of_mach(m):
        return PrandtlMeyer_Mach(m, gamma)
    return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)

def Mach_PrandtlMeyer(omega, gamma=1.4):
    def omega_of_mach(m):
        return PrandtlMeyer_Mach(m, gamma)-omega
    #return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)
    return fsolve(omega_of_mach, 2.+0.*omega)

def deflection_Mach_IsentropicPsratio(Mach, Pratio, gamma=1.4):
    m2 = Is.Mach_PiPs(Is.PiPs_Mach(Mach, gamma)/Pratio, gamma)
    return -PrandtlMeyer_Mach(Mach, gamma)+PrandtlMeyer_Mach(m2, gamma)

def IsentropicPsratio_Mach_deflection(Mach, dev, gamma=1.4):
    m2 = Mach_PrandtlMeyer(PrandtlMeyer_Mach(Mach, gamma)-dev, gamma)
    return Is.PiPs_Mach(Mach, gamma)/Is.PiPs_Mach(m2, gamma)

