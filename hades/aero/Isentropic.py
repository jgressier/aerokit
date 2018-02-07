"""@package Isentropic
  isentropic functions
"""

import math
import numpy as np
from . import IterativeSolve

# -- isentropic functions --

def TiTs_Mach(Mach, gamma=1.4):
    return 1.+.5*(gamma-1)*Mach**2

def PiPs_Mach(Mach, gamma=1.4):
    return (1.+.5*(gamma-1.)*Mach**2)**(gamma/(gamma-1.))

def Mach_TiTs(TiTs, gamma=1.4):
    return np.sqrt((TiTs-1.)*2./(gamma-1.))

def Mach_PiPs(PiPs, gamma=1.4):
    return np.sqrt((PiPs**((gamma-1.)/gamma)-1.)*2./(gamma-1.))

def Velocity_MachTi(Mach, Ti, r=287.1, gamma=1.4):
    return Mach*np.sqrt(gamma*r*Ti/TiTs_Mach(Mach, gamma))
