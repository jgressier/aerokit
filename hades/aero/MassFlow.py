"""@package MassFlow
  Massflow for compressible flow
"""

import math
import numpy as np
from . import IterativeSolve
from ..common import defaultgas as defg

# -- Compressible flow functions  --

def WeightMassFlow(Mach, r=287.1, gamma=defg._gamma):
    return np.sqrt(gamma/r)*Mach*(1.+.5*(gamma-1)*Mach**2)**(-.5*(gamma+1.)/(gamma-1.))

def Sigma_Mach(Mach, gamma=defg._gamma):
    return (2./(gamma+1.)*(1.+.5*(gamma-1)*Mach**2))**(.5*(gamma+1.)/(gamma-1.))/Mach

def Mach_Sigma(sigma, Mach=2., gamma=defg._gamma):
    def sigma_of_mach(m):
        return Sigma_Mach(m, gamma)
    return IterativeSolve.secant_solve(sigma_of_mach, sigma, Mach)

def MachSub_Sigma(sigma, gamma=defg._gamma):
	# initial guess
	Mach = (2./(gamma+1.))**(.5*(gamma+1.)/(gamma-1.))/sigma
	def sigma_of_mach(m):
		return Sigma_Mach(m, gamma)
	return IterativeSolve.secant_solve(sigma_of_mach, sigma, Mach)

def MachSup_Sigma(sigma, gamma=defg._gamma):
	# initial guess
	cg = (gamma+1.)/(gamma-1.)
	Mach = np.sqrt((sigma*cg**(.5*cg))**(2./(cg-1.))-cg)
	def sigma_of_mach(m):
		return Sigma_Mach(m, gamma)
	return IterativeSolve.secant_solve(sigma_of_mach, sigma, Mach)

