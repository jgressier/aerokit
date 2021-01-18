"""@package MassFlow
  Massflow for compressible flow
"""

import math
import numpy as np
#from . import IterativeSolve
from scipy.optimize import fsolve, newton
from ..common import defaultgas as defg

# internal initialization Mach for sigma computations
def __MachSub_sigma(sigma, gamma):
	return (2./(gamma+1.))**(.5*(gamma+1.)/(gamma-1.))/sigma

def __MachSup_sigma(sigma, gamma):
	return np.sqrt(1.+1./sigma)

# -- Compressible flow functions  --

def WeightMassFlow(Mach, r=287.1, gamma=defg._gamma):
    return np.sqrt(gamma)*Mach*(1.+.5*(gamma-1)*Mach**2)**(-.5*(gamma+1.)/(gamma-1.))

def Sigma_Mach(Mach, gamma=defg._gamma):
    return (2./(gamma+1.)*(1.+.5*(gamma-1)*Mach**2))**(.5*(gamma+1.)/(gamma-1.))/Mach

def Mach_Sigma(sigma, Mach=2., gamma=defg._gamma):
	if np.size(Mach)!=1:
		Minit = np.where(Mach<=1., __MachSub_sigma(sigma, gamma), __MachSup_sigma(sigma, gamma))
	else:
		Minit=Mach
	def sigma_of_mach(m):
		return Sigma_Mach(m, gamma)-sigma
	result = newton(sigma_of_mach, Minit)
	return result #if np.size(result)!=1 else np.asscalar(result)

def MachSub_Sigma(sigma, gamma=defg._gamma):
	# initial guess
	Mach = __MachSub_sigma(sigma, gamma)
	return Mach_Sigma(sigma, Mach, gamma)

def MachSup_Sigma(sigma, gamma=defg._gamma):
	# initial guess
	#cg = (gamma+1.)/(gamma-1.)
	#Mach = np.sqrt((sigma*cg**(.5*cg))**(2./(cg-1.))-cg)
	Mach = __MachSup_sigma(sigma, gamma)
	return Mach_Sigma(sigma, Mach, gamma)

