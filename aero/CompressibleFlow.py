# Python module : CompressibleFlow

import math
import numpy as np
import IterativeSolve

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

# -- Compressible flow functions  --

def WeightMassFlow(Mach, r=287.1, gamma=1.4):
    return np.sqrt(gamma/r)*Mach*(1.+.5*(gamma-1)*Mach**2)**(-.5*(gamma+1.)/(gamma-1.))

def Sigma_Mach(Mach, gamma=1.4):
    return (2./(gamma+1.)*(1.+.5*(gamma-1)*Mach**2))**(.5*(gamma+1.)/(gamma-1.))/Mach

def Mach_Sigma(sigma, Mach=2., gamma=1.4):
    def sigma_of_mach(m):
        return Sigma_Mach(m, gamma)
    return IterativeSolve.secant_solve(sigma_of_mach, sigma, Mach)

# -- 2D supersonic invariants --

def PrandtlMeyer_Mach(Mach, gamma=1.4):
    g    = gamma
    gm1  = gamma - 1.
    beta = np.sqrt(Mach**2-1.)
    return (np.sqrt((g+1.)/gm1)*np.arctan(np.sqrt(gm1/(g+1.))*beta) - np.arctan(beta))*180./math.pi

def Mach_PrandtlMeyer(omega, gamma=1.4):
    def omega_of_mach(m):
        return PrandtlMeyer_Mach(m, gamma)
    return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)
#    return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)

def deflection_Mach_IsentropicPsratio(Mach, Pratio, gamma=1.4):
    m2 = Mach_PiPs(PiPs_Mach(Mach, gamma)/Pratio, gamma)
    return -PrandtlMeyer_Mach(Mach, gamma)+PrandtlMeyer_Mach(m2, gamma)

def IsentropicPsratio_Mach_deflection(Mach, dev, gamma=1.4):
    m2 = Mach_PrandtlMeyer(PrandtlMeyer_Mach(Mach, gamma)-dev, gamma)
    return PiPs_Mach(Mach, gamma)/PiPs_Mach(m2, gamma)

# -- Thrust functions --

def ThrustFunction(Mach, gamma=1.4):
    return Sigma_Mach(Mach, gamma)/PiPs_Mach(Mach, gamma)*(1.+gamma*Mach**2)

def Mach_ThrustFunction(thrust, Mach=2., gamma=1.4):
    def thrust_of_mach(m):
        return ThrustFunction(m, gamma)
    return IterativeSolve.secant_solve(thrust_of_mach, thrust, Mach)
