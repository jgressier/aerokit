"""@package Propulsion
  thrust flow function
"""

import Isentropic
import MassFlow

def ThrustFunction(Mach, gamma=1.4):
    return MassFlow.Sigma_Mach(Mach, gamma)/Isentropic.PtPs_Mach(Mach, gamma)*(1.+gamma*Mach**2)

def Mach_ThrustFunction(thrust, Mach=2., gamma=1.4):
    def thrust_of_mach(m):
        return ThrustFunction(m, gamma)
    return IterativeSolve.secant_solve(thrust_of_mach, thrust, Mach)

