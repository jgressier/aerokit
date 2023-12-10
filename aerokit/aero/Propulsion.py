"""@package Propulsion
  thrust flow function
"""
import aerokit.aero.Isentropic as Is
import aerokit.aero.MassFlow as mf
import aerokit.IterativeSolve as ITS


def ThrustFunction(Mach, gamma=1.4):
    return mf.Sigma_Mach(Mach, gamma) / Is.PtPs_Mach(Mach, gamma) * (1.0 + gamma * Mach ** 2)


def Mach_ThrustFunction(thrust, Mach=2.0, gamma=1.4):
    def thrust_of_mach(m):
        return ThrustFunction(m, gamma)

    return ITS.secant_solve(thrust_of_mach, thrust, Mach)
