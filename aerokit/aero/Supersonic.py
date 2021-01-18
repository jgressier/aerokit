"""
    The ``Isentropic`` module
    =========================
 
    Provides supersonic flow functions
 
    :Example:
 
    >>> import hades.aero.Supersonic as sup
    >>> Is.TiTs_Mach(1.)
    1.2
    >>> Is.TiTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
 """

import math
import numpy as np
from hades.aero import IterativeSolve
from hades.aero import Isentropic as Is
from hades.aero import degree     as deg
from hades.common import defaultgas as defg # relative import is deprecated by doctest
from scipy.optimize import newton

# -- 2D supersonic invariants --

def PrandtlMeyer_Mach(Mach, gamma=defg._gamma):
    """
        Computes Prandtl-Meyer (or Busemann) function from Mach number

        Long comment
 
        :param Mach:  local Mach number
        :param gamma: specific heat ratio, default from hades.common.defaultgas
        :return:      result omega (degree)
 
        :Example:

        >>> PrandtlMeyer_Mach(2.) # with default gamma 1.4
        26.379760813416446

        .. seealso:: 
        .. note:: available for scalar or array (numpy) computations
    """     
    g    = gamma
    gm1  = gamma - 1.
    beta = np.sqrt(Mach**2-1.)
    cg   = np.sqrt((g+1.)/gm1)
    return (cg*np.arctan(beta/cg) - np.arctan(beta))*(180./math.pi)

def old_Mach_PrandtlMeyer(omega, gamma=defg._gamma):
    def omega_of_mach(m):
        return PrandtlMeyer_Mach(m, gamma)
    return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)

def Mach_PrandtlMeyer(omega, gamma=defg._gamma):
    def omega_of_mach(m):
        return PrandtlMeyer_Mach(m, gamma)-omega
    #return IterativeSolve.secant_solve(omega_of_mach, omega, 2.)
    result = newton(omega_of_mach, 2.+0.*omega)
    return result if np.size(result)!=1 else np.asscalar(result)

def Mach_PMFmmu(om_m_mun, gamma=defg._gamma):
    "return Mach number from omega-mu value (in degree)"
    g    = gamma
    gm1  = gamma - 1.
    cg   = np.sqrt((g+1.)/gm1)
    return np.sqrt((cg*deg.tan((om_m_mun + 90.)/cg))**2+1.)

def deflection_Mach_IsentropicPsratio(Mach, Pratio, gamma=defg._gamma):
    m2 = Is.Mach_PiPs(Is.PiPs_Mach(Mach, gamma)/Pratio, gamma)
    return -PrandtlMeyer_Mach(Mach, gamma)+PrandtlMeyer_Mach(m2, gamma)

def IsentropicPsratio_Mach_deflection(Mach, dev, gamma=defg._gamma):
    m2 = Mach_PrandtlMeyer(PrandtlMeyer_Mach(Mach, gamma)-dev, gamma)
    return Is.PiPs_Mach(Mach, gamma)/Is.PiPs_Mach(m2, gamma)



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()