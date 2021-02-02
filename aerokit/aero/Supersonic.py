"""
    The ``Isentropic`` module
    =========================
 
    Provides supersonic flow functions
 
    :Example:
 
    >>> import aerokit.aero.Supersonic as sup
    >>> Is.TtTs_Mach(1.)
    1.2
    >>> Is.TtTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
 """

import math
import numpy as np
from aerokit.aero import IterativeSolve
from aerokit.aero import Isentropic as Is
from aerokit.aero import degree     as deg
from aerokit.common import defaultgas as defg # relative import is deprecated by doctest
from scipy.optimize import newton

# -- 2D supersonic invariants --

def PrandtlMeyer_Mach(Mach, gamma=defg._gamma):
    """
        Computes Prandtl-Meyer (or Busemann) function from Mach number

        Long comment
 
        :param Mach:  local Mach number
        :param gamma: specific heat ratio, default from aerokit.common.defaultgas
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
    m2 = Is.Mach_PtPs(Is.PtPs_Mach(Mach, gamma)/Pratio, gamma)
    return -PrandtlMeyer_Mach(Mach, gamma)+PrandtlMeyer_Mach(m2, gamma)

def IsentropicPsratio_Mach_deflection(Mach, dev, gamma=defg._gamma):
    m2 = Mach_PrandtlMeyer(PrandtlMeyer_Mach(Mach, gamma)-dev, gamma)
    return Is.PtPs_Mach(Mach, gamma)/Is.PtPs_Mach(m2, gamma)



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()