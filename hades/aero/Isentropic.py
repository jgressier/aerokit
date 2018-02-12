"""
    The ``Isentropic`` module
    =========================
 
    Provides standart Stagnation to Static Ratios for Ideal gas
 
    :Example:
 
    >>> import hades.aero.Isentropic as Is
    >>> Is.TiTs_Mach(1.)
    1.2
    >>> Is.TiTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
    Provides Ti/Ts Pi/Ps ratios from Mach number and reversed functions.
    Specific heat ratio `gamma` is optionnal and can be specified in the functions itself
    or using hades.common.defaultgas module
 """

import math
import numpy as np
from hades.common import defaultgas as defg # relative import is deprecated by doctest

# ===============================================================
# implemented functions

def TiTs_Mach(Mach, gamma=defg._gamma):
    """
    	Computes Ti/Ts ratio from Mach number

		Long comment
 
		:param Mach:  local Mach number
        :param gamma: specific heat ratio, default from hades.common.defaultgas
		:return:      result Ti/Ts ratio
 
 		:Example:

		>>> TiTs_Mach(1.) # with default gamma 1.4
		1.2
		>>> TiTs_Mach(2., gamma=1.6)
		2.2

		.. seealso:: Mach_TiTs()
		.. note:: available for scalar or array (numpy) computations
    """	
    return 1.+.5*(gamma-1)*Mach**2

def PiPs_Mach(Mach, gamma=defg._gamma):
    return (1.+.5*(gamma-1.)*Mach**2)**(gamma/(gamma-1.))

def Mach_TiTs(TiTs, gamma=defg._gamma):
    return np.sqrt((TiTs-1.)*2./(gamma-1.))

def Mach_PiPs(PiPs, gamma=defg._gamma):
    return np.sqrt((PiPs**((gamma-1.)/gamma)-1.)*2./(gamma-1.))

def Velocity_MachTi(Mach, Ti, r=287.1, gamma=defg._gamma):
    return Mach*np.sqrt(gamma*r*Ti/TiTs_Mach(Mach, gamma))



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()