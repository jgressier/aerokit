"""
    The ``Isentropic`` module
    =========================
 
    Provides standart Stagnation to Static Ratios for Ideal gas
 
    :Example:
 
    >>> import aerokit.aero.Isentropic as Is
    >>> Is.TtTs_Mach(1.)
    1.2
    >>> Is.TtTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
    Provides Tt/Ts Pt/Ps ratios from Mach number and reversed functions.
    Specific heat ratio `gamma` is optional and can be specified in the functions itself
    or using aerokit.common.defaultgas module
 """

import math
import numpy as np
from aerokit.common import defaultgas as defg # relative import is deprecated by doctest

# ===============================================================
# implemented functions

def TtTs_Mach(Mach, gamma=defg._gamma):
    """Computes Tt/Ts ratio from Mach number

    Args:
      Mach: local Mach number (either scalar or numpy array)
      gamma: specific heat ratio, default from aerokit.common.defaultgas

    Returns:
      result Tt/Ts ratio
      
    Example:
      
        >>> TtTs_Mach(1.) # with default gamma 1.4
        1.2
    """
    return 1.+.5*(gamma-1)*Mach**2

def PtPs_Mach(Mach, gamma=defg._gamma):
    """Computes Pt/Ps ratio from Mach number
    
    		Long comment

    Args:
      Mach: local Mach number (either scalar or numpy array)
      gamma: specific heat ratio, default from aerokit.common.defaultgas

    Returns:
      result Pt/Ps ratio
      
    Example:

        >>> PtPs_Mach(1.) # with default gamma 1.4
        1.892929158737854
        >>> PtPs_Mach(2., gamma=1.6)
        8.187044460255244
    """
    return (1.+.5*(gamma-1.)*Mach**2)**(gamma/(gamma-1.))

def Mach_TtTs(TtTs, gamma=defg._gamma):
    """Computes Mach number from Tt/Ts ratio

    Args:
      TtTs: local Total over Static temperature ratio (either scalar or numpy array)
      gamma:  (Default value = defg._gamma)

    Returns:
        Mach number (same type as TtTs)

    """
    return np.sqrt((TtTs-1.)*2./(gamma-1.))

def Mach_PtPs(PtPs, gamma=defg._gamma):
    """Computes Mach number from Pt/Ps ratio

    Args:
      PtPs: local Total over Static pressure ratio (either scalar or numpy array)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return np.sqrt((PtPs**((gamma-1.)/gamma)-1.)*2./(gamma-1.))

def Velocity_MachTt(Mach, Tt, r=defg._r, gamma=defg._gamma):
    """Computes velocity from 

    Args:
      Mach: Mach number
      Tt: Total temperature
      r:  (Default value = defg._r)
      gamma:  (Default value = defg._gamma)

    Returns:
        Velocity

    """
    return Mach*np.sqrt(gamma*r*Tt/TtTs_Mach(Mach, gamma))

# backward compatibility
TiTs_Mach = TtTs_Mach
PiPs_Mach = PtPs_Mach
Mach_PiPs = Mach_PtPs
Mach_TiTs = Mach_TtTs
Velocity_MachTi = Velocity_MachTt

# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()