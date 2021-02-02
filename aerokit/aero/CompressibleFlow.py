"""@package  CompressibleFlow
  legacy package importing legacy functions
"""

from . import Isentropic as ist
from . import Supersonic as sup
from . import MassFlow   as mfl

# -- isentropic functions --

def TtTs_Mach(Mach, gamma=1.4):
    """ legacy interface to same function in Isentropic
    """
    return ist.TtTs_Mach(Mach, gamma)

def PtPs_Mach(Mach, gamma=1.4):
    """ legacy interface to same function in Isentropic
    """
    return ist.PtPs_Mach(Mach, gamma)

def Mach_TtTs(TtTs, gamma=1.4):
    """ legacy interface to same function in Isentropic
    """
    return ist.Mach_TtTs(TtTs, gamma)

def Mach_PtPs(PtPs, gamma=1.4):
    """ legacy interface to same function in Isentropic
    """
    return ist.Mach_PtPs(PtPs, gamma)

def Velocity_MachTi(Mach, Ti, r=287.1, gamma=1.4):
    """ legacy interface to same function in Isentropic
    """
    return ist.Velocity_MachTi(Mach, Ti, r, gamma)

# -- 2D supersonic invariants --

def PrandtlMeyer_Mach(Mach, gamma=1.4):
    """ legacy interface to same function in Supersonic
    """
    return sup.PrandtlMeyer_Mach(Mach, gamma)

def Mach_PrandtlMeyer(omega, gamma=1.4):
    """ legacy interface to same function in Supersonic
    """
    return sup.Mach_PrandtlMeyer(omega, gamma)

def deflection_Mach_IsentropicPsratio(Mach, Pratio, gamma=1.4):
    """ legacy interface to same function in Supersonic
    """
    return sup.deflection_Mach_IsentropicPsratio(Mach, Pratio, gamma)

def IsentropicPsratio_Mach_deflection(Mach, dev, gamma=1.4):
    """ legacy interface to same function in Supersonic
    """
    return sup.IsentropicPsratio_Mach_deflection(Mach, dev, gamma)

# -- Compressible flow functions  --

def WeightMassFlow(Mach, r=287.1, gamma=1.4):
    """ legacy interface to same function in MassFlow
    """
    return mfl.WeightMassFlow(Mach, r, gamma)

def Sigma_Mach(Mach, gamma=1.4):
    """ legacy interface to same function in MassFlow
    """
    return mfl.Sigma_Mach(Mach, gamma)

def Mach_Sigma(sigma, Mach=2., gamma=1.4):
    """ legacy interface to same function in MassFlow
    """
    return mfl.Mach_Sigma(sigma, Mach, gamma)

