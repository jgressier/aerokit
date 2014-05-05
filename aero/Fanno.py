"""@package Fanno
@brief functions to compute Fanno transformation in subsonic/supersonic ducts
"""

import numpy as np

def maxFparam_Mach(Mach, gamma=1.4):
    """ computes maximum Fanno parameter for a given Mach number state
    This maximum value makes the flow reach a sonic state
    \f$
      \left(\frac{1 - M^2}{\gamma M^2}\right) + \left(\frac{\gamma + 1}{2\gamma}\right)\ln\left[\frac{M^2}{\left(\frac{2}{\gamma + 1}\right)\left(1 + \frac{\gamma - 1}{2}M^2\right)}\right]
    \f$
    """
    m2 = np.square(Mach)
    return (1-m2)/(gamma*m2) + (gamma+1)/(2*gamma)*np.log(m2*(gamma+1)/2/(1+(gamma-1)/2*m2))

# --- ratio to critical/choking state ---

def Ps_Pscri(Mach, gamma=1.4):
    """ computes Ps/Ps_cri ratio from Mach number
    """
    return 1./(Mach*np.sqrt(2/(gamma+1)*(1+(gamma-1)/2*np.square(Mach))))

def Ts_Tscri(Mach, gamma=1.4):
    """ computes Ts/Ts_cri ratio from Mach number
    """
    return (gamma+1)/2/(1+(gamma-1)/2*np.square(Mach))
    
def Rho_Rhocri(Mach, gamma=1.4):
    """ computes rho/rho_cri ratio from Mach number
    """
    return 1./Mach*np.sqrt(2/(gamma+1)*(1+(gamma-1)/2*np.square(Mach)))
    
def Pi_Picri(Mach, gamma=1.4):
    """ computes Pi/Pi_cri ratio from Mach number
    """
    return 1./Mach*(2/(gamma+1)*(1+(gamma-1)/2*np.square(Mach)))**((gamma+1)/2/(gamma-1))

def V_Vcri(Mach, gamma=1.4):
    """ computes V/a_cri ratio from Mach number
    """
    return Mach / np.sqrt(2/(gamma+1)*(1+(gamma-1)/2*np.square(Mach)))

def NormdS(Mach, gamma=1.4):
    """ computes delta S / Cp up to critical point
    """
    return np.log(Mach**((gamma-1)/gamma) * (2/(gamma+1)*(1+(gamma-1)/2*np.square(Mach)))**(-(gamma+1)/2/gamma))

