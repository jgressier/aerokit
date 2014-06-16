"""@package Rayleigh
@brief functions to compute Rayleigh (heating) transformation in subsonic/supersonic ducts
"""

import numpy as np

def maxTiratio_Mach(Mach, gamma=1.4):
    """ computes maximum heating ratio parameter for a given Mach number state
    This maximum value makes the flow reach a sonic state
    """
    return 1./Ti_Ticri(Mach, gamma=1.4)

# --- ratio to critical/choking state ---

def Ps_Pscri(Mach, gamma=1.4):
    """ computes Ps/Ps_cri ratio from Mach number
    """
    return (gamma+1)/(1+gamma*np.square(Mach))

def Ts_Tscri(Mach, gamma=1.4):
    """ computes Ts/Ts_cri ratio from Mach number
    """
    m2 = np.square(Mach)
    return m2*np.square(((gamma+1)/(1+gamma*m2)))
    
def Ti_Ticri(Mach, gamma=1.4):
    """ computes Ti/Ti_cri ratio from Mach number
    """
    m2 = np.square(Mach)
    return 2*(gamma+1)*m2 / np.square(1+gamma*m2) * (1+(gamma-1)/2*m2)
    
def Rho_Rhocri(Mach, gamma=1.4):
    """ computes rho/rho_cri ratio from Mach number
    """
    m2 = np.square(Mach)
    return (1+gamma*m2)/((gamma+1)*m2)
    
def Pi_Picri(Mach, gamma=1.4):
    """ computes Pi/Pi_cri ratio from Mach number
    """
    m2 = np.square(Mach)
    return (gamma+1)/(1+gamma*m2)*(2/(gamma+1)*(1+(gamma-1)/2*m2))**(gamma/(gamma-1))

def V_Vcri(Mach, gamma=1.4):
    """ computes V/a_cri ratio from Mach number
    """
    m2 = np.square(Mach)
    return m2*(gamma+1)/(1+gamma*m2)

def NormdS(Mach, gamma=1.4):
    """ computes delta S / Cp up to critical point
    """
    m2 = np.square(Mach)
    return np.log(m2*((gamma+1)/(1+gamma*m2))**((gamma+1)/gamma))

def SubMach_TiTicri(Tiratio, gamma=1.4):
    """ computes Mach number from Ti/Ti_cri ratio (<1)
    """
    alpha = np.sqrt(1.-Tiratio)
    return np.sqrt((1.-alpha)/(alpha*gamma+1))

def SupMach_TiTicri(Tiratio, gamma=1.4):
    """ computes Mach number from Ti/Ti_cri ratio (<1)
    """
    alpha = np.sqrt(1.-Tiratio)
    return np.sqrt((1.+alpha)/(1.-alpha*gamma))

