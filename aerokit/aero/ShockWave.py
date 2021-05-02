"""@package ShockWave
  local Rankine Hugoniot equations for shock waves
"""

import math
import aerokit.aero.IterativeSolve as ITS
import numpy      as np
import aerokit.aero.Isentropic as Is
import aerokit.aero.degree as degree
from aerokit.common import defaultgas as defg # relative import is deprecated by doctest
from aerokit.common._ode import _rkf45

# --- NORMAL SHOCK WAVE ---

def Ps_ratio(Mn, gamma=defg._gamma):
    """ Computes static pressure ratio across a shockwave

    Args:
      Mn: upstream normal Mach number in shock reference frame
      gamma:  (Default value = defg._gamma)

    Returns:
        static pressure ratio

    """
    return 1.+ (2.*gamma/(gamma+1.))*(Mn**2-1.)

def Mn_Ps_ratio(Pratio: float, gamma=defg._gamma):
    """ Computes normal Mach number from static pressure ratio across a shockwave

    Args:
      Pratio: static pressure ratio
      gamma:  (Default value = defg._gamma)

    Returns:
        normal Mach number in shock reference frame

    """
    return np.sqrt(1+(Pratio-1.)*(gamma+1.)/(2.*gamma)) 

def Rho_ratio(Mn, gamma=defg._gamma):
    """

    Args:
      Mn: normal Mach number
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return ((gamma+1.)*Mn**2)/(2.+(gamma-1.)*Mn**2)

def Ts_ratio(Mn, gamma=defg._gamma):
    """

    Args:
      Mn: normal Mach number
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return Ps_ratio(Mn, gamma)/Rho_ratio(Mn, gamma)

def downstream_Mn(Mn, gamma=defg._gamma):
    """

    Args:
      Mn: param gamma:  (Default value = defg._gamma)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return np.sqrt((1.+.5*(gamma-1.)*Mn**2)/(gamma*Mn**2-.5*(gamma-1.)))

def Pt_ratio(Mn, gamma=defg._gamma):
    """Total pressure ration through shock wave

    Args:
      Mn: upstream relative normal Mach number to param 
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return Ps_ratio(Mn, gamma)*Is.PtPs_Mach(downstream_Mn(Mn, gamma))/Is.PtPs_Mach(Mn, gamma)

def Mn_Pt_ratio(ptratio, gamma=defg._gamma):
    """

    Args:
      ptratio: 
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    def ptratio_of_mach(m):
        """

        Args:
          m: 

        Returns:

        """
        return Pt_ratio(m, gamma)
    if ptratio > 1:
        print("!!! cannot find Mn for Pt_ratio > 1")
        return 1
    return ITS.secant_solve(ptratio_of_mach, ptratio, 1.5)

# --- LOCAL 2D SHOCK WAVE ---

def deflection_Mach_sigma(Mach, sigma, gamma=defg._gamma):
    """

    Args:
      Mach: param sigma:
      gamma: Default value = defg._gamma)
      sigma: 

    Returns:

    """
    sigrad = np.radians(sigma)
    Mn     = Mach*np.sin(sigrad)
    return sigma-np.degrees(np.arctan2(np.tan(sigrad), Rho_ratio(Mn, gamma)))

def deflection_Mach_ShockPsratio(Mach, Pratio, gamma=defg._gamma):
    """

    Args:
      Mach: param Pratio:
      gamma: Default value = defg._gamma)
      Pratio: 

    Returns:

    """
    return deflection_Mach_sigma(Mach, degree.asin(Mn_Ps_ratio(Pratio, gamma)/Mach), gamma)

def downstreamMach_Mach_ShockPsratio(Mach, Pratio, gamma=defg._gamma):
    """

    Args:
      Mach: param Pratio:
      gamma: Default value = defg._gamma)
      Pratio: 

    Returns:

    """
    Mn0 = Mn_Ps_ratio(Pratio, gamma)
    sig = degree.asin(Mn0/Mach)
    dev = deflection_Mach_sigma(Mach, sig, gamma)
    Mn1 = downstream_Mn(Mn0, gamma)
    return Mn1/degree.sin(sig-dev)

def weaksigma_Mach_deflection(Mach, deflection, gamma=defg._gamma):
    """

    Args:
      Mach: param deflection:
      gamma: Default value = defg._gamma)
      deflection: 

    Returns:

    """
    ka = (1. + .5*(gamma+1.)*Mach**2)*degree.tan(deflection)
    kb = 1. - Mach**2
    kc = (1. + .5*(gamma-1.)*Mach**2)*degree.tan(deflection)
    kd = (ka**2/3. - kb)/3.
    ke = 2.*ka**3/27. - ka*kb/3. + kc
    if ke**2 - 4.*kd**3 > 0:
        print("no weak shock wave solution")
        return degree.asin(1./Mach)
    else:
        phi = np.arccos(-.5*ke/np.sqrt(kd**3))
        kf  = 2.*np.sqrt(kd)*np.cos(phi/3.) - ka/3.
        return degree.atan(1./kf)

def strongsigma_Mach_deflection(Mach, deflection, gamma=defg._gamma):
    """

    Args:
      Mach: param deflection:
      gamma: Default value = defg._gamma)
      deflection: 

    Returns:

    """
    ka = (1. + .5*(gamma+1.)*Mach**2)*degree.tan(deflection)
    kb = 1. - Mach**2
    kc = (1. + .5*(gamma-1.)*Mach**2)*degree.tan(deflection)
    kd = (ka**2/3. - kb)/3.
    ke = 2.*ka**3/27. - ka*kb/3. + kc
    if ke**2 - 4.*kd**3 > 0:
        print("no strong shock wave solution")
        return 90.
    else:
        phi = np.arccos(-.5*ke/np.sqrt(kd**3)) + 4*math.pi
        kf  = 2.*np.sqrt(kd)*np.cos(phi/3.) - ka/3.
        return degree.atan(1./kf)

def sigma_Mach_deflection(Mach, deflection, init=None, gamma=defg._gamma):
    """computes the shock angle from upstream Mach number and deflection angle

    Args:
      Mach: param deflection:
      init: Default value = None)
      gamma: Default value = defg._gamma)
      deflection: 

    Returns:

    """
    if init is None:
        sig0 = degree.asin(1./Mach)+deflection
    else:
        sig0 = init
    def local_f(sig):
        """

        Args:
          sig: 

        Returns:

        """
        return deflection_Mach_sigma(Mach, sig, gamma)
    return ITS.secant_solve(local_f, np.abs(deflection), sig0)*np.sign(deflection) 

def dev_Max(Mach, gamma=defg._gamma):
    """computes the maximum deviation (always subsonic downstream flow), separation of weak/strong shock

    Args:
      Mach: param gamma:  (Default value = defg._gamma)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return deflection_Mach_sigma(Mach, sigma_DevMax(Mach, gamma=gamma), gamma=gamma)

def dev_Sonic(Mach, gamma=defg._gamma):
    """computes the deviation angle for a downstream SONIC Mach number

    Args:
      Mach: param gamma:  (Default value = defg._gamma)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    return deflection_Mach_sigma(Mach, sigma_Sonic(Mach, gamma=gamma), gamma=gamma)

def sigma_DevMax(Mach, gamma=defg._gamma):
    """computes the shock angle at maximum deviation (always subsonic downstream flow), separation of weak/strong shock

    Args:
      Mach: param gamma:  (Default value = defg._gamma)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    fogpu = 4./(gamma+1.)
    M2    = np.square(Mach)
    #ka = (M2-1.)*(1.+.5*(gamma-1.)*M2)
    #kb = .25*((gamma+1.)*M2-(3.-gamma))*M2 + 1.
    #return degree.atan(np.sqrt((kb+np.sqrt(kb**2+ka))/ka))
    return degree.asin(np.sqrt(1./fogpu/gamma/M2*(M2-fogpu+np.sqrt(M2*M2+2*(gamma-1.)*fogpu*M2+4.*fogpu))))

def sigma_Sonic(Mach, gamma=defg._gamma):
    """computes the shock angle for a downstream SONIC Mach number

    Args:
      Mach: param gamma:  (Default value = defg._gamma)
      gamma:  (Default value = defg._gamma)

    Returns:

    """
    M2 = np.square(Mach)
    ka = gamma-3+M2*(gamma+1)
    kb = (gamma+1)*(np.square(M2-3)+gamma*np.square(M2+1))
    return degree.asin(np.sqrt((ka+np.sqrt(kb))/(4*gamma*M2)))

# --- CONICAL SHOCK WAVE ---

def conical_deflection_Mach_sigma(Mach, sigma, gamma=defg._gamma, tol=1.0e-6):
    """

    Args:
      Mach: param sigma:
      gamma: Default value = defg._gamma)
      tol: Default value = 1.0e-6)
      sigma: 

    Returns:

    """
    def rhs(phi, data):
        """RHS for conical equations (2 equations)

        Args:
          phi: angle (coordinate of the ODE)
          data: theta (flow angle) and Mach number

        Returns: theta/mach variation rate (RHS)
        """
        th = data[0]
        ma = data[1]
        k  = 1.-(ma*degree.sin(phi-th))**2
        rhs=np.zeros(2)
        rhs[0] = -degree.sin(th)*degree.cos(phi-th)/degree.sin(phi)/k
        rhs[1] =  np.pi/180.*degree.sin(th)*degree.sin(phi-th)/degree.sin(phi)/k*ma*(1.+.5*(gamma-1)*ma*ma)
        return rhs

    # initial data for integration from downstream shock to wall
    th   = deflection_Mach_sigma(Mach, sigma, gamma)
    ma   = downstream_Mn(Mach*degree.sin(sigma), gamma)/degree.sin(sigma-th)
    phi  = sigma
    thma = np.array([th, ma])
    h    = -phi/10.
    conv = phi
    while (abs(conv) >= tol):
        dthma, err = _rkf45(rhs, phi, thma, h)
        # Accept integration step if error e is within tolerance
        if err <= tol:
            conv = thma[0] - phi
            if conv/(h-dthma[0]) < 1:
                h = h*conv/(h-dthma[0])
            else:
                thma = thma + dthma
                phi  = phi  + h
                #print(Mach, sigma, phi, thma)
        else:
            h = 0.9*h*(tol/err)**0.2
            #print("new h ",h)

    deflection = .5*(thma[0] + phi)
    return deflection

def conical_sigma_Mach_walldeflection(Mach, deflection, gamma=defg._gamma):
    """computes shock angle sigma from upstream Mach number and wall deflection

    Args:
      Mach: param deflection:
      gamma: Default value = defg._gamma)
      deflection: 

    Returns:

    """
    def local_def(sig):
        """internal wrapping function to iterative solve
        """
        return conical_deflection_Mach_sigma(Mach, sig, gamma)
    return ITS.secant_solve(local_def, deflection, degree.asin(1./Mach)+deflection)

def conical_Mach_walldeflection_sigma(deflection, sigma, gamma=defg._gamma):
    """computes upstream Mach number from wall deflection and shock angle sigma 

    Args:
      Mach: param deflection:
      gamma: Default value = defg._gamma)
      deflection: 

    Returns:

    """
    assert sigma > deflection
    def local_def(mach):
        """internal wrapping function to iterative solve
        """
        #print("local mach sigma", mach, sigma)
        return conical_deflection_Mach_sigma(mach, sigma, gamma)
    return ITS.secant_solve(local_def, deflection, 1./degree.sin(sigma-.5*deflection))

# backward compatibility
Pi_ratio = Pt_ratio
Mn_Pi_ratio = Mn_Pt_ratio