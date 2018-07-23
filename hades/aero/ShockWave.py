"""@package ShockWave
  locla Rankine Hugoniot equations for shock waves
"""

import math
import degree
import IterativeSolve
import numpy      as np
import Isentropic as Is
from hades.common import defaultgas as defg # relative import is deprecated by doctest

# --- NORMAL SHOCK WAVE ---

def Ps_ratio(Mn, gamma=defg._gamma):
    return 1.+ (2.*gamma/(gamma+1.))*(Mn**2-1.)

def Mn_Ps_ratio(Pratio, gamma=defg._gamma):
    return np.sqrt(1+(Pratio-1.)*(gamma+1.)/(2.*gamma)) 

def Rho_ratio(Mn, gamma=defg._gamma):
    return ((gamma+1.)*Mn**2)/(2.+(gamma-1.)*Mn**2)

def Ts_ratio(Mn, gamma=defg._gamma):
    return Ps_ratio(Mn, gamma)/Rho_ratio(Mn, gamma)

def downstream_Mn(Mn, gamma=defg._gamma):
    return np.sqrt((1.+.5*(gamma-1.)*Mn**2)/(gamma*Mn**2-.5*(gamma-1.)))

def Pi_ratio(Mn, gamma=defg._gamma):
    return Ps_ratio(Mn, gamma)*Is.PiPs_Mach(downstream_Mn(Mn, gamma))/Is.PiPs_Mach(Mn, gamma)

def Mn_Pi_ratio(piratio, gamma=defg._gamma):
    def piratio_of_mach(m):
        return Pi_ratio(m, gamma)
    if piratio > 1:
        print "!!! cannot find Mn for Pi_ratio > 1"
        return 1
    return IterativeSolve.secant_solve(piratio_of_mach, piratio, 1.5)

# --- LOCAL 2D SHOCK WAVE ---

def deflection_Mach_sigma(Mach, sigma, gamma=defg._gamma):
    sigrad = np.radians(sigma)
    Mn     = Mach*np.sin(sigrad)
    return sigma-np.degrees(np.arctan2(np.tan(sigrad), Rho_ratio(Mn, gamma)))

def deflection_Mach_ShockPsratio(Mach, Pratio, gamma=defg._gamma):
    return deflection_Mach_sigma(Mach, degree.asin(Mn_Ps_ratio(Pratio, gamma)/Mach), gamma)

def downstreamMach_Mach_ShockPsratio(Mach, Pratio, gamma=defg._gamma):
    Mn0 = Mn_Ps_ratio(Pratio, gamma)
    sig = degree.asin(Mn0/Mach)
    dev = deflection_Mach_sigma(Mach, sig, gamma)
    Mn1 = downstream_Mn(Mn0, gamma)
    return Mn1/degree.sin(sig-dev)

def weaksigma_Mach_deflection(Mach, deflection, gamma=defg._gamma):
    ka = (1. + .5*(gamma+1.)*Mach**2)*degree.tan(deflection)
    kb = 1. - Mach**2
    kc = (1. + .5*(gamma-1.)*Mach**2)*degree.tan(deflection)
    kd = (ka**2/3. - kb)/3.
    ke = 2.*ka**3/27. - ka*kb/3. + kc
    if ke**2 - 4.*kd**3 > 0:
        print "no weak shock wave solution"
        return degree.asin(1./Mach)
    else:
        phi = np.acos(-.5*ke/np.sqrt(kd**3))
        kf  = 2.*np.sqrt(kd)*np.cos(phi/3.) - ka/3.
        return degree.atan(1./kf)

def strongsigma_Mach_deflection(Mach, deflection, gamma=defg._gamma):
    ka = (1. + .5*(gamma+1.)*Mach**2)*degree.tan(deflection)
    kb = 1. - Mach**2
    kc = (1. + .5*(gamma-1.)*Mach**2)*degree.tan(deflection)
    kd = (ka**2/3. - kb)/3.
    ke = 2.*ka**3/27. - ka*kb/3. + kc
    if ke**2 - 4.*kd**3 > 0:
        print "no strong shock wave solution"
        return 90.
    else:
        phi = np.acos(-.5*ke/np.sqrt(kd**3)) + 4*math.pi
        kf  = 2.*np.sqrt(kd)*np.cos(phi/3.) - ka/3.
        return degree.atan(1./kf)

def sigma_Mach_deflection(Mach, deflection, gamma=defg._gamma):
    """
    computes the shock angle from upwstream Mach number and deflection angle
    """
    def local_f(sig):
        return deflection_Mach_sigma(Mach, sig, gamma)
    return IterativeSolve.secant_solve(local_f, deflection, degree.asin(1./Mach)+deflection)

def dev_Max(Mach, gamma=defg._gamma):
    """ computes the maximum deviation (always subsonic downstream flow), separation of weak/strong shock
    """
    return deflection_Mach_sigma(Mach, sigma_DevMax(Mach, gamma=gamma), gamma=gamma)

def dev_Sonic(Mach, gamma=defg._gamma):
    """ computes the deviation angle for a downstream SONIC Mach number
    """
    return deflection_Mach_sigma(Mach, sigma_Sonic(Mach, gamma=gamma), gamma=gamma)

def sigma_DevMax(Mach, gamma=defg._gamma):
    """ computes the shock angle at maximum deviation (always subsonic downstream flow), separation of weak/strong shock
    """
    fogpu = 4./(gamma+1.)
    M2    = np.square(Mach)
    #ka = (M2-1.)*(1.+.5*(gamma-1.)*M2)
    #kb = .25*((gamma+1.)*M2-(3.-gamma))*M2 + 1.
    #return degree.atan(np.sqrt((kb+np.sqrt(kb**2+ka))/ka))
    return degree.asin(np.sqrt(1./fogpu/gamma/M2*(M2-fogpu+np.sqrt(M2*M2+2*(gamma-1.)*fogpu*M2+4.*fogpu))))

def sigma_Sonic(Mach, gamma=defg._gamma):
    """ computes the shock angle for a downstream SONIC Mach number
    """
    M2 = np.square(Mach)
    ka = gamma-3+M2*(gamma+1)
    kb = (gamma+1)*(np.square(M2-3)+gamma*np.square(M2+1))
    return degree.asin(np.sqrt((ka+np.sqrt(kb))/(4*gamma*M2)))

# --- CONICAL SHOCK WAVE ---

def conical_deflection_Mach_sigma(Mach, sigma, gamma=defg._gamma, tol=1.0e-6):
    def rkf45(F, x, y, h):
        # Runge-Kutta-Fehlberg formulas
        C = [37./378, 0., 250./621, 125./594, 0., 512./1771]
        D = [2825./27648, 0., 18575./48384, 13525./55296, 277./14336, 1./4]
        n = len(y)
        K = numpy.zeros((6,n))
        K[0] = h*F(x,y)
        K[1] = h*F(x + 1./5*h, y + 1./5*K[0])
        K[2] = h*F(x + 3./10*h, y + 3./40*K[0] + 9./40*K[1])
        K[3] = h*F(x + 3./5*h, y + 3./10*K[0]- 9./10*K[1] + 6./5*K[2])
        K[4] = h*F(x + h, y - 11./54*K[0] + 5./2*K[1] - 70./27*K[2] + 35./27*K[3])
        K[5] = h*F(x + 7./8*h, y + 1631./55296*K[0] + 175./512*K[1] + 575./13824*K[2] + 44275./110592*K[3] + 253./4096*K[4])
        # Initialize arrays {dy} and {E}
        E  = numpy.zeros((n))
        dy = numpy.zeros((n))
        # Compute solution increment {dy} and per-step error {E}
        for i in range(6):
            dy = dy + C[i]*K[i]
            E  = E + (C[i] - D[i])*K[i]
        # Compute RMS error e
        e = np.sqrt(sum(E**2)/n)
        return dy, e

    def rhs(phi, data):
        th = data[0]
        ma = data[1]
        k  = 1.-(ma*degree.sin(phi-th))**2
        rhs=numpy.zeros(2)
        rhs[0] = -degree.sin(th)*degree.cos(phi-th)/degree.sin(phi)/k
        rhs[1] =  np.pi/180.*degree.sin(th)*degree.sin(phi-th)/degree.sin(phi)/k*ma*(1.+.5*(gamma-1)*ma*ma)
        return rhs

    th   = deflection_Mach_sigma(Mach, sigma, gamma)
    ma   = downstream_Mn(Mach*degree.sin(sigma), gamma)/degree.sin(sigma-th)
    phi  = sigma
    thma = numpy.array([th, ma])
    h    = -phi/10.
    conv = phi
    while (abs(conv) >= tol):
        dthma, err = rkf45(rhs, phi, thma, h)
        # Accept integration step if error e is within tolerance
        if err <= tol:
            conv = thma[0] - phi
            if conv/(h-dthma[0]) < 1:
                h = h*conv/(h-dthma[0])
            else:
                thma = thma + dthma
                phi  = phi  + h
                print phi, thma
        else:
            h = 0.9*h*(tol/err)**0.2
            print "new h ",h

    deflection = .5*(thma[0] + phi)
    return deflection

def conical_sigma_Mach_walldeflection(Mach, deflection, gamma=defg._gamma):
    def local_f(sig):
        return conical_deflection_Mach_sigma(Mach, sig, gamma)
    return IterativeSolve.secant_solve(local_f, deflection, degree.asin(1./Mach)+deflection)

