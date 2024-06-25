# python module: ThermoRealGas

import math

rconst = 287.04

def Gamma(T):
    a  = math.exp(3090./T)
    cp = ((3.5+(-2.8e-5+2.240e-8*T)*T)+(3090./T)**2 *a/(a-1.)**2)
    return cp/(cp-1.)

def Cp(T):
    a = math.exp(3090./T)
    return rconst*((3.5+(-2.8e-5+2.240e-8*T)*T)+(3090./T)**2 *a/(a-1.)**2)

def Enthalpy(T):
    a = math.exp(3090./T)
    return rconst*((3.5+(-1.4e-5+7.467e-9*T)*T)*T+3090./(a-1.)) 

def Phi(T):
    a = math.exp(3090./T)
    return (((3.5*math.log(T)+(-2.8e-5+1.12e-8*T)*T)+3090./T/(a-1.))-math.log((a-1.)/a))/math.log(10.)
