# python module: ThermoKerosen

#import math

r_air = 287.04

def Cp(T):
    return r_air*(4.47659+(8.01994e-3-1.8373e-6*T)*T)

def Enthalpy(T):
    return r_air*(-149.054+(4.47659+(4.00997e-3-6.12432e-7*T)*T)*T)

def Pc_inf():
    return 10300.*4184

def Pc_eff(T):
    return Pc_inf()+r_air*(1607.2+(-4.47659+(-4.00997e-3+6.12432e-7*T)*T)*T)
