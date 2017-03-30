# -*- coding: utf-8 -*-
"""
Description
-----------
Module containing the definition extended variables according to vtk calculator syntax

Contains
--------

Examples
--------

Notes
-----
"""
from dicovar import dicoVar as dico

# -----------------------------------------------------------------

def expressions_from_rhopv():
    extdico = {}
    # copy dicovar (same as deepcopy)
    for var in dico.keys():
        #print 'copy ',var, dico[var]
        extdico[var] = dico[var]
    # additional variables
    _g     = '1.4'
    _gsgmu = '3.5'
    _gmusd = '0.2'
    _a     = 'sqrt('+_g+'*'+dico['P']+'/'+dico['rho']+')'
    extdico['Pt']   = dico['P']+'*(1+'+_gmusd+'*(mag('+dico['V']+'))^2/('+_g+'*'+dico['P']+'/'+dico['rho']+'))^'+_gsgmu
    extdico['T']    = dico['P']+'/'+dico['rho']+'/287.1'
    extdico['S']    = 'ln('+dico['P']+'/'+dico['rho']+'^'+_g+')'
    extdico['V']    = 'mag('+dico['V']+')'
    extdico['a' ]   = _a
    extdico['Mach'] = 'mag('+dico['V']+')/'+_a
    extdico['Mx']   = 'U_X/'+_a
    extdico['alpha'] = 'atan('+dico['Vy']+'/'+dico['Vx']+')*180/acos(-1)'
    return extdico