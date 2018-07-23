
"""
    The ``nozzle`` module
    =========================
 
    Provides specific functions for computation of 1D defintion of nozzles via section law
 
    :Example:
 
    >>> import hades.aero.nozzle as nz
    >>> Is.TiTs_Mach(1.)
    1.2
    >>> Is.TiTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
	.. note:: Specific heat ratio `gamma` is defined only using hades.common.defaultgas module
"""


import math
import numpy as np
from hades.common import defaultgas as defg # relative import is deprecated by doctest
from hades.aero   import Isentropic as Is
from hades.aero   import MassFlow   as mf
from hades.aero   import ShockWave  as sw

# === NPR computation from As/Ac definition of nozzle ===

def NPR_choked_subsonic(AsAc):
	return Is.PiPs_Mach(mf.MachSub_Sigma(AsAc))

def NPR_choked_supersonic(AsAc):
	return Is.PiPs_Mach(mf.MachSup_Sigma(AsAc))

def NPR_shock_at_exit(AsAc):
	Msup  = mf.MachSup_Sigma(AsAc)
	Msh   = sw.downstream_Mn(Msup)
	return Is.PiPs_Mach(Msh) / sw.Pi_ratio(Msup)

def _NPR_Ms_list(AsAc):
	"""
    	Computes all NPR limits and associated exit Mach number

		internal function
 
		:param AsAc:  ratio of section at exit over throat
		:return:      result NPR and Mach numbers
 
 		:Example:

		>>> import hades.aero.MassFlow as mf ; mf.Sigma_Mach(Is.Mach_PiPs(np.array(_NPR_Ms_list(2.)[:3:2])))
		array([ 2.,  2.])

		.. seealso:: NPR_choked_subsonic(), NPR_choked_supersonic(), NPR_shock_at_exit()
		.. note:: available for scalar or array (numpy) computations
    """
	Msub  = mf.MachSub_Sigma(AsAc)
	NPR0  = Is.PiPs_Mach(Msub)
	Msup  = mf.MachSup_Sigma(AsAc)
	Msh   = sw.downstream_Mn(Msup)
	NPRsw = Is.PiPs_Mach(Msh) / sw.Pi_ratio(Msup)
	NPR1  = Is.PiPs_Mach(Msup)
	return NPR0, NPRsw, NPR1, Msub, Msh, Msup

def Ms_from_AsAc_NPR(AsAc, NPR):
	"""
    	Computes Mach number at exit of a nozzle given As/Ac and NPR

		This method checks the NPR to define regime and computes Mach number at exit
 
		:param AsAc:  ratio of section at exit over throat
        :param NPR:   ratio of total pressure at inlet over 'expected' static pressure at exit
		:return:      result Mach number at exit
 
 		:Example:

		>>> print round(Ms_from_AsAc_NPR(2.636, 1.5), 8) # case with shock in diffuser
		0.32586574

		.. seealso:: 
		.. note:: NOT available for array (numpy) computations
    """	
	NPR0, NPRsw, NPR1, Msub, Msh, Msup = _NPR_Ms_list(AsAc)
	if (NPR < NPR0):
		Ms = Is.Mach_PiPs(NPR)
	elif (NPR > NPRsw): 
		Ms = Msup
	else:
		gmu = defg._gamma-1.
		K   = NPR/AsAc/((defg._gamma+1.)/2)**((defg._gamma+1.)/2/gmu)
		Ms  = np.sqrt((np.sqrt(1.+2.*gmu*K*K)-1.)/gmu)
	return Ms

def Madapt_from_AsAc_NPR(AsAc, NPR):
	"""
    	Computes Mach number for pressure adapted flow of a nozzle given As/Ac and NPR

		This method checks the NPR to define regime and computes Mach number in jet. The switch between 
		overexpanded jet and underexpanded jet is 
 
		:param AsAc:  ratio of section at exit over throat
        :param NPR:   ratio of total pressure at inlet over 'expected' static pressure at exit
		:return:      result Mach number at exit
 
 		:Example:

		>>> print round(Ms_from_AsAc_NPR(2.636, 1.5), 8) # case with shock in diffuser
		0.32586574

		.. seealso:: 
		.. note:: NOT available for array (numpy) computations
    """	
	NPR0, NPRsw, NPR1, Msub, Msh, Msup = _NPR_Ms_list(AsAc)
	if (NPR < NPR0):
		Ms = Is.Mach_PiPs(NPR)
	elif (NPR > NPR1): # under expanded flow
		Ms = Is.Mach_PiPs(NPR)
	elif (NPR > NPRsw): # shock wave in jet
		Ms = sw.downstreamMach_Mach_ShockPsratio(Msup, NPR1/NPR)
	else:
		gmu = defg._gamma-1.
		K   = NPR/AsAc/((defg._gamma+1.)/2)**((defg._gamma+1.)/2/gmu)
		Ms  = np.sqrt((np.sqrt(1.+2.*gmu*K*K)-1.)/gmu)
	return Ms



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()