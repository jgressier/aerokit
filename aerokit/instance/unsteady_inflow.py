""" 
	The ``inflow`` module
	=========================
 
	Provides class and function for inflow solvers for 1D unsteady flows
 
	:Example:
 
    >>> import aerokit.aero.Isentropic as Is
    >>> Is.TtTs_Mach(1.)
    1.2
    >>> Is.TtTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
    Blabla
  
"""

import math
import numpy     as np
from scipy.optimize import fsolve
import aerokit.aero.ShockWave  as sw
import aerokit.aero.unsteady1D as uq

# ===============================================================
# implementation of inflow_PB class

class inflow_pb():
	"""
		Class inflow_pb for unsteady 1D Euler

		attributes:
			_qinit  : initial state before flow variation, is qR state for equivalent Riemann problem
			_pstar  : static pressure at inflow
			_ustar  : velocity at inflow
			_qstarL : q state on time axis
			_qstarR : q state between q*L and qR
	"""
	def __init__(self, qinit):
	 	self._qR    = qinit.copy()

	def __repr__(self):
		return "(rho, u, p)*_L : (%s, %s, %s)\n" % (self._qstarL.rho, self._qstarL.u, self._qstarL.p) + \
			   "(rho, u, p)*_R : (%s, %s, %s)\n" % (self._qstarR.rho, self._qstarR.u, self._qstarR.p) + \
			   "(rho, u, p)_R  : (%s, %s, %s)"   % (self._qR.rho, self._qR.u, self._qR.p)

	def left_fastest(self):
		return self._waves[0]

	def right_fastest(self):
		return self._waves[4]

	def ustar(self):
		return self._ustar

	def pstar(self):
		return self._pstar

	def qstarL(self):
		return self._qstarL

	def qstarR(self):
		return self._qstarR

	def _delta_uR(self, p):
		"""right contribution to delta u in inflow problem"""
		if (p > self._qR.p):
			_delta = self._qR.delta_u_shock(p)
		else:
			_delta = self._qR.delta_u_expansion(p)
		return _delta

	def _pt0_from_pstar_rtt0(self, p, rtt0):
		ustar = self._qR.u + self._delta_uR(p)
		gam   = self._qR._gamma
		gsgmu = gam/(gam-1.)
		#print 'it',p, ustar, p * ( 1. + .5*ustar**2 / (gsgmu*rtt0 - .5*ustar**2) )**gsgmu
		return p * ( 1. + .5*ustar**2 / (gsgmu*rtt0 - .5*ustar**2) )**gsgmu

	def _pstar_estimate_expansion(self, pt0):
		"""exact solution for 2 expansion waves (for uniform gamma only)

			..note: from Toro, page 128
		"""
		return np.sqrt(self._qR.p * pt0)

	def solve_with_pt0_rtt0(self, pt0, rtt0):
		def dpt(p):
			return self._pt0_from_pstar_rtt0(p, rtt0) - pt0
		self._pstar  = fsolve(dpt, self._pstar_estimate_expansion(pt0), xtol=1e-10)
		self._ustar  = self._qR.u + self._delta_uR(self._pstar)
		# finalize qstarR state
		if (self._pstar > self._qR.p):
			_rho = self._qR.rho_through_shock(self._pstar)
			self._qstarR = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qR._gamma)
		else:
			_rho = self._qR.rho_through_isentropic(self._pstar)
			self._qstarR = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qR._gamma)
		#
		self._qstarL = self._qR.copy()
		self._qstarL.compute_from_pt_rtt_u(pt0, rtt0, self._ustar)
		return
