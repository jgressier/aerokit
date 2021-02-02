""" 
	The ``riemann`` module
	=========================
 
	Provides class and function for Riemann solvers for 1D unsteady flows
 
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
# implementation of RIEMANN_PB class

class riemann_pb():
	"""
		Class riemann_pb for unsteady 1D Euler

		attributes:
			_qL, _qR : 
			_pstar
			_ustar
			_qstarL
			_qstarR
			_waves[5]

	"""
	def __init__(self, qL, qR):
		self._qL    = qL.copy()
		self._qR    = qR.copy()
		self._waves = 5*[0]
		# define _pstar
		self._solve()
		self._ustar    = self._ustar_from_pstar(self._pstar)
		self._waves[2] = self._ustar
		#
		# finalize qstarL state
		if (self._pstar > self._qL.p):
			_rho = self._qL.rho_through_shock(self._pstar)
			_wsh = (self._qL.massflow()-_rho*self._ustar)/(self._qL.rho-_rho)
			self._waves[0] = _wsh
			self._waves[1] = _wsh
			self._qstarL = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qL._gamma)
		else:
			_rho = self._qL.rho_through_isentropic(self._pstar)
			self._qstarL = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qL._gamma)
			self._waves[0] = self._qL.u  - self._qL.asound()
			self._waves[1] = self._ustar - self._qstarL.asound()
		#
		# finalize qstarR state
		if (self._pstar > self._qR.p):
			_rho = self._qR.rho_through_shock(self._pstar)
			_wsh = (self._qR.massflow()-_rho*self._ustar)/(self._qR.rho-_rho)
			self._waves[3] = _wsh
			self._waves[4] = _wsh
			self._qstarR = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qR._gamma)
		else:
			_rho = self._qR.rho_through_isentropic(self._pstar)
			self._qstarR = uq.unsteady_state(rho=_rho, u=self._ustar, p=self._pstar, gamma=self._qR._gamma)
			self._waves[3] = self._qR.u  + self._qR.asound()
			self._waves[4] = self._ustar + self._qstarR.asound()

	def __repr__(self):
		return "(rho, u, p)_L : (%s, %s, %s)\n" % (self._qL.rho, self._qL.u, self._qL.p) + \
			   "(rho, u, p)_R : (%s, %s, %s)"   % (self._qR.rho, self._qR.u, self._qR.p)

	def qsol(self, xot):
		"""
			xot: x/t numpy array
		"""
		# left uniform part (initialization)
		rho = self._qL.rho * np.ones(len(xot))
		u   = self._qL.u   * np.ones(len(xot))
		p   = self._qL.p   * np.ones(len(xot))
		# left star uniform part
		i = np.logical_and((xot >= self._waves[1]), (xot < self._waves[2]))
		rho[i] = self._qstarL.rho
		u[i]   = self._qstarL.u
		p[i]   = self._qstarL.p
		# right star uniform part
		i = np.logical_and((xot >= self._waves[2]), (xot < self._waves[3]))
		rho[i] = self._qstarR.rho
		u[i]   = self._qstarR.u
		p[i]   = self._qstarR.p
		# right uniform part
		i = (xot >= self._waves[4])
		rho[i] = self._qR.rho
		u[i]   = self._qR.u
		p[i]   = self._qR.p
		# if left expansion
		if (self._pstar < self._qL.p):
			i = np.logical_and((xot > self._waves[0]), (xot < self._waves[1]))
			gam = self._qL._gamma
			gp1 = gam+1.
			gm1 = gam-1.
			u[i]   = gm1/gp1*(self._qL.u + 2./gm1*(self._qL.asound()+xot[i]))
			p[i]   = self._qL.p * ( (u[i] - xot[i])/self._qL.asound() )**(2.*gam/gm1)
			rho[i] = self._qL.rho_through_isentropic(p[i])
		# if right expansion
		if (self._pstar < self._qR.p):
			i = np.logical_and((xot > self._waves[3]), (xot < self._waves[4]))
			gam = self._qR._gamma
			gp1 = gam+1.
			gm1 = gam-1.
			u[i]   = gm1/gp1*(self._qR.u + 2./gm1*(-self._qR.asound()+xot[i]))
			p[i]   = self._qR.p * ( (xot[i] - u[i] )/self._qR.asound() )**(2.*gam/gm1)
			rho[i] = self._qR.rho_through_isentropic(p[i])
		return uq.unsteady_state(rho, u, p, np.where(xot<0., self._qL._gamma, self._qR._gamma))

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

	def _delta_uL(self, p):
		"""left contribution to delta u in Riemann problem"""
		if (p > self._qL.p):
			_delta = self._qL.delta_u_shock(p)
		else:
			_delta = self._qL.delta_u_expansion(p)
		return _delta

	def _delta_uR(self, p):
		"""right contribution to delta u in Riemann problem"""
		if (p > self._qR.p):
			_delta = self._qR.delta_u_shock(p)
		else:
			_delta = self._qR.delta_u_expansion(p)
		return _delta

	def _ustar_from_pstar(self, p):
		return .5*(self._qR.u + self._qL.u + self._delta_uR(p) - self._delta_uL(p)) 

	def _pstar_estimate_expansion(self):
		"""exact solution for 2 expansion waves (for uniform gamma only)

			..note: from Toro, page 128
		"""
		gam   = .5*(self._qL._gamma+self._qR._gamma)
		gmusd = .5*(gam-1.)
		return ((self._qL.asound()+self._qR.asound()
			    -gmusd*(self._qR.u - self._qL.u)) / 
				( self._qL.asound()/self._qL.p**(gmusd/gam)
				+ self._qR.asound()/self._qR.p**(gmusd/gam) ))**(gam/gmusd)

	def _solve(self):
		def du(p):
			return self._delta_uL(p) + self._delta_uR(p) + self._qR.u - self._qL.u
		self._pstar = fsolve(du, self._pstar_estimate_expansion(), xtol=1e-10)
		return
