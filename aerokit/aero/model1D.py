"""@package model1D
  representation of 1D flows
"""

import math
import numpy     as np
#import ShockWave as sw
#import IterativeSolve
from ..common   import defaultgas as defg
from aerokit.aero import Isentropic as Is

# -- class --

class state():
	"""
	defines a one dimensional state class

	Attributes:
		_gamma
		rho
		u
		p
	"""
	def __init__(self, rho, u, p, gamma=defg._gamma):
	 	self._gamma = gamma
	 	self.rho    = rho
	 	self.u      = u
	 	self.p      = p

	def __repr__(self):
		return "state (rho, u, p) : (%s, %s, %s)" % (self.rho, self.u, self.p)

	def copy(self):
		return state(self.rho, self.u, self.p, self._gamma)

	def compute_from_pt_rtt_M(self, pt, rtt, M):
		ps  = pt /Is.PtPs_Mach(M, self._gamma)
		rts = rtt/Is.TtTs_Mach(M, self._gamma)
		self.__init__(rho=ps/rts, u=M*np.sqrt(self._gamma*rts), p=ps)
		
	def compute_from_pt_rtt_u(self, pt, rtt, u):
		gam   = self._gamma
		gsgmu = gam/(gam-1.)
		a2    = gam * rtt - .5*(gam-1.)*u**2
		ps    = pt / ( 1. + .5*u**2 / (gsgmu*rtt - .5*u**2) )**gsgmu
		self.__init__(rho=gam * ps / a2, u=u, p=ps)

	def compute_from_pt_rtt_p(self, pt, rtt, p):
		M   = Is.Mach_PtPs(pt/p, self._gamma)
		rts = rtt/Is.TtTs_Mach(M, self._gamma)
		self.__init__(rho=p/rts, u=M*np.sqrt(self._gamma*rts), p=p)
		
	def asound(self):
		"""returns speed of sound"""
		return math.sqrt(self._gamma*self.p/self.rho)

	def left_acoustic(self):
		"""returns left 'actual' speed of sound"""
		return self.u-math.sqrt(self._gamma*self.p/self.rho)

	def right_acoustic(self):
		"""returns left 'actual' speed of sound"""
		return self.u+math.sqrt(self._gamma*self.p/self.rho)

	def Mach(self):
		"""returns Mach number"""
		return self.u/self.asound()

	def massflow(self):
		"""returns massflow"""
		return self.rho * self.u

	def Ptot(self):
		"""returns Total pressure"""
		return self.p*Is.PtPs_Mach(self.Mach(), self._gamma)

	def rTtot(self):
		"""returns Total temperature"""
		return self.p/self.rho + .5*(self._gamma-1.)/self._gamma*self.u**2

	# def rTtot(self):
	# 	"""returns speed of sound"""
	# 	return self.u/self.asound()

