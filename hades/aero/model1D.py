"""@package model1D
  representation of 1D flows
"""

import math
import numpy     as np
#import ShockWave as sw
#import IterativeSolve
from ..common import defaultgas as defg

# -- class --

class state():
	"""
	defines a one dimensional state
	"""
	def __init__(self, rho, u, p, gamma=defg._gamma):
	 	self._gamma = gamma
	 	self.rho    = rho
	 	self.u      = u
	 	self.p      = p

	def __repr__(self):
		return "state (rho, u, p) : (%s, %s, %s)" % (self.rho, self.u, self.p)

	def asound(self):
		"""returns speed of sound"""
		return math.sqrt(self._gamma*self.p/self.rho)

	def mach(self):
		"""returns speed of sound"""
		return self.u/self.asound()

	def _rankinehugoniot(self, ushock):
		"""returns Rankine-Hugoniot state, given a shock velocity
		   this function is made private because there is no test about upstream/downstream consistency
		   nor the right ushock range (greater than u+a or lesser than u-a
		"""
		loc_Mn    = (self.u-ushock)/self.asound()          # this Mach number is signed
		rho_ratio = sw.Rho_ratio(loc_Mn, self._gamma)
		return state(self.rho * rho_ratio,
	 	             ushock + (self.u-ushock)/rho_ratio,
	 	             self.p * sw.Ps_ratio(loc_Mn, self._gamma))

	def shock_downstream(self, ushock):
		"""computes downstream state of genuine unsteady shock wave"""

		### TO BE IMPLEMENTED# -- functions --

