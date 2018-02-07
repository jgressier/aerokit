"""@package riemann
  Riemann solvers for 1D unsteady flows
"""

import math
import numpy     as np
import ShockWave as sw
import IterativeSolve

# -- class --

class riemann_pb():
	"""
	defines a Riemann problem
	"""
	def __init__(self, qL, qR):
	 	self._qL    = qL.copy()
	 	self._qR    = qR.copy()

	def __repr__(self):
		return "(rho, u, p)_L : (%s, %s, %s)\n" % (self._qL.rho, self._qL.u, self._qL.p) +
		       "(rho, u, p)_R : (%s, %s, %s)"   % (self._qR.rho, self._qR.u, self._qR.p)


