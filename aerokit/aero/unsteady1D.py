""" 
	The ``unsteady1D`` module
	=========================
 
	Specific functions to compute propagation and evolution of 1D state
 
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
import copy
import numpy     as np
import aerokit.aero.ShockWave as sw
import aerokit.aero.model1D   as model1D

# -- class --

class unsteady_state(model1D.state):
	"""
	defines a one dimensional state
	"""
	def copy(self):
		return copy.deepcopy(self)

	def copysymmetric(self):
		return unsteady_state(self.rho, -self.u, self.p, self._gamma)

	def _rankinehugoniot_from_ushock(self, ushock):
		"""
			returns Rankine-Hugoniot state, given a shock velocity
		
			..warning:: this function is made private because there is no test about upstream/downstream consistency
			nor the right ushock range (greater than u+a or lesser than u-a
		"""
		loc_Mn    = (self.u-ushock)/self.asound()          # this Mach number is signed
		rho_ratio = sw.Rho_ratio(loc_Mn, self._gamma)
		return unsteady_state(self.rho * rho_ratio,
		             ushock + (self.u-ushock)/rho_ratio,
		             self.p * sw.Ps_ratio(loc_Mn, self._gamma),
	 	             gamma=self._gamma)

	def delta_u_expansion(self, p):
		"""
			computes delta u through expansion given 'self' initial state and final pressure

			usage
			-----

			if C- expansion, from left state,  middle velocity ustar is (uL - delta_u_expansion)
			is C+ expansion, from right state, middle velocity ustar is (uR + delta_u_expansion)

			..warning:: There is no consistent test that (arg) p is lesser than self.p
			..note:: from Toro, page 119
		"""
		gmusd = .5*(self._gamma-1.)
		return self.asound()/gmusd*((p/self.p)**(gmusd/self._gamma)-1.)

	def delta_u_shock(self, p):
		"""
			computes delta u through shock given 'self' initial state and final pressure

			usage
			-----

			if C- shock, from left state,  middle velocity ustar is (uL - delta_u_shock)
			is C+ shock, from right state, middle velocity ustar is (uR + delta_u_shock)

			..warning:: There is no consistent test that (arg) p is greater than self.p
		"""
		return (p-self.p)*math.sqrt(2./self.rho/((self._gamma+1.)*p + (self._gamma-1.)*self.p))

	def rho_through_shock(self, p):
		"""
			computes rho increase through shock given 'self' initial state and final pressure

			..warning:: There is no consistent test that (arg) p is greater than self.p
		"""
		gmusgpu = (self._gamma-1.)/(self._gamma+1.)
		pratio  = p/self.p
		return self.rho * (pratio+gmusgpu)/(1.+gmusgpu*pratio)

	def rho_through_isentropic(self, p):
		"""
			computes rho increase through isentropic compression or expansion given 'self' initial state and final pressure
		"""
		return self.rho * (p/self.p) ** (1./self._gamma)

	def shock_downstream(self, ushock):
		"""computes downstream state of genuine unsteady shock wave"""

		### TO BE IMPLEMENTED# -- functions --

