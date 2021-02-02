"""
    The ``gasgenerator`` module
    =========================
 
    Provides standart Stagnation to Static Ratios for Ideal gas
 
    :Example:
 
    >>> import aerokit.aero.Isentropic as Is
    >>> Is.TtTs_Mach(1.)
    1.2
    >>> Is.TtTs_Mach(2., gamma=1.6)
    2.2
 
    Available functions
    -------------------
 
    Provides Ti/Ts Pi/Ps ratios from Mach number and reversed functions.
    Specific heat ratio `gamma` is optionnal and can be specified in the functions itself
    or using aerokit.common.defaultgas module
 """

import math
import numpy as np
from aerokit.common import defaultgas as defg # relative import is deprecated by doctest
import aerokit.aero.Isentropic as Is

# ===============================================================
# implemented functions

class base():
    """
    	virtual class for gas generator and other engine module

		Long comment
 
		:param Mach:  local Mach number
        :param gamma: specific heat ratio, default from aerokit.common.defaultgas
		:return:      result Ti/Ts ratio
 
 		:Example:

		>>> TtTs_Mach(1.) # with default gamma 1.4
		1.2

		.. seealso:: 
		.. note:: available for scalar or array (numpy) computations
    """	
    def __init__(self, OPR, Tt4, gam_cold=defg._gamma, gam_hot=defg._gamma, r_cold=defg._r, r_hot=defg._r,
                                 xi_inlet=1., xi_cc=1., eta_cc=1., Pci=42.8e6,
                                 etapolCHP=1., etapolTHP=1., eta_shaft=1.):
        self.OPR = OPR
        self.Tt4 = Tt4
        self.Pci = Pci
        self.gam_cold  = gam_cold
        self.gam_hot   = gam_hot
        self.r_cold    = r_cold
        self.r_hot     = r_hot
        self.xi_cc     = xi_cc
        self.xi_inlet  = xi_inlet
        self.etapolCHP = etapolCHP
        self.etapolTHP = etapolTHP
        self.eta_shaft = eta_shaft

    def flight_conditions(self, T0=1., P0=1., M0=0.):
        self.T0 = T0
        self.P0 = P0
        self.M0 = M0

    def update(self):
        gc  = self.gam_cold
        cpc = gc*self.r_cold/(gc-1.)
        self.Tt0 = self.T0*Is.TtTs_Mach(self.M0, gamma=self.gam_cold)
        self.Pt0 = self.P0*Is.PtPs_Mach(self.M0, gamma=self.gam_cold)
        self.V0  = self.M0*np.sqrt(gc*self.r_cold*self.T0)
        self.Pt2 = self.Pt0*self.xi_inlet
        self.Tt2 = self.Tt0
        self.Pt3 = self.Pt2*self.OPR
        self.Tt3 = self.Tt2*self.OPR**((gc-1.)/(gc*self.etapolCHP))
        #self.Tt4 = Ti_4
        self.Pt4 = self.Pt3*self.xi_cc
        gh  = self.gam_hot
        cph = gh*self.r_hot/(gh-1.)
        self.far = (cph*self.Tt4 - cpc*self.Tt3)/(self.xi_cc*self.Pci - cph*self.Tt4)
        self.Tt45 = self.Tt4 - cpc*(self.Tt3-self.Tt2)/(self.eta_shaft*cph*(1.+self.far))
        #print self.Tt3, self.Pt4, self.Tt45, self.Tt4, self.far
        self.Pt45 = self.Pt4*(self.Tt45/self.Tt4)**(gh/((gh-1.)*self.etapolTHP))

    def Wsp_cc(self):
        return self.far * self.Pci

    def Wsp_TBP_ideal(self, P19=None):
        if P19==None: P19=self.P0
        gogm1 = self.gam_hot/(self.gam_hot-1.)
        return (1.+self.far)*gogm1*self.r_hot/(self.gam_hot-1.)*self.Tt45*((P19/self.Pt45)**(1./gogm1)-1.)

    def thermal_efficiency(self):
        return self.Wsp_TBP_ideal()/self.Wsp_cc()

# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()