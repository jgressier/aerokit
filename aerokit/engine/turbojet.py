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
import aerokit.aero.Isentropic     as Is
import aerokit.engine.gasgenerator as gg

# ===============================================================
# implemented functions

class turbojet_opt(gg.base):
    """
    	Computes Ti/Ts ratio from Mach number

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
                                 xi_inlet=1., xi_nozzle=1., xi_cc=1., eta_cc=1., Pci=42.8e6,
                                 etapolCHP=1., etapolTHP=1., eta_shaft=1.):
        gg.base.__init__(self, OPR, Tt4, gam_cold=gam_cold, gam_hot=gam_hot, r_cold=r_cold, r_hot=r_hot,
                                 xi_inlet=xi_inlet, xi_cc=xi_cc, eta_cc=eta_cc, Pci=Pci,
                                 etapolCHP=etapolCHP, etapolTHP=etapolTHP, eta_shaft=eta_shaft)
        self.xi_nozzle = xi_nozzle


    def update(self):
        gg.base.update(self)
        gh  = self.gam_hot
        cph = gh*self.r_hot/(gh-1.)
        self.Pt9 = np.maximum(self.Pt45 * self.xi_nozzle, self.P0)
        self.M9  = Is.Mach_PtPs(self.Pt9/self.P0, gamma=self.gam_hot)
        self.V9  = Is.Velocity_MachTi(self.M9, self.Tt45, r=self.r_hot, gamma=self.gam_hot)

    def Wsp_kinEn(self):       
        return .5*((1.+self.far)*self.V9**2-self.V0**2)

    def norm_thrust(self):
        return (1.+self.far)*self.V9-self.V0

    def spec_thrust(self):
        return self.norm_thrust()

    def thermal_efficiency(self):
        return self.Wsp_kinEn()/self.Wsp_cc()

    def propulsive_efficiency(self):
        return self.norm_thrust()*self.V0/self.Wsp_kinEn()

    def thermoprop_efficiency(self):
        return self.norm_thrust()*self.V0/self.Wsp_cc()

    def Sfc(self):
        return self.far/self.norm_thrust()

# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()