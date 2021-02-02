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
import aerokit.engine.turbojet as tj

# ===============================================================
# implemented functions

class turbofan_adapt(tj.turbojet_opt):
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
    def __init__(self, OPR, Tt4, bpr, fanpower_ratio, gam_cold=defg._gamma, gam_hot=defg._gamma, r_cold=defg._r, r_hot=defg._r,
                                 xi_inlet=1., xi_nozzle=1., xi_cc=1., eta_cc=1., Pci=42.8e6,
                                 etapolCHP=1., etapolTHP=1., etapolTBP=1., etapolfan=1., eta_shaft=1.):
        tj.turbojet_opt.__init__(self, OPR, Tt4, gam_cold=gam_cold, gam_hot=gam_hot, r_cold=r_cold, r_hot=r_hot,
                                 xi_inlet=xi_inlet, xi_cc=xi_cc, eta_cc=eta_cc, Pci=Pci,
                                 etapolCHP=etapolCHP, etapolTHP=etapolTHP, eta_shaft=eta_shaft)
        self.bpr            = bpr
        self.etapolfan      = etapolfan
        self.etapolTBP      = etapolTBP
        self.fanpower_ratio = fanpower_ratio

    def update(self):
        tj.turbojet_opt.update(self)
        gh  = self.gam_hot
        cph = gh*self.r_hot/(gh-1.)
        Wsp_mono = (1.-(self.Pt45/self.P0*self.xi_nozzle)**(-(gh-1.)*self.etapolTBP/gh))*self.Tt45*cph*(1.+self.far)
        self.Tt5  = self.Tt45 - Wsp_mono*self.fanpower_ratio/cph/(1.+self.far)
        self.Pt5  = self.Pt45*(self.Tt5/self.Tt45)**(gh/((gh-1.)*self.etapolTBP))
        # core nozzle
        self.Pt9 = self.Pt5 * self.xi_nozzle
        self.M9  = Is.Mach_PtPs(self.Pt9/self.P0, gamma=self.gam_hot)
        self.V9  = Is.Velocity_MachTi(self.M9, self.Tt5, r=self.r_hot, gamma=self.gam_hot)
        # fan
        gc  = self.gam_cold
        cpc = gc*self.r_cold/(gc-1.)
        #print self.bpr, cpc
        self.Tt17 = self.Tt2 + self.eta_shaft*Wsp_mono*self.fanpower_ratio/(self.bpr*cpc)
        self.Pt17 = self.Pt2*(self.Tt17/self.Tt2)**((gc*self.etapolfan)/(gc-1.))
        # bypass nozzle    
        self.Tt19 = self.Tt17
        self.Pt19 = self.Pt17*self.xi_nozzle
        self.M19  = Is.Mach_PtPs(self.Pt19/self.P0, gamma=gc)
        self.V19  = Is.Velocity_MachTi(self.M19, self.Tt19, r=self.r_cold, gamma=gc)

    def Wsp_kinEn(self):       
        return .5*((1.+self.far)*self.V9**2-self.V0**2 + self.bpr*(self.V19**2-self.V0**2))

    def norm_thrust(self):
        return ((1.+self.far)*self.V9-self.V0 + self.bpr*(self.V19-self.V0))

    def spec_thrust(self):
        return self.norm_thrust()/(1.+self.bpr)


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()