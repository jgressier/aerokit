"""@package SW interaction
    2D planar Shock waves interaction using 
    local Rankine Hugoniot equations 
"""

import numpy as np
import aerokit.aero.Isentropic as Is
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw
from aerokit.common import defaultgas as defg  # relative import is deprecated by doctest
from scipy.optimize import newton


class ShockInteraction:
    """Class to define and compute shock interaction
    upstream flow is 0 ;

    Args:
        M0 (float): upstream Mach number
        sigma0 (float): C+ like shockwave (>0)
        sigma1 (float): C- like shockwave (<0)
        gamma (float, optional): ratio of specific quantities
    """

    def __init__(self, M0, sigma01, sigma02, gamma=defg._gamma):
        if M0 <= 1.0:
            raise ValueError("upstream Mach number M0 must be supersonic")
        self._M0 = M0
        self._gamma = gamma
        self.sigma01 = sigma01  # use setter
        self.sigma02 = sigma12  # use setter

    @property
    def sigma01(self):
        return self._sigma01

    @property
    def sigma02(self):
        return self._sigma02

    @sigma01.setter
    def sigma01(self, value):
        self._solved = False
        if value < deg.asin(1 / self._M0):
            print(value, deg.asin(1 / self._M0))
            raise ValueError("sigma01 angle must be more than µ0")
        self._sigma01 = value
        self._theta1 = sw.deflection_Mach_sigma(self._M0, self._sigma01, self._gamma)

    @sigma02.setter
    def sigma02(self, value):
        self._solved = False
        if value > -deg.asin(1 / self._M0):
            raise ValueError("sigma02 angle must be negative and (abs) more than µ0")
        self._sigma02 = value
        self._theta2 = sw.deflection_Mach_sigma(self._M0, self._sigma02, self._gamma)

    def solve(self):
        Mn01 = self._M0*deg.sin(self._sigma01)
        Mn02 = self._M0*deg.sin(self._sigma02)
        p1 = sw.Ps_ratio(Mn01, self._gamma)
        p2 = sw.Ps_ratio(Mn02, self._gamma)
        def delta_p(theta):
            
        th_init = (self._theta1 + self._theta2) / 2

        return th
