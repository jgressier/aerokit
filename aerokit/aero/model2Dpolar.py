"""@package model2Dpolar
  representation of 1D flows
"""

import numpy as np
from aerokit.common import defaultgas as defg
from aerokit.aero import Isentropic as Is
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw

# -- class --
class state2Dpolar():
    """
    defines a 2 dimensional state class

    Attributes:
            _gamma
            rho
            Mach, angle (deg)
            p
    """

    def __init__(self, Mach, angle=0., rho=1., p=1., gamma=defg._gamma):
        self._gamma = gamma
        self.rho = rho
        self.Mach = Mach
        self.angle = angle
        self.p = p

    def __repr__(self):
        return "state (rho, Mach, angle, p) : (%s, %s, %s, %s)" % (self.rho, self.Mach, self.angle, self.p)

    @property
    def size(self):
        return self.rho.size if isinstance(self.rho, np.ndarray) else 1

    def copy(self):
        return state2Dpolar(self.Mach, self.angle, self.rho, self.p, self._gamma)

    def rotate(self, deviation):
        self.angle += deviation

    def weakshock_deviation(self, deviation):
        sigma = sw.weaksigma_Mach_deflection(self.Mach, deviation, self._gamma)
        Mn0 = self.Mach * deg.sin(sigma)
        M1 = abs(sw.downstream_Mn(Mn0, self._gamma)/deg.sin(sigma-deviation))
        return state2Dpolar(M1, self.angle+deviation,
                            rho=self.rho*sw.Rho_ratio(Mn0, self._gamma),
                            p=self.p*sw.Ps_ratio(Mn0, self._gamma),)

    def shock_sigma(self, sigma):
        deviation = sw.deflection_Mach_sigma(self.Mach, sigma)
        Mn0 = self.Mach * deg.sin(sigma)
        M1 = abs(sw.downstream_Mn(Mn0, self._gamma)/deg.sin(sigma-deviation))
        return state2Dpolar(M1, self.angle+deviation,
                            rho=self.rho*sw.Rho_ratio(Mn0, self._gamma),
                            p=self.p*sw.Ps_ratio(Mn0, self._gamma),)
