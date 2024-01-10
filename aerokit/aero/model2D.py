"""@package model2D
  representation of 2D flows
"""

import numpy as np
from aerokit.aero import Isentropic, Supersonic
from aerokit.common import defaultgas as defg
from aerokit.aero.model1D import __state
import aerokit.aero.ShockWave as sw
import aerokit.aero.degree as deg
from typing import TypeVar

Afloat = TypeVar('Afloat', float, np.ndarray)

# -- class --

class State2d(__state):
    """
    defines a 2 dimensional state class

    Attributes:
            _gamma
            rho
            u, v
            p
    """

    def __init__(self, rho: Afloat, u: Afloat, v: Afloat, p: Afloat, gamma=defg._gamma):
        self._gamma = gamma
        self.rho = rho
        self.u = u
        self.v = v
        self.p = p

    def __repr__(self):
        return "state (rho, u, v, p) : (%s, %s, %s)" % (self.rho, self.u, self.v, self.p)

    @property
    def size(self):
        return max((q.size if isinstance(q, np.ndarray) else 1) for q in (self.rho, self.u, self.v, self.p))

    def copy(self):
        return State2d(self.rho, self.u, self.v, self.p, self._gamma)

    def KinE(self):
        return .5*(self.u**2+self.v**2)

    def Vmag(self):
        return np.sqrt(self.u**2+self.v**2)

    def Mach(self):
        """returns Mach number"""
        return self.Vmag() / self.asound()

    def angle(self):
        return deg.atan2(self.v, self.u)

    def omega(self):
        """returns Prandtl-Meyer (or Busemann) angle in degree"""
        return Supersonic.PrandtlMeyer_Mach(self.Mach(), self._gamma)
    
    # def state_RH(self):
    #     """return Rankine-Hugoniot jump state"""
    #     M = self.Mach()
    #     R = sw.Rho_ratio(M, self._gamma)
    #     return state(
    #         rho=self.rho * R, u=self.u / R, p=self.p * sw.Ps_ratio(M, self._gamma), gamma=self._gamma
    #     )

    # def state_isentropic_Mach(self, Mach):
    #     """return state defined by Mach number through isoenergetic isentropic transformation"""
    #     ps = self.rTtot() / Is.PtPs_Mach(Mach, self._gamma)
    #     rts = self.rTtot() / Is.TtTs_Mach(Mach, self._gamma)
    #     return state(rho=ps / rts, u=Mach * np.sqrt(self._gamma * rts), p=ps)

    # def compute_from_pt_rtt_M(self, pt, rtt, M):
    #     ps = pt / Is.PtPs_Mach(M, self._gamma)
    #     rts = rtt / Is.TtTs_Mach(M, self._gamma)
    #     self.__init__(rho=ps / rts, u=M * np.sqrt(self._gamma * rts), p=ps)

    # def compute_from_pt_rtt_u(self, pt, rtt, u):
    #     gam = self._gamma
    #     gsgmu = gam / (gam - 1.0)
    #     a2 = gam * rtt - 0.5 * (gam - 1.0) * u ** 2
    #     ps = pt / (1.0 + 0.5 * u ** 2 / (gsgmu * rtt - 0.5 * u ** 2)) ** gsgmu
    #     self.__init__(rho=gam * ps / a2, u=u, p=ps)

    # def compute_from_pt_rtt_p(self, pt, rtt, p):
    #     """Init state from Ptot r.Ttot and Ps (velocity sign is arbitrary and positive)

    #     Args:
    #             pt ([float]): [description]
    #             rtt ([float]): [description]
    #             p ([float]): [description]
    #     """
    #     M = Is.Mach_PtPs(pt / p, self._gamma)
    #     rts = rtt / Is.TtTs_Mach(M, self._gamma)
    #     self.__init__(rho=p / rts, u=M * np.sqrt(self._gamma * rts), p=p)

    def __getitem__(self, i):
        assert isinstance(self.rho, np.ndarray)
        return State2d(self.rho[i], self.u[i], self.v[i], self.p[i])
