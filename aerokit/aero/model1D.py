"""@package model1D
  representation of 1D flows
"""

import numpy as np
from aerokit.common import defaultgas as defg
from aerokit.aero import Isentropic as Is
import aerokit.aero.ShockWave as sw
from typing import TypeVar

Afloat = TypeVar('Afloat', float, np.ndarray)

# -- class --


class __state(): # parent class not specific to 1D
    __init_attributes = ('rho', 'u', 'p', 'gamma')

    def __init__(self):
        self.rho = None
        self.p = None
        self._gamma = None

    def asound(self):
        """returns speed of sound"""
        return np.sqrt(self._gamma * self.p / self.rho) #type: ignore

    def Vmag(self):
        pass

    def KinE(self):
        pass

    def Mach(self):
        pass

    def Ptot(self):
        """returns Total pressure"""
        return self.p * Is.PtPs_Mach(self.Mach(), self._gamma)

    def rTtot(self):
        """returns Total temperature"""
        return self.p / self.rho + (self._gamma - 1.0) / self._gamma * self.KinE() #type: ignore


class State1d(__state):
    """
    defines a one dimensional state class

    Attributes:
            _gamma
            rho
            u
            p
    """

    def __init__(self, rho: Afloat, u: Afloat, p: Afloat, gamma=defg._gamma):
        self._gamma = gamma
        self.rho = rho
        self.u = u
        self.p = p

    def __repr__(self):
        return "state (rho, u, p) : (%s, %s, %s)" % (self.rho, self.u, self.p)

    @property
    def size(self):
        return max((q.size if isinstance(q, np.ndarray) else 1) for q in (self.rho, self.u, self.p))

    def copy(self):
        return state(self.rho, self.u, self.p, self._gamma)

    def KinE(self):
        return .5*self.u**2

    def Vmag(self):
        return np.abs(self.u)

    def Mach(self):
        """returns Mach number"""
        return self.u / self.asound()

    def massflow(self):
        """returns massflow"""
        return self.rho * self.u

    def state_RH(self):
        """return steady Rankine-Hugoniot jump state"""
        M = self.Mach()
        R = sw.Rho_ratio(M, self._gamma)
        return state(
            rho=self.rho * R, u=self.u / R, p=self.p * sw.Ps_ratio(M, self._gamma), gamma=self._gamma
        )

    def state_isentropic_Mach(self, Mach):
        """return state defined by Mach number through isoenergetic isentropic transformation"""
        ps = self.rTtot() / Is.PtPs_Mach(Mach, self._gamma)
        rts = self.rTtot() / Is.TtTs_Mach(Mach, self._gamma)
        return State1d(rho=ps / rts, u=Mach * np.sqrt(self._gamma * rts), p=ps)

    def compute_from_pt_rtt_M(self, pt, rtt, M):
        ps = pt / Is.PtPs_Mach(M, self._gamma)
        rts = rtt / Is.TtTs_Mach(M, self._gamma)
        self.__init__(rho=ps / rts, u=M * np.sqrt(self._gamma * rts), p=ps)

    def compute_from_pt_rtt_u(self, pt, rtt, u):
        gam = self._gamma
        gsgmu = gam / (gam - 1.0)
        a2 = gam * rtt - 0.5 * (gam - 1.0) * u ** 2
        ps = pt / (1.0 + 0.5 * u ** 2 / (gsgmu * rtt - 0.5 * u ** 2)) ** gsgmu
        self.__init__(rho=gam * ps / a2, u=u, p=ps)

    def compute_from_pt_rtt_p(self, pt, rtt, p):
        """Init state from Ptot r.Ttot and Ps (velocity sign is arbitrary and positive)

        Args:
                pt ([float]): [description]
                rtt ([float]): [description]
                p ([float]): [description]
        """
        M = Is.Mach_PtPs(pt / p, self._gamma)
        rts = rtt / Is.TtTs_Mach(M, self._gamma)
        self.__init__(rho=p / rts, u=M * np.sqrt(self._gamma * rts), p=p)

    def __getitem__(self, i):
        assert isinstance(self.rho, np.ndarray)
        assert isinstance(self.u, np.ndarray)
        assert isinstance(self.p, np.ndarray)
        return State1d(self.rho[i], self.u[i], self.p[i])

    def left_acoustic(self):
        """returns left 'actual' speed of sound"""
        return self.u - np.sqrt(self._gamma * self.p / self.rho)

    def right_acoustic(self):
        """returns left 'actual' speed of sound"""
        return self.u + np.sqrt(self._gamma * self.p / self.rho)


class state(State1d): # for backward compatibility
    pass