"""@package SW interaction
    2D planar Shock waves interaction using 
    local Rankine Hugoniot equations 
"""

import numpy as np
from aerokit.common import defaultgas as defg  # relative import is deprecated by doctest
import aerokit.aero.Isentropic as Is
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw
import aerokit.aero.model2Dpolar as M2D
import aerokit.aero.plot.shockpolar as plotsw
from scipy import optimize


class ShockInteraction:
    """Class to define and compute shock interaction
    upstream flow is 0 ; 1 and 2 are respective bottom and top downstream states

    Args:
        M0 (float): upstream Mach number
        sigma0 (float): C+ like shockwave (>0)
        sigma1 (float): C- like shockwave (<0)
        gamma (float, optional): ratio of specific quantities
    """

    def __init__(self, M0, sigma01, sigma02, gamma=defg._gamma):
        if M0 <= 1.0:
            raise ValueError("upstream Mach number M0 must be supersonic")
        self._state = dict()
        self._state[0] = M2D.state2Dpolar(M0, 0.0)  # default angle=0., rho and p normalized
        self._gamma = gamma
        self.sigma01 = sigma01  # use setter
        self.sigma02 = sigma02  # use setter

    @property
    def M0(self):
        return self._state[0].Mach

    @property
    def sigma01(self):
        return self._sigma01

    @property
    def sigma02(self):
        return self._sigma02

    @sigma01.setter
    def sigma01(self, value):
        self._solved = False
        if value < deg.asin(1 / self.M0):
            raise ValueError("sigma01 angle must be more than µ0")
        self._sigma01 = value
        self._state[1] = self[0].shock_sigma(value)

    @sigma02.setter
    def sigma02(self, value):
        self._solved = False
        if value > -deg.asin(1 / self.M0):
            raise ValueError("sigma02 angle must be negative and (abs) more than µ0")
        self._sigma02 = value
        self._state[2] = self[0].shock_sigma(value)

    def __getitem__(self, value) -> M2D.state2Dpolar:
        return self._state[value]

    def check34balanced(self, tol=1.0e-6):
        return (abs(self[3].p - self[4].p) / self[0].p < tol) and (abs(self[3].angle - self[4].angle) < tol)

    def solve(self, verbose=False):
        """solve the interaction of shocks using intersection points
        in the p/angle diagram

        Returns:
            _type_: _description_
        """
        print(f"SWI solve init:\n  Q1:{self[1]}\n  Q2:{self[2]}")
        Q1i = self[1] if self[1].Mach > 1 else self[0]
        Q2i = self[2] if self[2].Mach > 1 else self[0]

        def delta_p(theta):
            if self[1].Mach > 1:
                Q3 = Q1i.weakshock_deviation(theta - Q1i.angle)
            else:
                Q3 = Q1i.strongshock_deviation(theta - self[0].angle)
            if self[2].Mach > 1:
                Q4 = Q2i.weakshock_deviation(theta - Q2i.angle)
            else:
                Q4 = Q2i.strongshock_deviation(theta - self[0].angle)
            if verbose: print(f"SWI iterations:\n  {Q3}\n  {Q4}")
            return Q4.p - Q3.p

        th_init = -1.
        sol = optimize.root_scalar(delta_p, x0=th_init, method='newton')
        theta = sol.root
        self.solved = sol.converged
        if self[1].Mach > 1:
            self._state[3] = self[1].weakshock_deviation(theta - self[1].angle)
        else:
            self._state[3] = self[0].strongshock_deviation(theta - self[0].angle)
            self._state[1] = self[3].copy()
            self._sigma01 = sw.deflection_Mach_sigma(self.M0, theta-self[0].angle)
        if self[2].Mach > 1:
            print('weak test')
            self._state[4] = self[2].weakshock_deviation(theta - self[2].angle)
        else:
            print('strong test')
            self._state[4] = self[0].strongshock_deviation(theta - self[0].angle)
            self._state[2] = self[4].copy()
            self._sigma02 = sw.strongsigma_Mach_deflection(self.M0, theta-self[0].angle)
        print(f"SWI solve end:\n  Q1:{self[1]}\n  Q2:{self[2]}\n  Q3:{self[3]}\n  Q4:{self[4]}")
        return sol

    def plot_angle_pressure(self, ax=plotsw.plt):
        # plot polar curves
        plotsw.plot_theta_pressure(self.M0, ax=ax)
        if self[1].Mach > 1:
            plotsw.plot_theta_pressure(self[1].Mach, thet_init=self[1].angle, p_init=self[1].p, color='red', ax=ax)
        if self[2].Mach > 1:
            plotsw.plot_theta_pressure(self[2].Mach, thet_init=self[2].angle, p_init=self[2].p, color='red', ax=ax)

        # # plot symbols for flow regions
        # plt.plot(0,    1.,   'bo')
        # plt.plot(wdev, p1p0, 'wo')
        # plt.plot(0.,   p2p0, 'go')