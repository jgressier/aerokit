"""@package SW interaction
    2D planar Shock waves interaction using 
    local Rankine Hugoniot equations 
"""

from dataclasses import dataclass

import numpy as np
from aerokit.common import defaultgas as defg  # relative import is deprecated by doctest
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw
import aerokit.aero.model2D as M2D
import aerokit.aero.plot.shockpolar as plotsw
from scipy import optimize


@dataclass(frozen=True)
class TriplePointSolution:
    """Balanced states and shock angles at a three-shock intersection."""

    theta: float
    reflected_state: M2D.State2DMach
    stem_state: M2D.State2DMach
    sigma_reflected: float
    sigma_stem: float
    root_result: object


class ShockInteraction:
    """Class to define and compute shock interaction
    upstream flow is 0 ; 1 and 2 are respective bottom and top downstream states

    Args:
        M0 (float): upstream Mach number
        sigma0 (float): C+ like shockwave (>0)
        sigma1 (float): C- like shockwave (<0)
        gamma (float, optional): ratio of specific quantities
    """

    def __init__(self, M0, sigma01, sigma02):
        if M0 <= 1.0:
            raise ValueError("upstream Mach number M0 must be supersonic")
        self._state = dict()
        self._state[0] = M2D.State2DMach(M0, 0.0)  # default angle=0., rho and p normalized
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

    def __getitem__(self, value) -> M2D.State2DMach:
        return self._state[value]

    def check34balanced(self, tol=1.0e-6):
        return (abs(self[3].p - self[4].p) / self[0].p < tol) and (abs(self[3].angle - self[4].angle) < tol)

    @staticmethod
    def solve_triple_point(upstream_state, incident_state, epsilon=1.0e-5):
        """Intersect a reflected weak polar with an upstream strong polar.

        The root is bracketed on the physically admissible reflected-shock
        branch, between the incident-flow angle and its maximum attached
        turning.  The returned reflected and stem states have equal pressure
        and flow angle but may have different density and Mach number across
        the slip line.

        Args:
            upstream_state: State ahead of the incident shock and Mach stem.
            incident_state: State immediately behind the incident shock.
            epsilon: Angular distance in degrees kept from polar endpoints.

        Returns:
            A :class:`TriplePointSolution`.

        Raises:
            ValueError: If the states have no bracketed triple-point solution.
        """
        incident_deviation = incident_state.angle - upstream_state.angle
        if incident_deviation == 0.0:
            raise ValueError("incident and upstream flow angles must differ")

        turning_limit = incident_state.devmax()
        if incident_deviation > 0.0:
            theta_left = incident_state.angle - turning_limit + epsilon
            theta_right = incident_state.angle - epsilon
        else:
            theta_left = incident_state.angle + epsilon
            theta_right = incident_state.angle + turning_limit - epsilon

        def pressure_mismatch(theta):
            reflected = incident_state.weakshock_deviation(theta - incident_state.angle)
            stem = upstream_state.strongshock_deviation(theta - upstream_state.angle)
            return reflected.p - stem.p

        mismatch_left = pressure_mismatch(theta_left)
        mismatch_right = pressure_mismatch(theta_right)
        if mismatch_left * mismatch_right > 0.0:
            raise ValueError("reflected and strong polars have no bracketed intersection")

        root_result = optimize.root_scalar(
            pressure_mismatch,
            bracket=(theta_left, theta_right),
            method="brentq",
        )
        theta = root_result.root
        reflected_state = incident_state.weakshock_deviation(theta - incident_state.angle)
        stem_state = upstream_state.strongshock_deviation(theta - upstream_state.angle)
        sigma_reflected = sw.weaksigma_Mach_deflection(
            incident_state.Mach, theta - incident_state.angle
        )
        sigma_stem = sw.strongsigma_Mach_deflection(
            upstream_state.Mach, theta - upstream_state.angle
        )
        return TriplePointSolution(
            theta=theta,
            reflected_state=reflected_state,
            stem_state=stem_state,
            sigma_reflected=sigma_reflected,
            sigma_stem=sigma_stem,
            root_result=root_result,
        )

    def solve(self, verbose=False):
        """solve the interaction of shocks using intersection points
        in the p/angle diagram

        Returns:
            _type_: _description_
        """
        zvar = True # change variable to optimize newton
        #print(f"SWI solve init:\n  Q1:{self[1]}\n  Q2:{self[2]}")
        Q1i = self[1] if self[1].Mach > 1 else self[0]
        Q2i = self[2] if self[2].Mach > 1 else self[0]
        deltarange = 1.1*max([abs(Q.angle)+Q.devmax() for Q in (Q1i, Q2i) ])

        def ang2z(a):
            return deg.tan(90./deltarange*a) if zvar else a
        
        def z2ang(z):
            return deg.atan(z)*deltarange/90. if zvar else z

        def delta_p(ztheta):
            theta = z2ang(ztheta)
            if self[1].Mach > 1:
                dev = theta - Q1i.angle
                Q3 = Q1i.weakshock_deviation(dev)
                if abs(dev) > Q1i.devmax(): # penalty
                    Q3.p += 10*deg.tan(-dev - Q1i.devmax())*self[0].p
            else:
                Q3 = Q1i.strongshock_deviation(theta - self[0].angle)
            if self[2].Mach > 1:
                dev = theta - Q2i.angle
                Q4 = Q2i.weakshock_deviation(dev)
                if abs(dev) > Q2i.devmax(): # penalty
                    Q4.p += 10*deg.tan(dev - Q2i.devmax())*self[0].p
            else:
                Q4 = Q2i.strongshock_deviation(theta - self[0].angle)
            if verbose: print(f"SWI iterations:\n  {Q3}\n  {Q4}")
            return Q4.p - Q3.p

        z0 = ang2z(.8*Q2i.angle if self[2].Mach > 1 else Q1i.angle-.5*Q1i.devmax())
        z1 = ang2z(.8*Q1i.angle if self[1].Mach > 1 else Q2i.angle+.5*Q2i.devmax())
        #z0 = ang2z(self[1].angle-.5*self[1].devmax())
        sol = optimize.root_scalar(delta_p, x0=z0, x1=z1, method='secant')#, bracket=[-30., 30.])
        #sol = optimize.root_scalar(delta_p, x0=z0, x1=z1, method='newton')#, bracket=[-30., 30.])
        sol.root = z2ang(sol.root)
        theta = sol.root
        self.solved = sol.converged
        if self[1].Mach > 1:
            self._state[3] = self[1].weakshock_deviation(theta - self[1].angle)
        else:
            self._state[3] = self[0].strongshock_deviation(theta - self[0].angle)
            self._state[1] = self[3].copy()
            self._sigma01 = sw.strongsigma_Mach_deflection(self.M0, theta-self[0].angle)
        if self[2].Mach > 1:
            self._state[4] = self[2].weakshock_deviation(theta - self[2].angle)
        else:
            self._state[4] = self[0].strongshock_deviation(theta - self[0].angle)
            self._state[2] = self[4].copy()
            self._sigma02 = sw.strongsigma_Mach_deflection(self.M0, theta-self[0].angle)
        #print(f"SWI solve end:\n  Q1:{self[1]}\n  Q2:{self[2]}\n  Q3:{self[3]}\n  Q4:{self[4]}")
        return sol

    def plot_angle_pressure(self, ax=plotsw.plt):
        # plot polar curves
        plotsw.plot_theta_pressure(self.M0, ax=ax)
        if self[1].Mach > 1:
            plotsw.plot_theta_pressure(self[1].Mach, thet_init=self[1].angle, p_init=self[1].p, color='red', ax=ax)
        if self[2].Mach > 1:
            plotsw.plot_theta_pressure(self[2].Mach, thet_init=self[2].angle, p_init=self[2].p, color='red', ax=ax)

        # plot symbols for flow regions
        for i, sty in zip([1, 2, 3, 4], ['ro', 'ro', 'rx', 'b+']):
            ax.plot(self[i].angle, self[i].p, sty)
