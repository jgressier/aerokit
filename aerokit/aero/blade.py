"""Two-dimensional, inviscid blade-row models.

The :class:`BladeRow` class works in an ``(x, y)`` plane.  The second
component may equivalently be the tangential component of an ``(x, theta)``
description.  Angles are expressed in degrees, as in :mod:`aerokit.aero.model2D`.

The blade velocity is ``U = (0, Uy)``.  It can describe a translating row
directly, or a rotating row locally with ``Uy = omega * radius``.  Absolute
and relative velocities obey ``C = W + U``.
"""

from __future__ import annotations

import numpy as np

from aerokit.aero import Isentropic
from aerokit.aero import MassFlow
from aerokit.aero import degree as deg
from aerokit.aero.model2D import State2d


class BladeRow:
    """Ideal two-dimensional blade row.

    Parameters
    ----------
    Uy : float, optional
        Tangential velocity of the blades.  Positive ``Uy`` is in the positive
        second-coordinate direction.
    omega, radius : float, optional
        Local rotating-row definition.  When supplied, ``Uy = omega * radius``.

    Notes
    -----
    ``downstream`` is an inviscid relative-frame model: relative total
    pressure and relative total temperature are conserved.  It is deliberately
    a blade-row building block, rather than a loss or spanwise mixing model.
    """

    def __init__(self, Uy=None, *, omega=None, radius=None):
        if Uy is not None and (omega is not None or radius is not None):
            raise ValueError("specify either Uy or omega and radius")
        if Uy is None:
            if (omega is None) != (radius is None):
                raise ValueError("omega and radius must be specified together")
            Uy = 0.0 if omega is None else omega * radius
        self.Uy = Uy

    @property
    def velocity(self):
        """Blade velocity vector ``(0, Uy)``."""
        return (0.0, self.Uy)

    def relative(self, state: State2d) -> State2d:
        """Return *state* expressed in the blade-relative reference frame."""
        return State2d(state.rho, state.u, state.v - self.Uy, state.p, state._gamma)

    def absolute(self, state: State2d) -> State2d:
        """Return a relative *state* expressed in the absolute reference frame."""
        return State2d(state.rho, state.u, state.v + self.Uy, state.p, state._gamma)

    def upstream(self, state: State2d, *, frame="absolute") -> State2d:
        """Return an independent copy of the inlet state in the requested frame."""
        self._check_frame(frame)
        return self.relative(state) if frame == "relative" else state.copy()

    def downstream(self, inlet: State2d, beta, *, frame="absolute",
                   conservation="vx", area_ratio=1.0, branch="inlet") -> State2d:
        """Turn an inlet flow to relative outlet angle *beta*.

        Parameters
        ----------
        inlet : :class:`~aerokit.aero.model2D.State2d`
            Inlet state, expressed in *frame*.
        beta : float
            Outlet relative-flow angle from the positive axial direction.
        frame : {"absolute", "relative"}, optional
            Frame of both the supplied inlet and returned outlet state.
        conservation : {"vx", "massflow"}, optional
            ``"vx"`` preserves axial velocity. ``"massflow"`` preserves the
            compressible mass flux ``rho * Vx * area`` and uses *area_ratio*
            (outlet area / inlet area).
        area_ratio : float, optional
            Outlet-to-inlet streamtube area ratio; used by ``"massflow"``.
        branch : {"inlet", "subsonic", "supersonic"}, optional
            Isentropic Mach branch used by ``"massflow"``.

        Returns
        -------
        :class:`~aerokit.aero.model2D.State2d`
            Outlet state in *frame*.
        """
        self._check_frame(frame)
        if conservation not in ("vx", "massflow"):
            raise ValueError("conservation must be 'vx' or 'massflow'")
        if area_ratio <= 0.0:
            raise ValueError("area_ratio must be positive")

        win = self.relative(inlet) if frame == "absolute" else inlet.copy()
        if np.ndim(win.rho) != 0:
            raise TypeError("BladeRow currently operates on scalar states")
        if win.u <= 0.0:
            raise ValueError("the inlet axial velocity must be positive")
        if abs(deg.cos(beta)) <= np.finfo(float).eps:
            raise ValueError("the outlet relative flow cannot be normal to x")

        if conservation == "vx":
            wx = win.u
            Mout = None
        else:
            target_flux = win.rho * win.u / area_ratio
            Mout = self._mach_for_massflux(target_flux, beta, win, branch)
            wx = None

        pt = win.Ptot()
        rtt = win.rTtot()
        gamma = win._gamma
        if Mout is None:
            # W_x specifies the velocity magnitude; the relative total
            # temperature determines the static temperature and hence Mach.
            wmag = wx / deg.cos(beta)
            a2 = gamma * rtt - 0.5 * (gamma - 1.0) * wmag**2
            if a2 <= 0.0:
                raise ValueError("outlet velocity exceeds available relative total enthalpy")
            Mout = wmag / np.sqrt(a2)

        p = pt / Isentropic.PtPs_Mach(Mout, gamma)
        rt = rtt / Isentropic.TtTs_Mach(Mout, gamma)
        wmag = Mout * np.sqrt(gamma * rt)
        outlet_relative = State2d(p / rt, wmag * deg.cos(beta),
                                  wmag * deg.sin(beta), p, gamma)
        return self.absolute(outlet_relative) if frame == "absolute" else outlet_relative

    def specific_work(self, inlet: State2d, outlet: State2d) -> float:
        """Return Euler work ``Uy * (Ctheta_out - Ctheta_in)`` per unit mass."""
        return self.Uy * (outlet.v - inlet.v)

    @staticmethod
    def _check_frame(frame):
        if frame not in ("absolute", "relative"):
            raise ValueError("frame must be 'absolute' or 'relative'")

    @staticmethod
    def _mach_for_massflux(target_flux, beta, inlet, branch):
        if branch not in ("inlet", "subsonic", "supersonic"):
            raise ValueError("branch must be 'inlet', 'subsonic', or 'supersonic'")
        if branch == "inlet":
            branch = "subsonic" if inlet.Mach() <= 1.0 else "supersonic"
        gamma = inlet._gamma
        cbeta = deg.cos(beta)

        # A positive target flux needs an outlet directed downstream.
        if target_flux * cbeta <= 0.0:
            raise ValueError("outlet angle is incompatible with the inlet mass-flow direction")
        # rho * W_x = Pt / sqrt(rTt) * WeightMassFlow(M) * cos(beta).
        # Sigma_Mach is the reciprocal of the normalized mass-flow function:
        # Sigma(M) = WeightMassFlow(1) / WeightMassFlow(M).
        sonic_flux = (inlet.Ptot() / np.sqrt(inlet.rTtot())
                      * MassFlow.WeightMassFlow(1.0, gamma=gamma) * cbeta)
        if target_flux > sonic_flux * (1.0 + 1.e-12):
            raise ValueError("requested mass flow exceeds the choked outlet capacity")
        if np.isclose(target_flux, sonic_flux):
            return 1.0
        sigma = sonic_flux / target_flux
        if branch == "subsonic":
            return MassFlow.MachSub_Sigma(sigma, gamma)
        return MassFlow.MachSup_Sigma(sigma, gamma)


# A short name is convenient in compressor-map scripts.
Blade = BladeRow
