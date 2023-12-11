"""
    The ``Euler`` module
    =========================
 
    Provides 1D compressible stability equations 
 
    :Example:
 
    >>> import aerokit.stability.Euler as Euler
 
    Available functions
    -------------------
 
	.. note:: 
"""

from typing import Union
import numpy as np
from aerokit.stability._base import LinOperator
import aerokit.aero.model1D as m1d
import aerokit.aero.ShockWave as sw


class Euler1D(LinOperator):
    def __init__(self, n, xmin=None, xmax=None, basestate=None) -> None:
        super().__init__(n, xmin, xmax)
        self._BC_dict = {
            "sym": self.setBC_sym,
            "per": self.setBC_per,
            "RH": self.setBC_RH,
        }
        if basestate is not None:
            assert isinstance(basestate, m1d.state)
            self.set_basestate(basestate)

    def set_basestate(self, state):
        super().set_basestate(state)
        self.nvar = 3

    def compute_operators(self):
        """compute operators for linearized
        given primitive variables P, linearized operator is
        At dP/dt + Bx dP/dt + B0 P = 0
        """
        D = self._diffop.matder(1)
        n = self.dim
        N = self.dim * self.nvar
        q = self._basestate
        assert isinstance(q, m1d.state)
        # order is rho, u, p
        self._At = np.eye(N)
        # Bx
        ip = 2 * n
        self._Bx = np.diag(np.tile(q.u, self.nvar)) @ np.kron(np.eye(self.nvar), D)
        self._Bx[0:n, n:ip] = np.diag(q.rho) @ D
        self._Bx[n:ip, ip:] = np.diag(1.0 / q.rho) @ D
        self._Bx[ip:, n:ip] = np.diag(q._gamma * q.p) @ D
        # derivatives of basestate
        dq = m1d.state(rho=D @ q.rho, u=D @ q.u, p=D @ q.p)
        # B0
        self._B0 = np.zeros((N, N))
        self._B0[:ip, :ip] = np.diag(np.tile(dq.u, 2))
        self._B0[ip:, ip:] = np.diag(q._gamma * dq.u)
        self._B0[:n, n:ip] = np.diag(dq.rho)
        self._B0[n:ip, :n] = np.diag(-dq.p / q.rho ** 2)
        self._B0[ip:, n:ip] = np.diag(dq.p)
        self._B = -self._Bx - self._B0
        self.compute_BC()

    def compute_BC(self):
        n = self.dim
        if self._BC_type is None:
            raise ValueError("BC have not been set: use self.set_BC(Ltype, Rtype)")
        if self._BC_type[0]["type"] == "per":
            assert self._BC_type[1]["type"] == "per"
            self.setBC_per()
        else:
            bcdict = self._BC_type[0]
            bctype = bcdict.pop('type')
            self._BC_dict[bctype](0, 0, **bcdict)
            bcdict = self._BC_type[1]
            bctype = bcdict.pop('type')
            self._BC_dict[bctype](n - 1, n - 1, **bcdict)

    def setBC_sym(self, istate, irow):
        n = self.dim
        i0 = istate if irow is None else irow
        irho = 0
        iu = n
        ip = 2 * n
        il = 3 * n  # must specify last in A and B augmented
        D = self._diffop.matder(1)
        for i in (irho, iu, ip):
            self._At[i0 + i, :] = 0.0
            self._B[i0 + i, :] = 0.0
        # drho = 0
        self._B[i0, :iu] = D[istate, :]
        # u = 0
        self._B[i0 + iu, iu + istate] = 1.0
        # dp = 0
        self._B[i0 + ip, ip:il] = D[istate, :]

    def setBC_per(self):
        n = self.dim
        D = self._diffop.matder(1)
        for i in (0, n - 1, n, 2 * n - 1, 2 * n, 3 * n - 1):
            self._At[i, :] = 0.0
            self._B[i, :] = 0.0
        for i in range(self.nvar):
            self._B[i * n, i * n] = 1.0
            self._B[i * n, i * n + n - 1] = -1.0
            self._B[i * n + n - 1, i * n : (i + 1) * n] = D[0, :]
            self._B[i * n + n - 1, i * n : (i + 1) * n] -= D[n - 1, :]

    def setBC_RH(self, istate: int, irow: Union[int, None] = None, **kwargs):
        """Apply linerized Rankine Hugoniot equations from constant upstream value.

        There are 3 equations but 4 unknows (3 downstream perturbations and shock velocity or position).
        Substituting the shock velocity s reduces the boundary equations to 2.
        At immediate downstream state 1s, the two RH BC can be written
            k_rho[0:1]*rho1s' + k_u[0:1]*u1s' + k_p[0:1]*p1s' = 0
        Then q1s' perturbations are replaced par q1m' perturbation at mean position using
            q1s' = q1m' + dq/dx|1m * dxs
            with s=dxs/dt

        Args:
            istate (int): downstream state of RH equations
            irow (int): local row corresponding of where equation is written
        """
        # sizes
        n = self.dim
        D = self._diffop.matder(1)
        # parse additional argument
        dM0dx = kwargs.pop('dM0dx') if 'dM0dx' in kwargs else None
        if kwargs:
            raise ValueError(f"unknown extra parameters for RH BC in Euler model: {kwargs.keys()}")
        # increase matrix size to add "dx" shock oscillation amplitude
        self._At = np.block([[self._At, np.zeros((3 * n, 1))], [np.zeros((1, 3 * n)), np.zeros((1, 1))]])
        self._B = np.block([[self._B, np.zeros((3 * n, 1))], [np.zeros((1, 3 * n)), np.zeros((1, 1))]])
        # --- base state ---
        q = self._basestate
        assert isinstance(q, m1d.state)
        dq = m1d.state(rho=D @ q.rho, u=D @ q.u, p=D @ q.p)
        q1 = q[istate]
        dq1 = dq[istate]
        g = q1._gamma
        M1 = q1.Mach()
        # M0 is upstream but equations are symmetric
        q0 = q1.state_RH()
        print(q0, q1)
        M0 = sw.downstream_Mn(M1, g)
        rhoratio = sw.Rho_ratio(M0, g)
        pratio = sw.Ps_ratio(M0, g)
        a0 = q0.asound()
        # parse additional argument
        print("dM0dx:", dM0dx)
        if dM0dx is None:  # if not defined, compute upstream gradient from downstream
            print("upstream gradient not set, computed from RH and base state")
            raise ValueError("not yet implemented")
        # BC coefs bck0[0:2,0:2] in 3 BC equations for rho', u', p' at 0i
        # BC coefs bck1[0:2,0:2] in 3 BC equations for rho', u', p' at 1i
        # BC coefs bckSh[0:2,0:1] in 3 BC equations for dxs and s
        bck0 = np.zeros((3, 3))
        bck1 = np.zeros((3, 3))
        bckSh = np.zeros((3, 2))
        # BC for rho1' (RH rho ratio from M0)
        drhoratio = 4 / (g + 1.0) * rhoratio ** 2 / M0 ** 3  # coefficient for rho1'/rho0 = kdrho * MO'
        bck0[0, 0] = -rhoratio / q0.rho
        bck1[0, 0] = 1.0 / q0.rho
        bckSh[0, 0] = -drhoratio * dM0dx
        bckSh[0, 1] = drhoratio / a0
        # BC for u1' (RH mass equation)
        bck0[1, 0] = -q0.u / q0.rho
        bck0[1, 1] = -1.0
        bck1[1, 0] = q1.u / q0.rho
        bck1[1, 1] = rhoratio
        bckSh[1, 1] = 1 - rhoratio
        # BC for p1' (RH p ratio from M0)
        dpratio = 4 * g / (g + 1.0) * M0
        bck0[2, 2] = -pratio / q0.p
        bck1[2, 2] = 1.0 / q0.p
        bckSh[2, 0] = -dpratio * dM0dx
        bckSh[2, 1] = dpratio / a0
        # change BC at 0i (immediate downstream) to 0m (mean 0 upstream)
        #   q0i' = dq/dx|0m * dxs
        # bck0 coeffs feed bckSh[:,0] dxs coef
        dx = self.x[1] - self.x[0]
        q0p = q0.state_isentropic_Mach(M0 + dM0dx * dx)
        drho0 = (q0p.rho - q0.rho) / dx
        du0 = (q0p.u - q0.u) / dx
        dp0 = (q0p.p - q0.p) / dx
        for ibc in range(3):
            for ivar, idq0 in enumerate((drho0, du0, dp0)):
                bckSh[ibc, 0] += bck0[ibc, ivar] * idq0
        # change BC at 1i (immediate downstream) to 1m (mean 1 downstream)
        #   q1i' = q1m' + dq/dx|1m * dxs
        # bck1 coeffs are the same for 1i et 1m AND feed bckSh[:,0] dxs coef
        for ibc in range(3):
            for ivar, dq0 in enumerate((dq1.rho, dq1.u, dq1.p)):
                bckSh[ibc, 0] += bck1[ibc, ivar] * dq0
        #
        i0 = istate if irow is None else irow
        irho = 0
        #iu = n
        ip = 2 * n
        idxs = 3 * n
        # Build  At dq/dx = B.q
        for ibc, ieq in enumerate((irho, idxs, ip)):
            self._At[ieq + i0, :] = 0.0
            self._B[ieq + i0, :] = 0.0
            # At contains only dependency on s = d(dxs)/dt : bckSh[:,1] coefs
            self._At[ieq + i0, idxs] = -bckSh[ibc, 1]
            # B contains dependency on dxs
            self._B[ieq, idxs] = bckSh[ibc, 0]
            # B contains dependencies on rho1m', u1m', p1m'
            for ivar in range(3):
                self._B[ieq, ivar * n + istate] = bck1[ibc, ivar]


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()
