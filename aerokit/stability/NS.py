"""
The ``NS`` module
=========================

Provides local compressible stability equations

:Example:

>>> import aerokit.stability.NS as NS

Available functions
-------------------

    .. note::
"""

from typing import Union
import numpy as np
from aerokit.stability._base import LinOperator
from aerokit.common.mapping import Mapping


class NSaxi(LinOperator):

    req_keys = ["kx", "m", "rho", "P", "Ux", "gamma"]

    def __init__(self, n, rmin=0.0, rmax=None, basestate={}, mapping=None) -> None:
        super().__init__(n, xmin=rmin, xmax=rmax)
        self._mapping = mapping
        if mapping is not None:
            self._r = mapping.xi_to_x(self._diffop.xi)[::-1]
        if basestate is not None:
            self.set_basestate(basestate)
        # if not provided, may be initialized later

    def set_basestate(self, state):
        super().set_basestate(state)
        self.nvar = 5
        # vars are perturbations of rho, ux, ur, utheta, rT

    @property
    def r(self):
        return self._diffop.x if self._mapping is None else self._r

    @property
    def x(self):
        return self.r

    def _radial_matder(self, order):
        if self._mapping is None:
            return self._diffop.matder(order)
        self._diffop.compute_matder(order)
        derivatives = self._mapping.transform_derivative_matrices(self._diffop._matder)
        return derivatives[::-1, ::-1, order - 1]

    def compute_operators(self):
        """
        compute operators for linearized
        given primitive variables P, linearized operator is
        At dP/dt + B1 dP/dr + B0 P = 0
        """
        D = self._radial_matder(1)
        n = self.dim
        N = self.dim * self.nvar
        irho = slice(0, n)
        iux = slice(n, 2 * n)
        iur = slice(2 * n, 3 * n)
        iut = slice(3 * n, 4 * n)
        irT = slice(4 * n, N)
        if not np.isclose(self.r[0], 0.0):
            raise ValueError("NSaxi requires the radial grid to start at the axis (r = 0)")
        rcorr = np.maximum(
            self.r, self.r[1] / 10.0
        )  # alias for cylindrical coordinates and corrected to avoid zero
        # parameters
        kx, m = self._basestate['kx'], self._basestate['m']
        gam = self._basestate['gamma']
        # base state and derivatives
        rho, Ux, P = (self._basestate[var] for var in ['rho', 'Ux', 'P'])
        # Compute base-flow gradients with the operator differentiation matrix.
        drho = D @ rho
        dUx = D @ Ux
        dT = D @ (P / rho)
        # order is rho, ux, ur, utheta, rT
        self._At = np.eye(N)
        # B0
        # self._B0 = np.zeros((N, N))
        self._B0 = 1j * kx * np.diag(np.tile(Ux, self.nvar))  # i*kx*Ux on diagonal
        # continuity
        self._B0[irho, iux] += 1j * kx * np.diag(rho)
        self._B0[irho, iur] += np.diag(drho + rho / rcorr)
        self._B0[irho, iut] += 1j * m * np.diag(rho / rcorr)
        # momentum ux
        self._B0[iux, irho] += 1j * kx * np.diag(P / rho**2)
        self._B0[iux, iur] += np.diag(dUx)
        self._B0[iux, irT] += 1j * kx * np.eye(n)
        # momentum ur
        self._B0[iur, irho] = np.diag(dT / rho)
        self._B0[iur, irT] = np.diag(drho / rho)
        # momentum utheta
        self._B0[iut, irho] = 1j * m * np.diag(P / rho**2 / rcorr)
        self._B0[iut, irT] = 1j * m * np.diag(1.0 / rcorr)
        # temperature
        self._B0[irT, iux] += 1j * kx * (gam - 1) * np.diag(P / rho)
        self._B0[irT, iur] += np.diag(dT + (gam - 1) * P / rho / rcorr)
        self._B0[irT, iut] += 1j * m * (gam - 1) * np.diag(P / rho / rcorr)
        # B1
        # self._B1 = np.diag(np.tile(q.u, self.nvar)) @ np.kron(np.eye(self.nvar), D)
        self._B1 = np.zeros((N, N))
        self._B1[irho, iur] = np.diag(rho) @ D
        self._B1[iur, irho] = np.diag(P / rho**2) @ D
        self._B1[iur, irT] = D
        self._B1[irT, iur] = np.diag((gam - 1) * P / rho) @ D
        #
        self._B = -self._B1 - self._B0
        self.compute_BC()

    def compute_BC(self):
        self.setBC_axis(self._basestate['m'])
        self.setBC_far()

    def _clear_row(self, irow):
        self._At[irow, :] = 0.0
        self._B[irow, :] = 0.0

    def setBC_axis(self, m):
        n = self.dim
        D = self._radial_matder(1)
        irho, iux, iur, iut, irT = (ivar * n for ivar in range(self.nvar))
        if m == 0:
            # Scalars are even; both transverse velocity components are odd.
            for irow, icol in ((irho, irho), (iux, iux), (irT, irT)):
                self._clear_row(irow)
                self._B[irow, icol : icol + n] = D[0, :]
            for irow in (iur, iut):
                self._clear_row(irow)
                self._B[irow, irow] = 1.0
        elif abs(m) == 1:
            # Scalars vanish, while regular Cartesian transverse velocity
            # requires utheta - i*sign(m)*ur = 0 and d(utheta + i*sign(m)*ur)/dr = 0.
            for irow in (irho, iux, irT):
                self._clear_row(irow)
                self._B[irow, irow] = 1.0
            sign_m = np.sign(m)
            self._clear_row(iur)
            self._B[iur, iur] = -1j * sign_m
            self._B[iur, iut] = 1.0
            self._clear_row(iut)
            self._B[iut, iur : iur + n] = 1j * sign_m * D[0, :]
            self._B[iut, iut : iut + n] = D[0, :]
        else:  # m >= 2
            # All perturbation components vanish as r**(|m|-1) or faster.
            for irow in (irho, iux, iur, iut, irT):
                self._clear_row(irow)
                self._B[irow, irow] = 1.0

    def setBC_far(self):
        n = self.dim
        Dn = self._radial_matder(1)[-1, :]
        for ivar in range(self.nvar):
            irow = (ivar + 1) * n - 1
            self._clear_row(irow)
            self._B[irow, irow] = 1.0


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()
