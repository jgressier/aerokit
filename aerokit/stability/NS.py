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

    def __init__(self, n, rmin=0., rmax=None, basestate={}, mapping=None) -> None:
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
        """compute operators for linearized
        given primitive variables P, linearized operator is
        At dP/dt + B1 dP/dr + B0 P = 0
        """
        D = self._radial_matder(1)
        n = self.dim
        N = self.dim * self.nvar
        assert np.isclose(self.r[0], 0)
        rcorr =  np.maximum(self.r, self.r[1]/10.)  # alias for cylindrical coordinates and corrected to avoid zero
        # parameters
        kx, m = self._basestate['kx'], self._basestate['m']
        gam = self._basestate['gamma']
        # base state and derivatives
        rho, Ux, P = self._basestate['rho'], self._basestate['Ux'], self._basestate['P']
        drho, dUx, dT = (D @ rho, D @ Ux, D @ (P/rho))
        # order is rho, ux, ur, utheta, rT
        self._At = np.eye(N)
        # B0
        #self._B0 = np.zeros((N, N))
        self._B0 = 1j*kx*np.diag(np.tile(Ux, self.nvar)) # i*kx*Ux on diagonal
        self._B0[0:n, n:2*n] += 1j*kx*np.diag(rho)
        self._B0[0:n, 2*n:3*n] += np.diag(drho + 1j*kx*rho/rcorr)
        self._B0[0:n, 3*n:4*n] += 1j*m*np.diag(rho/rcorr)
        self._B0[n:2*n, 2*n:3*n] += 1j*kx*np.diag(P/rho**2)
        self._B0[n:2*n, 3*n:4*n] += np.diag(dUx)
        self._B0[n:2*n, 4*n:N] += 1j*kx*np.eye(n)
        self._B0[2*n:3*n, 0:n] = np.diag(dT) 
        self._B0[2*n:3*n, 4*n:N] = np.diag(drho/rho)
        self._B0[3*n:4*n, 0:n] = 1j*m*np.diag(P/rho**2/rcorr) 
        self._B0[3*n:4*n, 4*n:N] = 1j*m*np.diag(1./rcorr)
        self._B0[4*n:N, n:2*n] += 1j*kx*(gam-1)*np.diag(P/rho)
        self._B0[4*n:N, 2*n:3*n] += np.diag(dT + (gam-1)*P/rho/rcorr)
        self._B0[4*n:N, 3*n:4*n] += 1j*m*np.diag(rho/rcorr)
        # B1
        #self._B1 = np.diag(np.tile(q.u, self.nvar)) @ np.kron(np.eye(self.nvar), D)
        self._B1 = np.zeros((N, N))
        self._B1[0:n, 2*n:3*n] = np.diag(rho) @ D
        self._B1[2*n:3*n, 0:n] = np.diag(P/rho**2) @ D
        self._B1[2*n:3*n, 4*n:N] = D
        self._B1[4*n:N, 2*n:3*n] = np.diag((gam-1)*P/rho) @ D
        #
        self._B = -self._B1 - self._B0
        self.compute_BC()

    def compute_BC(self):
        self.setBC_axis(self._basestate['m'])
        self.setBC_far()

    def setBC_axis(self, m):
        n = self.dim
        Dn = self._radial_matder(1)[-1,:]
        for ivar in range(self.nvar):
            irow = ivar*n
            self._At[irow, :] = 0.0
            self._B[irow, :] = 0.0
            self._B[irow, irow] = 1.0

    def setBC_far(self):
        n = self.dim
        Dn = self._radial_matder(1)[-1,:]
        for ivar in range(self.nvar):
            irow = (ivar + 1)*n - 1
            self._At[irow, :] = 0.0
            self._B[irow, :] = 0.0
            self._B[irow, irow] = 1.0


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()
