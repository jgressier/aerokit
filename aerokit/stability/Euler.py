
"""
    The ``Euler`` module
    =========================
 
    Provides 1D compressible stability equations 
 
    :Example:
 
    >>> import aerokit.stability.OrrSommerfeld as OS
 
    Available functions
    -------------------
 
	.. note:: 
"""

import numpy as np
import aerokit.common.numspectral as ns
from aerokit.stability import LinOperator
import aerokit.aero.model1D as m1d 
import aerokit.aero.ShockWave as sw


class Euler1D(LinOperator):

    def __init__(self, n, basestate: m1d.state = None) -> None:
        super().__init__(n)
        if basestate is not None:
            self.set_basestate(basestate)

    def set_basestate(self, state):
        super().set_basestate(state)
        assert state.size == self.dim
        self.size = self.dim
        self.nvar = 3

    def compute_operators(self):
        """compute operators for linearized
          given primitive variables P, linearized operator is
          At dP/dt + Bx dP/dt + B0 P = 0
        """
        D = self._diffop.matder(1)
        n = self.size
        N = self.size * self.nvar
        q = self._basestate
        # order is rho, u, p
        self._At = np.eye(N)
        # Bx
        self._Bx = np.diag(np.tile(q.u, self.nvar))
        self._Bx[0:n,n:2*n] = np.diag(q.rho)
        self._Bx[n:2*n,2*n:] = np.diag(1./q.rho)
        self._Bx[2*n:,n:2*n] = np.diag(q._gamma * q.p)
        self._Bx = self._Bx @ np.kron(np.eye(self.nvar), D)
        # derivatives of basestate
        dq = m1d.state(rho = D @ q.rho, u = D @ q.u, p= D @ q.p)
        # B0
        self._B0 = np.zeros((N, N))
        self._B0[:2*n,:2*n] = np.diag(np.tile(dq.u, 2))
        self._B0[2*n:,2*n:] = np.diag(q._gamma*dq.u)
        self._B0[:n:,n:2*n] = np.diag(dq.rho)
        self._B0[n:2*n:,:n] = np.diag(-dq.p/q.rho**2)
        self._B0[2*n:,n:2*n] = np.diag(dq.p)
        self._B = self._Bx + self._B0

    def setBC(self, Ltype: str, Rtype: str):
        n = self.size
        if Ltype == 'per':
            assert Rtype == 'per'
            self.setBC_per()
        BCfunc = { 'sym': self.setBC_sym, 'RH': self.setBC_RH}
        BCfunc[Ltype](0, 0)
        BCfunc[Rtype](n-1, n-1)

    def setBC_sym(self, istate, irow):
        n = self.size
        i0 = istate if irow is None else irow
        D = self._diffop.matder(1)
        for i in (0, n, 2*n):
            self._At[i0+i,:] = 0.
            self._B[i0+i,:] = 0.
        # drho = 0
        self._B[i0,:n] = D[istate,:]
        # u = 0
        self._B[i0+n,istate] = 1.
        # dp = 0
        self._B[i0+2*n,2*n:] = D[istate,:]

    def setBC_per(self):
        n = self.size
        D = self._diffop.matder(1)
        for i in (0, n-1, n, 2*n-1, 2*n, 3*n-1):
            self._At[i,:] = 0.
            self._B[i,:] = 0.
        for i in range(self.nvar):
            self._B[i*n,i*n] = 1.
            self._B[i*n,i*n+n-1] = -1.
            self._B[i*n+n-1,i*n:(i+1)*n] = D[0,:]
            self._B[i*n+n-1,i*n:(i+1)*n] -= D[n-1,:]

    def setBC_RH(self, istate: int, irow: int = None):
        """Apply linerized Rankine Hugoniot equations from constant upstream value.

        There are 3 equations but 4 unknows (3 downstream perturbations and shock velocity). Substituting
        the shock velocity s reduces the boundary equations to 2.
        At immediate downstream state 1i, the two RH BC can be written
            k_rho[0:1]*rho1i' + k_u[0:1]*u1i' + k_p[0:1]*p1i' = 0
        Then q1i' perturbations are replaced par q1m' perturbation at mean position using
            q1i' = q1m' + dq/dx|1m * dxs 
            with s=dxs/dt

        Args:
            istate (int): downstream state of RH equations
            irow (int): local row corresponding of where equation is written
        """
        n = self.size
        D = self._diffop.matder(1)
        q = self._basestate
        dq = m1d.state(rho = D @ q.rho, u = D @ q.u, p= D @ q.p) 
        g = q._gamma
        M1 = q.Mach()[istate]
        # M0 is upstream but equations are symmetric
        M0 = sw.downstream_Mn(M1, g)
        rhoratio = sw.Rho_ratio(M0, g)
        pratio = sw.Ps_ratio(M0, g)
        a0 = q.asound()[istate] / np.sqrt(pratio/rhoratio)
        rho0 = q.rho[istate]/rhoratio
        # BC coefs k[0:1,0:2] stands for coefficients in 1st and 2nd BC for rho', u', p'
        # RH continuity equation: a0*M0' = srho*rho1' + su*u1'
        srho = q.u[istate]/(rho0*(1.-rhoratio))
        su = rhoratio*(1.-rhoratio)
        k = np.zeros((2,3))
        # BC for rho1'
        kdrho = 4/(g+1.)*rhoratio**2 / M0**3  # coefficient for rho1'/rho0 = kdrho * MO'
        k[0,0] = 1./rho0 - kdrho/a0*srho
        k[0,1] = -kdrho/a0*su
        # BC for p1'
        kdp = 4*g/(g+1.)*M0
        k[1,0] = -kdp/a0*srho
        k[1,1] = -kdp/a0*su
        k[1,2] = 1./q.p[istate]
        # change BC at 1i (immediate downstream) to 1m (mean 1 downstream) and write BC to At and B
        #   d/dt(q1i') = d/dt(q1m') + dq/dx|1m * s
        i0 = istate if irow is None else irow
        irho = 0
        ip = 2*n
        # RH terms on temporal matrix
        for ibc, ieq in enumerate((irho, ip)):
            self._At[ieq+i0,:] = 0.
            self._B[ieq+i0,:] = 0.
            for j in range(3):
                self._At[ieq+i0,j*n+istate] = k[ibc,j]
        # shock shift terms on B: dq/dx * a0 * MO'
        dq1 = np.array([dq.rho[istate], dq.u[istate], dq.p[istate]])
        for ibc, ieq in enumerate((irho, ip)):
            for iterm in range(3):
                self._B[ieq,istate] += k[ibc,iterm]*dq1[iterm]*srho
                self._B[ieq,2*n+istate] += k[ibc,iterm]*dq1[iterm]*su



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()