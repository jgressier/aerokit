"""
    The ``OrrSommerfeld`` module
    =========================
 
    Provides 1D incompressible stability equations named Orr-Sommerfeld

    Orr-Sommerfeld equation by pseudo-spectral collocation method.

	\begin{equation*}
		\Bigg[ \frac{1}{i Re} (D^2 - \alpha^2)^2 - (\alpha U - \omega)(D^2 - \alpha^2) + \alpha U'' \Bigg] \hat{v} = 0 \qquad \text{Orr–Sommerfeld equation}
	\end{equation*}

	*   $D = \frac{d}{d y}$
	*   BCs: $\hat{v}(1) = \hat{v'}(1) = \hat{v}(-1) = \hat{v'}(-1) = 0$
	*   temporal theory: $\alpha \in \mathbb{R}$, $\omega \in \mathbb{C}$
	*   spatial theory: $\alpha \in \mathbb{C}$, $\omega \in \mathbb{R}$.
 
    :Example:
 
    >>> import aerokit.stability.OrrSommerfeld as OS
 
    Available functions
    -------------------
 
	.. note::
"""

import numpy as np
from aerokit.stability import LinOperator


class OrrSommerfeldModel(LinOperator):
    def __init__(self, n, xmin=None, xmax=None, basestate=None) -> None:
        """Initialization of OrrSommerfeld model

        Args:
            n (_type_): _description_
            xmin (_type_, optional): _description_. Defaults to None.
            xmax (_type_, optional): _description_. Defaults to None.
            basestate (_type_, optional): _description_. Defaults to None.
        """
        super().__init__(n, xmin, xmax)
        self._BC_dict = {"wall": self.setBC_wall}
        if basestate is not None:
            self.set_basestate(basestate)

    def set_basestate(self, state: dict):
        """set basestate for Orr-Sommerfeld model"""
        for k in ["alpha", "Reynolds", "uprofile"]:
            assert k in state.keys(), f"key {k} missing in params"
        # must check u profile ?
        super().set_basestate(state)

    def compute_operators(self):
        """compute operators for linearized
        given primitive variables P, linearized operator is
        At dv/dt = B v
        """
        assert self.check_basestate()
        alpha = self._basestate["alpha"]
        Rey = self._basestate["Reynolds"]
        D4 = self._diffop.matder(4)
        D2 = self._diffop.matder(2)
        n = self.dim
        # compute u and ddu
        u = self._basestate["uprofile"](self.x)
        ddu = D2 @ u
        # At = alpha**2 * I - D**2
        self._At = D2 - alpha**2 * np.eye(n)
        # B = 1/Re * [ D**4 - 2*alpha**2 * D**2 + alpha**4 * I ] +
        #   + alpha*j* [ u*(alpha**2-D**2) + ddu*I ]
        self._B = np.diag(-1j * alpha * u) @ self._At + np.diag(1j * alpha * ddu)
        self._B += (D4 - (2 * alpha**2) * D2 + alpha**4 * np.eye(n)) / Rey
        self.compute_BC()

    def compute_BC(self):
        n = self.dim
        self._BC_dict[self._BC_type[0]["type"]](0, 0)
        self._BC_dict[self._BC_type[1]["type"]](
            n - 1, n - 2
        )  # caution: function must write 2 BC at n-2 and n-1

    def setBC_wall(self, istate, irow):
        n = self.dim
        i0 = istate if irow is None else irow
        D = self._diffop.matder(1)
        for i in (i0, i0 + 1):
            self._At[i, :] = 0.0
            self._B[i, :] = 0.0
        # v = 0
        self._B[i0, istate] = 1.0
        # dv = 0
        self._B[i0 + 1, :n] = D[istate, :]


class Poiseuille(OrrSommerfeldModel):
    def __init__(self, n, alpha, Reynolds) -> None:
        super().__init__(n, xmin=-1.0, xmax=1.0)
        self.set_basestate(
            {"Reynolds": Reynolds, "alpha": alpha, "uprofile": lambda x: 1 - x**2}
        )
        self.set_BC({"type": "wall"}, {"type": "wall"})


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()

# def resol(alpha, omega, Rey, DiffOp):
#     """
#     We write the Orr-Sommerfeld problem
#     in the following form:    mat.X = b

#     u   : base flow velocity
#     ddu : second derivative of u
#     e2  : second derivative matrix e2 = e*e
#     """

#     npts = DiffOp.npts
#     id_mat = np.eye(npts)
#     b   = np.zeros(npts, dtype=complex)
#     mat = np.zeros((npts, npts), dtype=complex)
#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     for i in range(npts):
#       for j in range(npts):
#         mat[i,j] = alpha * (-u[i] * DiffOp.matder(2)[i,j] + (u[i] * alpha**2 + ddu[i]) * id_mat[i,j])

#     z1 = 1 / (ci * Rey)
#     z2 = -2 * alpha**2 / (ci * Rey) + omega
#     z3 = alpha**4 / (ci * Rey) - omega * alpha**2
#     mat = mat + z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_mat

#     mat[0,:] = DiffOp.matder(2)[0,:]
#     mat[1,:] = DiffOp.matder(1)[0,:]
#     mat[npts-2,:] = DiffOp.matder(1)[npts-1,:]
#     mat[npts-1,:] = 0
#     b[:] = 0
#     mat[-1,-1] = 1
#     b[0] = 1

#     sol_v = np.linalg.solve(mat, b)
#     csol = sol_v[0]
#     return csol, sol_v

# def spectre(alpha, Rey, DiffOp, plot=True):
#     n = DiffOp.npts
#     mat1 = np.zeros((n, n), dtype=complex)
#     mat2 = np.zeros((n, n), dtype=complex)

#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     id_matrix = np.eye(n) # identity matrix

#     for i in range(n):
#       for j in range(n):
#         mat1[i, j] = alpha * (-u[i] * DiffOp.matder(2)[i, j] + (u[i] * alpha**2 + ddu[i]) * id_matrix[i, j])

#     z1 = 1 / (ci * Rey)
#     z2 = -2 * alpha**2 / (ci * Rey)
#     z3 = alpha**4 / (ci * Rey)
#     mat1 += z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_matrix

#     z2 = -1
#     z3 = alpha**2
#     mat2 += z2 * DiffOp.matder(2) + z3 * id_matrix

#     # Boundary conditions
#     mat1[0, :] = 0
#     mat1[0, 0] = 1
#     mat2[0, :] = 0

#     mat1[1, :] = DiffOp.matder(1)[0, :]
#     mat2[1, :] = 0

#     mat1[n-2, :] = DiffOp.matder(1)[n-1, :]
#     mat2[n-2, :] = 0

#     mat1[n-1, :] = 0
#     mat1[n-1, n-1] = 1
#     mat2[n-1, :] = 0

#     l,v = eig(mat1, mat2)

#     if plot:
#         import matplotlib.pyplot as plt
#         plt.plot(l.real, l.imag, 'ob', markerdim=6)
#         plt.hlines(y=0.0, xmin=-2, xmax=2, color='k', linestyle='--')
#         plt.xlabel('$\omega_r$')
#         plt.ylabel('$\omega_i$')
#         plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
#         plt.title('Spectrum')
#         # plt.axis('square')
#         plt.axis([0.1, 1, -2, 0.5])
#         plt.show()
#     # else:
#     #     vp[:4] = l[:4]
#     return l, v


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()
