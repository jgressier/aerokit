r"""
    The ``OrrSommerfeld`` module
    =========================

    Provides 1D incompressible stability equations named Orr-Sommerfeld

    Orr-Sommerfeld equation by pseudo-spectral collocation method.

    \begin{equation*}
        \Bigg[    \frac{1}{i Re}   (D^2 - \alpha^2)^2
             - (\alpha U - \omega) (D^2 - \alpha^2)
             + \alpha U''
        \Bigg] \hat{v} = 0 \qquad \text{Orr–Sommerfeld equation}
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


class DictKeyError(Exception): pass

class OrrSommerfeldModel(LinOperator):

    req_keys = ["alpha", "Reynolds", "uprofile"]

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
        req_keys = self.__class__.req_keys
        mis_keys = [k for k in req_keys if k not in state.keys()]
        if mis_keys == req_keys:
            raise DictKeyError(f"basestate misses all required keys {req_keys}")
        if mis_keys != []:
            raise DictKeyError(f"basestate misses keys {mis_keys} (keys {req_keys} required)")
        # must check u profile? (condition on uprofile? callable for 0 or 1 or other?)
        # # # if required type for uprofile is function:
        # # check uprofile.__call__ exists or:
        # # _ = uprofile(0) # should not raise TypeError: <type> object is not callable
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
        self.set_basestate({"Reynolds": Reynolds, "alpha": alpha, "uprofile": lambda x: 1 - x**2})
        self.set_BC({"type": "wall"}, {"type": "wall"})


# ===============================================================
# historic (to be deleted ?)

# def resol(alpha, omega, Rey, DiffOp):
#     """
#     We write the Orr-Sommerfeld problem
#     in the following form:    mat.X = vec

#     u   : base flow velocity
#     ddu : second derivative of u
#     e2  : second derivative matrix e2 = e*e
#     """

#     npts = DiffOp.npts
#     vec = np.zeros(npts, dtype=complex)
#     mat = np.zeros((npts, npts), dtype=complex)

#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     id_mat = np.eye(npts) # identity matrix

#     for i in range(npts):
#       for j in range(npts):
#         mat[i, j] = alpha * (-u[i] * DiffOp.matder(2)[i, j]
#                            + (u[i] * alpha**2 + ddu[i]) * id_mat[i, j])

#     z1 = 1 / (ci * Rey)
#     z2 = -2 * alpha**2 / (ci * Rey) + omega
#     z3 = alpha**4 / (ci * Rey) - omega * alpha**2
#     mat += z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_mat

#     # Boundary conditions

#     # Line 0
#     mat[0, :] = DiffOp.matder(2)[0, :]
#     vec[0] = 1

#     # Line 1
#     mat[1, :] = DiffOp.matder(1)[0, :]

#     # Line npts-2
#     mat[-2, :] = DiffOp.matder(1)[-1, :]

#     # Line npts-1
#     mat[-1, :] = 0
#     mat[-1, -1] = 1

#     sol_v = np.linalg.solve(mat, vec)
#     csol = sol_v[0]

#     return csol, sol_v

# def spectre(alpha, Rey, DiffOp, plot=True):
#     """
#     """

#     npts = DiffOp.npts
#     mat1 = np.zeros((npts, npts), dtype=complex)
#     mat2 = np.zeros((npts, npts), dtype=complex)

#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     id_mat = np.eye(npts) # identity matrix

#     for i in range(npts):
#       for j in range(npts):
#         mat1[i, j] = alpha * (-u[i] * DiffOp.matder(2)[i, j]
#                             + (u[i] * alpha**2 + ddu[i]) * id_mat[i, j])

#     z1 = 1 / (ci * Rey)
#     z2 = -2 * alpha**2 / (ci * Rey)
#     z3 = alpha**4 / (ci * Rey)
#     mat1 += z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_mat

#     z2 = -1
#     z3 = alpha**2
#     mat2 += z2 * DiffOp.matder(2) + z3 * id_mat

#     # Boundary conditions

#     # Line 0
#     mat1[0, :] = 0
#     mat1[0, 0] = 1
#     mat2[0, :] = 0

#     # Line 1
#     mat1[1, :] = DiffOp.matder(1)[0, :]
#     mat2[1, :] = 0

#     # Line npts-2
#     mat1[-2, :] = DiffOp.matder(1)[-1, :]
#     mat2[-2, :] = 0

#     # Line npts-1
#     mat1[-1, :] = 0
#     mat1[-1, -1] = 1
#     mat2[-1, :] = 0

#     l, v = eig(mat1, mat2)

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
