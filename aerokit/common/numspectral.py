"""definition of classes to derive spectral based differentiation operators

"""
import numpy as np
from numpy.polynomial import chebyshev
from scipy.linalg import toeplitz
from aerokit.common._dev import lazyprop
from aerokit.common.mapping import AffineMapping


class ChebCollocation:
    """Chebyshev collocation"""

    def __init__(self, npts: int, xmin=None, xmax=None, mapping=None):
        self._npts = npts
        self._max_Dorder = 0
        # default is -1 to 1 (reverse original xi distribution)
        self._xmin = -1.0 if xmin is None else xmin
        self._xmax = 1.0 if xmax is None else xmax
        self._mapping = AffineMapping(self._xmin, self._xmax) if mapping is None else mapping

    @property
    def npts(self):
        return self._npts

    @lazyprop
    def xi(self):
        th = np.arange(self.npts) * np.pi / (self.npts - 1)
        return np.cos(th)

    @lazyprop
    def x(self):
        return self._mapping.xi_to_x(self.xi)

    def extrapol(self, fk, x):
        """
        Compute the polynomial interpolant of the data (xk, fk) to new x,
        where xk are the Chebyshev nodes.

        Parameters:
        - fk: Vector of y-coordinates of data, at Chebyshev points.
        - x: Vector of x-values where polynomial interpolant is to be evaluated.

        Returns:
        - p: Vector of interpolated values.
        """
        VdMxi = chebyshev.chebvander(self.xi, self._npts-1)
        coefs = np.linalg.solve(VdMxi, fk)
        fx = chebyshev.chebval(self._mapping.x_to_xi(x), coefs)
        return fx

    def fit_to_gauss(self, x, f):
        """compute value at collocation points that best fit (x,f) pairs

        Args:
            x (_type_): _description_
            f (_type_): _description_
        """
        VdM = chebyshev.chebvander(self._mapping.x_to_xi(x), self._npts-1)
        coefs = np.linalg.lstsq(VdM, f, rcond=None)[0]
        fxi = chebyshev.chebval(self.xi, coefs)
        return fxi


    def matder(self, order: int):
        assert order >= 1
        if order > self._max_Dorder:
            self.compute_matder(order)
        return self._mapping.transform_derivative_matrices(self._matder)[:, :, order - 1]

    def compute_matder(self, maxorder):
        """
        The function computes the differentiation
        matrices D1, D2, ..., DM on Chebyshev nodes.
        """
        N = self.npts
        M = maxorder
        I = np.eye(N)
        n1 = N // 2
        n2 = -(-N // 2)  # equivalent of ceil
        irange = np.arange(N)
        Th = np.tile(irange * np.pi / (N - 1) / 2, (N, 1))
        DX = 2 * np.sin(Th + Th.T) * np.sin(Th - Th.T)
        DX = np.vstack([DX[:n1, :], -np.flipud(np.fliplr(DX[:n2, :]))])

        np.fill_diagonal(DX, 1)

        first_col = (-1.0) ** irange
        first_row = first_col.copy()  # if you want it to be symmetric
        T = toeplitz(first_col, first_row)

        C = T

        C[0, :] *= 2
        C[-1, :] *= 2
        C[:, 0] /= 2
        C[:, -1] /= 2

        Z = 1.0 / DX
        # Z[L] = 0
        np.fill_diagonal(Z, 0)
        D = I

        DM = np.zeros((N, N, M))

        for ell in range(M):
            D = (ell + 1) * Z * (C * D.diagonal()[:, np.newaxis] - D)
            # D[L] = -D.sum(axis=1)
            np.fill_diagonal(D, -D.sum(axis=1))
            DM[:, :, ell] = D
        self._matder = DM
        self._max_Dorder = maxorder


class MultiDomainChebCollocation:
    """Two Chebyshev subdomains joined at an internal interface.

    The interface node is shared, giving a grid on ``[xmin, interface]`` and
    ``[interface, xmax]`` that is clustered at both outer boundaries and at
    the interface.
    """

    def __init__(self, npts, xmin=-1.0, interface=0.0, xmax=1.0):
        if npts < 4:
            raise ValueError("each subdomain requires at least four points")
        if not xmin < interface < xmax:
            raise ValueError("interface must lie strictly between xmin and xmax")
        self._left = ChebCollocation(npts, xmin=xmin, xmax=interface)
        self._right = ChebCollocation(npts, xmin=interface, xmax=xmax)
        self._npts_per_domain = npts
        self._matder = {}

    @property
    def npts(self):
        return 2 * self._npts_per_domain - 1

    @property
    def interface_index(self):
        return self._npts_per_domain - 1

    @property
    def x(self):
        return np.r_[self._left.x, self._right.x[1:]]

    def _embedded_interface_row(self, derivative, side):
        row = np.zeros(self.npts)
        if side == "left":
            row[: self._npts_per_domain] = derivative[-1, :]
        else:
            row[self.interface_index :] = derivative[0, :]
        return row

    def matder(self, order):
        """Return a global derivative matrix with an averaged interface row."""
        if order not in self._matder:
            left = self._left.matder(order)
            right = self._right.matder(order)
            derivative = np.zeros((self.npts, self.npts))
            derivative[: self._npts_per_domain, : self._npts_per_domain] = left
            derivative[self.interface_index :, self.interface_index :] = right
            derivative[self.interface_index, :] = 0.5 * (
                self._embedded_interface_row(left, "left") + self._embedded_interface_row(right, "right")
            )
            self._matder[order] = derivative
        return self._matder[order]

    def interface_derivative_jumps(self, maxorder=3):
        """Return rows representing left-minus-right derivative jumps."""
        jumps = []
        for order in range(1, maxorder + 1):
            left = self._left.matder(order)
            right = self._right.matder(order)
            jumps.append(
                self._embedded_interface_row(left, "left")
                - self._embedded_interface_row(right, "right")
            )
        return jumps
