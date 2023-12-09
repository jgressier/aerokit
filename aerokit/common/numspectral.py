"""definition of classes to derive spectral based differentiation operators

"""
import numpy as np
from numpy.polynomial import chebyshev
from scipy.linalg import toeplitz
from aerokit.common._dev import lazyprop


class ChebCollocation:
    """Chebyshev collocation"""

    def __init__(self, npts: int, xmin=None, xmax=None):
        self._npts = npts
        self._max_Dorder = 0
        # default is -1 to 1 (reverse original xi distribution)
        self._xmin = -1.0 if xmin is None else xmin
        self._xmax = 1.0 if xmax is None else xmax

    @property
    def npts(self):
        return self._npts

    @lazyprop
    def xi(self):
        th = np.arange(self.npts) * np.pi / (self.npts - 1)
        return np.cos(th)

    @lazyprop
    def x(self):
        return self._xi_to_x(self.xi)

    def _x_to_xi(self, x):
        return 1.-2.*(x-self._xmin)/ (self._xmax - self._xmin)

    def _xi_to_x(self, xi):
        """compute x for any xi"""
        return self._xmin + (xi - 1.0) / (-2.0) * (self._xmax - self._xmin)

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
        fx = chebyshev.chebval(self._x_to_xi(x), coefs)
        return fx

    def fit_to_gauss(self, x, f):
        """compute value at collocation points that best fit (x,f) pairs

        Args:
            x (_type_): _description_
            f (_type_): _description_
        """
        VdM = chebyshev.chebvander(self._x_to_xi(x), self._npts-1)
        coefs = np.linalg.lstsq(VdM, f)[0]
        fxi = chebyshev.chebval(self.xi, coefs)
        return fxi


    def matder(self, order: int):
        assert order >= 1
        if order > self._max_Dorder:
            self.compute_matder(order)
        return self._matder[:, :, order - 1]

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
            DM[:, :, ell] = D * (-2.0 / (self._xmax - self._xmin)) ** (ell + 1)
        self._matder = DM
        self._max_Dorder = maxorder
