"""definition of classes to derive spectral based differentiation operators

"""
import numpy as np
from scipy.linalg import toeplitz
import matplotlib.pyplot as plt
from aerokit.common._dev import lazyprop


class ChebCollocation():
    """Chebyshev collocation
    """
    def __init__(self, npts: int, xmin = None, xmax = None):
        self._npts = npts
        self._max_Dorder = 0
        # default is -1 to 1 (reverse original xi distribution)
        self._xmin = -1. if xmin is None else xmin
        self._xmax = 1. if xmax is None else xmax

    @property
    def npts(self):
        return self._npts

    @lazyprop
    def xi(self):
        th = np.arange(self.npts) * np.pi / (self.npts - 1)
        return np.cos(th)

    @lazyprop
    def x(self):
        return self._xmin + (self.xi-1.)/(-2.)*(self._xmax - self._xmin)

    def extrapol(self, fk, x):
        """
        Compute the polynomial interpolant of the data (xk, fk),
        where xk are the Chebyshev nodes.

        Parameters:
        - fk: Vector of y-coordinates of data, at Chebyshev points.
        - x: Vector of x-values where polynomial interpolant is to be evaluated.

        Returns:
        - p: Vector of interpolated values.
        """
        fk = np.array(fk).reshape(-1, 1)
        x = np.array(x).reshape(-1, 1)

        N = len(fk)
        M = len(x)

        xk = np.sin(
            np.pi
            * np.array(list(range(N - 1, 0, -2)) + list(range(0, -N, -2)))
            / (2 * (N - 1))
        ).reshape(-1, 1)

        # Compute weights for Chebyshev formula
        w = np.ones(N) * (-1) ** np.array(range(N))
        w[0] = w[0] / 2
        w[-1] = w[-1] / 2

        # Compute quantities x-x(k) and their reciprocals
        D = x - xk.T
        eps = np.finfo(float).eps
        D = 1.0 / (D + eps * (D == 0))

        # Evaluate interpolant as matrix-vector products
        p = np.sum(D * (w * fk.T), axis=1) / np.sum(D * w, axis=1)
        return p.reshape(-1, 1)

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
        x = self.x
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
            DM[:, :, ell] = D * ((self._xmax-self._xmin)/(-2.))**(ell+1)
        self._matder = DM
        self._max_Dorder = maxorder
