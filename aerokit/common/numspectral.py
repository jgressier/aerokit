"""Chebyshev collocation grids and spectral differentiation matrices.

The operators differentiate values sampled at Chebyshev--Lobatto points.  A
mapping converts the reference coordinate to the requested physical coordinate
and applies the chain rule to the resulting matrices.  The module provides a
single-domain grid (:class:`ChebCollocation`) and a two-domain grid
(:class:`MultiDomainChebCollocation`) with a shared interface point.
"""
import numpy as np
from numpy.polynomial import chebyshev
from scipy.linalg import toeplitz
from aerokit.common._dev import lazyprop
from aerokit.common.mapping import AffineMapping


class ChebCollocation:
    """Chebyshev--Lobatto collocation grid on a mapped one-dimensional domain.

    Parameters
    ----------
    npts:
        Number of collocation points, including both endpoints.
    xmin, xmax:
        Physical endpoints for the default affine mapping.  They are ignored
        when ``mapping`` is supplied.
    mapping:
        Object implementing the :class:`aerokit.common.mapping.Mapping`
        interface.  It maps the reference coordinate ``xi`` to ``x`` and
        transforms derivative matrices with the chain rule.

    Notes
    -----
    ``xi`` is the usual decreasing ``cos(k*pi/(npts - 1))`` sequence.  The
    default :class:`~aerokit.common.mapping.AffineMapping` reverses it, so
    ``x`` is increasing from ``xmin`` to ``xmax``.  ``matder(order)`` acts on
    samples ordered as ``x`` and returns a physical-coordinate derivative.
    """

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

        Parameters
        ----------
        fk:
            Values at this object's collocation points, in ``self.x`` order.
        x:
            Physical coordinates where the interpolating polynomial is
            evaluated.

        Returns
        -------
        numpy.ndarray
            Interpolated values at ``x``.  Evaluation outside the grid range
            is polynomial extrapolation and can be ill-conditioned.
        """
        VdMxi = chebyshev.chebvander(self.xi, self._npts-1)
        coefs = np.linalg.solve(VdMxi, fk)
        fx = chebyshev.chebval(self._mapping.x_to_xi(x), coefs)
        return fx

    def fit_to_gauss(self, x, f):
        """Fit data in a Chebyshev basis and evaluate it at collocation points.

        The least-squares fit is useful when data are available away from the
        collocation grid.  It is not a Gaussian quadrature operation despite
        the historic method name.

        Parameters
        ----------
        x, f:
            Coordinates and values of the data to fit.

        Returns
        -------
        numpy.ndarray
            Fitted values in ``self.x`` order.
        """
        VdM = chebyshev.chebvander(self._mapping.x_to_xi(x), self._npts-1)
        coefs = np.linalg.lstsq(VdM, f, rcond=None)[0]
        fxi = chebyshev.chebval(self.xi, coefs)
        return fxi


    def matder(self, order: int):
        """Return the ``order``-th physical-coordinate differentiation matrix.

        The matrix is cached after its first construction.  For a vector
        ``f`` sampled at ``self.x``, ``matder(order) @ f`` approximates
        ``d**order f / dx**order`` at the same points.
        """
        assert order >= 1
        if order > self._max_Dorder:
            self.compute_matder(order)
        return self._mapping.transform_derivative_matrices(self._matder)[:, :, order - 1]

    def compute_matder(self, maxorder):
        """Compute and cache reference differentiation matrices up to an order.

        The internally stored matrices differentiate with respect to the
        Chebyshev reference coordinate.  :meth:`matder` applies the mapping's
        chain-rule transformation before returning a matrix to callers.
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
    the interface.  If each subdomain has ``npts`` points, the global grid has
    ``2 * npts - 1`` degrees of freedom and ``interface_index == npts - 1``.

    The common nodal value makes the represented field continuous.  Each
    subdomain nevertheless has its own one-sided derivative at the interface.
    :meth:`matder` uses their arithmetic mean in its interface row, which is
    appropriate for differentiating a smooth field.  It does *not* impose
    derivative continuity.  For collocation eigenvalue or boundary-value
    systems, use :meth:`interface_derivative_jumps` to replace equations with
    explicit left-minus-right continuity constraints.

    This class assumes a smooth coordinate mapping inside each affine
    subdomain.  Problems with discontinuous coefficients, fluxes, or solution
    values need their own interface flux/jump treatment.
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
        """Embed a local one-sided interface row in global-column ordering.

        ``derivative`` is an ``npts x npts`` local derivative matrix.  Since
        the local grids are ordered from left to right, row ``-1`` is the
        interface row of the left subdomain and row ``0`` is that of the
        right subdomain.  The global shared interface value occurs only once:
        at :attr:`interface_index`.

        Parameters
        ----------
        derivative:
            Local derivative matrix for one subdomain.
        side:
            ``"left"`` or ``"right"``.  It selects the respective one-sided
            interface derivative and its global block of columns.

        Returns
        -------
        numpy.ndarray
            A row of length :attr:`npts`, ready for a global matrix product.
            It represents either ``D^-_Gamma`` or ``D^+_Gamma``; it is not an
            averaged row and it imposes no condition by itself.
        """
        row = np.zeros(self.npts)
        if side == "left":
            row[: self._npts_per_domain] = derivative[-1, :]
        else:
            row[self.interface_index :] = derivative[0, :]
        return row

    def matder(self, order):
        """Return a global derivative matrix with an averaged interface row.

        Away from the interface, the matrix contains independent left and
        right Chebyshev blocks.  At the shared interface it contains

        ``0.5 * (D^-_Gamma + D^+_Gamma)``,

        where the two one-sided rows have first been embedded into the global
        column ordering.  This convention is exact for polynomials represented
        identically on both subdomains and is convenient for smooth base-state
        differentiation.  It is not a substitute for interface conditions in
        a coupled differential problem.
        """
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
        """Return left-minus-right derivative-jump rows at the interface.

        Element ``order - 1`` of the returned list is a global row ``J_m``
        such that ``J_m @ f`` equals
        ``d**m f/dx**m|left - d**m f/dx**m|right``.  A smooth field therefore
        satisfies ``J_m @ f == 0``.  For example, the Orr--Sommerfeld
        assembly can replace rows at the interface by the first three returned
        rows to impose continuity of its derivatives.
        """
        jumps = []
        for order in range(1, maxorder + 1):
            left = self._left.matder(order)
            right = self._right.matder(order)
            jumps.append(
                self._embedded_interface_row(left, "left")
                - self._embedded_interface_row(right, "right")
            )
        return jumps
