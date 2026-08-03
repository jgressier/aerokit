"""Coordinate mappings."""
from abc import ABC, abstractmethod
import numpy as np


class Mapping(ABC):
    """Base class for coordinate mappings."""

    @abstractmethod
    def xi_to_x(self, xi):
        pass

    @abstractmethod
    def x_to_xi(self, x):
        pass

    @abstractmethod
    def scale_derivative(self, derivative, order):
        pass

    @abstractmethod
    def dxi_dx(self):
        pass

    @abstractmethod
    def d2xi_dx2(self):
        pass

    @abstractmethod
    def d3xi_dx3(self):
        pass

    @abstractmethod
    def d4xi_dx4(self):
        pass

    def transform_derivative_matrices(self, derivatives):
        """Transform reference-coordinate derivative matrices to physical ones."""
        if derivatives.shape[2] > 4:
            raise NotImplementedError("derivatives above fourth order are not supported")

        physical = derivatives.copy()
        dxi_dx = np.asarray(self.dxi_dx())[..., np.newaxis]
        if derivatives.shape[2] >= 1:
            physical[:, :, 0] = dxi_dx * derivatives[:, :, 0]
        if derivatives.shape[2] >= 2:
            d2xi_dx2 = np.asarray(self.d2xi_dx2())[..., np.newaxis]
            physical[:, :, 1] = dxi_dx**2 * derivatives[:, :, 1] + d2xi_dx2 * derivatives[:, :, 0]
        if derivatives.shape[2] >= 3:
            d3xi_dx3 = np.asarray(self.d3xi_dx3())[..., np.newaxis]
            physical[:, :, 2] = (
                dxi_dx**3 * derivatives[:, :, 2]
                + 3 * dxi_dx * d2xi_dx2 * derivatives[:, :, 1]
                + d3xi_dx3 * derivatives[:, :, 0]
            )
        if derivatives.shape[2] >= 4:
            d4xi_dx4 = np.asarray(self.d4xi_dx4())[..., np.newaxis]
            physical[:, :, 3] = (
                dxi_dx**4 * derivatives[:, :, 3]
                + 6 * dxi_dx**2 * d2xi_dx2 * derivatives[:, :, 2]
                + (3 * d2xi_dx2**2 + 4 * dxi_dx * d3xi_dx3) * derivatives[:, :, 1]
                + d4xi_dx4 * derivatives[:, :, 0]
            )
        return physical


class AffineMapping(Mapping):
    """Affine mapping between [-1, 1] and [xmin, xmax]."""

    def __init__(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax

    def xi_to_x(self, xi):
        return self.xmin + (xi - 1.0) / (-2.0) * (self.xmax - self.xmin)

    def x_to_xi(self, x):
        return 1. - 2. * (x - self.xmin) / (self.xmax - self.xmin)

    def scale_derivative(self, derivative, order):
        return derivative * (-2.0 / (self.xmax - self.xmin)) ** order

    def dxi_dx(self):
        return -2.0 / (self.xmax - self.xmin)

    def d2xi_dx2(self):
        return 0.0

    def d3xi_dx3(self):
        return 0.0

    def d4xi_dx4(self):
        return 0.0

    def transform_derivative_matrices(self, derivatives):
        if derivatives.shape[2] <= 4:
            return super().transform_derivative_matrices(derivatives)
        return np.stack(
            [self.scale_derivative(derivatives[:, :, order], order + 1) for order in range(derivatives.shape[2])],
            axis=2,
        )
