"""Coordinate mappings."""
from abc import ABC, abstractmethod


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
