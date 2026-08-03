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


class AffineMapping(Mapping):
    """Affine mapping between [-1, 1] and [xmin, xmax]."""

    def __init__(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax

    def xi_to_x(self, xi):
        return self.xmin + (xi - 1.0) / (-2.0) * (self.xmax - self.xmin)

    def x_to_xi(self, x):
        return 1. - 2. * (x - self.xmin) / (self.xmax - self.xmin)
