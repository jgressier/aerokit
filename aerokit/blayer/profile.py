"""Boundary-layer profile definitions and integral thicknesses."""

import numpy as np


def _trapezoid(values, coordinates):
    """Integrate with the NumPy 1.x and 2.x trapezoidal-rule APIs."""
    if hasattr(np, "trapezoid"):
        return np.trapezoid(values, coordinates)
    return np.trapz(values, coordinates)


class Profile:
    """Generic wall-normal boundary-layer profile.

    Subclasses add the thermodynamic fields and integral definitions needed for
    a particular flow model, such as incompressible or compressible flow.
    """

    def __init__(self, y, **fields):
        self.y = np.asarray(y, dtype=float)
        if self.y.ndim != 1 or self.y.size < 2 or np.any(np.diff(self.y) <= 0.0):
            raise ValueError("y must contain at least two strictly increasing points")
        self.fields = {}
        for name, values in fields.items():
            values = np.asarray(values, dtype=float)
            if values.ndim != 1 or values.size != self.y.size:
                raise ValueError(f"{name} must be a one-dimensional array with the same size as y")
            self.fields[name] = values
            setattr(self, name, values)


class IncompressibleProfile(Profile):
    """Boundary-layer profile with an incompressible streamwise velocity."""

    def __init__(self, y, velocity):
        super().__init__(y, velocity=velocity)

    def displacement_thickness(self, velocity_edge=None):
        r"""Compute \(\delta^*=\int(1-u/U_e)\,dy\)."""
        velocity_edge = self.velocity[-1] if velocity_edge is None else velocity_edge
        if velocity_edge == 0.0:
            raise ValueError("velocity_edge must be non-zero")
        return _trapezoid(1.0 - self.velocity / velocity_edge, self.y)

    def momentum_thickness(self, velocity_edge=None):
        r"""Compute \(\theta=\int(u/U_e)(1-u/U_e)\,dy\)."""
        velocity_edge = self.velocity[-1] if velocity_edge is None else velocity_edge
        if velocity_edge == 0.0:
            raise ValueError("velocity_edge must be non-zero")
        velocity_ratio = self.velocity / velocity_edge
        return _trapezoid(velocity_ratio * (1.0 - velocity_ratio), self.y)

    def thicknesses(self, velocity_edge=None):
        """Return displacement and momentum thicknesses as a tuple."""
        return (
            self.displacement_thickness(velocity_edge),
            self.momentum_thickness(velocity_edge),
        )
