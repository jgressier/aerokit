"""
    The ``refstate`` module
    =========================
 
    Provides reference state module for normalization and similarity computations
  
    Available functions
    -------------------
 
"""

class refstate:
    """Dimensional reference quantities used for nondimensionalization."""

    def __init__(self, density=1.0, velocity=1.0, viscosity=1.0, length=1.0, **kwargs):
        self.density = density
        self.velocity = velocity
        self.viscosity = viscosity
        self.length = length
        self._dict = {
            "density": density,
            "velocity": velocity,
            "viscosity": viscosity,
            "length": length,
            **kwargs,
        }

    @property
    def Reynolds(self):
        """Return ``rho * U * L / mu`` using dynamic viscosity ``mu``."""
        if self.viscosity == 0.0:
            raise ValueError("viscosity must be non-zero")
        return self.density * self.velocity * self.length / self.viscosity

    @property
    def reynolds(self):
        """Lower-case alias for :attr:`Reynolds`."""
        return self.Reynolds
