"""Backward-compatible imports for :mod:`aerokit.instance.riemann`."""

import warnings

from aerokit.instance.riemann import riemann_pb

warnings.warn(
    "aerokit.aero.riemann is deprecated; import aerokit.instance.riemann instead",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["riemann_pb"]
