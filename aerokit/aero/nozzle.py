"""Backward-compatible imports for :mod:`aerokit.instance.nozzle`."""

import warnings

from aerokit.instance.nozzle import (
    Madapt_from_AsAc_NPR,
    Ms_from_AsAc_NPR,
    NPR_choked_subsonic,
    NPR_choked_supersonic,
    NPR_shock_at_exit,
    nozzle,
)

warnings.warn(
    "aerokit.aero.nozzle is deprecated; import aerokit.instance.nozzle instead",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "Madapt_from_AsAc_NPR",
    "Ms_from_AsAc_NPR",
    "NPR_choked_subsonic",
    "NPR_choked_supersonic",
    "NPR_shock_at_exit",
    "nozzle",
]
