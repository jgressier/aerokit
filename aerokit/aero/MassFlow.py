"""
@package MassFlow
  Isentropic mass-flow and area-ratio relations for a perfect gas.

  ``WeightMassFlow`` is the dimensionless mass flux associated with a
  stagnation state.  ``Sigma_Mach`` is the corresponding area ratio ``A/A*``;
  it is one at the sonic section and has one subsonic and one supersonic Mach
  solution for every value greater than one.
"""

import warnings

import numpy as np

from scipy.optimize import brentq, newton
from ..common import defaultgas as defg

# -- Compressible flow functions  --


def WeightMassFlow(Mach, *legacy_args, gamma=defg._gamma, **legacy_kwargs):
    """Return the dimensionless isentropic mass-flow function.

    For a flow with total pressure ``Pt`` and total temperature ``Tt``, the
    dimensional mass flux is ``rho * u = Pt / sqrt(r * Tt)`` times this
    function.  ``r`` is not part of this dimensionless function. Passing it
    is accepted temporarily for backward compatibility and emits a
    ``DeprecationWarning``.
    """
    r = legacy_kwargs.pop("r", None)
    if legacy_kwargs:
        unexpected = next(iter(legacy_kwargs))
        raise TypeError("WeightMassFlow() got an unexpected keyword argument {!r}".format(unexpected))
    if len(legacy_args) > 2:
        raise TypeError("WeightMassFlow() accepts at most three positional arguments")
    if legacy_args:
        if r is not None:
            raise TypeError("WeightMassFlow() got multiple values for argument 'r'")
        r = legacy_args[0]
        if len(legacy_args) == 2:
            gamma = legacy_args[1]
    if r is not None:
        warnings.warn(
            "The 'r' argument is unused and deprecated; remove it from WeightMassFlow().",
            DeprecationWarning,
            stacklevel=2,
        )
    return (
        np.sqrt(gamma)
        * Mach
        * (1.0 + 0.5 * (gamma - 1) * Mach ** 2) ** (-0.5 * (gamma + 1.0) / (gamma - 1.0))
    )


def Sigma_Mach(Mach, gamma=defg._gamma):
    """Return the isentropic area ratio ``A/A*`` for a Mach number.

    Conservation of mass gives ``Sigma_Mach(M) =
    WeightMassFlow(1) / WeightMassFlow(M)``.  The throat is therefore the
    minimum-area, sonic section (``Sigma_Mach(1) == 1``).
    """
    return (2.0 / (gamma + 1.0) * (1.0 + 0.5 * (gamma - 1) * Mach ** 2)) ** (
        0.5 * (gamma + 1.0) / (gamma - 1.0)
    ) / Mach


def _validate_sigma(sigma):
    """Reject area ratios for which no isentropic Mach solution exists."""
    if np.any(~np.isfinite(sigma)) or np.any(np.asarray(sigma) < 1.0):
        raise ValueError("sigma must be finite and greater than or equal to one")


def __MachSub_sigma(sigma, gamma):
    """Return the historical subsonic initial estimate for ``A/A*``."""
    return 1.0 - np.sqrt((1.0 - 1.0 / sigma) * (4.0 * (gamma - 1.0) / (3.0 - gamma)))


def __MachSup_sigma(sigma, gamma):
    """Return the historical supersonic initial estimate for ``A/A*``."""
    return 1.0 + np.sqrt((sigma - 1.0) * (4.0 * (gamma - 1.0) / (3.0 - gamma)))


def _Mach_Sigma_branch(sigma, initial, subsonic, gamma):
    """Invert ``Sigma_Mach`` on one branch with Brent's method.

    ``Sigma_Mach`` decreases from infinity to one below Mach one and increases
    from one to infinity above it, so each bracket contains exactly one root.
    The initial estimate determines the first bracket; it is expanded only if
    necessary to enclose the root.
    """
    def residual(mach):
        return Sigma_Mach(mach, gamma) - sigma

    if sigma == 1.0:
        return 1.0

    if subsonic:
        initial = np.clip(initial, np.finfo(float).tiny, 1.0)
        if residual(initial) > 0.0:
            return brentq(residual, initial, 1.0)
        return brentq(residual, np.finfo(float).tiny, initial)

    upper = max(initial, 1.0)
    while residual(upper) < 0.0:
        upper *= 2.0
    return brentq(residual, 1.0, upper)


def _Mach_Sigma_solve(sigma, initial, subsonic, gamma, method):
    """Apply a selected scalar inversion while preserving array inputs."""
    _validate_sigma(sigma)
    sigma_array, initial_array, branch_array = np.broadcast_arrays(sigma, initial, subsonic)
    result = np.empty(sigma_array.shape, dtype=float)
    for index in np.ndindex(sigma_array.shape):
        if method == "brentq":
            result[index] = _Mach_Sigma_branch(
                sigma_array[index], initial_array[index], branch_array[index], gamma
            )
        elif method == "newton":
            result[index] = newton(
                lambda mach: Sigma_Mach(mach, gamma) - sigma_array[index], initial_array[index]
            )
        else:
            raise ValueError("method must be 'brentq' or 'newton'")
    return result.item() if result.ndim == 0 else result


def Mach_Sigma(sigma, Mach=2.0, gamma=defg._gamma, method="brentq"):
    """Return Mach number for area ratio *sigma* on the branch selected by *Mach*.

    Values of ``Mach <= 1`` select the subsonic branch; other values select the
    supersonic branch. ``Mach`` is also Newton's initial estimate. For array
    branch selectors, the historical branch-specific estimates are used.
    Select ``method='brentq'`` (the default) for a branch-preserving bracketed
    solve, or ``method='newton'`` for the historical Newton iteration.
    """
    _validate_sigma(sigma)
    branch = np.asarray(Mach) <= 1.0
    if np.size(Mach) != 1:
        initial = np.where(branch, __MachSub_sigma(sigma, gamma), __MachSup_sigma(sigma, gamma))
    else:
        initial = Mach
    return _Mach_Sigma_solve(sigma, initial, branch, gamma, method)


def MachSub_Sigma(sigma, gamma=defg._gamma, method="brentq"):
    """Return the subsonic Mach number associated with ``sigma = A/A*``."""
    _validate_sigma(sigma)
    return _Mach_Sigma_solve(sigma, __MachSub_sigma(sigma, gamma), True, gamma, method)


def MachSup_Sigma(sigma, gamma=defg._gamma, method="brentq"):
    """Return the supersonic Mach number associated with ``sigma = A/A*``."""
    _validate_sigma(sigma)
    return _Mach_Sigma_solve(sigma, __MachSup_sigma(sigma, gamma), False, gamma, method)
