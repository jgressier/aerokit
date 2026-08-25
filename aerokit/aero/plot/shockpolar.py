"""@package shockpolar
  Plot shock polar curves
"""

import numpy as np
import aerokit.common.defaultgas as defg  # relative import is deprecated by doctest
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw

import matplotlib.pyplot as plt
from aerokit.aero.plot.defaultstyle import figure_theta_pressure, figure_theta_sigma, set_grid  # noqa: F401


def plot_theta_sigma(
    mach,
    gamma=defg._gamma,
    npts=100,
    curve='both',
    devmax=False,
    sonic=False,
    color='k',
    linestyle='-',
    ax=plt,
    **kwargs
):
    """Plot a shock polar in deviation--shock-angle coordinates.

    Args:
        mach: Upstream Mach number.
        gamma: Specific-heat ratio.
        npts: Number of points on each branch.
        curve: Branch to plot.
    """

    sig = np.linspace(deg.asin(1.0 / mach), 90.0, npts + 1)
    dev = sw.deflection_Mach_sigma(mach, sig, gamma)
    if curve in ['right', 'both']:
        ax.plot(dev, sig, color=color, linestyle=linestyle, **kwargs)
    if curve in ['left', 'both']:
        ax.plot(-dev, sig, color=color, linestyle=linestyle, **kwargs)
    if devmax:
        thet = sw.dev_Max(mach, gamma=gamma)
        sig = sw.sigma_DevMax(mach, gamma=gamma)
        if curve in ['right', 'both']:
            ax.plot(thet, sig, 'ro', alpha=0.9)
        if curve in ['left', 'both']:
            ax.plot(-thet, sig, 'ro', alpha=0.9)
    if sonic:
        thet = sw.dev_Sonic(mach, gamma=gamma)
        sig = sw.sigma_Sonic(mach, gamma=gamma)
        if curve in ['right', 'both']:
            ax.plot(thet, sig, 'ko', markerfacecolor='white')
        if curve in ['left', 'both']:
            ax.plot(-thet, sig, 'ko', markerfacecolor='white')


def plot_theta_pressure(
    mach,
    gamma=defg._gamma,
    npts=100,
    thet_init=0.0,
    p_init=1.0,
    curve='both',
    devmax=False,
    sonic=False,
    color='k',
    linestyle='-',
    ax=plt,
    sigma_range=None,
    **kwargs
):
    """Plot a shock polar in deviation--pressure-ratio coordinates.

    Args:
        mach: Upstream Mach number.
        gamma: Specific-heat ratio.
        npts: Number of points on each branch.
        thet_init: Upstream angle offset.
        p_init: Pressure-ratio offset.
        curve: Branch to plot.
        sigma_range: Optional ``(sigma_start, sigma_end)`` shock-angle interval
            in degrees. By default, the complete polar from the Mach angle to
            the normal shock is drawn.
    """

    if sigma_range is None:
        sigma_start, sigma_end = deg.asin(1.0 / mach), 90.0
    else:
        if len(sigma_range) != 2:
            raise ValueError("sigma_range must contain exactly two angles")
        sigma_start, sigma_end = sigma_range
    sig = np.linspace(sigma_start, sigma_end, npts + 1)
    dev = sw.deflection_Mach_sigma(mach, sig, gamma)
    ps = p_init * sw.Ps_ratio(mach * deg.sin(sig), gamma)  # pressure ratio only depends on normal Mach number
    if curve in ['right', 'both']:
        ax.plot(thet_init + dev, ps, color=color, linestyle=linestyle, **kwargs)
    if curve in ['left', 'both']:
        ax.plot(thet_init - dev, ps, color=color, linestyle=linestyle, **kwargs)
    if devmax:
        thet = sw.dev_Max(mach, gamma=gamma)
        sig = sw.sigma_DevMax(mach, gamma=gamma)
        ps = sw.Ps_ratio(mach * deg.sin(sig), gamma=gamma)
        if curve in ['right', 'both']:
            ax.plot(thet_init + thet, p_init * ps, 'ro', alpha=0.9)
        if curve in ['left', 'both']:
            ax.plot(thet_init - thet, p_init * ps, 'ro', alpha=0.9)
    if sonic:
        thet = sw.dev_Sonic(mach, gamma=gamma)
        sig = sw.sigma_Sonic(mach, gamma=gamma)
        ps = sw.Ps_ratio(mach * deg.sin(sig), gamma=gamma)
        if curve in ['right', 'both']:
            ax.plot(thet_init + thet, p_init * ps, 'ko', markerfacecolor='white')
        if curve in ['left', 'both']:
            ax.plot(thet_init - thet, p_init * ps, 'ko', markerfacecolor='white')


def plot_state(
    state,
    label=None,
    ax=plt,
    marker='o',
    color='k',
    markersize=8,
    offset=(6, 5),
    annotation_kwargs=None,
    **plot_kwargs
):
    """Plot and optionally annotate a flow state on a pressure polar.

    Args:
        state: Object exposing ``angle`` and ``p`` attributes.
        label: Optional annotation text.
        ax: Matplotlib axes or pyplot-compatible object.
        marker: State marker.
        color: Marker and default annotation color.
        markersize: Marker size.
        offset: Annotation offset in display points.
        annotation_kwargs: Optional overrides passed to ``annotate``.
        **plot_kwargs: Additional keyword arguments passed to ``plot``.

    Returns:
        ``(line, annotation)``; annotation is ``None`` when no label is given.
    """
    line, = ax.plot(
        state.angle,
        state.p,
        marker=marker,
        color=color,
        markersize=markersize,
        **plot_kwargs
    )
    annotation = None
    if label is not None:
        annotation_style = {
            'xytext': offset,
            'textcoords': 'offset points',
            'color': color,
            'fontweight': 'bold',
            'bbox': {'boxstyle': 'round,pad=0.2', 'fc': 'white', 'ec': color, 'alpha': 0.9},
            'arrowprops': {'arrowstyle': '-', 'color': color, 'lw': 0.8},
        }
        if annotation_kwargs is not None:
            annotation_style.update(annotation_kwargs)
        annotation = ax.annotate(label, (state.angle, state.p), **annotation_style)
    return line, annotation
