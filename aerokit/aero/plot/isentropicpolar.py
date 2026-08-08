"""@package shockpolar
  Plot shock polar curves
"""

import numpy as np
import aerokit.common.defaultgas as defg  # relative import is deprecated by doctest
import aerokit.aero.Supersonic as sup

import matplotlib.pyplot as plt
from aerokit.aero.plot.defaultstyle import *


def plot_theta_pressure(
    mach,
    dev_range,
    gamma=defg._gamma,
    npts=100,
    thet_init=0.0,
    p_init=1.0,
    curve='C+',
    devmax=False,
    sonic=False,
    color='k',
    linestyle='-',
    **kwargs
):
    """Plot an isentropic polar in deviation--pressure-ratio coordinates.

    Args:
        mach: Upstream Mach number.
        dev_range: Deviation-angle range.
        gamma: Specific-heat ratio.
        npts: Number of points on each branch.
        thet_init: Upstream angle offset.
        p_init: Pressure-ratio offset.
        curve: Branch to plot.
    """

    dev = np.linspace(dev_range[0], dev_range[1], npts + 1)
    ps = p_init * sup.IsentropicPsratio_Mach_deflection(
        mach, dev, gamma=gamma
    )  # pressure ratio only depends on normal Mach number
    if curve in ['c+', 'C+', 'both']:
        plt.plot(thet_init + dev, ps, color=color, linestyle=linestyle, **kwargs)
    ps = p_init * sup.IsentropicPsratio_Mach_deflection(
        mach, -dev, gamma=gamma
    )  # pressure ratio only depends on normal Mach number
    if curve in ['c-', 'C-', 'both']:
        plt.plot(thet_init + dev, ps, color=color, linestyle=linestyle, **kwargs)
