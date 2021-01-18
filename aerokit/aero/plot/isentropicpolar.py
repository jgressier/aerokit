"""@package shockpolar
  Plot shock polar curves
"""

import numpy as np
import aerokit.common.defaultgas as defg # relative import is deprecated by doctest
import aerokit.aero.degree      as deg
import aerokit.aero.Supersonic  as sup

import matplotlib.pyplot as plt
from aerokit.aero.plot.defaultstyle import *


def plot_theta_pressure(mach, dev_range, gamma=defg._gamma, npts=100, thet_init=0., p_init=1., curve='C+', devmax=False, sonic=False, color='k', linestyle='-', **kwargs):
    """
    	Plot isentropic polar curve in deviation / pressure ratio axes

		Long comment
 
		:param mach:       upstream Mach number
        :param gamma:      specific heat ratio, default from aerokit.common.defaultgas
        :param npts:       number of computed points, curve accuracy
        :param thet_init:  upstream angle (shift the curve by this angle), default 0.
        :param p_init:     reference pressure (shift the curve by this ratio), default 1.
        :param curve:      choose which curve to plot (left, right or both)
		:return:     
 
 		:Example:

		.. seealso:: 
		.. note:: 
    """	
 
    dev = np.linspace(dev_range[0], dev_range[1], npts+1)
    ps  = p_init * sup.IsentropicPsratio_Mach_deflection(mach, dev, gamma=gamma)     # pressure ratio only depends on normal Mach number
    if curve in ['c+', 'C+', 'both']: plt.plot(thet_init+dev, ps, color=color, linestyle=linestyle, **kwargs)
    ps  = p_init * sup.IsentropicPsratio_Mach_deflection(mach, -dev, gamma=gamma)     # pressure ratio only depends on normal Mach number
    if curve in ['c-', 'C-', 'both']: plt.plot(thet_init+dev, ps, color=color, linestyle=linestyle, **kwargs)
