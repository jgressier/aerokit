"""@package shockpolar
  Plot shock polar curves
"""

import numpy as np
import hades.common.defaultgas as defg # relative import is deprecated by doctest
import hades.aero.degree      as deg
import hades.aero.ShockWave   as sw

import matplotlib.pyplot as plt
import hades.aero.plot.defaultstyle as pltdef


def figure_theta_sigma(**kwargs):
	fig = plt.figure(**kwargs)
	pltdef.set_grid()
	return fig

def figure_theta_pressure(**kwargs):
	fig = plt.figure(**kwargs)
	pltdef.set_grid()
	plt.yscale('log')
	return fig


def plot_theta_sigma(mach, gamma=defg._gamma, npts=100, curve='both', devmax=False, sonic=False, color='k', linestyle='-', **kwargs):
    """
    	Plot shock polar curve in deviation / shock angle axes

		Long comment
 
		:param mach:       upstream Mach number
        :param gamma:      specific heat ratio, default from hades.common.defaultgas
        :param npts:       number of computed points, curve accuracy
        :param curve:      choose which curve to plot (left, right or both)
		:return:     
 
 		:Example:

		.. seealso:: 
		.. note:: 
    """	
 
    sig = np.linspace(deg.asin(1./mach), 90., npts+1)
    dev = sw.deflection_Mach_sigma(mach, sig, gamma)
    if curve in ['right', 'both']: plt.plot( dev, sig, color=color, linestyle=linestyle, **kwargs)
    if curve in ['left',  'both']: plt.plot(-dev, sig, color=color, linestyle=linestyle, **kwargs)
    if devmax:
    	thet = sw.dev_Max(mach, gamma=gamma)
    	sig  = sw.sigma_DevMax(mach, gamma=gamma)
    	if curve in ['right', 'both']: plt.plot( thet, sig, 'ro', alpha=0.9)
    	if curve in ['left', 'both']:  plt.plot(-thet, sig, 'ro', alpha=0.9)
    if sonic:
    	thet = sw.dev_Sonic(mach, gamma=gamma)
    	sig  = sw.sigma_Sonic(mach, gamma=gamma)
    	if curve in ['right', 'both']: plt.plot( thet, sig, 'wo')
    	if curve in ['left', 'both']:  plt.plot(-thet, sig, 'wo')


def plot_theta_pressure(mach, gamma=defg._gamma, npts=100, thet_init=0., p_init=1., curve='both', devmax=False, sonic=False, color='k', linestyle='-', **kwargs):
    """
    	Plot shock polar curve in deviation / pressure ratio axes

		Long comment
 
		:param mach:       upstream Mach number
        :param gamma:      specific heat ratio, default from hades.common.defaultgas
        :param npts:       number of computed points, curve accuracy
        :param thet_init:  upstream angle (shift the curve by this angle), default 0.
        :param p_init:     reference pressure (shift the curve by this ratio), default 1.
        :param curve:      choose which curve to plot (left, right or both)
		:return:     
 
 		:Example:

		.. seealso:: 
		.. note:: 
    """	
 
    sig = np.linspace(deg.asin(1./mach), 90., npts+1)
    dev = sw.deflection_Mach_sigma(mach, sig, gamma)
    ps  = p_init * sw.Ps_ratio(mach*deg.sin(sig), gamma)     # pressure ratio only depends on normal Mach number
    if curve in ['right', 'both']: plt.plot(thet_init+dev, ps, color=color, linestyle=linestyle, **kwargs)
    if curve in ['left',  'both']: plt.plot(thet_init-dev, ps, color=color, linestyle=linestyle, **kwargs)
    if devmax:
    	thet = sw.dev_Max(mach, gamma=gamma)
    	sig  = sw.sigma_DevMax(mach, gamma=gamma)
    	ps   = sw.Ps_ratio(mach*deg.sin(sig), gamma=gamma)
    	if curve in ['right', 'both']: plt.plot(thet_init+thet, p_init*ps, 'ro', alpha=0.9)
    	if curve in ['left', 'both']:  plt.plot(thet_init-thet, p_init*ps, 'ro', alpha=0.9)
    if sonic:
    	thet = sw.dev_Sonic(mach, gamma=gamma)
    	sig  = sw.sigma_Sonic(mach, gamma=gamma)
    	ps   = sw.Ps_ratio(mach*deg.sin(sig), gamma=gamma)
    	if curve in ['right', 'both']: plt.plot(thet_init+thet, p_init*ps, 'wo')
    	if curve in ['left', 'both']:  plt.plot(thet_init-thet, p_init*ps, 'wo')


