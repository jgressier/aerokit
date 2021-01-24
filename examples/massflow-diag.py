# -*- coding: utf-8 -*-
"""
@author: j.gressier
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import aerokit.aero.MassFlow   as mf

plt.rcParams['font.size'] = 14 ; plt.rcParams['lines.linewidth'] = 1.5
# for slides
#plt.rcParams['font.size'] = 24 ; plt.rcParams['grid.linewidth'] = 2 ; plt.rcParams['lines.linewidth'] = 4

npoints = 100
gam     = 1.4

Mmin = 0.1
Mmax = 4.

Mach = np.log10(np.logspace(Mmin, Mmax, npoints+1))

debr = mf.WeightMassFlow(Mach, gam)
sig  = mf.Sigma_Mach(Mach, gam)

fig=plt.figure(1, figsize=(10,8))
#fig.suptitle('Ratio to critical state, $\gamma = %.1f$'%gam, y=0.93)
plt.plot(Mach, debr, 'k-')
#plt.axis([Mmin, Mmax, 0., 4.])
plt.xlabel('Mach number')
plt.ylabel('weighted mass flow')
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.8)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('fct-weightedmf.pdf', bbox_inches='tight')
plt.show()


fig=plt.figure(2, figsize=(10,8))
#fig.suptitle('Ratio to critical state, $\gamma = %.1f$'%gam, y=0.93)
plt.plot(Mach, sig, 'k-')
plt.axis([Mmin, 3.5, 0.5, 6.])
plt.xlabel('Mach number')
plt.ylabel('$\Sigma(M)$')
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.8)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('fct-sigma.pdf', bbox_inches='tight')
plt.show()
