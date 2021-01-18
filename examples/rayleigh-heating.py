# -*- coding: utf-8 -*-
"""
Plot of Rayleigh curves in several diagrams
@author: j.gressier
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import aero.Rayleigh         as ray

npoints = 100
gam     = 1.4

fig=plt.figure(1, figsize=(10,8))
fig.suptitle('Mach number evolution when heating, $\gamma = %.1f$'%gam, fontsize=12, y=0.93)

for M0 in [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.2, 1.5, 2., 3., 4., 6., 10., 20.]:
#for M0 in [.1, .5]:
    Timax = ray.maxTiratio_Mach(M0, gam)
    Ti    = np.log10(np.logspace(1., Timax, npoints+1))
    Mach  = ray.SubMach_TiTicri(Ti/Timax, gam) if M0 < 1 else ray.SupMach_TiTicri(Ti/Timax, gam)
    plt.plot(Ti, Mach, '-', color='#ff0000')

plt.axis([1., 10, .1, 10.])
plt.ylabel('Mach number', fontsize=10)
plt.xlabel('$T_i/T_{i0}$', fontsize=10)
plt.loglog()
#plt.ticklabel_format(axis='both', style='plain')
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Rayleigh-heating.pdf', bbox_inches='tight')
plt.show()
