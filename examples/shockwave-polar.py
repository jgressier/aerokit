# -*- coding: utf-8 -*-
"""
Plot of local Rankine-Hugoniot equations (2D shock waves)
@author: j.gressier
"""

import aero.degree           as deg
#import aero.CompressibleFlow as aerof
import aero.ShockWave        as aerosw
import numpy                 as np
import matplotlib.pyplot     as plt

npoints = 100
gam     = 1.4

fig=plt.figure(1, figsize=(10,8))
fig.suptitle('Polar of Shock-Waves, $\gamma = %.1f$'%gam, fontsize=12, y=0.93)

macharray = [ 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.5, 3., 3.5, 4., 5., 10., 100. ]

for m in macharray:
    sig = np.linspace(deg.asin(1./m), 90., npoints+1)
    dev = aerosw.deflection_Mach_sigma(m, sig, gam)
    plt.plot(dev, sig, 'k-')
    plt.text(dev[npoints/2], sig[npoints/2], '%.3g'%m, horizontalalignment='center', verticalalignment='top',
             fontsize=8, bbox=dict(facecolor='white', alpha=0.8),
             rotation='60')
             
    #labels.append(legends[i]+", t=%.1f"%results[i][t].time)
#legend(labels, loc='upper left',prop={'size':10})  
plt.axis([0., 50., 0., 90.])
plt.xlabel('deviation $\Delta\\theta$', fontsize=10)
plt.ylabel('shock angle $\sigma$', fontsize=10)
plt.minorticks_on()
plt.grid(which='major', linestyle='-', alpha=0.8)
plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('polar.png', bbox_inches='tight')
#show()
