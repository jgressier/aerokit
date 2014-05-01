# -*- coding: utf-8 -*-
"""
Plot of local Rankine-Hugoniot equations (2D shock waves)
@author: j.gressier
"""

import aero.degree           as deg
import aero.CompressibleFlow as aerof
import aero.ShockWave        as aerosw
import numpy                 as np
import matplotlib.pyplot     as plt

npoints = 100
gam     = 1.4

fig=plt.figure(1, figsize=(10,8))
fig.suptitle('Polar of isentropic compression and Shock-Waves, $\gamma = %.1f$'%gam, fontsize=12, y=0.93)

macharray = [ 1.2, 1.5, 2., 3., 5.]
rdev      = np.linspace(0., 50., npoints+1)

for m in macharray:
    sig  = np.linspace(deg.asin(1./m), 90., npoints+1)
    dev  = aerosw.deflection_Mach_sigma(m, sig, gam)
    kpsw = aerosw.Ps_ratio(m*deg.sin(sig), gam)     # pressure ratio only depends on normal Mach number
    rdev = np.linspace(0., .99*aerof.PrandtlMeyer_Mach(m, gam), npoints+1)
    kpis = aerof.IsentropicPsratio_Mach_deflection(m, rdev, gam) 
    plt.plot(dev,  kpsw, 'k-')
    plt.plot(rdev, kpis, 'k--', alpha=0.6)
    plt.plot(rdev[-1], kpis[-1], 'ro', alpha=0.5)
    plt.text(rdev[-1], kpis[-1], 'M=1', fontsize=7, clip_on='true', horizontalalignment='right') 
    plt.text(dev[npoints/3], kpsw[npoints/3], '%.3g'%m, horizontalalignment='left', verticalalignment='top',
             fontsize=8, bbox=dict(facecolor='white', alpha=0.8),
             rotation='0')
             
    #labels.append(legends[i]+", t=%.1f"%results[i][t].time)
#legend(labels, loc='upper left',prop={'size':10})  
plt.xlim(0., 50.)
plt.ylim(1., 50.)
plt.xlabel('deviation $\Delta\\theta$', fontsize=10)
plt.ylabel('pressure ratio through deviation', fontsize=10)
plt.yscale('log')
plt.minorticks_on()
plt.grid(which='major', linestyle='-', alpha=0.7)
plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('polar-p.png', bbox_inches='tight')
#show()
