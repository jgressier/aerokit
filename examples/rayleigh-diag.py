# -*- coding: utf-8 -*-
"""
Plot of Rayleigh curves in several diagrams
@author: j.gressier
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import aero.Rayleigh         as ray

npoints = 200
gam     = 1.4

Mmin = 0.1
Mmax = 4.

Mach = np.log10(np.logspace(Mmin, Mmax, npoints+1))
Tparam = ray.maxTiratio_Mach(Mach, gam)
Ts     = ray.Ts_Tscri(Mach, gam)
Ti     = ray.Ti_Ticri(Mach, gam)
Ps     = ray.Ps_Pscri(Mach, gam)
Pi     = ray.Pi_Picri(Mach, gam)
V      = ray.V_Vcri(Mach, gam)
dS     = ray.NormdS(Mach, gam)

fig=plt.figure(1, figsize=(10,8))
fig.suptitle('Ratio to critical state, $\gamma = %.1f$'%gam, fontsize=12, y=0.93)
#plt.plot(Mach, Tparam, 'k--')
plt.plot(Mach, Ti, '-', color='#ff0000')
plt.plot(Mach, Ts, '-', color='#882222')
plt.plot(Mach, Pi, '-', color='#0000ff')
plt.plot(Mach, Ps, '-', color='#222288')
plt.plot(Mach, V,  '-', color='#009999')
plt.legend(['$T_i/T_i^\star$', '$T_s/T^\star$', '$p_i/p_i^\star$', '$p_s/p^\star$', '$V/V^\star$'],
           loc='upper left',prop={'size':10})  
plt.axis([Mmin, Mmax, 0., 3.])
plt.xlabel('Mach number', fontsize=10)
#plt.ylabel('shock angle $\sigma$', fontsize=10)
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Rayleigh-ratio.pdf', bbox_inches='tight')
#show()

fig=plt.figure(2, figsize=(10,8))
fig.suptitle('Rayleigh curve in T/S diagram, $\gamma = %.1f$'%gam, fontsize=12, y=0.93)
plt.plot(dS, Ts, 'k')
plt.axis([-2, 0.1, .1, 1.1])
plt.xlabel('$\Delta S/C_p$', fontsize=10)
plt.ylabel('$h/C_p$', fontsize=10)
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Rayleigh-TS.pdf', bbox_inches='tight')
#show()
