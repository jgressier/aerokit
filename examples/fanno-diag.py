# -*- coding: utf-8 -*-
"""
Plot of Fanno curves in several diagrams
@author: j.gressier
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import hades.aero.Fanno      as fanno

plt.rcParams['font.size'] = 14 ; plt.rcParams['lines.linewidth'] = 1.5
# for slides
#plt.rcParams['font.size'] = 24 ; plt.rcParams['grid.linewidth'] = 2 ; plt.rcParams['lines.linewidth'] = 4

npoints = 100
gam     = 1.4

Mmin = 0.1
Mmax = 4.

Mach = np.log10(np.logspace(Mmin, Mmax, npoints+1))
Fparam = fanno.maxFparam_Mach(Mach, gam)
Ts     = fanno.Ts_Tscri(Mach, gam)
Ps     = fanno.Ps_Pscri(Mach, gam)
Pi     = fanno.Pi_Picri(Mach, gam)
V      = fanno.V_Vcri(Mach, gam)
dS     = fanno.NormdS(Mach, gam)

fig=plt.figure(1, figsize=(10,8))
fig.suptitle('Ratio to critical state, $\gamma = %.1f$'%gam, y=0.93)
plt.plot(Mach, Fparam, 'k--')
plt.plot(Mach, Ts, 'r-')
plt.plot(Mach, Ps, '-', color='#000088')
plt.plot(Mach, Pi, '-', color='#0000ff')
plt.plot(Mach, V, 'k-', color='#009999')
plt.legend(['Fanno parameter', '$T_s/T^\star$', '$p_s/p^\star$', '$p_i/p_i^\star$', '$V/V^\star$'], loc='best')  
plt.axis([Mmin, Mmax, 0., 4.])
plt.xlabel('Mach number')
#plt.ylabel('shock angle $\sigma$', fontsize=10)
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Fanno-ratio.pdf', bbox_inches='tight')
plt.show()

fig=plt.figure(2, figsize=(10,8))
fig.suptitle('Fanno curve in T/S diagram, $\gamma = %.1f$'%gam, y=0.93)
plt.plot(dS, Ts, 'k')
#==============================================================================
# Mlab   = np.append(np.linspace(0.1, 1, 10), np.linspace(1, 4, 7))
# Tslab  = fanno.Ts_Tscri(Mlab, gam)
# dSlab  = fanno.NormdS(Mlab, gam)
# plt.plot(dSlab, Tslab, 'ro')
# for i in range(len(Mlab)):
#     plt.text(dSlab[i], Tslab[i], '  M=%.2g'%Mlab[i], horizontalalignment='left', verticalalignment='center',
#              fontsize=8, rotation=-30*np.sign(Mlab[i]-1))
#==============================================================================
#plt.axis([Mmin, Mmax, 0., 4.])
plt.xlabel('$\Delta S/C_p$', fontsize=10)
plt.ylabel('$h/C_p$', fontsize=10)
#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Fanno-TS.pdf', bbox_inches='tight')
#show()
