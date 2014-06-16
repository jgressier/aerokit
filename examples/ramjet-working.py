# -*- coding: utf-8 -*-
"""
Design of Ramjet
@author: j.gressier
all sections are normalized by A0, section of upstream flow
- 1 is flow section in inlet (possibly decelerated by external compression)
- 2 is the inlet throat
- 21 or 22 are the upstream and downstream sections at the normal shock in the inlet
- 3 is the beginning of the combustor
- 4 is the end of the combustor (A3=A4)
- 8 is the throat of the nozzle
- 9 is the nozzle exit
"""

import numpy                 as np
import matplotlib.pyplot     as plt
import aero.Isentropic       as isent
import aero.MassFlow         as mf
import aero.ShockWave        as sw
import aero.Rayleigh         as ray

gam  = 1.4   # assumed constant
npts = 50
M0   = 3.
M2   = 1.5
A2A0 = mf.Sigma_Mach(M2, gam)/mf.Sigma_Mach(M0, gam)
print "A0/A2 for isentropic compression from Mach %4.2f to %4.2f : %6.3f"%(M0, M2, 1./A2A0)
A3A2 = 1./A2A0
A8A2 = .6*A3A2

fig=plt.figure(1, figsize=(10,12))
fig.suptitle('Flow features on RAMJET, $M_2=%.2f$, $A_8/A_2=%.2f$, $\gamma = %.1f$'
    %(M2, A8A2, gam), fontsize=12, y=0.93)

plt.subplot(311)
plt.ylabel('$p_{i3}/p_{i0}$', fontsize=10)
plt.grid(which='major', linestyle=':', alpha=0.5)

plt.subplot(312)
plt.ylabel('$p_{i4}/p_{i0}$', fontsize=10)
plt.grid(which='major', linestyle=':', alpha=0.5)

plt.subplot(313)
plt.ylabel('$M_{4}$', fontsize=10)
plt.grid(which='major', linestyle=':', alpha=0.5)

# --- (FS) fully supersonic flow (not expected on design)

M3sup    = mf.Mach_Sigma(A3A2*mf.Sigma_Mach(M2, gam), gam)
M4max    = mf.Mach_Sigma(A3A2/A8A2, 2., gam)  # look for supersonic value
alphamax = ray.Ti_Ticri(M4max, gam)/ray.Ti_Ticri(M3sup, gam)
print "unchoking of fully supersonic flow for Ti4/Ti0 = %6.3f"%(alphamax)

FSalpha = np.log10(np.logspace(1., alphamax, npts+1))
FSm4    = ray.SupMach_TiTicri(FSalpha/alphamax*ray.Ti_Ticri(M4max, gam), gam)
FSpi4   = ray.Pi_Picri(FSm4, gam)/ray.Pi_Picri(M3sup, gam)
FSpi3   = np.ones(npts+1)
plt.subplot(311)
plt.ylabel('$p_{i3}/p_{i0}$', fontsize=10)
plt.plot(FSalpha, FSpi3, '-', color='#ff0000')
plt.subplot(312)
plt.ylabel('$p_{i4}/p_{i0}$', fontsize=10)
plt.plot(FSalpha, FSpi4, '-', color='#ff0000')
plt.subplot(313)
plt.ylabel('$M_{4}$', fontsize=10)
plt.plot(FSalpha, FSm4, '-', color='#ff0000')

# --- (CW) conventional working state

# M8 is sonic so M4 is known
CWm4 = mf.Mach_Sigma(A3A2/A8A2, .1, gam)  # look for subsonic value

CWm3low  = sw.downstream_Mn(M3sup, gam)
alphamin = ray.Ti_Ticri(CWm4, gam)/ray.Ti_Ticri(CWm3low, gam)
CWm3high = mf.Mach_Sigma(A3A2*mf.Sigma_Mach(sw.downstream_Mn(M2, gam), gam), .1, gam)
alphamax = ray.Ti_Ticri(CWm4, gam)/ray.Ti_Ticri(CWm3high, gam)
print "         unchoking of inlet throat for Ti4/Ti0 = %6.3f"%(alphamax)
print "             fully supersonic flow for Ti4/Ti0 = %6.3f"%(alphamin)

CWalpha = np.log10(np.logspace(alphamin, alphamax, npts+1))
CWm3    = ray.SubMach_TiTicri(ray.Ti_Ticri(CWm4, gam)/CWalpha, gam)
CWpi3   = mf.Sigma_Mach(CWm3, gam)/mf.Sigma_Mach(M2, gam)/A3A2
CWpi4   = CWpi3*ray.Pi_Picri(CWm4, gam)/ray.Pi_Picri(CWm3, gam)
CWm4    = np.ones(npts+1)*CWm4

plt.subplot(311)
plt.plot(CWalpha, CWpi3, '-', color='#bb0000')
plt.subplot(312)
plt.plot(CWalpha, CWpi4, '-', color='#bb0000')
plt.subplot(313)
plt.plot(CWalpha, CWm4, '-', color='#bb0000')
plt.show

# --- unchoked channel


#plt.minorticks_on()
plt.grid(which='major', linestyle=':', alpha=0.5)
#plt.grid(which='minor', linestyle=':', alpha=0.5)
fig.savefig('Rayleigh-TS.pdf', bbox_inches='tight')
#show()
