import numpy as np
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw
import matplotlib.pyplot as plt

# check integration regularity

M0=2.
sigmin=deg.asin(1./M0)
sigmax=sw.sigma_Sonic(M0)

sig = np.linspace(sigmin*1.01, sigmax, 200)
dev = 0.*sig
dev2 = dev
for i,s in enumerate(sig):
    dev[i], dev2[i] = sw.conical_deflection_Mach_sigma(M0, s)
plt.plot(sig, (dev2-dev)*1e20)
plt.show()