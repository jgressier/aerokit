import matplotlib.pyplot as plt
import aerokit.aero.ShockWave as sw
import aerokit.aero.plot.shockpolar as shp
import numpy as np

M0 = 2.

dev = np.linspace(-10, 10, 100)

fig, ax = plt.subplots(1, 1, figsize=(4,3))

sigw = [ sw.weaksigma_Mach_deflection(M0, d) for d in dev ]
sigs = [ sw.strongsigma_Mach_deflection(M0, d) for d in dev ]

shp.plot_theta_sigma(M0, ax=ax)
ax.plot(dev,sigw)
ax.plot(dev,sigs)
plt.show()