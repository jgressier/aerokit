import numpy as np
import aerokit.aero.MassFlow as mf
import matplotlib.pyplot as plt
#
g=1.4
mmax = 5.
n = 1000
Mach = np.linspace(.05, mmax, n)
Sigma = mf.Sigma_Mach(Mach)
Mres = mf.Mach_Sigma(Sigma, Mach)
fig, ax = plt.subplots(1,2)
ax[0].plot(Mach, Mres, Mach, mf.__MachSub_sigma(Sigma, g), Mach, mf.__MachSup_sigma(Sigma, g))
ax[0].set_xlim(0, mmax)
ax[0].set_ylim(0, mmax)
# Taylor expansions 
m = np.linspace(0.05, 5, 1000)
gpuogmu=(g+1.)/(g-1.)
sig1 = 1./m /((g+1)/2.)**(gpuogmu/2.)
sig2 = (1.+(3-g)/4./(g-1.)*(m-1.)**2)
sig3 = ( m**2/gpuogmu )**(gpuogmu/2.) / m
ax[1].plot(m, sig1, m, sig2, m, sig3, m, mf.Sigma_Mach(m))

plt.show()