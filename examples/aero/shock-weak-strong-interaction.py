import matplotlib.pyplot as plt
import aerokit.aero.degree as deg
import aerokit.aero.ShockWave as sw
import aerokit.aero.plot.shockpolar as shp
import aerokit.instance.SWinteraction as SWI
import numpy as np

counter = 0

for M0 in [2., 3., 4., 6, 10.]:
    # define a range from 10% to 80% max deviation
    dev = np.linspace(.1*sw.dev_Sonic(M0), .8*sw.dev_Sonic(M0), 50)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    th = np.zeros_like(dev)
    sig1 = np.zeros_like(dev)
    sig2 = np.zeros_like(dev)

    for i, d in enumerate(dev):
        sig01 = sw.sigma_Mach_deflection(M0, d)
        swi = SWI.ShockInteraction(M0, sig01, -90.)
        sol = swi.solve()
        if not sol.converged: counter += 1
        th[i] = swi[4].angle if sol.converged else 0.
        sig1[i] = swi.sigma01
        sig2[i] = swi.sigma02
    ax1.plot(dev, th)
    ax1.plot(dev, sig1)
    ax1.plot(dev, deg.acos(deg.cos(sig2)))
    ax1.grid()
    #ax1.set_ylim(-95., 95.)
    ax1.set_title(f"M0={M0}")

    id = 40
    swi = SWI.ShockInteraction(M0, sig1[id], -90.)
    sol = swi.solve(verbose=True)
    swi.plot_angle_pressure(ax=ax2)
    ax2.set_title(f"$\\theta$={dev[id]:.1f}° - $\sigma_1$={sig1[id]:.1f}°")

    print(f"save weakstrong-M0={M0}.png")
    plt.savefig(f"weakstrong-M0={M0}.png")

print(f"{counter} shock interactions not solved")