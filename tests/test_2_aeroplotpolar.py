import matplotlib.pyplot as plt
import aerokit.aero.ShockWave as sw
import aerokit.aero.model2Dpolar as m2d
import aerokit.aero.plot.shockpolar as shp
import numpy as np


def test_sigmapolar():
    M0 = 6.
    dev = np.linspace(0.1, sw.dev_Max(M0), 40)

    fig = shp.figure_theta_sigma()

    sigw = [ sw.weaksigma_Mach_deflection(M0, d) for d in dev ]
    sigs = [ sw.strongsigma_Mach_deflection(M0, d) for d in dev ]

    shp.plot_theta_sigma(M0)
    plt.plot(dev,sigw, 'o', markersize=3)
    plt.plot(dev,sigs, 'x', markersize=3)
    #plt.show()

def test_pressurepolar():
    M0 = 6.
    dev = np.linspace(-30, sw.dev_Max(M0), 60)

    fig = shp.figure_theta_pressure()

    Q = m2d.state2Dpolar(M0)

    pw = [ Q.weakshock_deviation(d).p for d in dev ]
    ps = [ Q.strongshock_deviation(d).p for d in dev ]

    shp.plot_theta_pressure(M0)
    plt.plot(dev,pw, 'o', markersize=3)
    plt.plot(dev,ps, 'x', markersize=3)
    #plt.show()