import matplotlib.pyplot as plt
import aerokit.aero.ShockWave as sw
import aerokit.aero.model2D as m2d
import aerokit.aero.plot.isentropicpolar as isp
import aerokit.aero.plot.shockpolar as shp
import numpy as np
import pytest


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
    plt.close(fig)

def test_pressurepolar():
    M0 = 6.
    dev = np.linspace(-30, sw.dev_Max(M0), 60)

    fig = shp.figure_theta_pressure()

    Q = m2d.State2DMach(M0)

    pw = [ Q.weakshock_deviation(d).p for d in dev ]
    ps = [ Q.strongshock_deviation(d).p for d in dev ]

    shp.plot_theta_pressure(M0)
    plt.plot(dev,pw, 'o', markersize=3)
    plt.plot(dev,ps, 'x', markersize=3)
    #plt.show()
    plt.close(fig)


@pytest.mark.parametrize(("curve", "nlines"), [("right", 1), ("left", 1), ("both", 2)])
def test_shockpolar_curve_selection(curve, nlines):
    fig, axis = plt.subplots()
    shp.plot_theta_sigma(2.0, npts=4, curve=curve, ax=axis)
    assert len(axis.lines) == nlines
    plt.close(fig)


def test_shockpolar_extrema_and_sonic_markers():
    fig, axis = plt.subplots()
    shp.plot_theta_pressure(2.0, thet_init=3.0, p_init=2.0, npts=4,
                            curve="both", devmax=True, sonic=True, ax=axis)
    assert len(axis.lines) == 6
    np.testing.assert_allclose(axis.lines[0].get_xdata()[0], 3.0)
    np.testing.assert_allclose(axis.lines[0].get_ydata()[0], 2.0)
    plt.close(fig)


def test_shockpolar_pressure_sigma_range():
    fig, axis = plt.subplots()
    shp.plot_theta_pressure(2.0, sigma_range=(70.0, 90.0), npts=4, ax=axis)
    assert len(axis.lines[0].get_xdata()) == 5
    np.testing.assert_allclose(axis.lines[0].get_xdata()[-1], 0.0, atol=1.0e-12)
    np.testing.assert_allclose(axis.lines[0].get_ydata()[-1], sw.Ps_ratio(2.0))
    plt.close(fig)


def test_shockpolar_sigma_extrema_and_sonic_markers():
    fig, axis = plt.subplots()
    shp.plot_theta_sigma(2.0, npts=4, curve="both", devmax=True, sonic=True, ax=axis)
    assert len(axis.lines) == 6
    plt.close(fig)


@pytest.mark.parametrize(("curve", "nlines"), [("C+", 1), ("C-", 1), ("both", 2)])
def test_isentropicpolar_curve_selection(curve, nlines):
    fig, axis = plt.subplots()
    plt.sca(axis)
    isp.plot_theta_pressure(2.0, (-5.0, 5.0), npts=4, thet_init=3.0,
                            p_init=2.0, curve=curve)
    assert len(axis.lines) == nlines
    np.testing.assert_allclose(axis.lines[0].get_xdata()[0], -2.0)
    plt.close(fig)
