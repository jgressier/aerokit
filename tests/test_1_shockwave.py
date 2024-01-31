from aerokit.aero import degree
import aerokit.aero.ShockWave as sw
import numpy as np
import pytest


def test_shockwave_Ps():
    assert sw.Ps_ratio(1.0) == 1.0
    assert sw.Ps_ratio(1.1) == pytest.approx(1.245)
    assert sw.Ps_ratio(2.0) == pytest.approx(4.5)
    assert sw.Ps_ratio(1.0, gamma=1.3) == 1.0
    assert sw.Ps_ratio(1.1, gamma=1.3) == pytest.approx(1.2373913043478264)
    assert sw.Ps_ratio(2.0, gamma=1.3) == pytest.approx(4.391304347826088)


def test_shockwave_Rho():
    assert sw.Rho_ratio(1.0) == 1.0
    assert sw.Rho_ratio(1.1) == pytest.approx(1.1690821256038648)
    assert sw.Rho_ratio(2.0) == pytest.approx(2.666666666666667)
    assert sw.Rho_ratio(1.0, gamma=1.3) == 1.0
    assert sw.Rho_ratio(1.1, gamma=1.3) == pytest.approx(1.1777401608125266)
    assert sw.Rho_ratio(2.0, gamma=1.3) == pytest.approx(2.875)


def test_shockwave_Ts():
    assert sw.Ts_ratio(1.0) == 1.0
    assert sw.Ts_ratio(1.1) == pytest.approx(1.0649380165289257)
    assert sw.Ts_ratio(2.0) == pytest.approx(1.6875)
    assert sw.Ts_ratio(1.0, gamma=1.3) == 1.0
    assert sw.Ts_ratio(1.1, gamma=1.3) == pytest.approx(1.0506488150103892)
    assert sw.Ts_ratio(2.0, gamma=1.3) == pytest.approx(1.527410207939509)


def test_shockwave_Mn1():
    assert sw.downstream_Mn(1.0) == 1.0
    assert sw.downstream_Mn(1.1) == pytest.approx(0.9117704213259055)
    assert sw.downstream_Mn(2.0) == pytest.approx(0.5773502691896257)
    assert sw.downstream_Mn(1.0, gamma=1.3) == 1.0
    assert sw.downstream_Mn(1.1, gamma=1.3) == pytest.approx(0.911201472607656)
    assert sw.downstream_Mn(2.0, gamma=1.3) == pytest.approx(0.5628780357842335)


# @pytest.mark.parametrize("m", [.001, .01, .1, .5, 1., 2., 5.])
# def test_sigma_def_massflow(m):
#     assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(mf.WeightMassFlow(1., gamma=1.3)/mf.WeightMassFlow(m, gamma=1.3))


def test_shockwave_Mn1_involutive_numpy():
    m = np.linspace(1.0, 10.0, 30)
    np.testing.assert_allclose(m, sw.downstream_Mn(sw.downstream_Mn(m)))


def test_Mn_Ps():
    assert sw.Ps_ratio(2.0) == pytest.approx(4.5)
    assert sw.Mn_Ps_ratio(sw.Ps_ratio(3.0)) == pytest.approx(3.0)


@pytest.mark.parametrize("M0, dev", [(1.5, 5.0), (3.0, 10.0), (4.0, 30.0)])
def test_polar_iterative_vs_cubic_weak(M0, dev):
    sig_it = sw.sigma_Mach_deflection(M0, dev)  # default is weak shock
    sig_cub = sw.weaksigma_Mach_deflection(M0, dev)
    assert sig_it == pytest.approx(sig_cub)


@pytest.mark.parametrize("M0, dev", [(1.5, 5.0), (3.0, 10.0), (4.0, 30.0)])
def test_polar_iterative_vs_cubic_strong(M0, dev):
    sig_it = sw.sigma_Mach_deflection(M0, dev, init=80.0)
    sig_cub = sw.strongsigma_Mach_deflection(M0, dev)
    assert sig_it == pytest.approx(sig_cub)


def test_conical_shock_deflection():
    dev = sw.conical_deflection_Mach_sigma(2.0, 35.0)
    assert dev == pytest.approx(16.5322, rel=1.0e-4)  # solution of iterative process


def test_conical_shock_sigma():
    sig = sw.conical_sigma_Mach_walldeflection(2.0, 30.0)
    assert sig == pytest.approx(48.079078, rel=1.0e-4)  # solution of iterative process


def test_conical_shock_mach():
    mach = sw.conical_Mach_walldeflection_sigma(30.0, 45.0)
    assert mach == pytest.approx(2.2376, rel=1.0e-4)  # solution of iterative process


@pytest.mark.parametrize("mach", [1.1, 2., 10.])
def test_devmax_consistency(mach):
    sig = sw.sigma_DevMax(mach)
    devmax = sw.dev_Max(mach)
    devm, devp = (sw.deflection_Mach_sigma(mach, coef*sig) for coef in (.99, 1.01))
    assert devm < devmax
    assert devp < devmax


@pytest.mark.parametrize("gam", [1.2, 1.4, 1.6])
def test_devmaxMinf(gam):
    Minf = 1e10
    assert sw.dev_MaxMinf(gam) == pytest.approx(sw.dev_Max(Minf, gam))


@pytest.mark.parametrize("mach", [1.1, 2., 10.])
def test_M1sonic(mach):
    sig = sw.sigma_Sonic(mach)
    dev = sw.dev_Sonic(mach)
    Mn1 = sw.downstream_Mn(mach*degree.sin(sig))
    assert Mn1/degree.sin(sig-dev) == pytest.approx(1.)


@pytest.mark.parametrize("dev", [1., 5., 10., 40.])
def test_mach_devmax(dev):
    M0max = sw.Mach_DevMax(dev)
    assert M0max > 1.
    newdev = sw.dev_Max(M0max)
    assert newdev == pytest.approx(dev)


@pytest.mark.parametrize("dev", [1., 5., 10., 40.])
def test_mach_sonic(dev):
    M0 = sw.Mach_DevSonic(dev)
    sig = sw.sigma_Sonic(M0)
    Mn1 = sw.downstream_Mn(M0*degree.sin(sig))
    assert Mn1/degree.sin(sig-dev) == pytest.approx(1.)
