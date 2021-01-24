import aerokit.aero.ShockWave as sw
import numpy as np
import pytest

def test_shockwave_Ps():
    assert sw.Ps_ratio(1.)  == 1.
    assert sw.Ps_ratio(1.1) == pytest.approx(1.245)
    assert sw.Ps_ratio(2.)  == pytest.approx(4.5)
    assert sw.Ps_ratio(1., gamma=1.3)  == 1.
    assert sw.Ps_ratio(1.1, gamma=1.3) == pytest.approx(1.2373913043478264)
    assert sw.Ps_ratio(2., gamma=1.3)  == pytest.approx(4.391304347826088)

def test_shockwave_Rho():
    assert sw.Rho_ratio(1.)  == 1.
    assert sw.Rho_ratio(1.1) == pytest.approx(1.1690821256038648)
    assert sw.Rho_ratio(2.)  == pytest.approx(2.666666666666667)
    assert sw.Rho_ratio(1., gamma=1.3)  == 1.
    assert sw.Rho_ratio(1.1, gamma=1.3) == pytest.approx(1.1777401608125266)
    assert sw.Rho_ratio(2., gamma=1.3)  == pytest.approx(2.875)

def test_shockwave_Ts():
    assert sw.Ts_ratio(1.)  == 1.
    assert sw.Ts_ratio(1.1) == pytest.approx(1.0649380165289257)
    assert sw.Ts_ratio(2.)  == pytest.approx(1.6875)
    assert sw.Ts_ratio(1., gamma=1.3)  == 1.
    assert sw.Ts_ratio(1.1, gamma=1.3) == pytest.approx(1.0506488150103892)
    assert sw.Ts_ratio(2., gamma=1.3)  == pytest.approx(1.527410207939509)

def test_shockwave_Mn1():
    assert sw.downstream_Mn(1.)  == 1.
    assert sw.downstream_Mn(1.1) == pytest.approx(0.9117704213259055)
    assert sw.downstream_Mn(2.)  == pytest.approx(0.5773502691896257)
    assert sw.downstream_Mn(1., gamma=1.3)  == 1.
    assert sw.downstream_Mn(1.1, gamma=1.3) == pytest.approx(0.911201472607656)
    assert sw.downstream_Mn(2., gamma=1.3)  == pytest.approx(0.5628780357842335)

# @pytest.mark.parametrize("m", [.001, .01, .1, .5, 1., 2., 5.])
# def test_sigma_def_massflow(m):
#     assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(mf.WeightMassFlow(1., gamma=1.3)/mf.WeightMassFlow(m, gamma=1.3))

def test_shockwave_Mn1_involutive_numpy():
    m = np.linspace(1., 10., 30)
    np.testing.assert_allclose(m, sw.downstream_Mn(sw.downstream_Mn(m)))