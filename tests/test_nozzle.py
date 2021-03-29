import aerokit.instance.nozzle as noz
import aerokit.aero.MassFlow   as mf
import numpy as np
import pytest

# def test_sigma_scalar():
#     assert mf.Sigma_Mach(.5) == pytest.approx(1.33984375)
#     assert mf.Sigma_Mach(1.) == 1.
#     assert mf.Sigma_Mach(2.) == pytest.approx(1.687500)


# @pytest.mark.parametrize("m", [.001, .01, .1, .5, 1., 2., 5.])
# def test_sigma_def_massflow(m):
#     assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(mf.WeightMassFlow(1., gamma=1.3)/mf.WeightMassFlow(m, gamma=1.3))

def test_nozzle_class():
    target_AoAc = 6.
    length      = 8.
    Noz_x    = np.linspace(0., length, 200, endpoint=True)
    ma_max   = mf.Mach_Sigma(target_AoAc, Mach=2.)
    ma       = 1. + (ma_max-1.)*np.sin(.5*(Noz_x-1.)*np.pi/(length-1.))
    Noz_AoAc = mf.Sigma_Mach(ma)
    convdiv = noz.nozzle(Noz_x, Noz_AoAc, AsoAc=target_AoAc)
    assert True
    #np.testing.assert_allclose(m, mf.Mach_Sigma(mf.Sigma_Mach(m), m))