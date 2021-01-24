import aerokit.aero.MassFlow as mf
import numpy as np
import pytest

def test_sigma_scalar():
    assert mf.Sigma_Mach(.5) == pytest.approx(1.33984375)
    assert mf.Sigma_Mach(1.) == 1.
    assert mf.Sigma_Mach(2.) == pytest.approx(1.687500)


@pytest.mark.parametrize("m", [.001, .01, .1, .5, 1., 2., 5.])
def test_sigma_def_massflow(m):
    assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(mf.WeightMassFlow(1., gamma=1.3)/mf.WeightMassFlow(m, gamma=1.3))

def test_sigma_reverse_numpy():
    m = np.linspace(.01, 2., 30)
    np.testing.assert_allclose(m, mf.Mach_Sigma(mf.Sigma_Mach(m), m))