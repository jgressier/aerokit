import aerokit.aero.MassFlow as mf
import numpy as np
import pytest


def test_sigma_scalar():
    assert mf.Sigma_Mach(0.5) == pytest.approx(1.33984375)
    assert mf.Sigma_Mach(1.0) == 1.0
    assert mf.Sigma_Mach(2.0) == pytest.approx(1.687500)


@pytest.mark.parametrize("m", [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0])
def test_sigma_def_massflow(m):
    assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(
        mf.WeightMassFlow(1.0, gamma=1.3) / mf.WeightMassFlow(m, gamma=1.3)
    )


def test_sigma_reverse_numpy():
    m = np.linspace(0.01, 2.0, 30)
    np.testing.assert_allclose(m, mf.Mach_Sigma(mf.Sigma_Mach(m), m))


@pytest.mark.parametrize("AsAc", [1.1, 2.0, 5.0, 10.0])
def test_MachSub_Sigma(AsAc):
    mach = mf.MachSub_Sigma(AsAc)
    assert mach < 1
    assert mf.Sigma_Mach(mach) == pytest.approx(AsAc, rel=1.0e-6)


@pytest.mark.parametrize("AsAc", [1.1, 2.0, 5.0, 10.0])
def test_MachSup_Sigma(AsAc):
    mach = mf.MachSup_Sigma(AsAc)
    assert mach > 1
    assert mf.Sigma_Mach(mach) == pytest.approx(AsAc, rel=1.0e-6)
