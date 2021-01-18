import hades.aero.MassFlow as mf
import numpy as np
import pytest

def test_sigma_scalar():
    assert mf.Sigma_Mach(.5) == pytest.approx(1.33984375)
    assert mf.Sigma_Mach(1.) == 1.
    assert mf.Sigma_Mach(2.) == pytest.approx(1.687500)

def test_sigma_numpy():
    assert 1