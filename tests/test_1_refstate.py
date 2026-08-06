import pytest

from aerokit.aero.refstate import refstate


def test_reynolds_number_from_reference_quantities():
    state = refstate(density=1.2, velocity=50.0, viscosity=1.5e-5, length=0.4)
    assert state.Reynolds == pytest.approx(1.6e6)
    assert state.reynolds == pytest.approx(1.6e6)


def test_reynolds_requires_nonzero_viscosity():
    with pytest.raises(ValueError, match="viscosity"):
        _ = refstate(viscosity=0.0).Reynolds
