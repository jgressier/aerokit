import aerokit.aero.unsteady1D as uq
import aerokit.aero.Isentropic as Is
import numpy as np
import pytest


def test_init_q():
    q = uq.unsteady_state(rho=1.4, u=2.0, p=1.0)
    assert q.asound() == pytest.approx(1.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)


def test_init_qgam():
    q = uq.unsteady_state(rho=1.3, u=2.0, p=1.0, gamma=1.3)
    assert q.asound() == pytest.approx(1.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)

