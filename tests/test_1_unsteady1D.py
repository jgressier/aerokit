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


def test_init_Pt_rTt_M():
    q = uq.unsteady_state(1.0, 0.0, 1.0)
    q.compute_from_pt_rtt_M(pt=4.0, rtt=2.0, M=-0.3)
    assert q.rTtot() == pytest.approx(2.0, rel=1e-10)
    assert q.Ptot() == pytest.approx(4.0, rel=1e-10)
    assert q.u == pytest.approx(-0.3 * q.asound(), rel=1e-10)


def test_init_Pt_rTt_u():
    q = uq.unsteady_state(1.0, 0.0, 1.0)
    q.compute_from_pt_rtt_u(pt=4.0, rtt=2.0, u=-0.3)
    assert q.rTtot() == pytest.approx(2.0, rel=1e-10)
    assert q.Ptot() == pytest.approx(4.0, rel=1e-10)


def test_init_Pt_rTt_p():
    q = uq.unsteady_state(1.0, 0.0, 1.0)
    q.compute_from_pt_rtt_p(pt=1.5, rtt=2.0, p=1.0)
    assert q.rTtot() == pytest.approx(2.0, rel=1e-10)
    assert q.Ptot() == pytest.approx(1.5, rel=1e-10)
