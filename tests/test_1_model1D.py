import aerokit.aero.model1D as m1d
import numpy as np
import pytest


# def test_init_default():
#     q = m1d.State1d()
#     #assert q.asound() == pytest.approx(1.0, rel=1e-10)
#     #assert q.Mach() == pytest.approx(2.0, rel=1e-10)


def test_init_q():
    q = m1d.State1d(rho=1.4, u=2.0, p=1.0)
    assert q.asound() == pytest.approx(1.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)
    c = q.copy()
    assert q.u == c.u
    c.u += 2.
    assert q.u+2. == c.u


def test_init_q_backwardcomp():
    q = m1d.state(rho=1.4, u=2.0, p=1.0)
    assert q.asound() == pytest.approx(1.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)
    c = q.copy()
    assert q.u == c.u
    c.u += 2.
    assert q.u+2. == c.u


def test_init_qgam():
    q = m1d.State1d(rho=1.3, u=2.0, p=1.0, gamma=1.3)
    assert q.asound() == pytest.approx(1.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)


def test_item_requires_array_state():
    with pytest.raises(TypeError, match="numpy ndarray"):
        m1d.State1d(rho=1.0, u=0.0, p=1.0)[0]


def test_init_Mneg():
    q = m1d.State1d(rho=1.3, u=-20.0, p=100.0, gamma=1.3)
    assert q.asound() == pytest.approx(10.0, rel=1e-10)
    assert q.Mach() == pytest.approx(-2.0, rel=1e-10)


def test_init_Pt_rTt_M():
    n = 10
    q = m1d.State1d(1.0, 0.0, 1.0)
    m = np.linspace(-.4, 2.5, n)
    q.compute_from_pt_rtt_M(pt=4.0, rtt=2.0, M=m)
    assert q.size == n
    assert np.allclose(q.rTtot(), 2.0, rtol=1e-10)
    assert np.allclose(q.Ptot(), 4.0, rtol=1e-10)
    assert np.allclose(q.u,  m * q.asound(), rtol=1e-10)


def test_init_Pt_rTt_u():
    q = m1d.State1d(1.0, 0.0, 1.0)
    q.compute_from_pt_rtt_u(pt=4.0, rtt=2.0, u=-0.3)
    assert q.rTtot() == pytest.approx(2.0, rel=1e-10)
    assert q.Ptot() == pytest.approx(4.0, rel=1e-10)


def test_init_Pt_rTt_p():
    q = m1d.State1d(1.0, 0.0, 1.0)
    q.compute_from_pt_rtt_p(pt=1.5, rtt=2.0, p=1.0)
    assert q.rTtot() == pytest.approx(2.0, rel=1e-10)
    assert q.Ptot() == pytest.approx(1.5, rel=1e-10)


def test_new_isentropic_trans():
    n = 10
    q1 = m1d.State1d(1.0, 0.0, 100.0)
    m = np.linspace(-.4, 2.5, n)
    q2 = q1.state_isentropic_Mach(m)
    assert q2.size == n
    assert np.allclose(q1.rTtot(), q2.rTtot(), rtol=1e-10)
    assert np.allclose(q1.Ptot(), q2.Ptot(), rtol=1e-10)


def test_new_RH():
    q1 = m1d.State1d(1.4, 150., 1e4) # Mach 1.5
    q2 = q1.state_RH()
    q3 = q2.state_RH() # should be the same as q1
    assert q1.Mach() > 1.
    assert q2.Mach() < 1.
    assert q1.rho == pytest.approx(q3.rho, rel=1e-10)
    assert q1.u == pytest.approx(q3.u, rel=1e-10)
    assert q1.p == pytest.approx(q3.p, rel=1e-10)

