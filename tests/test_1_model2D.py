import aerokit.aero.model2D as m2d
import numpy as np
import pytest


def test_init_q():
    q = m2d.State2d(rho=1.4, u=20.0, v=0., p=100.0)
    assert q.asound() == pytest.approx(10.0, rel=1e-10)
    assert q.Mach() == pytest.approx(2.0, rel=1e-10)
    c = q.copy()
    assert q.u == c.u
    assert q.v == c.v
    assert q.p == c.p
    c.v += 20.
    assert c.angle() == pytest.approx(45., rel=1e-10)


def test_init_qarray():
    n = 10
    u = np.linspace(0., 20., n)
    q = m2d.State2d(rho=1.4, u=u, v=0., p=100.0)
    assert q.size == n
    assert q.asound() == pytest.approx(10.0, rel=1e-10) #
    assert np.allclose(q.Mach(), u/10., atol=1e-10)


def test_mach_angle_conversion():
    q = m2d.State2d(rho=1.4, u=20.0, v=20.0, p=100.0)
    qma = q.to_mach_angle()
    assert qma.Mach == pytest.approx(2.0 * np.sqrt(2.0))
    assert qma.angle == pytest.approx(45.0)
    q_components = qma.to_components()
    assert q_components.u == pytest.approx(q.u)
    assert q_components.v == pytest.approx(q.v)
    assert q_components.rho == q.rho
    assert q_components.p == q.p
