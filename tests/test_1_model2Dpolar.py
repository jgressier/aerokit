from aerokit.aero import degree
import aerokit.aero.model2Dpolar as M2P
import aerokit.aero.ShockWave as sw
import numpy as np
import pytest


def test_init():
    Q = M2P.state2Dpolar(Mach=2., angle=20, rho=.1, p=10.)
    assert Q.Mach == 2.
    assert Q.angle == 20.
    assert Q.rho == .1
    assert Q.p == 10.
    Q.rotate(10)
    assert Q.angle == 30.

def test_init_default():
    Q = M2P.state2Dpolar(Mach=2.)
    assert Q.Mach == 2.
    assert Q.angle == 0.
    assert Q.rho == 1.
    assert Q.p == 1.

def test_normalshock():
    Q0 = M2P.state2Dpolar(Mach=2.)
    Q1 = Q0.shock_sigma(90.)
    assert Q1.angle == 0.
    assert Q1.Mach < 1.

def test_weakshock():
    Q0 = M2P.state2Dpolar(Mach=2.)
    Q1p = Q0.shock_sigma(35.)
    assert Q1p.angle > 0.
    assert Q1p.Mach > 1.
    Q1m = Q0.shock_sigma(-35.)
    assert Q1m.angle == pytest.approx(-Q1p.angle)
    assert Q1m.Mach > 1.
    Q11p = Q0.weakshock_deviation(Q1p.angle)
    assert Q11p.p > 1
    Q11m = Q0.weakshock_deviation(Q1m.angle)
    assert Q11m.p > 1


