import aerokit.aero.model2D as M2P
import pytest


def test_init():
    Q = M2P.State2DMach(Mach=2., angle=20, rho=.1, p=10.)
    assert Q.Mach == 2.
    assert Q.angle == 20.
    assert Q.rho == .1
    assert Q.p == 10.
    Q.rotate(10)
    assert Q.angle == 30.

def test_init_default():
    Q = M2P.State2DMach(Mach=2.)
    assert Q.Mach == 2.
    assert Q.angle == 0.
    assert Q.rho == 1.
    assert Q.p == 1.

def test_normalshock():
    Q0 = M2P.State2DMach(Mach=2.)
    Q1 = Q0.shock_sigma(90.)
    assert Q1.angle == pytest.approx(0.0)
    assert Q1.Mach < 1.

def test_weakshock():
    Q0 = M2P.State2DMach(Mach=2.)
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

def test_weakshock_consistency():
    Q0 = M2P.State2DMach(Mach=2.)
    Q1p = Q0.weakshock_deviation(10.)
    Q1m = Q0.weakshock_deviation(-10., direction=-1)
    assert Q1p.angle == pytest.approx(-Q1m.angle)
    assert Q1p.Mach == pytest.approx(Q1m.Mach)
    assert Q1p.rho == pytest.approx(Q1m.rho)
    assert Q1p.p == pytest.approx(Q1m.p)
    #
    Q1 = Q0.shock_sigma(35.)
    Q2 = Q0.weakshock_deviation(Q1.angle)
    assert Q2.p == pytest.approx(Q1.p)

def test_stronghock():
    Q0 = M2P.State2DMach(Mach=2.)
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
