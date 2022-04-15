import aerokit.aero.Fanno as fan
import numpy as np
import pytest

def test_maxFparam_Mach():
    assert fan.maxFparam_Mach(1    ) == pytest.approx(0)
    assert fan.maxFparam_Mach(2,1.4) == pytest.approx(0.30499650258)
    assert fan.maxFparam_Mach(.1   ) == pytest.approx(66.9215600298)
    assert fan.maxFparam_Mach(1.5  ) == pytest.approx(0.13605021738)
    assert fan.maxFparam_Mach(2,2  ) == pytest.approx(0.14486038542)

def test_Ps_Pscri():
    assert fan.Ps_Pscri(1     ) == pytest.approx(1)
    assert fan.Ps_Pscri(2, 1.4) == pytest.approx(0.4082482904638631)
    assert fan.Ps_Pscri(.1    ) == pytest.approx(10.943513103291655)
    assert fan.Ps_Pscri(1.5   ) == pytest.approx(0.6064784348631228)
    assert fan.Ps_Pscri(2,2   ) == pytest.approx(0.3535533905932737)

def test_Ts_Tscri():
    assert fan.Ts_Tscri(1     ) == pytest.approx(1)
    assert fan.Ts_Tscri(2, 1.4) == pytest.approx(.6666667)
    assert fan.Ts_Tscri(.1    ) == pytest.approx(1.1976047904191616)
    assert fan.Ts_Tscri(1.5   ) == pytest.approx(0.8275862068965517)
    assert fan.Ts_Tscri(2,2   ) == pytest.approx(0.5)

def test_Rho_Rhocri():
    assert fan.Rho_Rhocri(1     ) == pytest.approx(1)
    assert fan.Rho_Rhocri(2, 1.4) == pytest.approx(0.6123724356957945)
    assert fan.Rho_Rhocri(.1    ) == pytest.approx(9.137833441248533)
    assert fan.Rho_Rhocri(1.5   ) == pytest.approx(0.7328281087929398)
    assert fan.Rho_Rhocri(2,2   ) == pytest.approx(0.7071067811865476)

def test_Pi_Picri():
    assert fan.Pi_Picri(1     ) == pytest.approx(1.)
    assert fan.Pi_Picri(2, 1.4) == pytest.approx(1.6875)
    assert fan.Pi_Picri(.1    ) == pytest.approx(5.82182875)
    assert fan.Pi_Picri(1.5   ) == pytest.approx(1.1761670524)
    assert fan.Pi_Picri(2,2   ) == pytest.approx(1.4142135624)

def test_V_Vcri():
    assert fan.V_Vcri(1     ) == pytest.approx(1.)
    assert fan.V_Vcri(2, 1.4) == pytest.approx(1.632993161855452)
    assert fan.V_Vcri(.1    ) == pytest.approx(0.109435131032916)
    assert fan.V_Vcri(1.5   ) == pytest.approx(1.364576478442026)
    assert fan.V_Vcri(2,2   ) == pytest.approx(1.414213562373095)

def test_NormdS():
    assert fan.NormdS(1     ) == pytest.approx(0)
    assert fan.NormdS(2, 1.4) == pytest.approx(-0.14949946964701377)
    assert fan.NormdS(.1    ) == pytest.approx(-0.5033184087429144)
    assert fan.NormdS(1.5   ) == pytest.approx(-0.046360254516405956)
    assert fan.NormdS(2,2   ) == pytest.approx(-0.17328679513998627)