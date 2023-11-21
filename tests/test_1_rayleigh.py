import aerokit.aero.Rayleigh as ray

import pytest


def test_maxTiratio_Mach():
    assert ray.maxTiratio_Mach(0.1, 1.4) == pytest.approx(21.377994011976046)
    assert ray.maxTiratio_Mach(0.2, 1.4) == pytest.approx(5.761904761904761)
    assert ray.maxTiratio_Mach(1.2, 1.4) == pytest.approx(1.0217463193006673)
    assert ray.maxTiratio_Mach(20.0, 1.4) == pytest.approx(2.0236689814814817)
    assert ray.maxTiratio_Mach(0.5, 1.3) == pytest.approx(1.4464285714285716)
    assert ray.maxTiratio_Mach(2, 1.3) == pytest.approx(1.2604166666666667)


def test_SubMach_TiTicri():
    assert ray.SubMach_TiTicri(0.04677707) == pytest.approx(0.1)
    assert ray.SubMach_TiTicri(0.264462809, 1.4) == pytest.approx(0.2543448)
    assert ray.SubMach_TiTicri(0.04677707, 1.5) == pytest.approx(0.09799919)
    assert ray.SubMach_TiTicri(0.04677707, 0.5) == pytest.approx(0.12611333)


def test_SupMach_TiTicri():
    assert ray.SupMach_TiTicri(0.50701675) == pytest.approx(10)
    assert ray.SupMach_TiTicri(0.58733289, 1.4) == pytest.approx(4.03950381)
    assert ray.SupMach_TiTicri(0.50701675, 0.5) == pytest.approx(1.61955266)
    assert ray.SupMach_TiTicri(0.50701675, 1) == pytest.approx(2.39045722)


def test_Ps_Pscri():
    assert ray.Ps_Pscri(1, 1.4) == pytest.approx(1.0)
    assert ray.Ps_Pscri(1.5, 1.4) == pytest.approx(0.5783132530120482)
    assert ray.Ps_Pscri(0.5, 1.4) == pytest.approx(1.777778)
    assert ray.Ps_Pscri(1.5, 1.3) == pytest.approx(0.5859872611464967)
    assert ray.Ps_Pscri(0.5, 1.3) == pytest.approx(1.7358490566037734)


def test_Ts_Tscri():
    assert ray.Ts_Tscri(1, 1.4) == pytest.approx(1.0)
    assert ray.Ts_Tscri(1.5, 1.4) == pytest.approx(0.7525039918710987)
    assert ray.Ts_Tscri(0.5, 1.4) == pytest.approx(0.7901234567901234)
    assert ray.Ts_Tscri(1.5, 1.3) == pytest.approx(0.7726074080084382)
    assert ray.Ts_Tscri(0.5, 1.3) == pytest.approx(0.7532929868280526)


def test_Rho_Rhocri():
    assert ray.Rho_Rhocri(1, 1.4) == pytest.approx(1.0)
    assert ray.Rho_Rhocri(1.5, 1.4) == pytest.approx(0.7685185185185187)
    assert ray.Rho_Rhocri(0.5, 1.4) == pytest.approx(2.25)
    assert ray.Rho_Rhocri(1.5, 1.3) == pytest.approx(0.7584541062801933)
    assert ray.Rho_Rhocri(0.5, 1.3) == pytest.approx(2.3043478260869565)


def test_Pi_Picri():
    assert ray.Pi_Picri(1, 1.4) == pytest.approx(1.0)
    assert ray.Pi_Picri(1.5, 1.4) == pytest.approx(1.1215452274944158)
    assert ray.Pi_Picri(0.5, 1.4) == pytest.approx(1.114052503180089)
    assert ray.Pi_Picri(1.5, 1.3) == pytest.approx(1.1275538522966038)
    assert ray.Pi_Picri(0.5, 1.3) == pytest.approx(1.111142534976025)


def test_V_Vcri():
    assert ray.V_Vcri(1, 1.4) == pytest.approx(1.0)
    assert ray.V_Vcri(1.5, 1.4) == pytest.approx(1.3012048192771082)
    assert ray.V_Vcri(0.5, 1.4) == pytest.approx(0.44444444)
    assert ray.V_Vcri(1.5, 1.3) == pytest.approx(1.3184713375796178)
    assert ray.V_Vcri(0.5, 1.3) == pytest.approx(0.43396226415094336)


def test_NormdS():
    assert ray.NormdS(1, 1.4) == pytest.approx(0.0)
    assert ray.NormdS(1.5, 1.4) == pytest.approx(-0.12788052130716918)
    assert ray.NormdS(0.5, 1.4) == pytest.approx(-0.3999558269994989)
    assert ray.NormdS(1.5, 1.3) == pytest.approx(-0.13464795692852982)
    assert ray.NormdS(0.5, 1.3) == pytest.approx(-0.4105694949330352)
