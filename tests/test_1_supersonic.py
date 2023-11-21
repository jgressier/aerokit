import aerokit.aero.Supersonic as sup
import numpy as np
import pytest


def test_PrandtlMeyer_Mach():
    assert sup.PrandtlMeyer_Mach(2.0) == pytest.approx(26.3797608)
    assert sup.PrandtlMeyer_Mach(2.0, 1.1) == pytest.approx(34.8813639)
    assert sup.PrandtlMeyer_Mach(1.5, 2) == pytest.approx(8.69455340)
    assert sup.PrandtlMeyer_Mach(1) == pytest.approx(0.0)


def test_Mach_PrandtlMeyer():
    assert sup.Mach_PrandtlMeyer(26.3797608) == pytest.approx(2.0)
    assert sup.Mach_PrandtlMeyer(34.8813639, 1.1) == pytest.approx(2.0)
    assert sup.Mach_PrandtlMeyer(8.69455340, 2) == pytest.approx(1.5)
    assert sup.Mach_PrandtlMeyer(0) == pytest.approx(1.0)


def test_old_Mach_PrandtlMeyer():
    assert sup.Mach_PrandtlMeyer(26) == pytest.approx(sup.old_Mach_PrandtlMeyer(26))
    assert sup.Mach_PrandtlMeyer(34, 1.1) == pytest.approx(
        sup.old_Mach_PrandtlMeyer(34, 1.1)
    )
    assert sup.Mach_PrandtlMeyer(8, 2) == pytest.approx(sup.old_Mach_PrandtlMeyer(8, 2))


def test_Mach_PrandtlMeyer_reverse():
    assert sup.PrandtlMeyer_Mach(sup.Mach_PrandtlMeyer(2)) == pytest.approx(
        sup.Mach_PrandtlMeyer(sup.PrandtlMeyer_Mach(2))
    )
    assert sup.PrandtlMeyer_Mach(sup.Mach_PrandtlMeyer(2.5)) == pytest.approx(
        sup.Mach_PrandtlMeyer(sup.PrandtlMeyer_Mach(2.5))
    )
    assert sup.PrandtlMeyer_Mach(sup.Mach_PrandtlMeyer(1)) == pytest.approx(
        sup.Mach_PrandtlMeyer(sup.PrandtlMeyer_Mach(1))
    )
    assert sup.PrandtlMeyer_Mach(sup.Mach_PrandtlMeyer(1.1, 1.4)) == pytest.approx(
        sup.Mach_PrandtlMeyer(sup.PrandtlMeyer_Mach(1.1, 1.4))
    )


def test_Mach_PMFmmu():
    assert sup.Mach_PMFmmu(1) == pytest.approx(2.1081987956)
    assert sup.Mach_PMFmmu(0) == pytest.approx(2.0841829663)
    assert sup.Mach_PMFmmu(-100) == pytest.approx(1.0151675434)
    assert sup.Mach_PMFmmu(100) == pytest.approx(11.155475552)
    assert sup.Mach_PMFmmu(35) == pytest.approx(3.1890623955)


def test_deflection_Mach_IsentropicPsratio():
    assert sup.deflection_Mach_IsentropicPsratio(4, 3) == pytest.approx(-12.029820531)
    assert sup.deflection_Mach_IsentropicPsratio(10, 5) == pytest.approx(-7.4094985182)
    assert sup.deflection_Mach_IsentropicPsratio(10, 1) == pytest.approx(0.0)


def test_IsentropicPsratio_Mach_deflection():
    assert sup.IsentropicPsratio_Mach_deflection(10, 1) == pytest.approx(1.27277769541)
    assert sup.IsentropicPsratio_Mach_deflection(12, 3) == pytest.approx(2.29313895245)
    assert sup.IsentropicPsratio_Mach_deflection(10, 10) == pytest.approx(8.09290197316)
