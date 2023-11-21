import aerokit.instance.nozzle as noz
import aerokit.aero.MassFlow as mf
import numpy as np
import pytest


def test_NPR_choked_subsonic():
    assert noz.NPR_choked_subsonic(1.8) == pytest.approx(1.085834)


def test_NPR_choked_supersonic():
    assert noz.NPR_choked_supersonic(1.8) == pytest.approx(8.81333)


def test_NPR_shock_at_exit():
    assert noz.NPR_shock_at_exit(1.8) == pytest.approx(1.812258)


@pytest.mark.parametrize(
    "npr, Ms", [(1.05, 0.2649307), (1.6, 0.51033), (7.0, 2.076365), (10.0, 2.076365)]
)
def test_Ms_from_AsAc_NPR(npr, Ms):
    assert noz.Ms_from_AsAc_NPR(1.8, npr) == pytest.approx(Ms)


@pytest.mark.parametrize(
    "npr, Ms", [(1.05, 0.2649307), (1.6, 0.51033), (7.0, 1.927457), (10.0, 2.1571946)]
)
def test_Madapt_from_AsAc_NPR(npr, Ms):
    assert noz.Madapt_from_AsAc_NPR(1.8, npr) == pytest.approx(Ms)


def test_nozzle_class():
    target_AoAc = 6.0
    length = 8.0
    Noz_x = np.linspace(0.0, length, 200, endpoint=True)
    ma_max = mf.Mach_Sigma(target_AoAc, Mach=2.0)
    ma = 1.0 + (ma_max - 1.0) * np.sin(0.5 * (Noz_x - 1.0) * np.pi / (length - 1.0))
    Noz_AoAc = mf.Sigma_Mach(ma)
    convdiv = noz.nozzle(Noz_x, Noz_AoAc, AsoAc=target_AoAc)
    assert convdiv.AsoAc == target_AoAc
    assert convdiv.ithroat == 25


def test_nozzle():
    target_AoAc = 6.0
    length = 8.0
    Noz_x = np.linspace(0.0, length, 200, endpoint=True)
    ma_max = mf.Mach_Sigma(target_AoAc, Mach=2.0)
    ma = 1.0 + (ma_max - 1.0) * np.sin(0.5 * (Noz_x - 1.0) * np.pi / (length - 1.0))
    Noz_AoAc = mf.Sigma_Mach(ma)
    convdiv = noz.nozzle(Noz_x, Noz_AoAc, AsoAc=target_AoAc, NPR=4.0)
    assert not np.any(np.isnan(convdiv.Mach()))
    assert not np.any(np.isnan(convdiv.Ps()))
    assert not np.any(np.isnan(convdiv.Ptot()))
