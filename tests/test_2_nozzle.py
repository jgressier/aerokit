import aerokit.instance.nozzle as noz
import aerokit.aero.MassFlow   as mf
import numpy as np
import pytest

# def test_sigma_scalar():
#     assert mf.Sigma_Mach(.5) == pytest.approx(1.33984375)
#     assert mf.Sigma_Mach(1.) == 1.
#     assert mf.Sigma_Mach(2.) == pytest.approx(1.687500)


# @pytest.mark.parametrize("m", [.001, .01, .1, .5, 1., 2., 5.])
# def test_sigma_def_massflow(m):
#     assert mf.Sigma_Mach(m, gamma=1.3) == pytest.approx(mf.WeightMassFlow(1., gamma=1.3)/mf.WeightMassFlow(m, gamma=1.3))

def test_NPR_choked_subsonic():
    assert noz.NPR_choked_subsonic(1.8) == pytest.approx(1.085834)

def test_NPR_choked_supersonic():
    assert noz.NPR_choked_supersonic(1.8) == pytest.approx(8.81333)

def test_NPR_shock_at_exit():
    assert noz.NPR_shock_at_exit(1.8) == pytest.approx(1.812258)

@pytest.mark.parametrize("npr, Ms", [(1.05, .2649307), (1.6, .51033), (7., 2.076365), (10., 2.076365)])
def test_Ms_from_AsAc_NPR(npr, Ms):
    assert noz.Ms_from_AsAc_NPR(1.8, npr) == pytest.approx(Ms)

@pytest.mark.parametrize("npr, Ms", [(1.05, .2649307), (1.6, .51033), (7., 1.927457), (10., 2.1571946)])
def test_Madapt_from_AsAc_NPR(npr, Ms):
    assert noz.Madapt_from_AsAc_NPR(1.8, npr) == pytest.approx(Ms)

def test_nozzle_class():
    target_AoAc = 6.
    length      = 8.
    Noz_x    = np.linspace(0., length, 200, endpoint=True)
    ma_max   = mf.Mach_Sigma(target_AoAc, Mach=2.)
    ma       = 1. + (ma_max-1.)*np.sin(.5*(Noz_x-1.)*np.pi/(length-1.))
    Noz_AoAc = mf.Sigma_Mach(ma)
    convdiv = noz.nozzle(Noz_x, Noz_AoAc, AsoAc=target_AoAc)
    assert convdiv.AsoAc == target_AoAc
    assert convdiv.ithroat == 25

def test_nozzle():
    target_AoAc = 6.
    length      = 8.
    Noz_x    = np.linspace(0., length, 200, endpoint=True)
    ma_max   = mf.Mach_Sigma(target_AoAc, Mach=2.)
    ma       = 1. + (ma_max-1.)*np.sin(.5*(Noz_x-1.)*np.pi/(length-1.))
    Noz_AoAc = mf.Sigma_Mach(ma)
    convdiv = noz.nozzle(Noz_x, Noz_AoAc, AsoAc=target_AoAc, NPR=4.)
    assert not np.any(np.isnan(convdiv.Mach()))
    assert not np.any(np.isnan(convdiv.Ps()))
    assert not np.any(np.isnan(convdiv.Ptot()))
