import aerokit.instance.SWinteraction as SWI
import pytest


def test_wronginit():
    with pytest.raises(ValueError) as e_info:
        P1 = SWI.ShockInteraction(.5, 40., -45.)
    with pytest.raises(ValueError) as e_info:
        P1 = SWI.ShockInteraction(2., 10., -45.)
    with pytest.raises(ValueError) as e_info:
        P1 = SWI.ShockInteraction(2., 30., -15.)
    with pytest.raises(ValueError) as e_info:
        P1 = SWI.ShockInteraction(2., -40., 40.)

def test_init():
    P1 = SWI.ShockInteraction(2., 40., -45.)
    P1.solve()
    assert P1.solved
    assert abs(P1[3].p - P1[4].p) / P1[0].p < 1.e-6