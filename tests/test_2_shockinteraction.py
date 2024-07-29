import aerokit.instance.SWinteraction as SWI
import pytest


def test_wronginit():
    with pytest.raises(ValueError) as e_info:
        Pb = SWI.ShockInteraction(.5, 40., -45.)
    with pytest.raises(ValueError) as e_info:
        Pb = SWI.ShockInteraction(2., 10., -45.)
    with pytest.raises(ValueError) as e_info:
        Pb = SWI.ShockInteraction(2., 30., -15.)
    with pytest.raises(ValueError) as e_info:
        Pb = SWI.ShockInteraction(2., -40., 40.)

def test_init():
    Pb = SWI.ShockInteraction(2., 40., -45.)
    Pb.solve()
    assert Pb.solved
    assert Pb.check34balanced()

def test_weakweak():
    Pb = SWI.ShockInteraction(2., 35., -40.)
    Pb.solve()
    assert Pb.check34balanced()
    assert Pb[3].angle == pytest.approx(-4.8651971)
    Pb.plot_angle_pressure()
    SWI.plotsw.plt.show()

def test_weakstrong():
    Pb = SWI.ShockInteraction(4., 35., -90.)
    print(Pb.solve(verbose=True))
    assert Pb.check34balanced()
    assert Pb[3].angle == pytest.approx(2.3707989)
    # 2 and 4 must be the same states and subsonic
    assert Pb[2].Mach < 1
    assert Pb[2].Mach == Pb[4].Mach
    assert Pb[2].angle == Pb[4].angle
    assert Pb[2].p == Pb[4].p
    # 0 to 2 shock angle must be corrected
    assert Pb.sigma01 == 35.
    assert Pb.sigma02 == pytest.approx(89.33539)
    Pb.plot_angle_pressure()
    SWI.plotsw.plt.show()
