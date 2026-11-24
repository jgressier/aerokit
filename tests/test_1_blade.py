import pytest

from aerokit.aero.blade import BladeRow
from aerokit.aero.model2D import State2d


def test_absolute_relative_round_trip():
    row = BladeRow(Uy=80.0)
    absolute = State2d(rho=1.2, u=100.0, v=30.0, p=1.e5, gamma=1.4)
    relative = row.relative(absolute)
    assert relative.u == absolute.u
    assert relative.v == -50.0
    recovered = row.absolute(relative)
    assert recovered.v == absolute.v
    assert recovered.p == absolute.p


def test_rotating_definition_and_euler_work():
    row = BladeRow(omega=100.0, radius=0.4)
    assert row.Uy == 40.0
    inlet = State2d(1.0, 100.0, 10.0, 1.e5)
    outlet = State2d(1.0, 100.0, 30.0, 1.e5)
    assert row.specific_work(inlet, outlet) == 800.0


def test_vx_model_preserves_axial_velocity_and_relative_totals():
    row = BladeRow(Uy=60.0)
    inlet = State2d(1.2, 100.0, 20.0, 1.e5)
    outlet = row.downstream(inlet, beta=20.0)
    assert outlet.u == pytest.approx(inlet.u)
    relative_in = row.relative(inlet)
    relative_out = row.relative(outlet)
    assert relative_out.Ptot() == pytest.approx(relative_in.Ptot())
    assert relative_out.rTtot() == pytest.approx(relative_in.rTtot())


def test_compressible_massflow_with_area_change():
    row = BladeRow(Uy=40.0)
    inlet = State2d(1.2, 90.0, 20.0, 1.e5, gamma=1.33)
    outlet = row.downstream(inlet, beta=15.0, conservation="massflow", area_ratio=0.8)
    assert outlet.rho * outlet.u * 0.8 == pytest.approx(inlet.rho * inlet.u)
    relative_in = row.relative(inlet)
    relative_out = row.relative(outlet)
    assert relative_out.Ptot() == pytest.approx(relative_in.Ptot())
    assert relative_out.rTtot() == pytest.approx(relative_in.rTtot())


def test_relative_frame_interface():
    row = BladeRow(Uy=50.0)
    inlet = State2d(1.2, 100.0, -25.0, 1.e5)
    outlet = row.downstream(inlet, beta=10.0, frame="relative")
    assert outlet.angle() == pytest.approx(10.0)
