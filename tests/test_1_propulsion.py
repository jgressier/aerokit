import pytest

from aerokit.aero import Propulsion as propulsion


@pytest.mark.parametrize("mach", [0.4, 2.0])
def test_thrust_function_inverse_on_selected_branch(mach):
    thrust = propulsion.ThrustFunction(mach)
    assert propulsion.Mach_ThrustFunction(thrust, Mach=mach) == pytest.approx(mach)
