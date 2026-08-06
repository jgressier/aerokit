import numpy as np
import pytest

from aerokit.blayer.profile import (
    IncompressibleProfile,
    Profile,
)


def test_linear_boundary_layer_thicknesses():
    y = np.linspace(0.0, 1.0, 10001)
    velocity = y
    profile = IncompressibleProfile(y, velocity)

    assert isinstance(profile, IncompressibleProfile)
    assert profile.y is y
    assert profile.velocity is velocity
    assert profile.displacement_thickness() == pytest.approx(0.5)
    assert profile.momentum_thickness() == pytest.approx(1.0 / 6.0, abs=1e-8)
    assert profile.thicknesses() == pytest.approx((0.5, 1.0 / 6.0), abs=1e-8)


def test_profile_requires_increasing_matching_coordinates():
    with pytest.raises(ValueError, match="same size"):
        IncompressibleProfile([0.0, 1.0], [0.0])
    with pytest.raises(ValueError, match="strictly increasing"):
        IncompressibleProfile([0.0, 0.0], [0.0, 1.0])


def test_generic_profile_accepts_additional_fields():
    profile = Profile([0.0, 1.0], density=[1.0, 0.9])
    np.testing.assert_allclose(profile.density, [1.0, 0.9])
