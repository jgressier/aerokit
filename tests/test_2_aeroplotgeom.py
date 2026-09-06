import matplotlib.pyplot as plt
import numpy as np
import pytest

from aerokit.aero.plot.geom import Geom, Wall


def test_geom_draws_multiple_boundaries_and_solid_side():
    geometry = Geom()
    geometry.add_wall((-1.0, 2.0), (1.0, 1.0))
    lower_wall = geometry.add_wall((-1.0, 0.0, 2.0), (0.0, 0.0, 1.0), location="bottom")
    fig, axis = geometry.subplots(xlim=(-1.0, 2.0), ylim=(-0.2, 1.2))

    lines = geometry.plot()

    assert len(lines) == 2
    assert len(axis.patches) == 1
    assert axis.patches[0].get_hatch() == "///"
    assert axis.patches[0].get_facecolor()[-1] == pytest.approx(0.55)
    assert lower_wall.angle(1.5) == pytest.approx(np.degrees(np.arctan(0.5)))
    plt.close(fig)


def test_wall_supports_callable_coordinates():
    wall = Wall((-1.0, 1.0), lambda x: x**2, npts=5)
    np.testing.assert_allclose(wall.x, np.linspace(-1.0, 1.0, 5))
    np.testing.assert_allclose(wall.y, wall.x**2)

    vertical_wall = Wall(lambda y: y**2, (0.0, 1.0), npts=5)
    np.testing.assert_allclose(vertical_wall.x, vertical_wall.y**2)


def test_wall_rejects_ambiguous_or_invalid_definitions():
    with pytest.raises(ValueError, match="only one"):
        Wall(lambda x: x, lambda y: y)
    with pytest.raises(ValueError, match="location"):
        Wall((0.0, 1.0), (0.0, 1.0), location="inside")


def test_wall_accepts_fill_style_overrides():
    fig, axis = plt.subplots()
    axis.set(xlim=(0.0, 1.0), ylim=(-1.0, 1.0))
    wall = Wall((0.0, 1.0), (0.0, 0.0), location="bottom",
                fill_style={"facecolor": "red", "hatch": None})
    wall.plot(ax=axis)

    assert axis.patches[0].get_hatch() is None
    np.testing.assert_allclose(axis.patches[0].get_facecolor()[:3], (1.0, 0.0, 0.0))
    plt.close(fig)
