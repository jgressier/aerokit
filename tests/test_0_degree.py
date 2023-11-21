import aerokit.aero.degree as deg
import numpy as np
import pytest


def test_cos():
    assert deg.cos(0.0) == pytest.approx(1.0)
    assert deg.cos(90.0) == pytest.approx(0.0)


def test_sin():
    assert deg.cos(0.0) == pytest.approx(1.0)
    assert deg.cos(90.0) == pytest.approx(0.0)


def test_tan():
    assert deg.tan(0.0) == pytest.approx(0.0)
    assert deg.tan(45.0) == pytest.approx(1.0)


def test_acos():
    assert deg.acos(0.0) == pytest.approx(90.0)
    assert deg.acos(1.0) == pytest.approx(0.0)


def test_asin():
    assert deg.asin(0.0) == pytest.approx(0.0)
    assert deg.asin(1.0) == pytest.approx(90.0)


def test_atan():
    assert deg.atan(0.0) == pytest.approx(0.0)
    assert deg.atan(1.0) == pytest.approx(45.0)


def test_degree_numpy():
    a = np.linspace(0.0, 360, 50)
    np.testing.assert_allclose(deg.cos(a), deg.sin(a + 90.0), rtol=1.0e-12)
