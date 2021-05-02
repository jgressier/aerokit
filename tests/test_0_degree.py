import aerokit.aero.degree as deg
import numpy as np
import pytest

def test_cos():
	assert deg.cos(0.)  == pytest.approx(1.)
	assert deg.cos(90.) == pytest.approx(0.)
	
def test_sin():
	assert deg.cos(0.)  == pytest.approx(1.)
	assert deg.cos(90.) == pytest.approx(0.)
	
def test_tan():
	assert deg.tan(0.)  == pytest.approx(0.)
	assert deg.tan(45.) == pytest.approx(1.)
	
def test_acos():
	assert deg.acos(0.) == pytest.approx(90.)
	assert deg.acos(1.) == pytest.approx(0.)
	
def test_asin():
	assert deg.asin(0.) == pytest.approx(0.)
	assert deg.asin(1.) == pytest.approx(90.)

def test_atan():
	assert deg.atan(0.) == pytest.approx(0.)
	assert deg.atan(1.) == pytest.approx(45.)

def test_degree_numpy():
    a = np.linspace(0., 360, 50)
    np.testing.assert_allclose(deg.cos(a), deg.sin(a+90.), rtol=1.e-12)
