#import unittest
import aerokit.common.defaultgas as defgas
#import aerokit

def test_default_gamma():
	"""test default value of gamma in defaultgas module"""
	assert defgas._gamma == 1.4
	
def test_gamma():
	"""test set value of gamma in defaultgas module"""
	defgas.set_gamma(1.3)
	assert defgas.gamma() == 1.3
	
def test_r_ideal():
	"""test set value of r in defaultgas module"""
	defgas.set_r_ideal(518.)
	assert defgas.r_ideal() == 518.
	
