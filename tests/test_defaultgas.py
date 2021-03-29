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
	
def test_save_restore():
	g1 = 1.3
	r1 = 300.
	g2 = 1.8
	r2 = 600.
	defgas.set_gamma(g1)
	defgas.set_r_ideal(r1)
	defgas.save_default()
	# properties are still there
	assert defgas.gamma() == g1
	assert defgas.r_ideal() == r1
	defgas.set_gamma(g2)
	defgas.set_r_ideal(r2)
	# properties have changed
	assert defgas.gamma() == g2
	assert defgas.r_ideal() == r2
	defgas.restore_default()
	# properties are back
	assert defgas.gamma() == g1
	assert defgas.r_ideal() == r1
