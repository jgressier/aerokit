#import unittest
import hades.common.defaultgas
import hades

def test_default_gamma():
	"""test default value of gamma in defaultgas module"""
	assert hades.common.defaultgas._gamma == 1.4
	
