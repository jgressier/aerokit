#import unittest
import aerokit.common.defaultgas
import hades

def test_default_gamma():
	"""test default value of gamma in defaultgas module"""
	assert aerokit.common.defaultgas._gamma == 1.4
	
