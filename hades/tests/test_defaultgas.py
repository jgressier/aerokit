import unittest
#import common.defaultgas
import hades

class test_defaultgas(unittest.TestCase):

	def test_gamma(self):
		"""test default value of gamma in defaultgas module"""
		self.assertEqual(hades.common.defaultgas._gamma, 1.4)
	
