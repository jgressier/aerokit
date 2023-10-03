import aerokit.common.numspectral as ns
import numpy as np
import pytest


def test_spectral_diff1_gauss():
	def f(x):
		return np.exp(-x**2)

	def df(x):
		return -2*x*np.exp(-x**2)

	n = 20
	SpOp = ns.ChebCollocation(n)
	x = SpOp.x
	df_th = df(x)
	df_num = SpOp.matder(1) @ f(x)
	assert np.sqrt(np.sum((df_th-df_num)**2)/n ) < 1.e-10
