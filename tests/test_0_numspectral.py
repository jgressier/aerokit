import aerokit.common.numspectral as ns
import numpy as np
import pytest


def test_cheb_initdefault():
    n = 20
    SpOp = ns.ChebCollocation(n)
    assert SpOp.x[0] == -1.0
    assert SpOp.x[-1] == 1.0


def test_cheb_initx():
    n = 20
    SpOp = ns.ChebCollocation(n, 0, 10.0)
    assert SpOp.x[0] == 0.0
    assert SpOp.x[-1] == 10.0


def test_cheb_diff():
    def f(x):
        return np.exp(-((x - 0.1) ** 2))

    def df(x):
        return -2 * (x - 0.1) * f(x)

    def df2(x):
        return -2 * (1 - 2.0 * (x - 0.1) ** 2) * f(x)

    n = 20
    SpOp = ns.ChebCollocation(n)
    x = SpOp.x
    # check f' on x distribution (1 to -1.)
    df_th = df(x)
    df_num = SpOp.matder(1) @ f(x)
    assert np.sqrt(np.sum((df_th - df_num) ** 2) / n) < 1.0e-10
    # check f" on x distribution (1 to -1.)
    d2f_th = df2(x)
    d2f_num = SpOp.matder(1) @ df_num
    assert np.sqrt(np.sum((d2f_th - d2f_num) ** 2) / n) < 1.0e-8
    d2f_num = SpOp.matder(2) @ f(x)
    assert np.sqrt(np.sum((d2f_th - d2f_num) ** 2) / n) < 1.0e-8


def test_cheb_diff_map():
    def f(x):
        return np.exp(-((x - 0.1) ** 2))

    def df(x):
        return -2 * (x - 0.1) * f(x)

    def df2(x):
        return -2 * (1 - 2.0 * (x - 0.1) ** 2) * f(x)

    n = 20
    SpOp = ns.ChebCollocation(n, xmin=-0.5, xmax=0.8)
    x = SpOp.x
    # check f' on x distribution (1 to -1.)
    df_th = df(x)
    df_num = SpOp.matder(1) @ f(x)
    assert np.sqrt(np.sum((df_th - df_num) ** 2) / n) < 1.0e-10
    # check f" on x distribution (1 to -1.)
    d2f_th = df2(x)
    d2f_num = SpOp.matder(1) @ df_num
    assert np.sqrt(np.sum((d2f_th - d2f_num) ** 2) / n) < 1.0e-8
    d2f_num = SpOp.matder(2) @ f(x)
    assert np.sqrt(np.sum((d2f_th - d2f_num) ** 2) / n) < 1.0e-8
