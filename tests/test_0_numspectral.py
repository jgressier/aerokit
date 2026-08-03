import aerokit.common.numspectral as ns
from aerokit.common.mapping import AffineMapping, SemiInfiniteAlgebraicMapping, SemiInfiniteRationalMapping
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


def test_affine_mapping_derivatives():
    mapping = AffineMapping(0.0, 10.0)
    assert mapping.dxi_dx() == -0.2
    assert mapping.d2xi_dx2() == 0.0
    assert mapping.d3xi_dx3() == 0.0
    assert mapping.d4xi_dx4() == 0.0


@pytest.mark.parametrize(
    "mapping, dxi_dx, d2xi_dx2, d3xi_dx3, d4xi_dx4",
    [
        (
            SemiInfiniteRationalMapping(2.0),
            lambda xi: (1.0 - xi) ** 2 / 4.0,
            lambda xi: -(1.0 - xi) ** 3 / 8.0,
            lambda xi: 3.0 * (1.0 - xi) ** 4 / 32.0,
            lambda xi: -3.0 * (1.0 - xi) ** 5 / 32.0,
        ),
        (
            SemiInfiniteAlgebraicMapping(2.0),
            lambda xi: (1.0 - xi) * np.sqrt(1.0 - xi**2) / 2.0,
            lambda xi: -(1.0 - xi) ** 2 * (2.0 * xi + 1.0) / 4.0,
            lambda xi: 3.0 * xi * (1.0 - xi) ** 2 * np.sqrt(1.0 - xi**2) / 4.0,
            lambda xi: 3.0 * (1.0 - xi) ** 3 * (1.0 - 2.0 * xi - 4.0 * xi**2) / 8.0,
        ),
    ],
)
def test_semi_infinite_mappings(mapping, dxi_dx, d2xi_dx2, d3xi_dx3, d4xi_dx4):
    xi = np.array([-0.8, -0.2, 0.5])
    assert np.allclose(mapping.x_to_xi(mapping.xi_to_x(xi)), xi)
    assert mapping.xi_to_x(-1.0) == 0.0
    assert np.isinf(mapping.xi_to_x(1.0))
    assert np.allclose(mapping.dxi_dx(xi), dxi_dx(xi))
    assert np.allclose(mapping.d2xi_dx2(xi), d2xi_dx2(xi))
    assert np.allclose(mapping.d3xi_dx3(xi), d3xi_dx3(xi))
    assert np.allclose(mapping.d4xi_dx4(xi), d4xi_dx4(xi))


def test_cheb_extrapol_exact():
    n = 5
    SpOp = ns.ChebCollocation(n)
    f = lambda x: (x-.5)**2
    x = np.linspace(-1., 1., 11, endpoint=True)
    fx = SpOp.extrapol(f(SpOp.x), x)
    assert np.allclose(fx, f(x), rtol=1e-12)


def test_cheb_extrapol_approx():
    n = 12
    SpOp = ns.ChebCollocation(n)
    f = lambda x: 1/(1.+x**2)
    x = np.linspace(-1., 1., 11, endpoint=True)
    fx = SpOp.extrapol(f(SpOp.x), x)
    assert np.allclose(fx, f(x), rtol=1e-4)


def test_cheb_fit_exact():
    n = 5
    SpOp = ns.ChebCollocation(n, 0., 10.)
    f = lambda x: (x-1.)*(x-3.)*(x-7.)
    x = np.linspace(0., 10., 51, endpoint=True)
    fxi = SpOp.fit_to_gauss(x, f(x))
    assert np.allclose(fxi, f(SpOp.x), rtol=1e-12)


def test_cheb_fit_approx():
    n = 20
    SpOp = ns.ChebCollocation(n, 0., 10.)
    f = lambda x: 2+np.sin(2*x)
    x = np.linspace(0., 10., 51, endpoint=True)
    fxi = SpOp.fit_to_gauss(x, f(x))
    #print(np.abs(fxi-f(SpOp.x)).max())
    assert np.allclose(fxi, f(SpOp.x), rtol=1e-4)


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


def test_cheb_diff_map_scaling():
    n = 20
    reference = ns.ChebCollocation(n)
    mapped = ns.ChebCollocation(n, xmin=-0.5, xmax=0.8)
    reference.compute_matder(4)
    mapped.compute_matder(4)
    scale = 2.0 / 1.3

    assert np.allclose(mapped._matder, reference._matder)
    for order in range(1, 5):
        assert np.allclose(mapped.matder(order), reference.matder(order) * scale**order)
