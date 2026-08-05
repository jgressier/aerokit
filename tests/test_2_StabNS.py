import numpy as np

from aerokit.common.mapping import SemiInfiniteAlgebraicMapping
from aerokit.stability.NS import NSaxi


def test_nsaxi_semi_infinite_mapping():
    operator = NSaxi(20, basestate=None, mapping=SemiInfiniteAlgebraicMapping(10.0))
    assert operator.r[0] == 0.0
    assert np.isinf(operator.r[-1])
    assert np.all(np.isfinite(operator._radial_matder(1)))


def _operator_with_axis_bc(m):
    n = 10
    state = {
        "kx": 1.0,
        "m": m,
        "rho": np.ones(n),
        "P": np.ones(n),
        "Ux": np.zeros(n),
        "gamma": 1.4,
    }
    operator = NSaxi(n, rmin=0.0, rmax=1.0, basestate=state)
    operator.compute_operators()
    return operator


def test_nsaxi_axis_bc_m0():
    operator = _operator_with_axis_bc(0)
    n = operator.dim
    D = operator._radial_matder(1)

    assert np.allclose(operator._B[0, :n], D[0, :])
    assert np.allclose(operator._B[n, n : 2 * n], D[0, :])
    assert operator._B[2 * n, 2 * n] == 1.0
    assert operator._B[3 * n, 3 * n] == 1.0
    assert np.allclose(operator._B[4 * n, 4 * n :], D[0, :])


def test_nsaxi_axis_bc_m1():
    operator = _operator_with_axis_bc(1)
    n = operator.dim
    D = operator._radial_matder(1)

    assert operator._B[0, 0] == 1.0
    assert operator._B[n, n] == 1.0
    assert operator._B[4 * n, 4 * n] == 1.0
    assert operator._B[2 * n, 2 * n] == -1j
    assert operator._B[2 * n, 3 * n] == 1.0
    assert np.allclose(operator._B[3 * n, 2 * n : 3 * n], 1j * D[0, :])
    assert np.allclose(operator._B[3 * n, 3 * n : 4 * n], D[0, :])
