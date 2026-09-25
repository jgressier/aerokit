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


def test_nsaxi_requires_axis_on_radial_grid():
    n = 10
    state = {
        "kx": 1.0,
        "m": 0,
        "rho": np.ones(n),
        "P": np.ones(n),
        "Ux": np.zeros(n),
        "gamma": 1.4,
    }
    operator = NSaxi(n, rmin=0.1, rmax=1.0, basestate=state)

    with np.testing.assert_raises_regex(ValueError, "start at the axis"):
        operator.compute_operators()


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


def test_nsaxi_computes_baseflow_gradients_internally():
    n = 10
    rho = np.linspace(1.0, 2.0, n)
    ux = np.linspace(0.2, 0.8, n) ** 2
    state = {
        "kx": 1.0,
        "m": 0,
        "rho": rho,
        "P": np.ones(n),
        "Ux": ux,
        "gamma": 1.4,
    }
    operator = NSaxi(n, rmin=0.0, rmax=1.0, basestate=state)
    operator.compute_operators()

    derivative = operator._radial_matder(1)
    density_gradient = derivative @ rho
    velocity_gradient = derivative @ ux
    temperature_gradient = derivative @ (state["P"] / rho)
    corrected_radius = np.maximum(operator.r, operator.r[1] / 10.0)

    continuity_radial_block = operator._B0[:n, 2 * n : 3 * n]
    axial_radial_block = operator._B0[n : 2 * n, 2 * n : 3 * n]
    radial_density_block = operator._B0[2 * n : 3 * n, :n]
    radial_temperature_block = operator._B0[2 * n : 3 * n, 4 * n :]

    assert np.allclose(np.diag(continuity_radial_block), density_gradient + rho / corrected_radius)
    assert np.allclose(np.diag(axial_radial_block), velocity_gradient)
    assert np.allclose(np.diag(radial_density_block), temperature_gradient / rho)
    assert np.allclose(np.diag(radial_temperature_block), density_gradient / rho)
