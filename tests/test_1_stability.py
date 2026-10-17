import aerokit.stability as stab
import numpy as np
import pytest


def test_base_initdefault():
    SpOp = stab.LinOperator()
    assert SpOp.x[0] == -1.0
    assert SpOp.x[-1] == 1.0
    assert SpOp.dim == 100
    assert SpOp._basestate is None


def test_base_init():
    n = 20
    u = 100.0
    SpOp = stab.LinOperator(n, 0, 10.0)
    assert SpOp.dim == n
    assert SpOp.x[0] == 0.0
    assert SpOp.x[-1] == 10.0
    SpOp.set_basestate(u)
    SpOp.check_basestate()


class Test_mat33:
    def init(self) -> None:
        n = 3
        self.Op = stab.LinOperator(n)
        # should be defined by a class "compute_operators"
        self.Op._B = np.array([1, 2, 3, 1, 2, 1, 3, 2, 1]).reshape((n, n))
        self.Op._At = np.eye(n) / self.Op._harmonic_time_coef

    def test_eig(self):
        self.init()
        B, _ = self.Op.get_RHS_time_matrices()  # check only RHS
        vals, _ = np.linalg.eig(B)
        assert np.max(np.abs(vals)) == pytest.approx(3 + np.sqrt(5))

    def test_rayleigh(self):
        self.init()
        vpmax, _ = self.Op.converge_eigenpair(200.0, [1.0, 1.0, 1.0])
        assert vpmax == pytest.approx(3 + np.sqrt(5))


def test_rayleigh_exact_eigenpair_exits_on_singular_correction():
    """An exact eigenpair makes the Rayleigh correction matrix singular."""
    operator = stab.LinOperator(2)
    omega_exact = 3.0
    operator._At = np.eye(2)
    operator._B = -1j * omega_exact * operator._At
    initial_mode = np.array([3.0, 4.0])

    omega, mode = operator.converge_eigenpair(omega_exact, initial_mode)

    assert omega == omega_exact
    assert np.allclose(mode, initial_mode / np.linalg.norm(initial_mode))
