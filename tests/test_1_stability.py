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
