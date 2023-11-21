import aerokit.stability.OrrSommerfeld as OS
import numpy as np
import pytest


def test_initdefault():
    Op = OS.OrrSommerfeldModel(100)
    assert Op.x[0] == -1.0
    assert Op.x[-1] == 1.0
    assert Op.dim == 100
    assert Op._basestate is None


class Test_Poiseuille:
    def test_eig(self):
        model = OS.Poiseuille(n=101, alpha=1.0, Reynolds=6000.0)
        model.solve_eig()
        vals, _, order = model.select_and_sort(sort="imag")  # default is real order
        assert vals[order[0]] == pytest.approx(complex(0.2598, 0.00032), rel=1e-3)
        assert vals[order[1]] == pytest.approx(complex(0.95433, -0.04532), rel=1e-3)
