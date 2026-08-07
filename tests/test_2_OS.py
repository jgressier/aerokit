import aerokit.stability.OrrSommerfeld as OS
from aerokit.common.numspectral import MultiDomainChebCollocation
import numpy as np
import pytest


def test_initdefault():
    with pytest.raises(OS.DictKeyError):
        Op = OS.OrrSommerfeldModel(100, basestate={})
    with pytest.raises(OS.DictKeyError):
        Op = OS.OrrSommerfeldModel(100, basestate={'alpha':1.0})
    #
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

    def test_finite_eig(self):
        model = OS.Poiseuille(n=101, alpha=1.0, Reynolds=6000.0)
        vals, _ = model.solve_eig(finite_only=True)
        assert np.all(np.isfinite(vals))
        assert vals.size < model.dim

    def test_multidomain_connection(self):
        diffop = MultiDomainChebCollocation(31)
        model = OS.OrrSommerfeldModel(diffop.npts, basestate=None, diffop=diffop)
        model.set_basestate({"Reynolds": 1.0e4, "alpha": 1.0, "uprofile": lambda y: np.tanh(y / 0.1)})
        model.set_BC({"type": "wall"}, {"type": "wall"})
        model.compute_operators()
        interface = diffop.interface_index
        for irow, jump in zip(range(interface - 1, interface + 2), diffop.interface_derivative_jumps()):
            assert np.allclose(model._At[irow], 0.0)
            assert np.allclose(model._B[irow], jump)
