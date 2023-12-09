import aerokit.stability.Euler as Euler
from aerokit.aero.model1D import state
import numpy as np
import pytest


# def test_base_initdefault():
#     Op = Euler.LinOperator()
#     assert Op.x[0] == -1.
#     assert Op.x[-1] == 1.
#     assert Op.dim == 100
#     assert Op._basestate is None


# def test_base_init():
#     n = 20
#     u = 100.
#     Op = Euler.LinOperator(n, 0, 10.)
#     assert Op.dim == n
#     assert Op.x[0] == 0.
#     assert Op.x[-1] == 10.
#     Op.set_basestate(u)
#     Op.check_basestate()


class Test_ChannelPer:
    def init(self, n, Mach) -> Euler.Euler1D:
        ones = np.ones(n)
        asound = 100.0
        Q = state(rho=1.4 * ones, u=(Mach * asound) * ones, p=asound**2.0 * ones)
        Op = Euler.Euler1D(n, xmin=0.0, xmax=1.0, basestate=Q)
        Op.set_BC("per", "per")
        return Op

    @pytest.mark.parametrize("Mach", [-0.2, 0.3, 0.8])
    def test_Channel_Mach(self, Mach):
        Op = self.init(n=51, Mach=Mach)
        Op.solve_eig()
        Q = Op._basestate
        assert isinstance(Q, state)
        um, am, L = Q.u.mean(), Q.asound().mean(), Op.x.max() - Op.x.min()
        omega, _, order = Op.select_and_sort(0.01, 5000.0, -1.0, 1.0, sort="real")
        puls = []
        for u in [am - um, um, um + am]:
            puls += [(2 * np.pi * np.abs(u) / L * (i + 1)) for i in range(10)]
        puls = np.sort(puls)
        for i in range(10):  # check first 10 pulsations
            assert omega[order[i]].real == pytest.approx(puls[i])
            assert omega[order[i]].imag < 1.0e-7
