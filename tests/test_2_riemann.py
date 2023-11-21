import aerokit.instance.riemann as riem
import aerokit.aero.unsteady1D as uq
import numpy as np
import pytest


def test_expansion_shock():
    qL = uq.unsteady_state(rho=2.0, u=0.0, p=10.0)
    qR = uq.unsteady_state(rho=1.0, u=0.0, p=1.0)
    pb = riem.riemann_pb(qL, qR)
    assert pb.left_fastest() == pytest.approx(
        qL.left_acoustic(), rel=1e-10
    )  # expansion has C- from qL
    assert pb.right_fastest() > qR.right_acoustic()  # shock wave is faster than qR C+


def test_expansion_shock_symmetry():
    qL = uq.unsteady_state(rho=2.0, u=0.0, p=10.0)
    qR = uq.unsteady_state(rho=1.0, u=0.0, p=1.0)
    pb1 = riem.riemann_pb(qL, qR)
    pb2 = riem.riemann_pb(qR, qL)
    assert pb2.pstar() == pytest.approx(pb1.pstar(), rel=1e-10)
    assert pb2.left_fastest() == pytest.approx(-pb1.right_fastest(), rel=1e-10)
    assert pb2.right_fastest() == pytest.approx(-pb1.left_fastest(), rel=1e-10)


def test_expansion_expansion():
    qL = uq.unsteady_state(rho=1.0, u=-0.3, p=1.2)
    qR = uq.unsteady_state(rho=1.0, u=0.2, p=1.0)
    pb = riem.riemann_pb(qL, qR)
    assert pb.left_fastest() == pytest.approx(
        qL.left_acoustic(), rel=1e-10
    )  # expansion has C- from qL
    assert pb.right_fastest() == pytest.approx(
        qR.right_acoustic(), rel=1e-10
    )  # expansion has C+ from qR


def test_expansion_expansion_shift():
    qL1 = uq.unsteady_state(rho=2.0, u=-1.2, p=1.2)
    qR1 = uq.unsteady_state(rho=1.0, u=0.3, p=1.3)
    pb1 = riem.riemann_pb(qL1, qR1)
    qL2 = uq.unsteady_state(rho=2.0, u=-0.2, p=1.2)
    qR2 = uq.unsteady_state(rho=1.0, u=1.3, p=1.3)
    pb2 = riem.riemann_pb(qL2, qR2)
    assert pb2.pstar() == pytest.approx(pb1.pstar(), rel=1e-10)
    assert pb2.left_fastest() - pb1.left_fastest() == pytest.approx(
        1.0, rel=1e-10
    )  # shift u = 1
    assert pb2.right_fastest() - pb1.right_fastest() == pytest.approx(1.0, rel=1e-10)
    assert pb2.ustar() - pb1.ustar() == pytest.approx(1.0, rel=1e-10)


def test_shock_shock():
    qL = uq.unsteady_state(rho=1.0, u=0.3, p=1.0)
    qR = uq.unsteady_state(rho=1.0, u=-0.2, p=1.1)
    pb = riem.riemann_pb(qL, qR)
    assert pb.left_fastest() < qL.left_acoustic()  # shock wave is faster than qL C-
    assert pb.right_fastest() > qR.right_acoustic()  # shock wave is faster than qR C+
    for Ws, q0, q1 in (
        (pb.left_fastest(), qL, pb.qstarL()),
        (pb.right_fastest(), qR, pb.qstarR()),
    ):
        # mass balance
        assert q1.massflow() - q0.massflow() == pytest.approx(Ws * (q1.rho - q0.rho))
        # momentum balance
        assert (
            q1.u * q1.massflow() - q0.u * q0.massflow() + q1.p - q0.p
            == pytest.approx(Ws * (q1.massflow() - q0.massflow()))
        )
