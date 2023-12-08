import numpy as np
import numpy.linalg as nplin
import scipy.linalg as scilin
import aerokit.common.numspectral as ns
from typing import Any


class LinOperator:
    # this definition is related to arbitratry definition of omega: q(t) ~ exp(-j * omega * t)
    _harmonic_time_coef = -1j
    _BC_dict = {}

    def __init__(self, n=100, xmin=None, xmax=None) -> None:
        self._diffop = ns.ChebCollocation(n, xmin, xmax)
        self._basestate = None
        self._B = None
        self._BC_type = None

    @property
    def dim(self):
        return self._diffop.npts

    @property
    def x(self):
        return self._diffop.x

    def set_basestate(self, state: Any):
        """set base state as it will be used by derived class"""
        self._basestate = state

    def _check_BC(self, bc):
        """BC should be a dict with at least a 'type' key"""
        xbc = {"type": bc} if isinstance(bc, str) else bc
        if xbc["type"] not in self._BC_dict.keys():
            raise ValueError(f"{xbc['type']} key not found in available BC keys: {self._BC_dict.keys()}")
        return xbc

    def set_BC(self, Ltype: dict, Rtype: dict):
        self._BC_type = [self._check_BC(Ltype), self._check_BC(Rtype)]

    def check_basestate(self):
        return self._basestate is not None

    def compute_operators(self):
        raise NotImplementedError("this class must be overcharged and operators defined")

    def compute_BC(self):
        pass

    def get_RHS_time_matrices(self):
        if self._B is None:
            self.compute_operators()
        return self._B, self._At

    def solve_eig(self):
        B, Id = self.get_RHS_time_matrices()
        self._vals, self._vects = scilin.eig(B, self._harmonic_time_coef * Id)
        return self._vals, self._vects

    def select_and_sort(self, realmin=0.0, realmax=1.0e99, imagmin=-1e10, imagmax=1e10, sort="real"):
        """select and sort"""
        vp = self._vals
        condition = (vp.real < realmax) & (vp.real >= realmin) & (vp.imag > imagmin) & (vp.imag < imagmax)
        vals = np.compress(condition, self._vals)
        vects = np.compress(condition, self._vects, axis=1)
        np.set_printoptions(formatter={"float_kind": "{:.4f}".format})
        print(f"{vals.size}/{vp.size} remaining eigenvalues")
        sortmethod = {
            "real": np.real,
            "imag": lambda x: -np.imag(x),
        }
        order = np.argsort(sortmethod[sort](vals))
        return vals, vects, order

    def converge_eigenpair(self, omega, eigenmode, niter=10):
        """apply a Rayleigh quotient iteration convergence process to an estimate of eigenvalue, eigenmode
        of the (B- M*j*omega)*eigenmode = 0 stability problem

        Args:
            omega (_type_): _description_
            eigenmode (_type_): _description_
        """
        B, M = self.get_RHS_time_matrices()
        M = self._harmonic_time_coef * M  # must use explicit multiplication instead of *= (kind exception)
        new_eigm = eigenmode / nplin.norm(eigenmode)
        # print(omega, np.vdot(new_eigm, B @ new_eigm)/np.vdot(new_eigm, M @ new_eigm))
        new_omega = omega
        for i in range(niter):
            # # since matrix becomes singular with convergence, lstsqr is prefered
            # w, *_ = nplin.lstsq(B-new_omega*M, new_eigm, rcond=None)
            try:
                w = nplin.solve(B - new_omega * M, new_eigm)
                new_eigm = w / nplin.norm(w)
                new_omega = np.vdot(new_eigm, B @ new_eigm) / np.vdot(new_eigm, M @ new_eigm)
                # print(i, nplin.norm((B-new_omega*M)@ new_eigm), new_omega)
            except nplin.LinAlgError as err:
                if "Singular matrix" in str(err):
                    exit
                else:
                    raise
        return new_omega, new_eigm


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest

    doctest.testmod()
