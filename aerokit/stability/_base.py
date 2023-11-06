import numpy as np
import numpy.linalg as nplin
import aerokit.common.numspectral as ns

class LinOperator():

    # this definition is related to arbitratry definition of omega q(t) ~ exp(-j * omega * t)
    _harmonic_time_coef = -1j

    def __init__(self, n=100, xmin=None, xmax=None) -> None:
        self._diffop = ns.ChebCollocation(n, xmin, xmax)

    @property
    def dim(self):
        return self._diffop.npts
    
    @property
    def x(self):
        return self._diffop.x
    
    def set_basestate(self, state):
        """set base state as it will be used by derived class"""
        self._basestate = state

    def get_RHS_time_matrices(self):
        return self._B, self._At
    
    def converge_eigenpair(self, omega, eigenmode, niter=10):
        """apply a Rayleigh quotient iteration convergence process to an estimate of eigenvalue, eigenmode
        of the (B- M*i*omega)*eigenmode = 0 stability problem

        Args:
            omega (_type_): _description_
            eigenmode (_type_): _description_
        """
        B, M = self.get_RHS_time_matrices()
        new_eigm = eigenmode/nplin.norm(eigenmode)
        print(omega, (new_eigm.T @ (B @ new_eigm))/self._harmonic_time_coef)
        new_omega = omega
        for i in range(niter):
            w = nplin.solve(B-(self._harmonic_time_coef*new_omega)*M, new_eigm)
            new_eigm = w / nplin.norm(w)
            new_omega = np.vdot(new_eigm.T, B @ new_eigm)/self._harmonic_time_coef
            print(i, new_omega)
        return new_omega, new_eigm


# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()