import numpy as np
from scipy.linalg import eig, norm
from scipy.sparse.linalg import eigs
from aerokit.aero.model1D import state
from aerokit.stability.Euler import Euler1D

class Tube(Euler1D):
    def __init__(self, n, Mach) -> None:
        ones = np.ones(n)
        asound = 100.
        self.Q = state(rho = 1.4*ones, u = (Mach*asound)*ones, p = asound**2.*ones)
        super().__init__(n, xmin=0., xmax=1., basestate=self.Q)
        self.set_BC('per', 'per')


# --- MAIN ---
np.set_printoptions(formatter={'float_kind':"{:.4f}".format})

n = 101
model = Tube(n, Mach=.6)
model.solve_eig()
Aomega, Avects, order = model.select_and_sort(0.01, 5000., -1., 1.)

um, am, L = model.Q.u.mean(), model.Q.asound().mean(), model.x.max()-model.x.min()
for label, a  in zip(["u-a","u  ", "u+a"], [am-um, um, um+am]):
    print(f"theoretical eigenvalue (k=1) {label}: {2*np.pi*a/L:.2f}")

nsel = 20
print(f"{nsel} first eigenvalues")
omega0 = Aomega[order[:nsel]]
vects0 = Avects[:,order[:nsel]]
omega = np.empty_like(omega0)
vects = np.empty_like(vects0)
for i in range(nsel):
    omega[i], vects[:,i] = model.converge_eigenpair(omega0[i], vects0[:,i], niter=5)
    print(omega[i], np.abs(omega[i]-omega[0])/np.abs(omega[i]))
#print([norm(B @ v[:,i]-M @ v[:,i]) for i in order[:10]])

# --- PLOT ---
import matplotlib.pyplot as plt
plt.plot(Aomega.real, Aomega.imag, 'ob', markersize=10, alpha=.2)
plt.plot(omega.real, omega.imag, 'or', markersize=4)
plt.hlines(y=0.0, xmin=-2, xmax=2, color='k', linestyle='--')
plt.xlabel('$\omega_r$')
plt.ylabel('$\omega_i$')
#plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
plt.title('Spectrum')
# plt.axis('square')
#plt.axis([-1000, 1000, -2, 0.5])
plt.show()
icv = order[5]
#omega[icv], vects[:,icv] = model.converge_eigenpair(omega[icv], vects[:,icv], niter=50)
nvec = 8
fig, ax = plt.subplots(nvec, model.nvar)
for i in range(nvec):
    v0 = vects0[:,i]
    v1 = vects[:,i]
    #shift = np.vdot(v0,v1)/np.vdot(v0,v0)
    phase = np.average(np.mod(np.angle(v1)-np.angle(v0), 2*np.pi), weights=np.abs(v1))
    shift = np.exp(1j*phase)
    print(i, np.angle(shift, deg=True))
    for j in range(model.nvar):
        ax[i,j].plot(model.x, (shift*v0[j*n:(j+1)*n]).real, 'b')
        ax[i,j].plot(model.x, v1[j*n:(j+1)*n].real, 'r')
        #ax[i,j].plot(model.x, np.angle(v1[j*n:(j+1)*n])-np.angle(v0[j*n:(j+1)*n]), 'g')
        #ax[i,j].plot(model.x, np.angle(v0[j*n:(j+1)*n]), 'g', alpha=.5)
plt.show()
