import numpy as np
from scipy.linalg import eig
from aerokit.common.defaultgas import set_gamma
from aerokit.aero.model1D import state
from aerokit.stability.Euler import Euler1D
import aerokit.aero.ShockWave as sw
import aerokit.aero.Isentropic as Is

gam = 1.4
set_gamma(gam)

print("> check linearized")
eps = .001
M0 = 2.
print("d(p1/p0)", 
    (sw.Ps_ratio(M0*(1+eps))-sw.Ps_ratio(M0))/(eps*M0),
    4*gam/(gam+1.)*M0
    )
print("d(rho1/rho0)", 
    (sw.Rho_ratio(M0*(1+eps))-sw.Rho_ratio(M0))/(eps*M0),
    4/(gam+1.)*sw.Rho_ratio(M0)**2 / M0**3
    )

class ShockModel(Euler1D):
    def __init__(self, n, length, M1=None, M0=None) -> None:
        super().__init__(n, xmin=0., xmax=length)
        if M1 is None:
            assert M0 is not None, "at least one Mach number should be given"
            M1 = sw.downstream_Mn(M0)
        print(f"Initialize Shock impact model with M1={M1:.4f}")
        # scale u1 so that du1/dx=-1
        rttot = length**2*(.5+1./(gam-1.)/M1**2)
        u = -Is.Velocity_MachTt(M1, rttot, r=1.)*(self.x-self.x[-1])/(self.x[0]-self.x[-1])
        Q = state(gamma=gam)
        Q.compute_from_pt_rtt_u(1e5, rttot, u)
        self.set_basestate(Q)
        self.compute_operators()
        self.setBC('RH', 'sym')
        B, Id = self.get_RHS_time_matrices()
        self._vals, self._vects = eig(B, self._harmonic_time_coef*Id)

    def get_select(self, realmin=0., realmax=1.e99, imagmin=-1e10, imagmax=1e10, sort='real'):
        """select and sort"""
        vp = self._vals
        condition = (vp.real < realmax) & (vp.real >= realmin) & (vp.imag > imagmin) & (vp.imag < imagmax)
        vals = np.compress(condition, self._vals) 
        vects = np.compress(condition, self._vects, axis=1) 
        np.set_printoptions(formatter={'float_kind':"{:.4f}".format})
        print(f"{vals.size}/{vp.size} remaining eigenvalues")
        sortmethod = { 
            'real': np.real, 
            'imag': lambda x: -np.imag(x),
            }
        order = np.argsort(sortmethod[sort](vals))
        return vals, vects, order

# --- MAIN ---

n = 101

model = ShockModel(n, M0=2.13, length=.56)
vals0, vects0, order = model.get_select(realmin=0., realmax=4., imagmin=-10., imagmax=100.)
print("first eigenvalues")
print(vals0[order][:20])

n = 2*(n-1)+1

model = ShockModel(n, M0=2.13, length=.56)
vals, vects, order = model.get_select(realmin=0., realmax=10., imagmin=-10., imagmax=-.15)
print("first eigenvalues")
print(vals[order][:20])

#--- PLOT ---
import matplotlib.pyplot as plt
plt.plot(vals0.real, vals0.imag, 'or', markersize=3, alpha=.5)
plt.plot(vals.real, vals.imag, 'ob', markersize=5, alpha=.5)
for i,v in enumerate(vals[order[:20]]):
    plt.text(v.real, v.imag, str(i))
plt.xlabel('$\omega_r$')
plt.ylabel('$\omega_i$')
plt.grid()
#plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
plt.title('Spectrum')
# plt.axis('square')
#plt.axis([-1000, 1000, -2, 0.5])
plt.show()
nv = 10
fig, ax = plt.subplots(nv,model.nvar)
for i in range(nv):
    for j in (0,2):
        print(j, model._diffop.matder(1)[n-1,:] @ vects[j*n:(j+1)*n,order[i]])
    for j in range(model.nvar):
        ax[i,j].plot(model._diffop.x, vects[j*n:(j+1)*n,order[i]].real)
plt.show()