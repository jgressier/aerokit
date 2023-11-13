import numpy as np
from scipy.linalg import eig
from aerokit.stability.OrrSommerfeldt import OrrSommerfeldModel


class Poiseuille(OrrSommerfeldModel):
    def __init__(self, n, alpha, Reynolds) -> None:
        super().__init__(n, xmin=-1., xmax=1.)
        self.set_basestate({
            'Reynolds' : Reynolds,
            'alpha': alpha,
            'uprofile': lambda x: 1-x**2
        })
        self.compute_operators()
        self.setBC('wall', 'wall')
        B, Id = self.get_RHS_time_matrices()
        self._vals, self._vects = eig(B, -self._harmonic_time_coef*Id)

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

model = Poiseuille(n, 1., 6000.)
vals0, vects0, order = model.get_select(realmin=0., realmax=10., imagmin=-1., imagmax=5.)
print("first eigenvalues")
print(vals0[order][:20])

n = 2*(n-1)+1

model = Poiseuille(n, 1., 6000.)
vals, vects, order = model.get_select(realmin=0., realmax=10., imagmin=-1., imagmax=5., sort='imag')
print("first eigenvalues")
print(vals[order][:20])

#--- PLOT ---
import matplotlib.pyplot as plt
plt.plot(vals0.real, vals0.imag, 'or', markersize=3, alpha=.9)
plt.plot(vals.real, vals.imag, 'ob', markersize=5, alpha=.3)
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
nv = 6
fig, ax = plt.subplots(nv,2)
for i in range(nv):
    ax[i,1].plot(model.x, vects[:,order[i]].real)
plt.show()