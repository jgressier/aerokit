import numpy as np
from scipy.linalg import eig
from aerokit.aero.model1D import state
from aerokit.stability.Euler import Euler1D

n = 101
ones = np.ones(n)
Q = state(rho = 1.4*ones, u = 20*ones, p = 10000.*ones)

model = Euler1D(n, Q)
model.compute_operators()
model.setBC_per()

# debug
#np.set_printoptions(precision=2)
# np.set_printoptions(formatter={'float_kind':"{:.2f}".format})
# print(model._diffop.matder(1))
# print(model._B)
#print(model._B0)

Id = 1j*model._At
B = model._B
l, v = eig(B, Id)

# --- TRI ---
condition = np.logical_and(l.real >= 0., np.abs(l.imag) < 10.)
condition = (l.real >0) & (np.abs(l) >= 10) & (np.abs(l.real) < 5000.) & (np.abs(l.imag) < 1.)
#condition = np.logical_and(l.real >= 0., l.real < 500.)
l = np.compress(condition, l) 
v = np.compress(condition, v, axis=1) 
np.set_printoptions(formatter={'float_kind':"{:.4f}".format})
print(f"{l.size} remaining vp")
#order = np.argsort(-l.imag)
order = np.argsort(np.abs(l.real))
#order = np.argsort(np.abs(l))
#order = np.argsort(np.std(v, axis=0))
um, am = Q.u.mean(), Q.asound().mean()
for label, a  in zip(["u-a","u  ", "u+a"], [am-um, um, um+am]):
    print(f"theoretical eigenvalue (k=1) {label}: {2*np.pi*a/2.}")
print("first eigenvalues")
print(l[order][:10])

# --- PLOT ---
import matplotlib.pyplot as plt
plt.plot(l.real, l.imag, 'ob', markersize=6)
plt.hlines(y=0.0, xmin=-2, xmax=2, color='k', linestyle='--')
plt.xlabel('$\omega_r$')
plt.ylabel('$\omega_i$')
#plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
plt.title('Spectrum')
# plt.axis('square')
#plt.axis([-1000, 1000, -2, 0.5])
plt.show()
nv = 8
fig, ax = plt.subplots(nv,model.nvar)
for i in range(nv):
    for j in range(model.nvar):
        ax[i,j].plot(model._diffop.x, v[j*n:(j+1)*n,order[i]].real)
plt.show()