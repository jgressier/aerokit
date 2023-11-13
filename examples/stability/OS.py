import numpy as np
from scipy.linalg import eig
import aerokit.stability.OrrSommerfeld as OS

# --- MAIN ---

n = 101

model = OS.Poiseuille(n, 1., 6000.)
model.solve_eig()
vals0, vects0, order = model.select_and_sort(realmin=0., realmax=10., imagmin=-1., imagmax=5.)
print("first eigenvalues")
print(vals0[order][:20])

n = 2*(n-1)+1

model = OS.Poiseuille(n, 1., 6000.)
model.solve_eig()
vals, vects, order = model.select_and_sort(realmin=0., realmax=10., imagmin=-1., imagmax=5., sort='imag')
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
    ax[i,0].plot(model.x, model._diffop.matder(1) @ vects[:,order[i]].real)
    ax[i,1].plot(model.x, vects[:,order[i]].real)
plt.show()