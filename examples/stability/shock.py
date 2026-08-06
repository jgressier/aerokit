import numpy as np
from aerokit.common.defaultgas import set_gamma
from aerokit.aero.model1D import state
from aerokit.stability.Euler import Euler1D
import aerokit.aero.ShockWave as sw
import aerokit.aero.Isentropic as Is

gam = 1.4
set_gamma(gam)

print("> check linearized")
eps = 0.001
M0 = 2.0
print(
    "d(p1/p0)",
    (sw.Ps_ratio(M0 * (1 + eps)) - sw.Ps_ratio(M0)) / (eps * M0),
    4 * gam / (gam + 1.0) * M0,
)
print(
    "d(rho1/rho0)",
    (sw.Rho_ratio(M0 * (1 + eps)) - sw.Rho_ratio(M0)) / (eps * M0),
    4 / (gam + 1.0) * sw.Rho_ratio(M0) ** 2 / M0**3,
)


class ShockModel(Euler1D):
    def __init__(self, n, length, M1=None, M0=None, dM0dx=None) -> None:
        super().__init__(n, xmin=0.0, xmax=length)
        if M1 is None:
            assert M0 is not None, "at least one Mach number should be given"
            M1 = sw.downstream_Mn(M0)
        print(f"Initialize Shock impact model with M1={M1:.4f}")
        rttot = 8.35e4
        ptot = 1.8e5
        u = Is.Velocity_MachTt(M1, rttot, r=1.) * (self.x - self.x[-1]) / (self.x[0] - self.x[-1])
        print(f"u1: {u[0]}")
        Q = state(rho=1., u=0., p=1., gamma=gam)
        Q.compute_from_pt_rtt_u(ptot, rttot, u)
        self.set_basestate(Q)
        self.set_BC(
            {'type': 'RH', 'dM0dx': dM0dx}, 
            {'type': 'sym'}
        )

    # def get_select(
    #     self, realmin=0.0, realmax=1.0e99, imagmin=-1e10, imagmax=1e10, sort="real"
    # ):
    #     """select and sort"""
    #     vp = self._vals
    #     condition = (
    #         (vp.real < realmax)
    #         & (vp.real >= realmin)
    #         & (vp.imag > imagmin)
    #         & (vp.imag < imagmax)
    #     )
    #     vals = np.compress(condition, self._vals)
    #     vects = np.compress(condition, self._vects, axis=1)
    #     np.set_printoptions(formatter={"float_kind": "{:.4f}".format})
    #     print(f"{vals.size}/{vp.size} remaining eigenvalues")
    #     sortmethod = {
    #         "real": np.real,
    #         "imag": lambda x: -np.imag(x),
    #     }
    #     order = np.argsort(sortmethod[sort](vals))
    #     return vals, vects, order


# --- MAIN ---

M0 = 2.495
length=2.28e-3
dM0dx=850.

n = 201

model = ShockModel(n, M0=M0, length=length, dM0dx=dM0dx)

import matplotlib.pyplot as plt
if False:
    fig, ax = plt.subplots(4,1)
    Q: state = model._basestate
    for iax, ival in zip(ax, (Q.rho, Q.u, Q.p, Q.Mach())):
        iax.plot(model.x, ival)
    plt.show()

model.solve_eig()
vals0, vects0, order = model.select_and_sort(
    realmin=0., realmax=1200e3,
    imagmin=-1e5, imagmax=1e6
)
print("first eigenvalues")
print(vals0[order][:20])

n = 2 * (n - 1) + 1

model = ShockModel(n, M0=M0, length=length, dM0dx=dM0dx)
vals, vects = model.solve_eig()
vals, vects, order = model.select_and_sort(
    realmin=0., realmax=1000e3,
    imagmin=-1e5, imagmax=1e6,
    sort='imag'
)
print("first eigenvalues")
print(vals[order][:20])

# --- PLOT ---
plt.plot(vals0.real, vals0.imag, "ob", markersize=5, alpha=0.3)
plt.plot(vals.real, vals.imag, "or", markersize=3, alpha=0.6)
for i, v in enumerate(vals[order[:20]]):
    plt.text(v.real, v.imag, str(i))
plt.xlabel("$\omega_r$")
plt.ylabel("$\omega_i$")
plt.grid()
# plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
plt.title("Spectrum")
# plt.axis('square')
# plt.axis([-1000, 1000, -2, 0.5])
plt.show()
nv = 10
fig, ax = plt.subplots(nv, model.nvar)
for i in range(nv):
    # for j in (0, 2): # ? .real
    #     print(
    #         j, model._diffop.matder(1)[n - 1, :] @ vects[j * n : (j + 1) * n, order[i]]
    #     )
    for j in range(model.nvar):
        ax[i, j].plot(model._diffop.x, np.abs(vects[j * n : (j + 1) * n, order[i]]), 'b-')
        ax[i, j].plot(model._diffop.x, vects[j * n : (j + 1) * n, order[i]].real, 'b--')
plt.show()
