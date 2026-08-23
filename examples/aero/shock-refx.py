"""Regular reflection of an oblique shock, shown in a pressure polar.

This is the script version of ``notebooks/aero/Shock-21-reflexion.ipynb``.
Pressures are normalized by the upstream static pressure p0 and angles are
expressed in degrees.
"""

import matplotlib.pyplot as plt

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero.plot import shockpolar


# Problem parameters
M0 = 2.8
wall_deviation = 18.0

# Incident shock: state 0 -> state 1.
sigma01 = sw.weaksigma_Mach_deflection(M0, wall_deviation)
Mn0 = M0 * deg.sin(sigma01)
p1_p0 = sw.Ps_ratio(Mn0)
Mn1 = sw.downstream_Mn(Mn0)
M1 = Mn1 / deg.sin(sigma01 - wall_deviation)

# Reflected shock: state 1 -> state 2.  It turns the flow back to zero degrees.
sigma12 = sw.weaksigma_Mach_deflection(M1, wall_deviation)
Mn1 = M1 * deg.sin(sigma12)
p2_p1 = sw.Ps_ratio(Mn1)
Mn2 = sw.downstream_Mn(Mn1)
M2 = Mn2 / deg.sin(sigma12 - wall_deviation)
p2_p0 = p1_p0 * p2_p1

print(f"incident shock:  sigma01 = {sigma01:.3f} deg, M1 = {M1:.3f}")
print(f"reflected shock: sigma12 = {sigma12:.3f} deg, M2 = {M2:.3f}")

fig, ax = plt.subplots(figsize=(7, 5))
shockpolar.set_grid(ax)

# Only the branches used by the incident and reflected shocks are highlighted.
shockpolar.plot_theta_pressure(
    M0, curve="right", devmax=True, sonic=True, label="polar from state 0", ax=ax
)
shockpolar.plot_theta_pressure(
    M1,
    thet_init=wall_deviation,
    p_init=p1_p0,
    curve="left",
    color="tab:red",
    label="polar from state 1",
    ax=ax,
)

states = {
    "0": (0.0, 1.0),
    "1": (wall_deviation, p1_p0),
    "2": (0.0, p2_p0),
}
for label, (theta, pressure) in states.items():
    ax.plot(theta, pressure, "o", markersize=7)
    ax.annotate(label, (theta, pressure), xytext=(6, 5), textcoords="offset points")

ax.set_xlabel(r"flow deviation $\theta$ (deg)")
ax.set_ylabel(r"normalized static pressure $p/p_0$")
ax.set_title(f"Regular shock reflection ($M_0={M0:g}$)")
ax.legend()
fig.tight_layout()
plt.show()
