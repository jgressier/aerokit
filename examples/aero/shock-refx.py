"""Regular reflection of an oblique shock, shown in a pressure polar.

This is the script version of ``notebooks/aero/Shock-21-reflexion.ipynb``.
Pressures are normalized by the upstream static pressure p0 and angles are
expressed in degrees.
"""

import matplotlib.pyplot as plt
import numpy as np

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero import model2D as m2d
from aerokit.aero.plot.geom import Geom
from aerokit.aero.plot import shockpolar


# Problem parameters
M0 = 2.8
wall_deviation = 15.0

# Incident shock: state 0 -> state 1.
sigma01 = sw.weaksigma_Mach_deflection(M0, wall_deviation)
state = {0: m2d.State2DMach(M0)}
state[1] = state[0].weakshock_deviation(wall_deviation)
M1 = state[1].Mach
p1_p0 = state[1].p

# Reflected shock: state 1 -> state 2.  It turns the flow back to zero degrees.
sigma12 = sw.weaksigma_Mach_deflection(M1, wall_deviation)
state[2] = state[1].weakshock_deviation(-wall_deviation)
M2 = state[2].Mach
p2_p0 = state[2].p

print(f"incident shock:  sigma01 = {sigma01:.3f} deg, M1 = {M1:.3f}")
print(f"reflected shock: sigma12 = {sigma12:.3f} deg, M2 = {M2:.3f}")

fig, (ax_geom, ax) = plt.subplots(1, 2, figsize=(13, 5.5))
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

label_offsets = {"0": (-6, -15), "1": (10, -8), "2": (8, 8)}
for number in range(3):
    label = str(number)
    shockpolar.plot_state(
        state[number], label=label, offset=label_offsets[label], markersize=7, ax=ax
    )

ax.set_xlabel(r"flow deviation $\theta$ (deg)")
ax.set_ylabel(r"normalized static pressure $p/p_0$")
ax.set_title(f"Regular shock reflection ($M_0={M0:g}$)")
ax.legend()

# Physical configuration: an incident shock from the compression corner is
# reflected by the horizontal upper wall.
x_left, x_right = -0.2, 2.
y_bottom, y_top = 0.0, 1.0
y_min, y_max = -0.15, 1.15
corner = np.array([0.0, y_bottom])
reflection = np.array([1.0 / deg.tan(sigma01), y_top])
reflected_angle = wall_deviation - sigma12
wall_slope = deg.tan(wall_deviation)
reflected_slope = deg.tan(reflected_angle)
x_reflected_end = (
    reflection[1] - reflected_slope * reflection[0]
) / (wall_slope - reflected_slope)
reflected_end = np.array([x_reflected_end, wall_slope * x_reflected_end])
bottom_wall_end = np.array([x_right, wall_slope * x_right])

ax_geom.set(xlim=(x_left, x_right), ylim=(y_min, y_max))
geometry = Geom()
geometry.add_wall((x_left, x_right), (y_top, y_top), location="top")
geometry.add_wall(
    (x_left, corner[0], bottom_wall_end[0]),
    (y_bottom, corner[1], bottom_wall_end[1]),
    location="bottom",
)
geometry.plot(ax=ax_geom)

shock_style = {"color": "tab:red", "lw": 2.5}
ax_geom.plot(*zip(corner, reflection), **shock_style)
ax_geom.plot(*zip(reflection, reflected_end), **shock_style)
ax_geom.plot(*reflection, "ko", ms=5)

region_labels = [
    ("0", (0.05, 0.72)),
    ("1", (1, 0.5)),
    ("2", reflection + np.array([0.4, -0.1])),
]
for label, position in region_labels:
    ax_geom.text(
        *position,
        label,
        ha="center",
        va="center",
        fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.25", "fc": "white", "ec": "0.3"},
    )

ax_geom.text(
    0.5,
    0.3,
    rf"$\sigma_{{01}}={sigma01:.1f}^\circ$",
    bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
)
ax_geom.text(
    *(reflection + np.array([0.1, -0.25])),
    rf"$\sigma_{{12}}={-sigma12:.1f}^\circ$",
    ha="left",
    va="center",
    bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
)
ax_geom.set_xlim(x_left, x_right)
ax_geom.set_ylim(y_min, y_max)
ax_geom.set_aspect("equal", adjustable="box")
ax_geom.set_xlabel("x")
ax_geom.set_ylabel("y")
ax_geom.set_title("Physical configuration")
ax_geom.grid(ls=":", alpha=0.4)
fig.tight_layout()
plt.show()
