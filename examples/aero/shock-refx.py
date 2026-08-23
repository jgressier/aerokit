"""Regular reflection of an oblique shock, shown in a pressure polar.

This is the script version of ``notebooks/aero/Shock-21-reflexion.ipynb``.
Pressures are normalized by the upstream static pressure p0 and angles are
expressed in degrees.
"""

import matplotlib.pyplot as plt
import numpy as np

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero.plot import shockpolar


# Problem parameters
M0 = 2.8
wall_deviation = 15.0

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

states = {
    "0": (0.0, 1.0),
    "1": (wall_deviation, p1_p0),
    "2": (0.0, p2_p0),
}
label_offsets = {"0": (-6, -15), "1": (10, -8), "2": (8, 8)}
for label, (theta, pressure) in states.items():
    ax.plot(theta, pressure, "o", markersize=7)
    ax.annotate(
        label,
        (theta, pressure),
        xytext=label_offsets[label],
        textcoords="offset points",
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "0.3", "alpha": 0.9},
        arrowprops={"arrowstyle": "-", "color": "0.3", "lw": 0.8},
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
reflection = np.array([1.0 / np.tan(np.deg2rad(sigma01)), y_top])
reflected_angle = wall_deviation - sigma12
wall_slope = np.tan(np.deg2rad(wall_deviation))
reflected_slope = np.tan(np.deg2rad(reflected_angle))
x_reflected_end = (
    reflection[1] - reflected_slope * reflection[0]
) / (wall_slope - reflected_slope)
reflected_end = np.array([x_reflected_end, wall_slope * x_reflected_end])
bottom_wall_end = np.array([x_right, wall_slope * x_right])

wall_fill = {"facecolor": "0.75", "edgecolor": "0.55", "hatch": "///", "alpha": 0.55, "zorder": 0}
ax_geom.fill(
    [x_left, x_right, x_right, x_left],
    [y_top, y_top, y_max, y_max],
    **wall_fill,
)
ax_geom.fill(
    [x_left, x_right, bottom_wall_end[0], corner[0], x_left],
    [y_min, y_min, bottom_wall_end[1], corner[1], y_bottom],
    **wall_fill,
)
wall_style = {"color": "black", "lw": 3.0, "solid_capstyle": "round"}
ax_geom.plot([x_left, x_right], [y_top, y_top], **wall_style)
ax_geom.plot([x_left, corner[0]], [y_bottom, y_bottom], **wall_style)
ax_geom.plot(*zip(corner, bottom_wall_end), **wall_style)

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
