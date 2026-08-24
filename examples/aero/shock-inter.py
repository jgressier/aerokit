"""Interaction of two oblique shocks with a common upstream Mach number.

States 1 and 2 are behind the lower and upper incident shocks, respectively.
The transmitted shocks produce states 3 and 4 at the common downstream flow
angle and pressure.  Pressures are normalized by the pressure in state 0.
"""

import matplotlib.pyplot as plt
import numpy as np

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero.plot import shockpolar
from aerokit.instance.SWinteraction import ShockInteraction


# Problem parameters.  By convention, the bottom deviation is positive and
# the top deviation is negative.
M0 = 3.0
bottom_deviation = 10.0
top_deviation = -15.0

# Convert both prescribed incident-shock deviations to signed shock angles.
sigma01 = sw.weaksigma_Mach_deflection(M0, bottom_deviation)
sigma02 = sw.weaksigma_Mach_deflection(M0, top_deviation)

interaction = ShockInteraction(M0, sigma01, sigma02)
solution = interaction.solve()
if not solution.converged or not interaction.check34balanced():
    raise RuntimeError("the shock-interaction calculation did not converge")

theta_downstream = interaction[3].angle
sigma13 = sw.weaksigma_Mach_deflection(
    interaction[1].Mach, theta_downstream - interaction[1].angle
)
sigma24 = sw.weaksigma_Mach_deflection(
    interaction[2].Mach, theta_downstream - interaction[2].angle
)
print(f"bottom incident shock: sigma01 = {sigma01:.3f} deg")
print(f"top incident shock:    sigma02 = {sigma02:.3f} deg")
print(f"bottom transmitted shock: sigma13 = {sigma13:.3f} deg")
print(f"top transmitted shock:    sigma24 = {sigma24:.3f} deg")
print(f"common downstream deviation = {theta_downstream:.3f} deg")
print(f"common downstream pressure  = {interaction[3].p:.5f} p0")

fig, (ax_geom, ax) = plt.subplots(1, 2, figsize=(13, 5.5))
shockpolar.set_grid(ax)

# Incident-shock polar and the two translated polars for transmitted shocks.
shockpolar.plot_theta_pressure(
    M0, curve="both", color="0.25", label="incident-shock polar (state 0)", ax=ax
)
shockpolar.plot_theta_pressure(
    interaction[1].Mach,
    thet_init=interaction[1].angle,
    p_init=interaction[1].p,
    curve="left" if theta_downstream < interaction[1].angle else "right",
    color="tab:blue",
    label="transmitted polar from state 1",
    ax=ax,
)
shockpolar.plot_theta_pressure(
    interaction[2].Mach,
    thet_init=interaction[2].angle,
    p_init=interaction[2].p,
    curve="left" if theta_downstream < interaction[2].angle else "right",
    color="tab:orange",
    label="transmitted polar from state 2",
    ax=ax,
)

colors = ["black", "tab:blue", "tab:orange", "tab:blue", "tab:orange"]
markers = ["o", "o", "o", "X", "+"]
label_offsets = [(-6, -12), (10, -7), (-15, -7), (0, 15), (0, -22)]
for number, color, marker in zip(range(5), colors, markers):
    state = interaction[number]
    ax.plot(state.angle, state.p, marker=marker, color=color, markersize=8)
    ax.annotate(
        str(number),
        (state.angle, state.p),
        xytext=label_offsets[number],
        textcoords="offset points",
        color=color,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": color, "alpha": 0.9},
        arrowprops={"arrowstyle": "-", "color": color, "lw": 0.8},
    )

ax.set_xlabel(r"flow deviation $\theta$ (deg)")
ax.set_ylabel(r"normalized static pressure $p/p_0$")
ax.set_title(f"Two-shock interaction ($M_0={M0:g}$)")
handles, labels = ax.get_legend_handles_labels()
unique_legend = dict(zip(labels, handles))
ax.legend(unique_legend.values(), unique_legend.keys())

# Physical configuration.  Each shock direction is its upstream-flow angle
# plus its signed local shock angle sigma.
interaction_point = np.array([1.4, 0.0])
y_bottom, y_top = -1.0, 1.0

def point_on_horizontal(angle, y):
    """Intersection of a line through the interaction point with a given y."""
    return np.array([
        interaction_point[0] + (y - interaction_point[1]) / deg.tan(angle),
        y,
    ])


bottom_start = point_on_horizontal(sigma01, y_bottom)
top_start = point_on_horizontal(sigma02, y_top)
x_end = 2.2
bottom_transmitted_angle = interaction[1].angle + sigma13
top_transmitted_angle = interaction[2].angle + sigma24
bottom_end = interaction_point + np.array([
    x_end - interaction_point[0],
    (x_end - interaction_point[0]) * deg.tan(bottom_transmitted_angle),
])
top_end = interaction_point + np.array([
    x_end - interaction_point[0],
    (x_end - interaction_point[0]) * deg.tan(top_transmitted_angle),
])
slip_end = interaction_point + np.array([
    x_end - interaction_point[0],
    (x_end - interaction_point[0]) * deg.tan(theta_downstream),
])

shock_style = {"color": "tab:red", "lw": 2.5}
ax_geom.plot(*zip(bottom_start, interaction_point), **shock_style)
ax_geom.plot(*zip(top_start, interaction_point), **shock_style)
ax_geom.plot(*zip(interaction_point, bottom_end), **shock_style)
ax_geom.plot(*zip(interaction_point, top_end), **shock_style)
ax_geom.plot(*zip(interaction_point, slip_end), color="tab:green", lw=1.8, ls="--")
ax_geom.plot(*interaction_point, "ko", ms=5)

# Compression-corner walls that generate the two incident shocks.
x_left = -0.75
x_wall_end = 3.
y_min, y_max = -1.15, 1.15
bottom_wall_end = bottom_start + np.array([
    x_wall_end - bottom_start[0],
    (x_wall_end - bottom_start[0]) * deg.tan(bottom_deviation),
])
top_wall_end = top_start + np.array([
    x_wall_end - top_start[0],
    (x_wall_end - top_start[0]) * deg.tan(top_deviation),
])
wall_style = {"color": "black", "lw": 3.0, "solid_capstyle": "round"}

# Shade the solid side of each wall; the unshaded region is the flow passage.
wall_fill = {"facecolor": "0.75", "edgecolor": "0.55", "hatch": "///", "alpha": 0.55, "zorder": 0}
ax_geom.fill(
    [x_left, x_wall_end, bottom_wall_end[0], bottom_start[0], x_left],
    [y_min, y_min, bottom_wall_end[1], bottom_start[1], y_bottom],
    **wall_fill,
)
ax_geom.fill(
    [x_left, x_wall_end, top_wall_end[0], top_start[0], x_left],
    [y_max, y_max, top_wall_end[1], top_start[1], y_top],
    **wall_fill,
)
ax_geom.plot([x_left, bottom_start[0]], [y_bottom, y_bottom], **wall_style)
ax_geom.plot(*zip(bottom_start, bottom_wall_end), **wall_style)
ax_geom.plot([x_left, top_start[0]], [y_top, y_top], **wall_style)
ax_geom.plot(*zip(top_start, top_wall_end), **wall_style)

# Region labels are deliberately placed away from shock lines.
region_labels = [
    ("0", (-0.05, 0.0)),
    ("1", (0.75, -0.62)),
    ("2", (0.75, 0.62)),
    ("3", (2.20, -0.28)),
    ("4", (2.20, 0.18)),
]
for label, position in region_labels:
    ax_geom.text(
        *position, label, ha="center", va="center", fontweight="bold",
        bbox={"boxstyle": "circle,pad=0.25", "fc": "white", "ec": "0.3"},
    )

angle_labels = [
    (rf"$\sigma_{{01}}={sigma01:.1f}^\circ$", (0.3, -0.35), "center", "center"),
    (rf"$\sigma_{{02}}={sigma02:.1f}^\circ$", (0.3, 0.40), "center", "center"),
    (
        rf"$\sigma_{{13}}={sigma13:.1f}^\circ$",
        bottom_end + np.array([0.08, 0.08 * deg.tan(bottom_transmitted_angle)]),
        "left",
        "top",
    ),
    (
        rf"$\sigma_{{24}}={sigma24:.1f}^\circ$",
        top_end + np.array([0.08, 0.08 * deg.tan(top_transmitted_angle)]),
        "left",
        "bottom",
    ),
]
for label, position, horizontal_alignment, vertical_alignment in angle_labels:
    ax_geom.text(
        *position, label, ha=horizontal_alignment, va=vertical_alignment,
        bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
    )

ax_geom.annotate(
    "slip line\n" + rf"$\theta={theta_downstream:.1f}^\circ$",
    slip_end,
    xytext=(5, 0),
    textcoords="offset points",
    color="tab:green",
)
ax_geom.set_xlim(x_left, 3.)
ax_geom.set_ylim(y_min, y_max)
ax_geom.set_aspect("equal", adjustable="box")
ax_geom.set_xlabel("x")
ax_geom.set_ylabel("y")
ax_geom.set_title("Physical configuration")
ax_geom.grid(ls=":", alpha=0.4)
fig.tight_layout()
plt.show()
