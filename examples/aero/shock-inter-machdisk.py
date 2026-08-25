"""Shock interaction for which the regular solution is detached.

At this lower Mach number and larger wall deflection, the two downstream
shock polars do not intersect.  Their detachment-limit points retain an
angular gap.  States 3 and 4 are computed independently at triple points 03
and 04 by intersecting each reflected polar with the strong state-0 polar.

The physical view is a slightly asymmetric Mach-disk schematic.  A complete Mach
reflection also contains curved shocks and slip layers near both triple
points; the straight segments drawn here use the locally computed angles.
"""

import matplotlib.pyplot as plt
import numpy as np

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero import model2D as m2d
from aerokit.aero.plot import shockpolar
from aerokit.aero.plot.geom import Geom
from aerokit.instance.SWinteraction import ShockInteraction


# A regular weak/weak interaction exists in shock-inter.py.  These stronger,
# unequal compressions deliberately exceed the post-incident turning range.
M0 = 2.3
bottom_deviation = 20.0
top_deviation = -15.0

state = {0: m2d.State2DMach(M0)}
sigma01 = sw.weaksigma_Mach_deflection(M0, bottom_deviation)
sigma02 = sw.weaksigma_Mach_deflection(M0, top_deviation)
state[1] = state[0].shock_sigma(sigma01)
state[2] = state[0].shock_sigma(sigma02)

# The limiting transmitted shocks are at the maximum-deflection points of
# the translated polars.  Their endpoint angles expose the detachment gap.
bottom_transmitted_devmax = sw.dev_Max(state[1].Mach)
top_transmitted_devmax = sw.dev_Max(state[2].Mach)
sigma13 = -sw.sigma_DevMax(state[1].Mach)
sigma24 = sw.sigma_DevMax(state[2].Mach)
detachment_state_bottom = state[1].shock_sigma(sigma13)
detachment_state_top = state[2].shock_sigma(sigma24)
polar_gap = detachment_state_bottom.angle - detachment_state_top.angle

if polar_gap <= 0.0:
    raise RuntimeError("this case still admits a regular shock-polar intersection")

triple_03 = ShockInteraction.solve_triple_point(state[0], state[1])
triple_04 = ShockInteraction.solve_triple_point(state[0], state[2])
state[3], stem_state_03 = triple_03.reflected_state, triple_03.stem_state
state[4], stem_state_04 = triple_04.reflected_state, triple_04.stem_state
sigma13, sigma03 = triple_03.sigma_reflected, triple_03.sigma_stem
sigma24, sigma04 = triple_04.sigma_reflected, triple_04.sigma_stem

print(f"bottom incident shock: sigma01 = {sigma01:.3f} deg")
print(f"top incident shock:    sigma02 = {sigma02:.3f} deg")
print(f"post-incident Mach numbers: M1 = {state[1].Mach:.4f}, M2 = {state[2].Mach:.4f}")
print(
    "available transmitted turning: "
    f"state 1 = {bottom_transmitted_devmax:.3f} deg, "
    f"state 2 = {top_transmitted_devmax:.3f} deg"
)
print(f"detachment gap between polars = {polar_gap:.3f} deg")
print(
    f"triple point 03: theta = {state[3].angle:.3f} deg, "
    f"p/p0 = {state[3].p:.4f}, sigma03 = {sigma03:.3f} deg"
)
print(
    f"triple point 04: theta = {state[4].angle:.3f} deg, "
    f"p/p0 = {state[4].p:.4f}, sigma04 = {sigma04:.3f} deg"
)

fig, (ax_geom, ax_polar) = plt.subplots(1, 2, figsize=(13, 5.5))

# ---------------------------------------------------------------------------
# Pressure/deviation polar
shockpolar.set_grid(ax_polar)
shockpolar.plot_theta_pressure(
    M0, curve="both", color="0.25", label="incident-shock polar (state 0)", ax=ax_polar
)
shockpolar.plot_theta_pressure(
    state[1].Mach,
    thet_init=state[1].angle,
    p_init=state[1].p,
    curve="left",
    color="tab:blue",
    label="polar from state 1",
    ax=ax_polar,
)
shockpolar.plot_theta_pressure(
    state[2].Mach,
    thet_init=state[2].angle,
    p_init=state[2].p,
    curve="right",
    color="tab:orange",
    label="polar from state 2",
    ax=ax_polar,
)
strong_interval_style = {
    "color": "tab:red",
    "linewidth": 3.5,
    "label": "strong-polar interval 03--04",
    "ax": ax_polar,
}
if state[3].angle * state[4].angle >= 0.0:
    # Both triple points are on the same side of the strong polar.
    shockpolar.plot_theta_pressure(
        M0,
        sigma_range=tuple(sorted((abs(sigma03), abs(sigma04)))),
        curve="right" if state[3].angle >= 0.0 else "left",
        **strong_interval_style,
    )
else:
    # Opposite branches require one stroke from each triple point to the
    # normal-shock point at sigma=90 degrees.
    shockpolar.plot_theta_pressure(
        M0,
        sigma_range=(abs(sigma03), 90.0),
        curve="right" if state[3].angle > 0.0 else "left",
        **strong_interval_style,
    )
    strong_interval_style.pop("label")
    shockpolar.plot_theta_pressure(
        M0,
        sigma_range=(abs(sigma04), 90.0),
        curve="right" if state[4].angle > 0.0 else "left",
        **strong_interval_style,
    )

point_colors = ["black", "tab:blue", "tab:orange", "tab:blue", "tab:orange"]
point_markers = ["o", "o", "o", "X", "+"]
label_offsets = [(-6, -15), (10, -8), (-17, -8), (10, 10), (-20, 10)]
for number, color, marker, offset in zip(range(5), point_colors, point_markers, label_offsets):
    shockpolar.plot_state(
        state[number], label=str(number), marker=marker, color=color, offset=offset, ax=ax_polar
    )

for label, triple_state, offset in (
    ("03", stem_state_03, (10, -28)),
    ("04", stem_state_04, (-28, -28)),
):
    shockpolar.plot_state(
        triple_state,
        label=label,
        marker="D",
        color="tab:red",
        markersize=6,
        offset=offset,
        ax=ax_polar,
    )

gap_pressure = 0.5 * (detachment_state_bottom.p + detachment_state_top.p)
ax_polar.annotate(
    "",
    xy=(detachment_state_bottom.angle, gap_pressure),
    xytext=(detachment_state_top.angle, gap_pressure),
    arrowprops={"arrowstyle": "<->", "color": "tab:red", "lw": 1.5},
)
ax_polar.text(
    0.0,
    gap_pressure + 0.35,
    rf"no intersection: $\Delta\theta={polar_gap:.1f}^\circ$",
    color="tab:red",
    ha="center",
    bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.9},
)
ax_polar.set_xlabel(r"flow deviation $\theta$ (deg)")
ax_polar.set_ylabel(r"normalized static pressure $p/p_0$")
ax_polar.set_title(rf"Detached interaction ($M_0={M0:g}$)")
handles, labels = ax_polar.get_legend_handles_labels()
unique_legend = dict(zip(labels, handles))
ax_polar.legend(unique_legend.values(), unique_legend.keys())

# ---------------------------------------------------------------------------
# Physical configuration
x_left, x_right = -0.5, 2.5
y_bottom, y_top = -1.0, 1.0
y_min, y_max = -1.15, 1.15
x_triple = 1.4
triple_bottom = np.array([x_triple, -0.20])
triple_top = np.array([x_triple, 0.16])


def shock_origin(point, shock_angle, wall_y):
    """Return the wall point joined to ``point`` by a shock at ``shock_angle``."""
    return np.array([
        point[0] + (wall_y - point[1]) / deg.tan(shock_angle),
        wall_y,
    ])


bottom_corner = shock_origin(triple_bottom, sigma01, y_bottom)
top_corner = shock_origin(triple_top, sigma02, y_top)

# Compression walls and shaded solid material.
bottom_wall_end = np.array([
    x_right,
    y_bottom + (x_right - bottom_corner[0]) * deg.tan(bottom_deviation),
])
top_wall_end = np.array([
    x_right,
    y_top + (x_right - top_corner[0]) * deg.tan(top_deviation),
])
ax_geom.set(xlim=(x_left, x_right), ylim=(y_min, y_max))
geometry = Geom()
geometry.add_wall(
    (x_left, bottom_corner[0], bottom_wall_end[0]),
    (y_bottom, bottom_corner[1], bottom_wall_end[1]),
    location="bottom",
)
geometry.add_wall(
    (x_left, top_corner[0], top_wall_end[0]),
    (y_top, top_corner[1], top_wall_end[1]),
    location="top",
)
geometry.plot(ax=ax_geom)

# Incident shocks, triple-point reflected shocks, and the Mach disk.
bottom_transmitted_angle = state[1].angle + sigma13
top_transmitted_angle = state[2].angle + sigma24
reflected_dx = 0.25
bottom_transmitted_end = triple_bottom + np.array([
    reflected_dx,
    reflected_dx * deg.tan(bottom_transmitted_angle),
])
top_transmitted_end = triple_top + np.array([
    reflected_dx,
    reflected_dx * deg.tan(top_transmitted_angle),
])
shock_style = {"color": "tab:red", "lw": 2.5}
ax_geom.plot(*zip(bottom_corner, triple_bottom), **shock_style)
ax_geom.plot(*zip(top_corner, triple_top), **shock_style)
ax_geom.plot(*zip(triple_bottom, bottom_transmitted_end), **shock_style)
ax_geom.plot(*zip(triple_top, top_transmitted_end), **shock_style)
ax_geom.plot(
    [x_triple, x_triple],
    [triple_bottom[1], triple_top[1]],
    color="tab:red",
    lw=6,
    solid_capstyle="round",
    label="Mach disk",
)
ax_geom.plot(*triple_bottom, "D", color="tab:red", markeredgecolor="black", ms=7)
ax_geom.plot(*triple_top, "D", color="tab:red", markeredgecolor="black", ms=7)
ax_geom.annotate(
    "03", triple_bottom, xytext=(5, 5), textcoords="offset points",
    color="tab:red", fontweight="bold",
)
ax_geom.annotate(
    "04", triple_top, xytext=(5, -8), textcoords="offset points",
    color="tab:red", fontweight="bold",
)

# The two limiting downstream directions bound the unresolved slip-layer
# region.  They demonstrate geometrically the same angular gap as the polar.
x_slip_end = 2.30
for triple, theta in ((triple_bottom, state[3].angle), (triple_top, state[4].angle)):
    slip_end = triple + np.array([
        x_slip_end - x_triple,
        (x_slip_end - x_triple) * deg.tan(theta),
    ])
    ax_geom.plot(*zip(triple, slip_end), color="tab:green", lw=1.7, ls="--")

bottom_incident_midpoint = 0.5 * (bottom_corner + triple_bottom)
top_incident_midpoint = 0.5 * (top_corner + triple_top)
bottom_reflected_midpoint = 0.5 * (triple_bottom + bottom_transmitted_end)
top_reflected_midpoint = 0.5 * (triple_top + top_transmitted_end)
region_labels = [
    ("0", (0.10, 0.0)),
    ("1", bottom_incident_midpoint + np.array([0.10, -0.10])),
    ("2", top_incident_midpoint + np.array([0.10, 0.10])),
    ("3", bottom_reflected_midpoint + np.array([0.15, 0.0])),
    ("4", top_reflected_midpoint + np.array([0.15, 0.0])),
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

angle_labels = [
    (
        rf"$\sigma_{{01}}={sigma01:.1f}^\circ$",
        bottom_incident_midpoint + np.array([-0.05, 0.05]),
        "right",
        "center",
    ),
    (
        rf"$\sigma_{{02}}={sigma02:.1f}^\circ$",
        top_incident_midpoint + np.array([-0.05, -0.05]),
        "right",
        "center",
    ),
    (rf"$\sigma_{{13}}={sigma13:.1f}^\circ$", bottom_transmitted_end + (0.05, -0.02), "left", "top"),
    (rf"$\sigma_{{24}}={sigma24:.1f}^\circ$", top_transmitted_end + (0.05, 0.02), "left", "bottom"),
]
for label, position, horizontal_alignment, vertical_alignment in angle_labels:
    ax_geom.text(
        *position,
        label,
        ha=horizontal_alignment,
        va=vertical_alignment,
        bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
    )

ax_geom.annotate(
    "Mach disk",
    (x_triple, 0.0),
    xytext=(-10, 0),
    textcoords="offset points",
    color="tab:red",
    ha="right",
    va="center",
    arrowprops={"arrowstyle": "->", "color": "tab:red", "lw": 1.0},
)
ax_geom.set_xlim(x_left, x_right)
ax_geom.set_ylim(y_min, y_max)
ax_geom.set_aspect("equal", adjustable="box")
ax_geom.set_xlabel("x")
ax_geom.set_ylabel("y")
ax_geom.set_title("Mach-disk configuration")
ax_geom.grid(ls=":", alpha=0.4)

fig.tight_layout()
plt.show()
