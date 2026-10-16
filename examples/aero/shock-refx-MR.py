"""Mach reflection when a regular weak reflected shock cannot exist.

The incident shock turns the flow through the ramp angle.  At this lower Mach
number, state 1 cannot be turned all the way back to the horizontal wall by a
second attached shock.  The Mach-reflection solution is obtained from the
intersection of the reflected-shock polar from state 1 with the strong branch
from state 0 (the Mach stem).
"""

import matplotlib.pyplot as plt
import numpy as np

from aerokit.aero import ShockWave as sw
from aerokit.aero import degree as deg
from aerokit.aero import model2D as m2d
from aerokit.aero.plot.geom import Geom
from aerokit.aero.plot import shockpolar
from aerokit.instance.SWinteraction import ShockInteraction


# Problem parameters: deliberately beyond the regular-reflection limit.
M0 = 2.3
wall_deviation = 20.0

state = {0: m2d.State2DMach(M0)}
sigma01 = sw.weaksigma_Mach_deflection(M0, wall_deviation)
state[1] = state[0].shock_sigma(sigma01)
reflected_devmax = sw.dev_Max(state[1].Mach)

if reflected_devmax >= wall_deviation:
    raise RuntimeError("this case still admits a regular weak reflection")

# Three-shock construction. Across the slip line, the reflected and Mach-stem
# states have equal pressure and flow angle.
triple = ShockInteraction.solve_triple_point(state[0], state[1])
theta_downstream = triple.theta
state[2] = triple.reflected_state
state[3] = triple.stem_state
sigma12 = triple.sigma_reflected
sigma03 = triple.sigma_stem

pressure_error = abs(state[2].p - state[3].p)
angle_error = abs(state[2].angle - state[3].angle)
if pressure_error > 1.0e-8 or angle_error > 1.0e-8:
    raise RuntimeError("the three-shock construction is not balanced")

print(f"incident shock: sigma01 = {sigma01:.3f} deg")
print(f"state 1: M1 = {state[1].Mach:.4f}, p1/p0 = {state[1].p:.4f}")
print(f"required regular-reflection turn = {wall_deviation:.3f} deg")
print(f"maximum attached turn from state 1 = {reflected_devmax:.3f} deg")
print(f"Mach-reflection direction = {theta_downstream:.3f} deg")
print(f"reflected shock: sigma12 = {sigma12:.3f} deg")
print(f"Mach stem: sigma03 = {sigma03:.3f} deg")
print(f"matched pressure: p2/p0 = p3/p0 = {state[2].p:.4f}")

fig, (ax_geom, ax_polar) = plt.subplots(1, 2, figsize=(13, 5.5))

# ---------------------------------------------------------------------------
# Three-shock pressure polar
shockpolar.set_grid(ax_polar)
shockpolar.plot_theta_pressure(
    M0,
    curve="right",
    color="0.25",
    label="polar from state 0 / Mach stem",
    ax=ax_polar,
)
shockpolar.plot_theta_pressure(
    state[1].Mach,
    thet_init=state[1].angle,
    p_init=state[1].p,
    curve="left",
    color="tab:blue",
    label="reflected-shock polar from state 1",
    ax=ax_polar,
)
shockpolar.plot_theta_pressure(
    M0,
    sigma_range=(sigma03, 90.0),
    curve="right",
    color="tab:red",
    linewidth=3.5,
    label="strong branch: state 3 to normal shock",
    ax=ax_polar,
)

colors = ["black", "tab:blue", "tab:blue", "tab:red"]
markers = ["o", "o", "X", "+"]
offsets = [(-6, -15), (10, -8), (12, 10), (12, -18)]
for number, color, marker, offset in zip(range(4), colors, markers, offsets):
    shockpolar.plot_state(
        state[number], label=str(number), marker=marker, color=color, offset=offset, ax=ax_polar
    )

ax_polar.axvline(0.0, color="tab:red", ls=":", lw=1.2)
ax_polar.text(
    0.6,
    state[1].p,
    "regular-reflection target\nis outside the reflected polar",
    color="tab:red",
    va="bottom",
    bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.9},
)
ax_polar.set_xlabel(r"flow deviation $\theta$ (deg)")
ax_polar.set_ylabel(r"normalized static pressure $p/p_0$")
ax_polar.set_title(rf"Mach-reflection polar ($M_0={M0:g}$)")
ax_polar.legend()

# ---------------------------------------------------------------------------
# Physical Mach-reflection configuration
x_left, x_right = -0.4, 1.5
y_bottom, y_top = 0.0, 1.0
y_min, y_max = -0.15, 1.15
corner = np.array([0.0, y_bottom])
triple_height = 0.78
triple_point = np.array([
    triple_height / deg.tan(sigma01),
    triple_height,
])

stem_global_angle = sigma03
stem_height = y_top - triple_point[1]

# Parabolic Mach stem x(y).  Its slope is prescribed analytically:
#   dx/dy = cot(sigma03) at the triple point,
#   dx/dy = 0            at the wall (a vertical tangent).
stem_y = np.linspace(triple_point[1], y_top, 80)
stem_dy = stem_y - triple_point[1]
stem_cotangent = 1.0 / deg.tan(stem_global_angle)
stem_x = (
    triple_point[0]
    + stem_cotangent * stem_dy
    - stem_cotangent * stem_dy**2 / (2.0 * stem_height)
)
stem_wall = np.array([stem_x[-1], stem_y[-1]])

wall_slope = deg.tan(wall_deviation)
bottom_wall_end = np.array([x_right, wall_slope * x_right])
reflected_global_angle = state[1].angle + sigma12
downstream_dx = 0.3
reflected_end = triple_point + np.array([
    downstream_dx,
    downstream_dx * deg.tan(reflected_global_angle),
])
slip_end = triple_point + np.array([
    downstream_dx,
    downstream_dx * deg.tan(theta_downstream),
])

# Walls and solid shading.
ax_geom.set(xlim=(x_left, x_right), ylim=(y_min, y_max))
geometry = Geom()
geometry.add_wall((x_left, x_right), (y_top, y_top), location="top")
geometry.add_wall(
    (x_left, corner[0], bottom_wall_end[0]),
    (y_bottom, corner[1], bottom_wall_end[1]),
    location="bottom",
)
geometry.plot(ax=ax_geom)

# Incident shock, reflected shock, Mach stem, and slip line.
shock_style = {"color": "tab:red", "lw": 2.5}
ax_geom.plot(*zip(corner, triple_point), **shock_style)
ax_geom.plot(*zip(triple_point, reflected_end), **shock_style)
ax_geom.plot(stem_x, stem_y, color="tab:red", lw=6, solid_capstyle="round")
ax_geom.plot(*zip(triple_point, slip_end), color="tab:green", lw=1.8, ls="--")
ax_geom.plot(*triple_point, "ko", ms=5)

region_labels = [
    ("0", (0.05, 0.6)),
    ("1", (0.50, 0.35)),
    ("2", (1.12, 0.67)),
    ("3", tuple(triple_point + np.array([0.06, 0.06]))),
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
    (rf"$\sigma_{{01}}={sigma01:.1f}^\circ$", (0.28, 0.22), "center", "center"),
    (rf"$\sigma_{{12}}={sigma12:.1f}^\circ$", reflected_end + (0.05, 0.0), "left", "center"),
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
    "Mach-stem tangent:\n"
    + rf"$\sigma_{{03}}={sigma03:.1f}^\circ$",
    (stem_x[0], stem_y[0]),
    xytext=(-10, 0),
    textcoords="offset points",
    color="tab:red",
    ha="right",
    va="center",
    arrowprops={"arrowstyle": "->", "color": "tab:red", "lw": 1.0},
)
ax_geom.annotate(
    "slip line\n" + rf"$\theta={theta_downstream:.1f}^\circ$",
    slip_end,
    xytext=(5, 4),
    textcoords="offset points",
    color="tab:green",
)
ax_geom.set_xlim(x_left, x_right)
ax_geom.set_ylim(y_min, y_max)
ax_geom.set_aspect("equal", adjustable="box")
ax_geom.set_xlabel("x")
ax_geom.set_ylabel("y")
ax_geom.set_title("Mach reflection")
ax_geom.grid(ls=":", alpha=0.4)

fig.tight_layout()
plt.show()
