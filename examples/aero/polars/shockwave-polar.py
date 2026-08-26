"""Shock-wave polar maps for a calorically perfect gas.

Six figures are produced: the complete deviation/shock-angle polars; a
deviation/shock-angle map of shock-to-isentropic compression; the weak shock
branches; separate weak and strong deviation/upstream-Mach maps; and an
upstream/downstream-Mach map.
Run with ``--debug`` to add the underlying data mesh beside every figure.
"""

import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

import aerokit.aero.ShockWave as sw
import aerokit.aero.Isentropic as Is
import aerokit.aero.Supersonic as sup
import aerokit.aero.degree as deg


npoints = 160
npoints_mach = 500
debug = "--debug" in sys.argv
debug_mesh_row_skip = 20
debug_mesh_column_skip = 10
figure6_mach_max = 6.0
gamma = 1.4
mach_values = np.array([
    1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
    2.2, 2.5, 3.0, 3.5, 4.0, 5.0, 10.0, 100.0,
])


def set_axes_style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.minorticks_on()
    ax.grid(which="major", linestyle="-", alpha=0.7)
    ax.grid(which="minor", linestyle=":", alpha=0.5)


def make_figure(title):
    """Create a normal plot, or a plot plus its data-mesh debug panel."""
    if debug:
        fig, (ax, mesh_ax) = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle(title)
        ax.set_title("result")
        return fig, ax, mesh_ax
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_title(title)
    return fig, ax, None


def plot_debug_mesh(mesh_ax, x, y, xlabel, ylabel, row_skip, column_skip):
    """Display selected row and column lines from a structured 2-D mesh."""
    if mesh_ax is None:
        return
    row_indices = np.unique(
        np.r_[np.arange(0, x.shape[0], row_skip), x.shape[0] - 1]
    )
    column_indices = np.unique(
        np.r_[np.arange(0, x.shape[1], column_skip), x.shape[1] - 1]
    )
    for index in row_indices:
        mesh_ax.plot(x[index, :], y[index, :], color="tab:blue", linewidth=0.7)
    for index in column_indices:
        mesh_ax.plot(x[:, index], y[:, index], color="tab:orange", linewidth=0.7)
    mesh_ax.set_title(
        f"data mesh (rows / {row_skip}, columns / {column_skip})"
    )
    set_axes_style(mesh_ax, xlabel, ylabel)


def shock_to_isentropic_ratio(mach_grid, sigma_grid, deviation_grid):
    """Return shock/isentropic static-pressure ratio for the same deviation."""
    shock_pressure = sw.Ps_ratio(mach_grid * deg.sin(sigma_grid), gamma)
    isentropic_pressure = np.full_like(shock_pressure, np.nan)
    for row, mach in enumerate(mach_grid[:, 0]):
        max_characteristic_deviation = sup.PrandtlMeyer_Mach(mach, gamma)
        valid = deviation_grid[row] <= max_characteristic_deviation + 1.0e-10
        remaining_pm_angle = np.maximum(
            max_characteristic_deviation - deviation_grid[row, valid], 0.0
        )
        isentropic_mach = np.ones_like(remaining_pm_angle)
        supersonic = remaining_pm_angle > 1.0e-10
        if np.any(supersonic):
            isentropic_mach[supersonic] = sup.Mach_PrandtlMeyer(
                remaining_pm_angle[supersonic], gamma
            )
        isentropic_pressure[row, valid] = (
            Is.PtPs_Mach(mach, gamma)
            / Is.PtPs_Mach(isentropic_mach, gamma)
        )
    return np.ma.masked_invalid(shock_pressure / isentropic_pressure)


# =============================================================================
# FIGURE 1 - COMPLETE WEAK AND STRONG SHOCK POLARS
# =============================================================================
# Each curve sweeps the shock angle from the Mach angle (vanishing shock) to a
# normal shock.  The two coloured loci separate weak/strong and
# supersonic/subsonic downstream solutions.
fig_full, ax_full, ax_full_mesh = make_figure(
    rf"Complete shock-wave polars, $\gamma={gamma:.1f}$"
)
full_deviation_rows = []
full_sigma_rows = []
for mach in mach_values:
    sigma = np.linspace(deg.asin(1.0 / mach), 90.0, npoints + 1)
    deviation = sw.deflection_Mach_sigma(mach, sigma, gamma)
    full_deviation_rows.append(deviation)
    full_sigma_rows.append(sigma)
    ax_full.plot(deviation, sigma, "k-")
    # Put Mach labels on the weak branch while the downstream flow is still
    # supersonic, i.e. below the M2=1 sonic locus.
    label_sigma = deg.asin(1.0 / mach) + 0.65 * (
        sw.sigma_Sonic(mach, gamma) - deg.asin(1.0 / mach)
    )
    label_index = np.abs(sigma - label_sigma).argmin()
    tangent_angle = deg.atan2(
        sigma[label_index + 1] - sigma[label_index - 1],
        deviation[label_index + 1] - deviation[label_index - 1],
    )
    ax_full.text(
        deviation[label_index], sigma[label_index], rf"$M_0={mach:g}$",
        ha="left", va="top", fontsize=7, rotation=tangent_angle,
        rotation_mode="anchor", transform_rotates_text=True,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.0},
    )

boundary_mach = np.geomspace(1.001, 100.0, 300)
max_sigma = sw.sigma_DevMax(boundary_mach, gamma)
max_deviation = sw.deflection_Mach_sigma(boundary_mach, max_sigma, gamma)
sonic_sigma = sw.sigma_Sonic(boundary_mach, gamma)
sonic_deviation = sw.deflection_Mach_sigma(boundary_mach, sonic_sigma, gamma)
ax_full.plot(
    max_deviation, max_sigma, color="tab:green", linewidth=2.2,
    label="maximum-deviation locus",
)
ax_full.plot(
    sonic_deviation, sonic_sigma, color="tab:blue", linestyle="--", linewidth=2.2,
    label=r"sonic locus, $M_2=1$",
)
ax_full.set_xlim(0.0, 50.0)
ax_full.set_ylim(0.0, 90.0)
set_axes_style(ax_full, r"deviation $\Delta\theta$ (deg)", r"shock angle $\sigma$ (deg)")
plot_debug_mesh(
    ax_full_mesh,
    np.asarray(full_deviation_rows),
    np.asarray(full_sigma_rows),
    r"deviation $\Delta\theta$ (deg)",
    r"shock angle $\sigma$ (deg)",
    row_skip=1,
    column_skip=debug_mesh_column_skip,
)
ax_full.legend(loc="lower right")
fig_full.tight_layout()
fig_full.savefig("polar.png", bbox_inches="tight")


# =============================================================================
# SHARED SHOCK-POLAR DATA FOR FIGURES 2 TO 5
# =============================================================================
# Figures 2 and 3 share the deviation/shock-angle coordinates.  The dashed
# curves in Figure 2 include weak and strong solutions up to a normal shock;
# the solid curves in Figure 3 stop at the maximum-deflection point.
fig_ratio, ax_ratio, ax_ratio_mesh = make_figure(
    rf"Shock-to-isentropic compression ratio, $\gamma={gamma:.1f}$"
)


fig_weak, ax_weak, ax_weak_mesh = make_figure(
    rf"Weak shock-wave polars, $\gamma={gamma:.1f}$"
)
for mach in mach_values:
    sigma = np.linspace(deg.asin(1.0 / mach), sw.sigma_DevMax(mach, gamma), npoints + 1)
    deviation = sw.deflection_Mach_sigma(mach, sigma, gamma)
    ax_weak.plot(deviation, sigma, "k-", linewidth=1.0)
    # Figure 2 continues past maximum deviation along the strong branch to a
    # normal shock, whereas Figure 3 intentionally stops at the weak branch.
    complete_sigma = np.linspace(deg.asin(1.0 / mach), 90.0, npoints + 1)
    complete_deviation = sw.deflection_Mach_sigma(
        mach, complete_sigma, gamma
    )
    ax_ratio.plot(complete_deviation, complete_sigma, "k--", linewidth=1.0)
    ax_weak.annotate(
        rf"$M_0={mach:g}$", (deviation[-1], sigma[-1]),
        xytext=(0, 3), textcoords="offset points",
        ha="center", va="bottom", fontsize=7, rotation=90,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.0},
    )
    ax_ratio.annotate(
        rf"$M_0={mach:g}$", (deviation[-1], sigma[-1]),
        xytext=(0, 3), textcoords="offset points",
        ha="center", va="bottom", fontsize=7, rotation=90,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.0},
    )

ax_weak.plot(
    max_deviation, max_sigma, color="tab:green", linewidth=2.2,
    label="maximum-deviation locus",
)
ax_ratio.plot(
    max_deviation, max_sigma, color="tab:green", linewidth=2.2,
    label="maximum-deviation locus",
)

# Build a curvilinear mesh spanning all weak polars.  Include the
# Mach numbers drawn and labelled above exactly: besides making the contours
# smoother, this lets those M0 curves remain direct mesh lines rather than M0
# contours reconstructed by interpolation.
contour_mach_limits = (1.01, 10.0)
labelled_mesh_mach = mach_values[
    (mach_values >= contour_mach_limits[0])
    & (mach_values <= contour_mach_limits[1])
]
mesh_mach = np.unique(
    np.concatenate(
        (
            np.geomspace(*contour_mach_limits, npoints_mach),
            labelled_mesh_mach,
        )
    )
)
branch_fraction = np.linspace(0.0, 1.0, npoints + 1)
mach_mesh, fraction_mesh = np.meshgrid(mesh_mach, branch_fraction, indexing="ij")
mach_angle = deg.asin(1.0 / mach_mesh)
sigma_at_devmax = sw.sigma_DevMax(mach_mesh, gamma)
sigma_mesh = mach_angle + fraction_mesh * (sigma_at_devmax - mach_angle)
deviation_mesh = sw.deflection_Mach_sigma(mach_mesh, sigma_mesh, gamma)
downstream_mach_mesh = sw.downstream_Mach(mach_mesh, sigma_mesh, gamma)

# Complete mesh for Figures 2, 4, and 5.  Extending sigma to 90 degrees adds
# the strong solution beyond the maximum-deviation locus.
complete_sigma_mesh = mach_angle + fraction_mesh * (90.0 - mach_angle)
complete_deviation_mesh = sw.deflection_Mach_sigma(
    mach_mesh, complete_sigma_mesh, gamma
)
complete_downstream_mach_mesh = sw.downstream_Mach(
    mach_mesh, complete_sigma_mesh, gamma
)

# Compare the static-pressure rise through a shock with the pressure rise from
# an isentropic compression through the same deviation.  Points outside the
# isentropic-compression domain are masked from the contour plots.
compression_ratio = shock_to_isentropic_ratio(
    mach_mesh, sigma_mesh, deviation_mesh
)
complete_compression_ratio = shock_to_isentropic_ratio(
    mach_mesh, complete_sigma_mesh, complete_deviation_mesh
)
ratio_levels = (0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0, 1.02, 1.05, 1.1, 1.2, 1.5, 2.)

plot_debug_mesh(
    ax_ratio_mesh,
    complete_deviation_mesh,
    complete_sigma_mesh,
    r"deviation $\Delta\theta$ (deg)",
    r"shock angle $\sigma$ (deg)",
    row_skip=debug_mesh_row_skip,
    column_skip=debug_mesh_column_skip,
)
plot_debug_mesh(
    ax_weak_mesh,
    deviation_mesh,
    sigma_mesh,
    r"deviation $\Delta\theta$ (deg)",
    r"shock angle $\sigma$ (deg)",
    row_skip=debug_mesh_row_skip,
    column_skip=debug_mesh_column_skip,
)


# =============================================================================
# FIGURE 2 - SHOCK-TO-ISENTROPIC PRESSURE-RATIO MAP
# =============================================================================
# Fill ratios above unity, then draw labelled ratio isolines over complete
# dashed M0 curves.  The 0.95 contour highlights the transition close to equal
# pressure rise without crowding the lower-ratio portion of the map.
ratio_above_one = np.ma.masked_less_equal(complete_compression_ratio, 1.0)
ax_ratio.contourf(
    complete_deviation_mesh, complete_sigma_mesh, ratio_above_one,
    levels=(1.0, float(ratio_above_one.max())), colors=("tab:red",), alpha=0.18,
)
ratio_sigma_contours = ax_ratio.contour(
    complete_deviation_mesh, complete_sigma_mesh, complete_compression_ratio,
    levels=ratio_levels, colors="tab:red", linewidths=1.2,
)
ax_ratio.clabel(
    ratio_sigma_contours,
    fmt=lambda value: rf"$\Pi_{{shock}}/\Pi_{{is}}={value:g}$",
    fontsize=7,
)
ax_ratio.set_xlim(0.0, 50.0)
ax_ratio.set_ylim(0.0, 90.0)
set_axes_style(
    ax_ratio, r"deviation $\Delta\theta$ (deg)", r"shock angle $\sigma$ (deg)"
)
ratio_legend_handles, _ = ax_ratio.get_legend_handles_labels()
ratio_legend_handles.append(
    Patch(facecolor="tab:red", alpha=0.18, label=r"$\Pi_{shock}/\Pi_{is}>1$")
)
ax_ratio.legend(handles=ratio_legend_handles, loc="lower right")
fig_ratio.tight_layout()
fig_ratio.savefig("polar-ratio.png", bbox_inches="tight")


# =============================================================================
# FIGURE 3 - WEAK SHOCK POLARS WITH DOWNSTREAM-MACH CONTOURS
# =============================================================================
# Red dashed isolines show the downstream Mach number on the same weak-branch
# mesh used to construct Figure 2 before its strong-branch extension.
downstream_levels = (0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0)
mach2_contours = ax_weak.contour(
    deviation_mesh, sigma_mesh, downstream_mach_mesh,
    levels=downstream_levels, colors="tab:red", linestyles="--", linewidths=1.0,
)
ax_weak.clabel(mach2_contours, fmt=lambda value: rf"$M_2={value:g}$", fontsize=7)
ax_weak.set_xlim(0.0, 50.0)
ax_weak.set_ylim(0.0, 90.0)
set_axes_style(ax_weak, r"deviation $\Delta\theta$ (deg)", r"shock angle $\sigma$ (deg)")
ax_weak.legend(loc="lower right")
fig_weak.tight_layout()
fig_weak.savefig("polar-weak.png", bbox_inches="tight")


# =============================================================================
# FIGURES 4A/4B - SEPARATE WEAK AND STRONG SHOCK MAPS
# =============================================================================
# The weak mesh already spans the Mach angle to maximum deviation.  Construct
# a second structured mesh from maximum deviation to the normal-shock limit so
# neither plot contains the folded coordinates of the complete mesh.
strong_sigma_mesh = sigma_at_devmax + fraction_mesh * (
    90.0 - sigma_at_devmax
)
strong_deviation_mesh = sw.deflection_Mach_sigma(
    mach_mesh, strong_sigma_mesh, gamma
)
strong_downstream_mach_mesh = sw.downstream_Mach(
    mach_mesh, strong_sigma_mesh, gamma
)
strong_compression_ratio = shock_to_isentropic_ratio(
    mach_mesh, strong_sigma_mesh, strong_deviation_mesh
)

sigma_levels = (10.0, 15.0, 20.0, 25., 30.0, 35., 40., 45.0, 60.0, 75.0, 85.0, 90.0)


def plot_deviation_mach_map(
    title, deviation_grid, sigma_grid, downstream_mach_grid,
    pressure_ratio_grid, shade_ratio_above_one=False,
):
    """Plot one shock branch in deviation/upstream-Mach coordinates."""
    fig, ax, mesh_ax = make_figure(title)
    if shade_ratio_above_one:
        shaded_ratio = np.ma.masked_less_equal(pressure_ratio_grid, 1.0)
        ax.contourf(
            deviation_grid, mach_mesh, shaded_ratio,
            levels=(1.0, float(shaded_ratio.max())),
            colors=("tab:red",), alpha=0.18,
        )
    sigma_contours = ax.contour(
        deviation_grid, mach_mesh, sigma_grid,
        levels=sigma_levels, colors="black", linewidths=1.0,
    )
    ax.clabel(
        sigma_contours, fmt=lambda value: rf"$\sigma={value:g}^\circ$", fontsize=7
    )
    mach2_contours = ax.contour(
        deviation_grid, mach_mesh, downstream_mach_grid,
        levels=downstream_levels, colors="tab:red",
        linestyles="--", linewidths=1.0,
    )
    ax.clabel(mach2_contours, fmt=lambda value: rf"$M_2={value:g}$", fontsize=7)
    pressure_ratio_contours = ax.contour(
        deviation_grid, mach_mesh, pressure_ratio_grid,
        levels=ratio_levels, colors="tab:purple", linestyles=":", linewidths=1.1,
    )
    ax.clabel(
        pressure_ratio_contours,
        fmt=lambda value: rf"$\Pi_{{shock}}/\Pi_{{is}}={value:g}$",
        fontsize=7,
    )
    ax.plot(
        sw.dev_Max(mesh_mach, gamma), mesh_mach,
        color="tab:blue", linewidth=2.5, label="maximum-deviation boundary",
    )
    ax.plot(
        sup.PrandtlMeyer_Mach(mesh_mach, gamma), mesh_mach,
        color="tab:green", linestyle="-.", linewidth=2.0,
        label="isentropic-compression admissibility limit",
    )
    figure4_mach_max = 4.0
    ax.set_xlim(0.0, 40.0)
    ax.set_ylim(mesh_mach[0], figure4_mach_max)
    set_axes_style(
        ax, r"deviation $\Delta\theta$ (deg)", r"upstream Mach number $M_0$"
    )
    plot_debug_mesh(
        mesh_ax, deviation_grid, mach_mesh,
        r"deviation $\Delta\theta$ (deg)", r"upstream Mach number $M_0$",
        row_skip=debug_mesh_row_skip,
        column_skip=debug_mesh_column_skip,
    )
    if mesh_ax is not None:
        mesh_ax.set_ylim(mesh_mach[0], figure4_mach_max)
    legend_handles, _ = ax.get_legend_handles_labels()
    if shade_ratio_above_one:
        legend_handles.append(
            Patch(facecolor="tab:red", alpha=0.18,
                  label=r"$\Pi_{shock}/\Pi_{is}>1$")
        )
    ax.legend(handles=legend_handles)
    fig.tight_layout()
    return fig


# Figure 4a: weak branch, with the ratio-above-one domain highlighted.
fig_mach_weak = plot_deviation_mach_map(
    rf"Weak shock map, $\gamma={gamma:.1f}$",
    deviation_mesh, sigma_mesh, downstream_mach_mesh, compression_ratio,
    shade_ratio_above_one=True,
)
fig_mach_weak.savefig("polar-mach.png", bbox_inches="tight")

# Figure 4b: strong branch, left unshaded to facilitate isoline-crossing checks.
fig_mach_strong = plot_deviation_mach_map(
    rf"Strong shock map, $\gamma={gamma:.1f}$",
    strong_deviation_mesh, strong_sigma_mesh, strong_downstream_mach_mesh,
    strong_compression_ratio,
)
fig_mach_strong.savefig("polar-mach-strong.png", bbox_inches="tight")


# =============================================================================
# FIGURE 6 - UPSTREAM / DOWNSTREAM MACH MAP
# =============================================================================
# The final projection uses M0 and M2 as coordinates.  Shock-angle and
# deviation isolines describe the same complete weak/strong data used above.
fig_mach2, ax_mach2, ax_mach2_mesh = make_figure(
    rf"Upstream/downstream Mach map, $\gamma={gamma:.1f}$"
)
sigma_mach_contours = ax_mach2.contour(
    mach_mesh, complete_downstream_mach_mesh, complete_sigma_mesh,
    levels=sigma_levels, colors="black", linewidths=1.0,
)
ax_mach2.clabel(
    sigma_mach_contours, fmt=lambda value: rf"$\sigma={value:g}^\circ$", fontsize=7
)
deviation_levels = (2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0)
deviation_contours = ax_mach2.contour(
    mach_mesh, complete_downstream_mach_mesh, complete_deviation_mesh,
    levels=deviation_levels, colors="tab:red", linestyles="--", linewidths=1.0,
)
ax_mach2.clabel(
    deviation_contours,
    fmt=lambda value: rf"$\Delta\theta={value:g}^\circ$", fontsize=7,
)
ax_mach2.plot(mesh_mach, mesh_mach, color="0.5", linestyle=":", label="no shock")
ax_mach2.plot(
    mesh_mach,
    sw.downstream_Mach(mesh_mach, sw.sigma_DevMax(mesh_mach, gamma), gamma),
    color="tab:blue", linewidth=2.0, label="maximum-deviation boundary",
)
ax_mach2.set_xlim(mesh_mach[0], figure6_mach_max)
ax_mach2.set_ylim(
    max(.5, np.nanmin(complete_downstream_mach_mesh)),
    min(4.0, np.nanmax(complete_downstream_mach_mesh)),
)
#ax_mach2.set_xscale("log")
#ax_mach2.set_yscale("log")
set_axes_style(
    ax_mach2, r"upstream Mach number $M_0$", r"downstream Mach number $M_2$"
)
plot_debug_mesh(
    ax_mach2_mesh,
    mach_mesh,
    complete_downstream_mach_mesh,
    r"upstream Mach number $M_0$",
    r"downstream Mach number $M_2$",
    row_skip=debug_mesh_row_skip,
    column_skip=debug_mesh_column_skip,
)
if ax_mach2_mesh is not None:
    ax_mach2_mesh.set_xlim(mesh_mach[0], figure6_mach_max)
ax_mach2.legend(loc="upper left")
fig_mach2.tight_layout()
fig_mach2.savefig("polar-mach-mach.png", bbox_inches="tight")

plt.show()
