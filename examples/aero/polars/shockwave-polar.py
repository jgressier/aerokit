"""Shock-wave polar maps for a calorically perfect gas.

Five figures are produced: the complete deviation/shock-angle polars; a
deviation/shock-angle map of shock-to-isentropic compression; the weak shock
branches; a deviation/upstream-Mach map; and an upstream/downstream-Mach map.
"""

import matplotlib.pyplot as plt
import numpy as np

import aerokit.aero.ShockWave as sw
import aerokit.aero.Isentropic as Is
import aerokit.aero.Supersonic as sup
import aerokit.aero.degree as deg


npoints = 160
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


# Complete weak and strong shock polars.
fig_full, ax_full = plt.subplots(figsize=(10, 8))
ax_full.set_title(rf"Complete shock-wave polars, $\gamma={gamma:.1f}$")
for mach in mach_values:
    sigma = np.linspace(deg.asin(1.0 / mach), 90.0, npoints + 1)
    deviation = sw.deflection_Mach_sigma(mach, sigma, gamma)
    ax_full.plot(deviation, sigma, "k-")
    label_index = 3 * npoints // 4
    ax_full.text(
        deviation[label_index], sigma[label_index], rf"$M_0={mach:g}$",
        ha="left", va="bottom", fontsize=7, rotation=-25,
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
ax_full.legend(loc="lower right")
fig_full.tight_layout()
fig_full.savefig("polar.png", bbox_inches="tight")


# Compression-ratio map in the same coordinates as a shock polar.
fig_ratio, ax_ratio = plt.subplots(figsize=(10, 8))
ax_ratio.set_title(rf"Shock-to-isentropic compression ratio, $\gamma={gamma:.1f}$")


# Weak branches only: sigma <= sigma at maximum deviation.
fig_weak, ax_weak = plt.subplots(figsize=(10, 8))
ax_weak.set_title(rf"Weak shock-wave polars, $\gamma={gamma:.1f}$")
for mach in mach_values:
    sigma = np.linspace(deg.asin(1.0 / mach), sw.sigma_DevMax(mach, gamma), npoints + 1)
    deviation = sw.deflection_Mach_sigma(mach, sigma, gamma)
    ax_weak.plot(deviation, sigma, "k-", linewidth=1.0)
    ax_ratio.plot(deviation, sigma, "k-", linewidth=1.0)
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

# Curvilinear mesh spanning all weak polars, used for contouring.
mesh_mach = np.geomspace(1.01, 10.0, 180)
branch_fraction = np.linspace(0.0, 1.0, npoints + 1)
mach_mesh, fraction_mesh = np.meshgrid(mesh_mach, branch_fraction, indexing="ij")
mach_angle = deg.asin(1.0 / mach_mesh)
sigma_at_devmax = sw.sigma_DevMax(mach_mesh, gamma)
sigma_mesh = mach_angle + fraction_mesh * (sigma_at_devmax - mach_angle)
deviation_mesh = sw.deflection_Mach_sigma(mach_mesh, sigma_mesh, gamma)
downstream_mach_mesh = sw.downstream_Mach(mach_mesh, sigma_mesh, gamma)

shock_pressure_mesh = sw.Ps_ratio(mach_mesh * deg.sin(sigma_mesh), gamma)
isentropic_pressure_mesh = np.full_like(shock_pressure_mesh, np.nan)
for row, mach in enumerate(mesh_mach):
    max_characteristic_deviation = sup.PrandtlMeyer_Mach(mach, gamma)
    valid = deviation_mesh[row] <= max_characteristic_deviation + 1.0e-10
    remaining_pm_angle = np.maximum(
        max_characteristic_deviation - deviation_mesh[row, valid], 0.0
    )
    isentropic_mach = np.ones_like(remaining_pm_angle)
    supersonic = remaining_pm_angle > 1.0e-10
    if np.any(supersonic):
        isentropic_mach[supersonic] = sup.Mach_PrandtlMeyer(
            remaining_pm_angle[supersonic], gamma
        )
    isentropic_pressure_mesh[row, valid] = (
        Is.PtPs_Mach(mach, gamma) / Is.PtPs_Mach(isentropic_mach, gamma)
    )

compression_ratio = np.ma.masked_invalid(shock_pressure_mesh / isentropic_pressure_mesh)
ratio_levels = (0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.02, 1.05, 1.1)

ratio_sigma_contours = ax_ratio.contour(
    deviation_mesh, sigma_mesh, compression_ratio,
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
ax_ratio.legend(loc="lower right")
fig_ratio.tight_layout()
fig_ratio.savefig("polar-ratio.png", bbox_inches="tight")

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


# Weak-shock map in deviation/upstream-Mach coordinates.
fig_mach, ax_mach = plt.subplots(figsize=(10, 8))
ax_mach.set_title(rf"Weak shock map, $\gamma={gamma:.1f}$")
sigma_levels = (10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 75.0)
sigma_contours = ax_mach.contour(
    deviation_mesh, mach_mesh, sigma_mesh,
    levels=sigma_levels, colors="black", linewidths=1.0,
)
ax_mach.clabel(sigma_contours, fmt=lambda value: rf"$\sigma={value:g}^\circ$", fontsize=7)
mach2_map_contours = ax_mach.contour(
    deviation_mesh, mach_mesh, downstream_mach_mesh,
    levels=downstream_levels, colors="tab:red", linestyles="--", linewidths=1.0,
)
ax_mach.clabel(mach2_map_contours, fmt=lambda value: rf"$M_2={value:g}$", fontsize=7)
ax_mach.plot(
    sw.dev_Max(mesh_mach, gamma), mesh_mach,
    color="tab:blue", linewidth=2.5, label="maximum-deviation boundary",
)
ax_mach.set_xlim(0.0, 50.0)
ax_mach.set_ylim(mesh_mach[0], mesh_mach[-1])
#ax_mach.set_yscale("log")
set_axes_style(ax_mach, r"deviation $\Delta\theta$ (deg)", r"upstream Mach number $M_0$")
ax_mach.legend()
fig_mach.tight_layout()
fig_mach.savefig("polar-mach.png", bbox_inches="tight")


# Weak-shock map in upstream/downstream-Mach coordinates.
fig_mach2, ax_mach2 = plt.subplots(figsize=(10, 8))
ax_mach2.set_title(rf"Upstream/downstream Mach map, $\gamma={gamma:.1f}$")
sigma_mach_contours = ax_mach2.contour(
    mach_mesh, downstream_mach_mesh, sigma_mesh,
    levels=sigma_levels, colors="black", linewidths=1.0,
)
ax_mach2.clabel(
    sigma_mach_contours, fmt=lambda value: rf"$\sigma={value:g}^\circ$", fontsize=7
)
deviation_levels = (2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0)
deviation_contours = ax_mach2.contour(
    mach_mesh, downstream_mach_mesh, deviation_mesh,
    levels=deviation_levels, colors="tab:red", linestyles="--", linewidths=1.0,
)
ax_mach2.clabel(
    deviation_contours,
    fmt=lambda value: rf"$\Delta\theta={value:g}^\circ$", fontsize=7,
)
ax_mach2.plot(mesh_mach, mesh_mach, color="0.5", linestyle=":", label="no shock")
ax_mach2.set_xlim(mesh_mach[0], mesh_mach[-1])
ax_mach2.set_ylim(np.nanmin(downstream_mach_mesh), min(4., np.nanmax(downstream_mach_mesh)))
#ax_mach2.set_xscale("log")
#ax_mach2.set_yscale("log")
set_axes_style(
    ax_mach2, r"upstream Mach number $M_0$", r"downstream Mach number $M_2$"
)
ax_mach2.legend(loc="upper left")
fig_mach2.tight_layout()
fig_mach2.savefig("polar-mach-mach.png", bbox_inches="tight")

plt.show()
