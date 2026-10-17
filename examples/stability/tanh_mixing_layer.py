"""Temporal Orr--Sommerfeld spectrum of a centred hyperbolic-tangent shear layer."""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

# Allow ``python examples/stability/tanh_mixing_layer.py`` from a checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aerokit.stability.OrrSommerfeld import OrrSommerfeldModel
from aerokit.common.mapping import CenterClusteredMapping
from aerokit.common.numspectral import ChebCollocation, MultiDomainChebCollocation


def tanh_profile(y, uavg, delta_u, thickness):
    """Return ``uavg + delta_u / 2 * tanh(y / thickness)``."""
    return uavg + 0.5 * delta_u * np.tanh(y / thickness)


def spectrum(diffop, alpha, reynolds, uavg, delta_u, thickness):
    """Compute the finite temporal Orr--Sommerfeld eigenvalues."""
    def profile(y):
        return tanh_profile(y, uavg, delta_u, thickness)

    model = OrrSommerfeldModel(diffop.npts, diffop=diffop)
    model.set_basestate({"alpha": alpha, "Reynolds": reynolds, "uprofile": profile})
    model.set_BC({"type": "wall"}, {"type": "wall"})
    omega, modes = model.solve_eig(finite_only=True)
    return model, omega, modes, profile


def most_unstable(omega):
    """Return the temporal mode with the greatest growth rate."""
    return omega[np.argmax(omega.imag)]


def most_unstable_mode(omega, modes):
    """Return the most amplified eigenvalue and its normalized eigenfunction."""
    index = np.argmax(omega.imag)
    mode = modes[:, index]
    imax = np.argmax(np.abs(mode))
    return omega[index], mode * np.exp(-1j * np.angle(mode[imax])) / np.abs(mode[imax])


def richardson_analysis(values, sizes):
    """Return local orders and a Richardson-extrapolated finest-grid value."""
    values = np.asarray(values, dtype=float)
    steps = 1.0 / (np.asarray(sizes, dtype=float) - 1.0)
    refinement = steps[:-1] / steps[1:]
    if not np.allclose(refinement, refinement[0]):
        raise ValueError("Richardson analysis requires a constant refinement ratio")

    differences = np.abs(np.diff(values))
    orders = np.full(values.size - 2, np.nan)
    valid = (differences[:-1] > 0.0) & (differences[1:] > 0.0)
    orders[valid] = np.log(differences[:-1][valid] / differences[1:][valid]) / np.log(refinement[0])

    order = orders[-1]
    if not np.isfinite(order):
        raise ValueError("could not estimate the final local convergence order")
    factor = refinement[-1] ** order
    reference = values[-1] + (values[-1] - values[-2]) / (factor - 1.0)
    return orders, reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npts", type=int, default=121, help="number of mono-domain Chebyshev points")
    parser.add_argument("--alpha", type=float, default=4.446, help="real streamwise wavenumber")
    parser.add_argument("--reynolds", type=float, default=1.0e5, help="Reynolds number")
    parser.add_argument("--uavg", type=float, default=0.5, help="mean velocity")
    parser.add_argument("--delta-u", type=float, default=2.0, help="velocity difference across the layer")
    parser.add_argument("--thickness", type=float, default=0.1, help="tanh thickness")
    parser.add_argument(
        "--center-strength",
        type=float,
        default=0.95,
        help="strength of the centre-clustered mono-domain mapping",
    )
    parser.add_argument("--convergence", action="store_true", help="plot sensitivity to the number of points")
    parser.add_argument("--save", help="optional output figure path")
    parser.add_argument("--no-show", action="store_true", help="do not open the figure window")
    args = parser.parse_args()

    if args.thickness <= 0.0:
        parser.error("--thickness must be positive")
    if not 0.0 <= args.center_strength < 1.0:
        parser.error("--center-strength must be in [0, 1)")

    mono = ChebCollocation(args.npts, xmin=-2.0, xmax=2.0)
    centered = ChebCollocation(
        args.npts, mapping=CenterClusteredMapping(-2.0, 2.0, args.center_strength)
    )
    npts_per_domain = (args.npts + 1) // 2
    multi = MultiDomainChebCollocation(npts_per_domain, xmin=-2.0, interface=0.0, xmax=2.0)
    _, omega, modes, profile = spectrum(
        mono, args.alpha, args.reynolds, args.uavg, args.delta_u, args.thickness
    )
    _, omega_multi, modes_multi, _ = spectrum(
        multi, args.alpha, args.reynolds, args.uavg, args.delta_u, args.thickness
    )
    _, omega_centered, modes_centered, _ = spectrum(
        centered, args.alpha, args.reynolds, args.uavg, args.delta_u, args.thickness
    )
    omega_max, mode = most_unstable_mode(omega, modes)
    omega_multi_max, mode_multi = most_unstable_mode(omega_multi, modes_multi)
    omega_centered_max, mode_centered = most_unstable_mode(omega_centered, modes_centered)
    print(f"Computed {omega.size} finite mono-domain temporal eigenvalues.")
    print(f"Most unstable: omega = {omega_max.real:.6g} {omega_max.imag:+.6g}i")
    print(f"Computed {omega_multi.size} finite multi-domain temporal eigenvalues.")
    print(f"Most unstable: omega = {omega_multi_max.real:.6g} {omega_multi_max.imag:+.6g}i")
    print(f"Computed {omega_centered.size} finite centre-clustered temporal eigenvalues.")
    print(f"Most unstable: omega = {omega_centered_max.real:.6g} {omega_centered_max.imag:+.6g}i")

    y = np.linspace(-2.0, 2.0, 1000)
    fig, (ax_profile, ax_spectrum, ax_nodes) = plt.subplots(1, 3, figsize=(15, 4.5))
    ax_profile.plot(profile(y), y)
    ax_profile.set(
        xlabel=r"$U$",
        ylabel=r"$y$",
        title=rf"$U={args.uavg:g}+{args.delta_u / 2:g}\tanh(y/{args.thickness:g})$",
    )
    ax_profile.grid()

    ax_spectrum.plot(omega.real, omega.imag, ".", alpha=0.5, label="mono-domain")
    ax_spectrum.plot(omega_multi.real, omega_multi.imag, "+", alpha=0.5, label="two domains")
    ax_spectrum.plot(omega_centered.real, omega_centered.imag, "x", alpha=0.5, label="centre-clustered")
    ax_spectrum.plot(omega_max.real, omega_max.imag, "ro", label="most unstable")
    ax_spectrum.set(xlabel=r"$\Re(\omega)$", ylabel=r"$\Im(\omega)$", title="Temporal spectrum")
    ax_spectrum.grid()
    ax_spectrum.set_xlim(-2, 4)
    ax_spectrum.set_ylim(-2, 2)
    ax_spectrum.legend()

    ax_nodes.plot(mono.x, np.zeros_like(mono.x), "o", label="mono-domain")
    ax_nodes.plot(multi.x, np.ones_like(multi.x), "+", label="two domains")
    ax_nodes.plot(centered.x, 2.0 * np.ones_like(centered.x), "x", label="centre-clustered")
    ax_nodes.axvline(0.0, color="k", linewidth=0.8)
    ax_nodes.set(
        xlabel=r"$y$",
        yticks=[0, 1, 2],
        yticklabels=["mono", "multi", "centred"],
        title="Collocation points",
    )
    ax_nodes.grid()
    ax_nodes.legend()
    fig.tight_layout()

    fig_mode, (ax_mode_real, ax_mode_imag) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    ax_mode_real.plot(mode.real, mono.x, label="mono-domain")
    ax_mode_real.plot(mode_multi.real, multi.x, "--", label="two domains")
    ax_mode_real.plot(mode_centered.real, centered.x, ":", label="centre-clustered")
    ax_mode_real.axhline(0.0, color="k", linewidth=0.8)
    ax_mode_real.set(xlabel=r"$\Re(\hat v)/|\hat v|_{max}$", ylabel=r"$y$", title="Most unstable eigenfunction")
    ax_mode_real.grid()
    ax_mode_real.legend()
    ax_mode_imag.plot(mode.imag, mono.x, label="mono-domain")
    ax_mode_imag.plot(mode_multi.imag, multi.x, "--", label="two domains")
    ax_mode_imag.plot(mode_centered.imag, centered.x, ":", label="centre-clustered")
    ax_mode_imag.axhline(0.0, color="k", linewidth=0.8)
    ax_mode_imag.set(xlabel=r"$\Im(\hat v)/|\hat v|_{max}$", title="Phase-aligned eigenfunction")
    ax_mode_imag.grid()
    fig_mode.tight_layout()

    fig_spacing, ax_spacing = plt.subplots(figsize=(7, 4.5))
    for grid, label, marker in (
        (mono, "mono-domain", "o-"),
        (centered, "centre-clustered", "^-"),
        (multi, "two domains", "s-"),
    ):
        spacing = np.diff(grid.x)
        locations = 0.5 * (grid.x[:-1] + grid.x[1:])
        ax_spacing.semilogy(locations, spacing, marker, label=label)
    ax_spacing.axvline(0.0, color="k", linewidth=0.8)
    ax_spacing.set(
        xlabel=r"$y$",
        ylabel=r"local spacing $\Delta y$",
        title="Collocation-point spacing",
    )
    ax_spacing.grid()
    ax_spacing.legend()
    fig_spacing.tight_layout()

    if args.convergence:
        # n - 1 doubles at each refinement, as required by Richardson analysis.
        sizes = np.array([21, 41, 81, 161, 321])
        mono_growth = []
        centered_growth = []
        multi_growth = []
        multi_sizes = []
        for npts in sizes:
            _, values, _, _ = spectrum(
                ChebCollocation(npts, xmin=-2.0, xmax=2.0),
                args.alpha,
                args.reynolds,
                args.uavg,
                args.delta_u,
                args.thickness,
            )
            multidomain = MultiDomainChebCollocation((npts + 1) // 2, xmin=-2.0, interface=0.0, xmax=2.0)
            _, values_multi, _, _ = spectrum(
                multidomain, args.alpha, args.reynolds, args.uavg, args.delta_u, args.thickness
            )
            _, values_centered, _, _ = spectrum(
                ChebCollocation(npts, mapping=CenterClusteredMapping(-2.0, 2.0, args.center_strength)),
                args.alpha,
                args.reynolds,
                args.uavg,
                args.delta_u,
                args.thickness,
            )
            mono_growth.append(most_unstable(values).imag)
            multi_growth.append(most_unstable(values_multi).imag)
            centered_growth.append(most_unstable(values_centered).imag)
            multi_sizes.append(multidomain.npts)

        mono_orders, _ = richardson_analysis(mono_growth, sizes)
        centered_orders, _ = richardson_analysis(centered_growth, sizes)
        multi_orders, reference_growth = richardson_analysis(multi_growth, multi_sizes)
        print(f"Richardson reference amplification rate: {reference_growth:.12g}")
        print(
            "Final local order: "
            f"mono={mono_orders[-1]:.3g}, centre={centered_orders[-1]:.3g}, multi={multi_orders[-1]:.3g}"
        )

        fig_convergence, (ax_convergence, ax_error, ax_order) = plt.subplots(1, 3, figsize=(16, 4))
        ax_convergence.plot(sizes, mono_growth, "o-", label="mono-domain")
        ax_convergence.plot(sizes, centered_growth, "^-", label="centre-clustered")
        ax_convergence.plot(multi_sizes, multi_growth, "s-", label="two domains")
        ax_convergence.axhline(
            reference_growth, color="k", linestyle=":", label="Richardson reference (multi)"
        )
        ax_convergence.set(
            xlabel="total collocation points",
            ylabel=r"$\max\,\Im(\omega)$",
            title="Temporal amplification rate",
        )
        ax_convergence.grid()
        ax_convergence.legend()
        ax_error.semilogy(sizes, np.abs(np.asarray(mono_growth) - reference_growth), "o-", label="mono-domain")
        ax_error.semilogy(
            sizes, np.abs(np.asarray(centered_growth) - reference_growth), "^-", label="centre-clustered"
        )
        ax_error.semilogy(multi_sizes, np.abs(np.asarray(multi_growth) - reference_growth), "s-", label="two domains")
        ax_error.set(
            xlabel="total collocation points",
            ylabel=r"$|\Im(\omega)-\Im(\omega)_{ref}|$",
            title="Amplification-rate convergence",
        )
        ax_error.grid()
        ax_error.legend()
        ax_order.plot(sizes[1:-1], mono_orders, "o-", label="mono-domain")
        ax_order.plot(sizes[1:-1], centered_orders, "^-", label="centre-clustered")
        ax_order.plot(multi_sizes[1:-1], multi_orders, "s-", label="two domains")
        ax_order.set(
            xlabel="total collocation points",
            ylabel="local order",
            title="Estimated local convergence order",
        )
        ax_order.grid()
        ax_order.legend()
        fig_convergence.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150)
        output = Path(args.save)
        fig_mode.savefig(output.with_name(f"{output.stem}_mode{output.suffix}"), dpi=150)
        fig_spacing.savefig(output.with_name(f"{output.stem}_spacing{output.suffix}"), dpi=150)
        if args.convergence:
            fig_convergence.savefig(output.with_name(f"{output.stem}_convergence{output.suffix}"), dpi=150)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
