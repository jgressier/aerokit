#!/usr/bin/env python3
"""Temporal round-jet spectrum using :class:`aerokit.stability.NS.NSaxi`.

This is the temporal LSA where the axial
wavenumber ``kx`` is prescribed and the generalized eigenproblem is solved
for ``omega``. 
"""

#from __future__ import annotations

import argparse
from pathlib import Path
import sys

from matplotlib.pyplot import xlim
import numpy as np
from scipy.optimize import linear_sum_assignment

# Allow ``python examples/stability/temporal_round_jet_nsaxi.py`` from a checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aerokit.stability.NS import NSaxi
from aerokit.common.mapping import SemiInfiniteAlgebraicMapping


def round_jet_basestate(r, kx, m, mach, r_theta, gamma):
    """Return the Crocco--Busemann round-jet base state used by the MATLAB main."""
    velocity = 0.5 * (1.0 + np.tanh(r_theta / 4.0 * (0.5 / np.maximum(r, 1e-14) - r / 0.5)))
    velocity *= mach
    t_inf_tj = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    temperature = t_inf_tj + (1.0 - t_inf_tj) * velocity + 0.5 * (gamma - 1.0) * mach**2 * velocity * (1.0 - velocity)
    pressure = np.ones_like(r)
    density = pressure / temperature
    # print(f"r: {r.min()} - {r.max()}")
    return {
        "kx": kx,
        "m": m,
        "rho": density,
        "Ux": velocity,
        "P": pressure,
        "gamma": gamma,
    }


def solve_temporal_spectrum(npts=101, kx=1.0, m=1, mach=1.5, r_theta=5.0, gamma=1.4, mapping_scale=10.0):
    """Build the ``NSaxi`` operator and return its finite temporal spectrum."""
    mapping = SemiInfiniteAlgebraicMapping(mapping_scale)
    model = NSaxi(npts, basestate=None, mapping=mapping)
    model.set_basestate(round_jet_basestate(model.r, kx, m, mach, r_theta, gamma))
    omega, modes = model.solve_eig()
    finite = np.isfinite(omega)
    return model, omega[finite], modes[:, finite]


def select_modes(omega, modes, nmodes=1, mode="guided", reference=None):
    """Select positive-real-frequency eigenpairs.

    ``mode='growth'`` selects the modes with the largest growth rates.
    ``mode='continuation'`` matches modes to ``reference`` in the complex
    frequency plane.  Returned modes are always ordered by increasing real
    frequency.
    """
    if nmodes < 1:
        raise ValueError("nmodes must be positive")
    if mode == "growth":
        indices = np.argsort(omega.imag)[::-1][:nmodes]
    elif mode == "guided":
        positive = omega.real > 1.e-3
        if positive.sum() < nmodes:
            raise ValueError(
                f"requested {nmodes} modes but only {positive.sum()} have positive real frequency"
            )
        omega = omega[positive]
        modes = modes[:, positive]
        indices = np.argsort(omega.real)[:][:nmodes]
    elif mode == "continuation":
        if reference is None:
            raise ValueError("mode='continuation' requires reference frequencies")
        reference = np.asarray(reference)[:nmodes]
        cost = np.abs(reference[:, np.newaxis] - omega[np.newaxis, :])
        rows, columns = linear_sum_assignment(cost)
        indices = np.empty(nmodes, dtype=int)
        indices[rows] = columns
    else:
        raise ValueError("mode must be 'guided', 'growth' or 'continuation'")
    return omega[indices], modes[:, indices]


def plot_modes(model, omegas, modes):
    """Plot the selected modes as rows and primitive variables as columns."""
    import matplotlib.pyplot as plt

    finite_r = np.isfinite(model.r)
    labels = (r"$|\rho'|$", r"$|u_x'|$", r"$|u_r'|$", r"$|u_\theta'|$", r"$|rT'|$")
    fig, axes = plt.subplots(len(omegas), model.nvar, figsize=(16, 2.8 * len(omegas)), sharex=True, squeeze=False)
    for imode, omega in enumerate(omegas):
        for ivar, (axis, label) in enumerate(zip(axes[imode], labels)):
            component = modes[ivar * model.dim : (ivar + 1) * model.dim, imode]
            axis.plot(model.r[finite_r], np.abs(component[finite_r]))
            axis.set(xlabel=r"$r$", ylabel=label, xlim=(0, 4.0))
            axis.grid(True)
            if imode == 0:
                axis.set_title(label)
        axes[imode, 0].annotate(
            rf"$\omega={omega.real:.4g}{omega.imag:+.4g}i$",
            xy=(0, 0.95),
            xycoords="axes fraction",
            va="top",
        )
    fig.suptitle(f"Selected temporal modes $k_x={model._basestate['kx']}$")
    fig.tight_layout()


def animate_spectrum(k_values, spectra, selected_spectra):
    """Animate the full spectrum over a streamwise-wavenumber scan."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    all_omega = np.concatenate(spectra)
    real_min, real_max = np.min(all_omega.real), np.max(all_omega.real)
    imag_min, imag_max = np.min(all_omega.imag), np.max(all_omega.imag)
    real_pad = max(1.0, 0.05 * (real_max - real_min))
    imag_pad = max(1.0, 0.05 * (imag_max - imag_min))

    fig, axis = plt.subplots(figsize=(6, 5))
    spectrum, = axis.plot([], [], "o", alpha=0.35, label="Spectrum")
    selected, = axis.plot([], [], "ro", label="Selected modes")
    axis.axhline(0.0, color="k", linewidth=0.8)
    axis.set(
        # xlim=(real_min - real_pad, real_max + real_pad),
        # ylim=(imag_min - imag_pad, imag_max + imag_pad),
        xlim=(-20, 20),
        ylim=(-10, 10),
        xlabel=r"$\Re(\omega)$",
        ylabel=r"$\Im(\omega)$",
    )
    axis.legend()

    def update(frame):
        omega = spectra[frame]
        spectrum.set_data(omega.real, omega.imag)
        selected.set_data(selected_spectra[frame].real, selected_spectra[frame].imag)
        axis.set_title(rf"Temporal spectrum: $k_x={k_values[frame]:.4g}$")
        return spectrum, selected

    return FuncAnimation(fig, update, frames=len(k_values), interval=150, blit=False, repeat=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npts", type=int, default=101)
    parser.add_argument("--mapping-scale", type=float, default=1.0, help="L in the semi-infinite algebraic map")
    parser.add_argument("--kx", type=float, default=-4.0, help="prescribed axial wavenumber")
    parser.add_argument("--m", type=int, default=0, help="azimuthal mode number")
    parser.add_argument("--mach", type=float, default=1.5)
    parser.add_argument("--r-theta", type=float, default=5.0)
    parser.add_argument("--scan-k", action="store_true", help="scan kx from zero to --kx")
    parser.add_argument("--nk", type=int, default=41, help="number of kx values for --scan-k")
    parser.add_argument("--nmodes", type=int, default=2, help="number of most unstable modes to plot")
    parser.add_argument("--animation-output", help="save the scan animation (for example, spectrum.gif or spectrum.mp4)")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    if args.animation_output and not args.scan_k:
        parser.error("--animation-output requires --scan-k")

    model, omega, modes = solve_temporal_spectrum(
        npts=args.npts,
        kx=args.kx,
        m=args.m,
        mach=args.mach,
        r_theta=args.r_theta,
        mapping_scale=args.mapping_scale,
    )
    print(f"Computed {omega.size} finite temporal eigenvalues with NSaxi.")
    initial_omega_selected, initial_modes_selected = select_modes(omega, modes, args.nmodes)
    for rank, omega_mode in enumerate(initial_omega_selected, start=1):
        print(f"Mode {rank}: omega = {omega_mode.real:.6g} {omega_mode.imag:+.6g}i")

    k_values = None
    omega_selected = None
    spectra = None
    spectrum_animation = None
    if args.scan_k:
        k_values = np.linspace(0.0, args.kx, args.nk)
        omega_selected = np.full((k_values.size, args.nmodes), np.nan + 1j * np.nan)
        spectra = []
        reference = None
        for i, kx in enumerate(k_values):
            _, omega_k, modes_k = solve_temporal_spectrum(
                npts=args.npts,
                kx=kx,
                m=args.m,
                mach=args.mach,
                r_theta=args.r_theta,
                mapping_scale=args.mapping_scale,
            )
            spectra.append(omega_k)
            selection_mode = 'guided'
            omega_k_selected, _ = select_modes(
                omega_k, modes_k, args.nmodes, mode=selection_mode, reference=reference
            )
            omega_selected[i, :omega_k_selected.size] = omega_k_selected
            reference = omega_k_selected
        print(f"Maximum growth among selected branches: {omega_selected.imag.max():.6g}")

    if args.animation_output:
        spectrum_animation = animate_spectrum(k_values, spectra, omega_selected)
        print(f"Saving animation to {args.animation_output}")
        spectrum_animation.save(args.animation_output)

    if not args.no_plot:
        import matplotlib.pyplot as plt

        finite_r = np.isfinite(model.r)
        state = model._basestate
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), sharex=True)
        for axis, values, label in zip(
            axes,
            (state["Ux"], state["rho"], state["P"]),
            (r"$U_x$", r"$\rho$", r"$P$"),
        ):
            axis.plot(model.r[finite_r], values[finite_r], 'o', ms=2)
            axis.set(xlabel=r"$r$", ylabel=label)
            axis.set_xlim(0, 4.)
            axis.grid(True)
        fig.suptitle("Temporal round-jet base state")
        fig.tight_layout()

        plt.figure()
        plt.plot(omega.real, omega.imag, "o", alpha=0.5)
        plt.axhline(0.0, color="k", linewidth=0.8)
        plt.xlabel(r"$\Re(\omega)$")
        plt.ylabel(r"$\Im(\omega)$")
        plt.xlim(-10, 10.)
        plt.title(f"NSaxi temporal spectrum: $k_x={args.kx}$, $m={args.m}$")
        plot_modes(model, initial_omega_selected, initial_modes_selected)

        if args.scan_k:
            fig, (ax_real, ax_imag) = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
            for i in range(args.nmodes):
                ax_real.plot(k_values, omega_selected[:, i].real, "o", label=f"Mode {i+1}", ms=3)
                ax_imag.plot(k_values, omega_selected[:, i].imag, "o", label=f"Mode {i+1}", ms=3)
            ax_real.set(xlabel=r"$k_x$", ylabel=r"$\Re(\omega)$", title="Frequency branches")
            #ax_real.set_ylim(0, 10.)
            ax_imag.set(xlabel=r"$k_x$", ylabel=r"$\Im(\omega)$", title="Growth-rate branches")
            ax_real.grid(True)
            ax_imag.grid(True)
            ax_imag.legend()
            fig.tight_layout()
            if spectrum_animation is None:
                spectrum_animation = animate_spectrum(k_values, spectra, omega_selected)
        plt.show()


if __name__ == "__main__":
    main()
