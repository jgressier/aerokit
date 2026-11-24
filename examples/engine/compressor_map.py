"""Draw a preliminary compressor map with an ideal :class:`BladeRow`.

The map is parameterized by blade speed ``U`` and the imposed back-pressure
ratio ``Pt_in / p_out``.  For each pair, the axial velocity is adjusted until
the ideal blade-row outlet static pressure is ``p_out``.  The right panel is
the same operating points in the conventional loading/flow-coefficient plane:
``psi = Delta ht / U**2`` and ``phi = Vx / U``.

This is a single-radius, lossless streamtube example.  It is not a calibrated
compressor characteristic and does not represent stall, choke-line losses, or
radial stacking.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# Let the documented command work directly from a source checkout too; an
# installed aerokit package needs no adjustment.
if __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aerokit.aero.blade import BladeRow
from aerokit.aero.model2D import State2d


def operating_point(U, pt_in_over_p_out, *, p_out=1.e5, rtt_in=287.05 * 288.15,
                    beta_out=-50.0, gamma=1.4):
    """Solve one ideal operating point at imposed ``Pt_in / p_out``.

    ``beta_out`` is a relative outlet angle.  With a positive ``U`` and a
    modestly negative angle, the row produces positive exit absolute swirl and
    positive Euler work.
    """
    row = BladeRow(Uy=U)
    pt_in = pt_in_over_p_out * p_out
    max_vx = np.sqrt(2.0 * gamma * rtt_in / (gamma - 1.0)) * (1.0 - 1.e-10)

    def pressure_residual(vx):
        inlet = State2d.from_pt_rtt_uv(pt_in, rtt_in, vx, gamma=gamma)
        try:
            outlet = row.downstream(inlet, beta_out, conservation="massflow")
        except ValueError:
            return np.nan
        return outlet.p - p_out

    # The root is bracketed by scanning, avoiding assumptions about the
    # monotonicity of every velocity triangle.
    velocities = np.linspace(1.e-6, max_vx, 500)
    residuals = np.array([pressure_residual(vx) for vx in velocities])
    valid = np.isfinite(residuals[:-1]) & np.isfinite(residuals[1:])
    changes = np.flatnonzero(valid & (residuals[:-1] * residuals[1:] <= 0.0))
    if not len(changes):
        return None
    i = changes[0]
    vx = brentq(pressure_residual, velocities[i], velocities[i + 1])
    inlet = State2d.from_pt_rtt_uv(pt_in, rtt_in, vx, gamma=gamma)
    outlet = row.downstream(inlet, beta_out, conservation="massflow")
    dHt = row.specific_work(inlet, outlet)
    return {"U": U, "Pt_in/p_out": pt_in_over_p_out, "Vx": vx,
            "dHt": dHt, "phi": vx / U, "psi": dHt / U**2,
            "Pt_out/Pt_in": outlet.Ptot() / inlet.Ptot()}


def make_map(speeds, pressure_ratios, **kwargs):
    """Return one speed line of operating points for each value of ``U``."""
    return [[point for ratio in pressure_ratios
             if (point := operating_point(U, ratio, **kwargs)) is not None]
            for U in speeds]


def plot_map(curves):
    """Plot the physical and non-dimensional views of the same map."""
    fig, (ax_map, ax_coeff) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(curves)))
    for curve, color in zip(curves, colors):
        if not curve:
            continue
        vx = np.array([point["Vx"] for point in curve])
        pressure_ratio = np.array([point["Pt_in/p_out"] for point in curve])
        psi = np.array([point["psi"] for point in curve])
        phi = np.array([point["phi"] for point in curve])
        label = r"$U=%.0f$ m/s" % curve[0]["U"]
        ax_map.plot(vx, pressure_ratio, "o-", color=color, label=label)
        ax_coeff.plot(phi, psi, "o-", color=color, label=label)

    ax_map.set(xlabel=r"axial velocity $V_x$ (m/s)", ylabel=r"$P_{t,in}/p_{out}$",
               title="Ideal compressor map")
    ax_coeff.set(xlabel=r"flow coefficient $\phi=V_x/U$",
                 ylabel=r"loading coefficient $\psi=\Delta h_t/U^2$",
                 title=r"Same points in $\psi$--$\phi$ coordinates")
    for axis in (ax_map, ax_coeff):
        axis.grid(True)
        axis.legend(fontsize="small")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", help="write the figure to this path instead of displaying it")
    args = parser.parse_args()
    speeds = np.linspace(100.0, 200.0, 11)
    pressure_ratios = (1.02, 1.05, 1.10, 1.15, 1.20, 1.3, 1.4, 1.5)
    figure = plot_map(make_map(speeds, pressure_ratios))
    if args.save:
        figure.savefig(args.save, dpi=180)
    else:
        plt.show()


if __name__ == "__main__":
    main()
