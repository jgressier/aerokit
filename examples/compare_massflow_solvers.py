"""Compare the accuracy and CPU time of MassFlow's inversion methods.

Run from the repository root, for example:
    python scripts/compare_massflow_solvers.py --size 10000 --repeats 20
"""

import argparse
from pathlib import Path
import sys
import time

import numpy as np

# When invoked by path, Python otherwise searches ``scripts/`` before this
# checkout and may import an unrelated installed version of aerokit.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aerokit.aero import MassFlow


def benchmark(sigma, subsonic, method, gamma, repeats):
    """Return CPU seconds per call and the largest area-ratio residual."""
    solver = MassFlow.MachSub_Sigma if subsonic else MassFlow.MachSup_Sigma
    cpu_start = time.process_time()
    for _ in range(repeats):
        mach = solver(sigma, gamma=gamma, method=method)
    cpu_seconds = (time.process_time() - cpu_start) / repeats
    residual = np.max(np.abs(MassFlow.Sigma_Mach(mach, gamma=gamma) - sigma))
    return cpu_seconds, residual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=10000, help="number of area ratios per branch")
    parser.add_argument("--repeats", type=int, default=10, help="timed calls per method")
    parser.add_argument("--gamma", type=float, default=1.4, help="specific-heat ratio")
    args = parser.parse_args()

    # Include the difficult near-throat region and progressively larger areas.
    sigma = np.geomspace(1.0 + 1.0e-12, 100.0, args.size)
    print("method   branch       CPU / call (ms)   max |Sigma(M)-sigma|")
    for method in ("brentq", "newton"):
        for branch, subsonic in (("subsonic", True), ("supersonic", False)):
            cpu_seconds, residual = benchmark(sigma, subsonic, method, args.gamma, args.repeats)
            print("{:<8} {:<11} {:>14.3f}   {:.3e}".format(method, branch, cpu_seconds * 1.0e3, residual))


if __name__ == "__main__":
    main()
