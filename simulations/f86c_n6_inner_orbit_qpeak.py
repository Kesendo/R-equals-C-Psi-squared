"""The c = 2, N = 6 per-F71-orbit Q_peak of F86c on the earlier 21-point t-grid reading (a fixed-t slice at the winning
grid point). PerF71OrbitKTable now reads the envelope max_t, and its second route is
`python simulations/f86_envelope_second_route.py 6 2 4.0`; this script reproduces the earlier reading.

K_b(Q, t) from the F89 eigendecomposition route of `f89_to_f86_kbond_via_eigendecomp.py` (its compute_kbond builds the
same operator as the C# table; it sets J = Q * gamma in its own coupling book, so the F86 Q enters as Q/2), maximised
over the t-grid linspace(0.6, 1.6) * t_peak (21 points, t_peak = 1/(4 gamma_0)), then over Q with a three-point parabolic
refinement at the grid maximum. Grids: the C# default [0.2, 4] x 153 (ResonanceScan.DefaultQGrid), and the coarse
[0.2, 6] x 60 and x 40 grids on which the frozen tables of PerF71OrbitObservation and the earlier F86NewIdeasProbe read.

Run: python simulations/f86c_n6_inner_orbit_qpeak.py   (about a minute)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f89_to_f86_kbond_via_eigendecomp as route  # noqa: E402

GAMMA0, N = 0.05, 6
T_PEAK = 1 / (4 * GAMMA0)
T_GRID = np.linspace(0.6 * T_PEAK, 1.6 * T_PEAK, 21)


def q_peak(grid, bond):
    k = np.array([route.compute_kbond(q / 2, GAMMA0, N, bond, T_GRID)[0] for q in grid])
    i = int(np.argmax(k))
    c = np.polyfit(grid[i - 1:i + 2], k[i - 1:i + 2], 2)
    return -c[1] / (2 * c[0])


for name, grid in (("[0.2, 4] x 153", np.linspace(0.2, 4.0, 153)),
                   ("[0.2, 6] x 60", np.linspace(0.2, 6.0, 60)),
                   ("[0.2, 6] x 40", np.linspace(0.2, 6.0, 40))):
    print(f"{name:16s} endpoint b=0 {float(q_peak(grid, 0))!r:20s} flanking b=1 {float(q_peak(grid, 1))!r:20s} central b=2 {float(q_peak(grid, 2))!r}")
