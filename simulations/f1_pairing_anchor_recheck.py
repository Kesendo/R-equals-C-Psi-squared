"""Re-pair the three Python F1 anchors whose max pairing distance stands out
(ring_N6, ring_N4, star_N6; simulations/results/f1_n8_n9_metrics/*_python.json).

The anchor scripts (f1_topology_heisenberg_small_n_anchor.py, importing the pairing of
f1_chain_heisenberg_small_n_anchor.py) pair by a SORT on the key real*1e6 + imag. This
rebuilds each spectrum in the same convention (spin J/4, J = 1, gamma = 0.5) and pairs
it three ways: that sort, a min-sum assignment (scipy linear_sum_assignment on
|lambda_i - mirror_j|), and the per-eigenvalue nearest mirror (a lower bound for the
largest distance of any matching). Where the sort exceeds the min-sum assignment by
orders of magnitude, the recorded number is the sort's misordering, not the spectrum.

Run: python simulations/f1_pairing_anchor_recheck.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent))
from framework.lindblad import lindbladian_z_dephasing  # noqa: E402
from f1_chain_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    compute_palindromic_pairing_distances,
    heisenberg_graph_h,
)

GAMMA = 0.5


def star_bonds(n):
    return [(0, k) for k in range(1, n)]


def ring_bonds(n):
    return chain_bonds(n) + [(n - 1, 0)]


CASES = [("ring_N6", 6, ring_bonds(6)), ("ring_N4", 4, ring_bonds(4)), ("star_N6", 6, star_bonds(6))]

for name, n, bonds in CASES:
    sigma = n * GAMMA
    h = heisenberg_graph_h(n, bonds, J=1.0)
    ev = np.linalg.eigvals(lindbladian_z_dephasing(h, [GAMMA] * n))
    mirror = -2.0 * sigma - ev
    sort_max = compute_palindromic_pairing_distances(ev, sigma)["MaxPairingDistance"]
    cost = np.abs(ev[:, None] - mirror[None, :])
    r, c = linear_sum_assignment(cost)
    opt_max = cost[r, c].max()
    nearest_max = cost.min(axis=1).max()
    print(f"{name}: sort {sort_max:.6e}   min-sum {opt_max:.3e}   nearest {nearest_max:.3e}")
