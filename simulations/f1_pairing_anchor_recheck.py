"""Re-pair the three Python F1 anchors whose max pairing distance stands out
(ring_N6, ring_N4, star_N6; simulations/results/f1_n8_n9_metrics/*_python.json).

The anchor scripts (f1_topology_heisenberg_small_n_anchor.py, importing the pairing of
f1_chain_heisenberg_small_n_anchor.py) pair by a SORT on the key real*1e6 + imag. This
rebuilds each spectrum in the same convention (spin J/4, J = 1, gamma = 0.5) and pairs
it three ways: that sort, a min-sum assignment (scipy linear_sum_assignment on
|lambda_i - mirror_j|), and the per-eigenvalue nearest mirror (a lower bound for the
largest distance of any matching). Where the sort exceeds the min-sum assignment by
orders of magnitude, the recorded number is the sort's misordering, not the spectrum.

Where the sort and the min-sum assignment agree, the second block reads the eigenvalue cluster the largest nearest-mirror
distance sits in: its centre c, the nullities of L - c and (L - c)^2 by SVD (a Jordan-2 block
adds one to the second and not the first; the cluster size equal to the second nullity
leaves no block larger than 2), and the eigensolver's scale at a Jordan-2 eigenvalue,
sqrt(eps) * ||L||_2. The nullity read is valid only at the exact eigenvalue: every c here
is an integer, recovered by rounding to 1e-6. The K4 factor of the disconnected N = 8 anchor (K4 plus a
disjoint 4-chain, a C# block-spectrum run) is read at its centre -sigma, where its defect is.

Run: python simulations/f1_pairing_anchor_recheck.py   (a few minutes; the N = 6 SVDs dominate)
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


def jordan_reading(name, n, bonds, at_centre=False):
    sigma = n * GAMMA
    L = lindbladian_z_dephasing(heisenberg_graph_h(n, bonds, J=1.0), [GAMMA] * n)
    ev = np.linalg.eigvals(L)
    nearest = np.abs(ev[:, None] - (-2.0 * sigma - ev)[None, :]).min(axis=1)
    worst = ev[int(np.argmax(nearest))]
    c = complex(-sigma) if at_centre else complex(round(worst.real, 6), round(worst.imag, 6))
    a = L - c * np.eye(L.shape[0])
    norm = np.linalg.norm(L, 2)
    tol = 1e-9 * norm
    null1 = int(np.sum(np.linalg.svd(a, compute_uv=False) < tol))
    null2 = int(np.sum(np.linalg.svd(a @ a, compute_uv=False) < tol * norm))
    alg = int(np.sum(np.abs(ev - c) < 1e-5))
    scale = np.sqrt(np.finfo(float).eps) * norm
    eps_l = np.finfo(float).eps * norm
    print(f"{name}: c = {c} (centre -sigma = {-sigma}), cluster {alg}, "
          f"nullity(L - c) {null1}, nullity((L - c)^2) {null2}, Jordan-2 blocks {null2 - null1}; "
          f"nearest max {nearest.max():.3e} = {nearest.max() / scale:.2f} sqrt(eps) ||L|| "
          f"= {nearest.max() / eps_l:.0f} eps ||L||")


def complete_bonds(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


for name, n, bonds in CASES[1:] + [("K4 (factor of the N = 8 anchor)", 4, complete_bonds(4))]:
    jordan_reading(name, n, bonds, at_centre=name.startswith("K4"))
