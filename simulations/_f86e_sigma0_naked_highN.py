"""F86e: σ_0(N) − 2 (the naked change) over a long N-range, read as an angle.

Tom 2026-05-20: the naked change σ_0(N) − 2 "wanders like an angle";
needs more N to see it.

V_inter is built DIRECTLY from the popcount-block Hamiltonians H_p1 (N×N)
and H_p2 (C(N,2)²); the full block_dim² M_H is never materialised, so N
goes far past 20.

  V_inter[(i',j'),(i,j)] = δ_{j',j}·H_p1[i',i] − δ_{i',i}·H_p2[j',j]
  rows  = HD-1 pairs (i' ⊆ j',  popcount-1 / popcount-2)
  cols  = HD-3 pairs (i ∩ j = ∅)
σ_0 = top singular value of V_inter (commutator C=[H,·] is real; the −i of
M_H is a global phase, does not change singular values).
"""
from __future__ import annotations

import math
import sys
import time
from itertools import combinations

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def popcount(x: int) -> int:
    return bin(x).count("1")


def build_h_block(N: int, p: int):
    """XY-chain hop Hamiltonian restricted to the popcount-p subspace.
    H = (J/2)Σ(XX+YY), J=1 → hop amplitude 1 between bits-differ neighbours.

    States are generated directly from bit-position combinations, never
    iterate range(2**N), which is catastrophic past N≈25."""
    states = sorted(sum(1 << b for b in combo)
                     for combo in combinations(range(N), p))
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        for b in range(N - 1):
            if ((s >> b) & 1) != ((s >> (b + 1)) & 1):
                s2 = s ^ ((1 << b) | (1 << (b + 1)))
                H[idx[s2], idx[s]] += 1.0
    return H, states, idx


def sigma_0_direct(N: int) -> float:
    H_p1, states_p1, idx_p1 = build_h_block(N, 1)
    H_p2, states_p2, idx_p2 = build_h_block(N, 2)

    hd1 = [(ip, jp) for ip in states_p1 for jp in states_p2
           if (ip & jp) == ip]                      # i' ⊆ j'  → HD=1
    hd3 = [(ii, jj) for ii in states_p1 for jj in states_p2
           if (ii & jj) == 0]                       # disjoint → HD=3

    i1 = np.array([p[0] for p in hd1])
    j1 = np.array([p[1] for p in hd1])
    i3 = np.array([p[0] for p in hd3])
    j3 = np.array([p[1] for p in hd3])

    i1x = np.array([idx_p1[x] for x in i1])
    i3x = np.array([idx_p1[x] for x in i3])
    j1x = np.array([idx_p2[x] for x in j1])
    j3x = np.array([idx_p2[x] for x in j3])

    term1 = H_p1[i1x[:, None], i3x[None, :]] * (j1[:, None] == j3[None, :])
    term2 = -H_p2[j1x[:, None], j3x[None, :]] * (i1[:, None] == i3[None, :])
    V = term1 + term2
    return float(np.linalg.svd(V, compute_uv=False)[0])


def main() -> None:
    # Verify against the fw-based values
    print("Verification (direct builder vs known fw values):")
    for N, known in [(5, 2.7650951722), (7, 2.8284271247)]:
        s = sigma_0_direct(N)
        print(f"  N={N}: direct={s:.10f}  known={known:.10f}  "
              f"diff={s-known:+.2e}")
    print()

    print("=" * 84)
    print("σ_0(N) − 2 (the naked change) and read as angles")
    print("=" * 84)
    print(f"{'N':>3}  {'σ_0(N)':>15}  {'σ_0−2':>13}  "
          f"{'arccos(σ_0−2)°':>15}  {'arcsin°':>10}  {'Δθ_acos°':>10}")
    print("-" * 84)

    prev_acos = None
    rows = []
    for N in range(4, 45):
        t0 = time.time()
        s = sigma_0_direct(N)
        naked = s - 2.0
        acos_deg = math.degrees(math.acos(naked)) if -1 <= naked <= 1 else float("nan")
        asin_deg = math.degrees(math.asin(naked)) if -1 <= naked <= 1 else float("nan")
        dtheta = (acos_deg - prev_acos) if prev_acos is not None else float("nan")
        prev_acos = acos_deg
        elapsed = time.time() - t0
        flag = "" if elapsed < 3 else f"  ({elapsed:.1f}s)"
        print(f"{N:>3}  {s:>15.10f}  {naked:>13.9f}  "
              f"{acos_deg:>15.6f}  {asin_deg:>10.5f}  {dtheta:>10.5f}{flag}")
        rows.append((N, s, naked, acos_deg, asin_deg))

    print("-" * 84)
    print()
    # Reference angles
    print("Reference: σ_0(∞) − 2 ≈ 0.86287  →  arccos ≈ 30.39°,  arcsin ≈ 59.61°")
    print(f"  cos(30°) = √3/2 = {math.sqrt(3)/2:.6f}  (naked-limit 0.86287 is NEAR but not √3/2)")
    print(f"  N=7 naked = 2√2−2 = 2(√2−1) = {2*(math.sqrt(2)-1):.9f}  "
          f"→ arccos = {math.degrees(math.acos(2*(math.sqrt(2)-1))):.5f}°")


if __name__ == "__main__":
    main()
