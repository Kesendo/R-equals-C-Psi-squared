#!/usr/bin/env python3
"""F50 per-weight breakdown: where does the centralizer excess of high-symmetry
topologies actually live?

Companion to [`f50_topology_anomaly_sweep.py`](f50_topology_anomaly_sweep.py).
That sweep counts the weight-1 ker dim and finds K_3 N=3 as the unique
anomaly. This script asks: when the centralizer of H_G is larger than chain's
(true for many topologies, e.g. ring, star, K_N), where does the excess
live across weight sectors?

Findings (resolution of the open question from PROOF_WEIGHT1_DEGENERACY):

  Excess (vs chain) per weight at central topologies:

    K_3 N=3:  +2 at w=1, +2 at w=2  (palindromic pair around N/2 = 1.5)
    K_4 N=4:  +23 at w=2             (self-palindromic at N/2 = 2)
    K_5 N=5:  +40 at w=2, +40 at w=3 (palindromic pair around N/2 = 2.5)

  Other topologies at N=4 also show w=2 excess (chain baseline 13):
    star: +3, K_4-e: +9, ring: +10, K_4: +23.

The pattern: every graph with non-trivial Aut beyond chain has centralizer
excess at the central weights {floor(N/2), ceil(N/2)}, palindromic by F1.
K_N has the largest. K_3 N=3 only stands out for F50 (weight-1) because
N=3 puts central weight = 1.

`dim(ker[H, ·]|_w)` is computed via the rank of the constraint matrix
C|_{weight-w} where C[A] = [H, A] (matrix commutator). The pure-weight ker
is the # of pure-weight-w operators that commute with H.
"""
from __future__ import annotations

import sys
from itertools import product
from collections import Counter

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string_op(N, letters):
    op = np.array([[1]], dtype=complex)
    for L in letters:
        op = np.kron(op, PAULI[L])
    return op


def heisenberg(N, bonds, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (a, b) in bonds:
        for L in "XYZ":
            ops = ["I"] * N
            ops[a] = L
            ops[b] = L
            H += (J / 4) * pauli_string_op(N, ops)
    return H


def weight(letters):
    return sum(1 for c in letters if c in "XY")


def weight_w_ker_dim(H, N, target_w, tol=1e-8):
    """Number of pure-weight-w operators A with [H, A] = 0."""
    d = 2 ** N
    all_letters = list(product("IXYZ", repeat=N))
    weight_ops = [pauli_string_op(N, L) for L in all_letters if weight(L) == target_w]
    if not weight_ops:
        return 0
    M = np.column_stack([op.flatten("F") for op in weight_ops])
    I_d = np.eye(d, dtype=complex)
    C = np.kron(I_d, H) - np.kron(H.T, I_d)
    CM = C @ M
    return M.shape[1] - np.linalg.matrix_rank(CM, tol=tol)


def centralizer_dim_from_spectrum(H):
    eigs = np.linalg.eigvalsh(H)
    deg = Counter(np.round(eigs, 6))
    return sum(m ** 2 for m in deg.values())


def per_weight_breakdown(H, N):
    return [weight_w_ker_dim(H, N, w) for w in range(N + 1)]


def main():
    print("=" * 72)
    print("F50 per-weight centralizer breakdown")
    print("=" * 72)
    print()

    graphs = [
        ("N=3 chain", 3, [(0, 1), (1, 2)]),
        ("N=3 K_3", 3, [(0, 1), (1, 2), (2, 0)]),
        ("N=4 chain", 4, [(0, 1), (1, 2), (2, 3)]),
        ("N=4 ring", 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
        ("N=4 star", 4, [(0, 1), (0, 2), (0, 3)]),
        ("N=4 K_4-e", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]),
        ("N=4 K_4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
        ("N=5 chain", 5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
        ("N=5 ring", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
        ("N=5 K_5", 5, [(i, j) for i in range(5) for j in range(i + 1, 5)]),
    ]

    print(f"{'graph':>14} " + "".join(f"{w:>5}" for w in range(6))
          + f"{'sum_w':>7} {'Centr':>7}")
    print("-" * 75)
    results = {}
    for name, N, bonds in graphs:
        H = heisenberg(N, bonds)
        w_ker = per_weight_breakdown(H, N)
        centr = centralizer_dim_from_spectrum(H)
        results[name] = (N, w_ker, centr)
        cols = "".join(f"{w_ker[w]:>5}" if w < len(w_ker) else f"{'-':>5}"
                       for w in range(6))
        print(f"{name:>14} {cols} {sum(w_ker):>7} {centr:>7}")

    print()
    print("Excess (K_N - chain) per weight at each N:")
    for N_target in [3, 4, 5]:
        chain_name = f"N={N_target} chain"
        kn_name = f"N={N_target} K_{N_target}"
        if kn_name not in results:
            continue
        _, ch_w, ch_c = results[chain_name]
        _, kn_w, kn_c = results[kn_name]
        print(f"  N={N_target}: K_N vs chain centralizer = {kn_c} vs {ch_c} (Δ={kn_c - ch_c})")
        excess_str = " ".join(f"w={w}:Δ={kn_w[w]-ch_w[w]:+d}" for w in range(len(kn_w)))
        print(f"    weights: {excess_str}")

    print()
    print("Excess also for other high-symmetry topologies at N=4:")
    chain_n4 = results["N=4 chain"][1]
    for name in ["N=4 ring", "N=4 star", "N=4 K_4-e", "N=4 K_4"]:
        if name not in results:
            continue
        _, w_ker, _ = results[name]
        excess_at_w2 = w_ker[2] - chain_n4[2]
        print(f"  {name}: w=2 excess = {excess_at_w2:+d}  (chain N=4 w=2 = {chain_n4[2]})")

    print()
    print("Pattern: high-symmetry topologies have centralizer excess at CENTRAL weights")
    print("(w near N/2), palindromic by F1 Π-conjugation. K_N has the largest excess.")
    print("K_3 N=3 is special only because its central weight happens to be 1 = F50's tracked weight.")


if __name__ == "__main__":
    main()
