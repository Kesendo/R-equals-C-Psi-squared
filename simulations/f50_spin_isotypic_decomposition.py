#!/usr/bin/env python3
"""F50 spin-isotypic decomposition: closed-form attempt for K_N central-weight excess.

Computes the per-spin × per-weight pure-weight ker dimension for K_N, separating
the "single-block" contribution (operators supported entirely in one spin sector)
from "multi-block" contribution (operators that span multiple spin sectors but
remain pure-weight Pauli combinations).

Findings (2026-05-17):

  Single-block decomposition (per spin S, count of pure-weight-w Pauli
  combinations supported only in the S=eigenspace block):

    Max-spin (S = N/2, m=1, dim N+1): UNIVERSAL palindromic pattern
        (2, 4, 4, ..., 4, 2), sum = 4N for all N >= 3.

    Sub-max spin (S < N/2): concentrates pure-weight content at CENTRAL
    weights only. Width and parity of the central window depends on
    (m, d, N) in a non-trivial way.

  Empirical sub-max patterns:
    K_3 S=1/2  (m=2, d=2):           (0, 2, 2, 0)         w=1,2 (central window 2)
    K_4 S=1    (m=3, d=3):           (0, 0, 26, 0, 0)     w=2 only
    K_4 S=0    (m=2, d=1):           (0, 0, 1, 0, 0)      w=2 only
    K_5 S=3/2  (m=4, d=4):           (0, 0, 22, 22, 0, 0) w=2,3 (central window 2)
    K_5 S=1/2  (m=5, d=2):           (0, 0, 8, 8, 0, 0)   w=2,3 (central window 2)
    K_6 S=2    (m=5, d=5):           (0, 0, 38, 0, 38, 0, 0) w=2,4 (NOT w=3, parity!)
    K_6 S=1    (m=9, d=3):           (0, 0, 30, 124, 30, 0, 0) w=2,3,4
    K_6 S=0    (m=5, d=1):           (0, 0, 0, 0, 0, 0, 0)  ALL ZERO

  Multi-block contribution: weight-w operators that span multiple spin
  sectors block-diagonally. At K_3 N=3 w=1 multi-block = 2, identical to
  chain N=3 w=1 multi-block = 2. So the K_3 anomaly (+2 over chain) is
  ENTIRELY in single-block.

  Step toward closed-form:
    ker(K_N, w) = sum_S single_block(K_N, S, w) + multi_block(K_N, w)
    central_weight_excess(K_N) = (max-spin contributes uniformly across w)
                                + (sub-max contributes concentrated at central w)

  The max-spin palindromic 4N pattern is the universal "baseline" present
  at every K_N. The "anomaly" relative to chain comes from sub-max sectors
  concentrating their pure-weight mass at central weights.

  Open question remaining: closed-form for sub-max central pure-weight counts
  as f(m, d, N). The parity structure at K_6 (S=2 only even w, S=0 vanishes)
  suggests the spin-isotypic decomposition needs to be refined by an additional
  parity / symmetry label beyond just (m, d).
"""
from __future__ import annotations

import sys
import numpy as np
from itertools import product

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
    op = np.array([[1]], complex)
    for L in letters:
        op = np.kron(op, PAULI[L])
    return op


def heisenberg_KN(N, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for i in range(N):
        for j in range(i + 1, N):
            for L in "XYZ":
                ops = ["I"] * N
                ops[i] = L
                ops[j] = L
                H += (J / 4) * pauli_string_op(N, ops)
    return H


def weight(letters):
    return sum(1 for c in letters if c in "XY")


def spin_projectors(H, N, tol=1e-7):
    eigs, V = np.linalg.eigh(H)
    rounded = np.round(eigs, 7)
    unique = sorted(set(rounded))
    projs = {}
    for lam in unique:
        idx = [i for i, e in enumerate(rounded) if abs(e - lam) < tol]
        P = np.zeros((2 ** N, 2 ** N), complex)
        for i in idx:
            P += np.outer(V[:, i], V[:, i].conj())
        projs[lam] = (P, len(idx))
    return projs


def single_block_per_weight(N, target_w, P_S, tol=1e-7):
    d = 2 ** N
    all_letters = list(product("IXYZ", repeat=N))
    weight_ops = [pauli_string_op(N, L) for L in all_letters if weight(L) == target_w]
    M = np.column_stack([op.flatten("F") for op in weight_ops])
    n_w = M.shape[1]
    T = np.zeros((d * d, n_w), complex)
    for k in range(n_w):
        A = M[:, k].reshape(d, d, order="F")
        T[:, k] = (A - P_S @ A @ P_S).flatten("F")
    return n_w - np.linalg.matrix_rank(T, tol=tol)


def total_block_diagonal_per_weight(N, target_w, projs, tol=1e-7):
    """Dim of weight-w Pauli combinations A with P_lam A P_lam' = 0 for lam != lam'."""
    d = 2 ** N
    all_letters = list(product("IXYZ", repeat=N))
    weight_ops = [pauli_string_op(N, L) for L in all_letters if weight(L) == target_w]
    M = np.column_stack([op.flatten("F") for op in weight_ops])
    n_w = M.shape[1]
    constraints = []
    items = list(projs.items())
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            P_i, _ = items[i][1]
            P_j, _ = items[j][1]
            T1 = np.zeros((d * d, n_w), complex)
            T2 = np.zeros((d * d, n_w), complex)
            for k in range(n_w):
                A = M[:, k].reshape(d, d, order="F")
                T1[:, k] = (P_i @ A @ P_j).flatten("F")
                T2[:, k] = (P_j @ A @ P_i).flatten("F")
            constraints.append(T1)
            constraints.append(T2)
    if not constraints:
        return n_w
    big_T = np.vstack(constraints)
    return n_w - np.linalg.matrix_rank(big_T, tol=tol)


def main():
    print("=" * 78)
    print("F50 spin-isotypic decomposition for K_N (single + multi-block)")
    print("=" * 78)
    print()

    for N in [3, 4, 5, 6]:
        H = heisenberg_KN(N)
        projs = spin_projectors(H, N)
        spin_list = [(round(lam, 4), dim) for lam, (_, dim) in projs.items()]
        print(f"--- K_{N} (N={N}, eigenspaces {spin_list}) ---")
        header = f"  {'lam':>6} {'dim':>4} " + " ".join(f"{'w=' + str(w):>5}" for w in range(N + 1)) + f" {'sum':>6}"
        print(header)
        print("  " + "-" * (len(header) - 2))
        sb_per_w = {w: 0 for w in range(N + 1)}
        for lam, (P, dim) in projs.items():
            rows = [single_block_per_weight(N, w, P) for w in range(N + 1)]
            for w, r in enumerate(rows):
                sb_per_w[w] += r
            cols = " ".join(f"{r:>5}" for r in rows)
            print(f"  {lam:>+6.2f} {dim:>4} {cols} {sum(rows):>6}")
        # Total block-diagonal (single + multi)
        totals = [total_block_diagonal_per_weight(N, w, projs) for w in range(N + 1)]
        totals_cols = " ".join(f"{t:>5}" for t in totals)
        print(f"  {'total':>6} {'-':>4} {totals_cols} {sum(totals):>6}")
        mb = [totals[w] - sb_per_w[w] for w in range(N + 1)]
        mb_cols = " ".join(f"{m:>5}" for m in mb)
        print(f"  {'multi':>6} {'-':>4} {mb_cols} {sum(mb):>6}")
        print()


if __name__ == "__main__":
    main()
