#!/usr/bin/env python3
"""EQ-031 step (a): test the absolute closed-form for ‖M(N)‖_F².

The ratio formulas in framework.py are:
  main class         ratio² = 4·k / (k − 1)
  single-body class  ratio² = 4·(2k − 1) / (2k − 3)

If these come from a clean per-bond × Liouvillian-extension structure,
the absolute formulas should be (telescoping the ratios from N=2):
  main class         ‖M(N)‖² = c_H · (N − 1) · 4^(N−2)
  single-body class  ‖M(N)‖² = c_H · (2N − 3) · 4^(N−2)

where c_H = ‖M(2)‖² is a per-Hamiltonian anchor at N=2 (one bond).

This script computes ‖M(N)‖² at N=2, 3, 4, 5 for representative Hamiltonians
across the trichotomy and tests both forms by extracting c_H from N=2 and
predicting N=3, 4, 5.

Prediction passes if measured / predicted = 1.000 to 4-decimal precision.
"""
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


GAMMA = 0.1
J = 1.0


def m_norm_sq(N, terms):
    bonds = [(i, i + 1) for i in range(N - 1)]
    bilinear = [(t[0], t[1], J) for t in terms]
    H = fw._build_bilinear(N, bonds, bilinear)
    L = fw.lindbladian_z_dephasing(H, [GAMMA] * N)
    M = fw.palindrome_residual(L, N * GAMMA, N)
    return float(np.linalg.norm(M)) ** 2


# Representative Hamiltonians, named (label, terms, expected_class)
# main class: any Pauli pair with at least one genuine 2-body term, OR
# single-body pair where both terms are on the same side of the bond
# (e.g., IX+IY: both letters are I·X and I·Y, same site=2)
# single-body class: pairs of form (Iσ, σI) for σ ∈ {X,Y,Z}, which give
# bond-overlap doubling on inner sites
TEST_CASES = [
    # truly Hamiltonians (M=0 expected, can't extract c)
    ('XX+YY (truly)',  [('X','X'), ('Y','Y')], 'truly'),
    ('IX+XI (truly)',  [('I','X'), ('X','I')], 'truly'),
    # main soft
    ('IX+IY (soft)',   [('I','X'), ('I','Y')], 'main'),
    ('XY+YZ (hard)',   [('X','Y'), ('Y','Z')], 'main'),
    ('XX+YZ (hard)',   [('X','X'), ('Y','Z')], 'main'),
    ('IX+IZ (hard)',   [('I','X'), ('I','Z')], 'main'),
    ('XY+ZX (hard)',   [('X','Y'), ('Z','X')], 'main'),
    # single-body class
    ('IY+YI (soft)',   [('I','Y'), ('Y','I')], 'single_body'),
    ('IZ+ZI (hard)',   [('I','Z'), ('Z','I')], 'single_body'),
]


def F_main(N):
    """Universal scaling factor: (N-1) * 4^(N-2)."""
    if N < 2:
        return float('nan')
    return (N - 1) * (4 ** (N - 2))


def F_single_body(N):
    """Single-body scaling factor: (2N-3) * 4^(N-2)."""
    if N < 2:
        return float('nan')
    return (2 * N - 3) * (4 ** (N - 2))


def main():
    print("EQ-031 step (a): absolute ‖M(N)‖² scaling test")
    print(f"  γ={GAMMA}, J={J}, N ∈ {{2, 3, 4, 5}}")
    print(f"  Predicted main:        ‖M(N)‖² = c_H · (N − 1) · 4^(N−2)")
    print(f"  Predicted single-body: ‖M(N)‖² = c_H · (2N − 3) · 4^(N−2)")
    print()

    Ns = [2, 3, 4, 5]
    print(f"  {'Hamiltonian':>16s}  {'class':>11s}  "
          + "  ".join(f"N={N:<2d}{'measured':>13s}" for N in Ns))
    for label, terms, cls in TEST_CASES:
        norms = [m_norm_sq(N, terms) for N in Ns]
        cells = "  ".join(f"      {n:>13.4e}" for n in norms)
        print(f"  {label:>16s}  {cls:>11s}  {cells}")
    print()

    print("Anchor c_H from N=2 (single bond), then predict N=3, 4, 5:")
    print()
    print(f"  {'Hamiltonian':>16s}  {'class':>11s}  "
          f"{'c_H (N=2)':>11s}  "
          + "  ".join(f"N={N:<2d}{'meas/pred':>11s}" for N in Ns[1:]))

    all_main_ratios = []
    all_sb_ratios = []
    for label, terms, cls in TEST_CASES:
        if cls == 'truly':
            continue  # ‖M‖²=0 by definition, can't anchor
        norms = [m_norm_sq(N, terms) for N in Ns]
        c_H = norms[0]
        F = F_main if cls == 'main' else F_single_body
        ratios = [norms[i] / (c_H * F(Ns[i])) for i in range(1, len(Ns))]
        if cls == 'main':
            all_main_ratios.extend(ratios)
        else:
            all_sb_ratios.extend(ratios)
        cells = "  ".join(f"      {r:>11.6f}" for r in ratios)
        print(f"  {label:>16s}  {cls:>11s}  {c_H:>11.4e}  {cells}")

    print()
    print("Summary (meas/pred = 1.000000 means absolute formula is exact):")
    if all_main_ratios:
        arr = np.array(all_main_ratios)
        print(f"  main class:        n={len(arr)}, "
              f"min={arr.min():.6f}, max={arr.max():.6f}, "
              f"mean={arr.mean():.6f}, std={arr.std():.2e}")
    if all_sb_ratios:
        arr = np.array(all_sb_ratios)
        print(f"  single-body class: n={len(arr)}, "
              f"min={arr.min():.6f}, max={arr.max():.6f}, "
              f"mean={arr.mean():.6f}, std={arr.std():.2e}")


if __name__ == "__main__":
    main()
