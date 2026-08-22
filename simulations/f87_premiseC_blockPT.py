#!/usr/bin/env python3
"""DEFINITIVE attack #1: is "O(gamma) shift = M_omega block eigenvalues" correct?

The earlier slope check used nearest-neighbour matching to the gamma=0 spectrum, which is
ILL-POSED because the gamma=0 spectrum -i*omega is MASSIVELY degenerate (many coherences
|E_a><E_b| share one omega). That is precisely the degenerate case. The correct test groups by
omega and compares SETS within each omega-group:

  For small gamma, eigenvalues of L cluster near each -i*omega. Within the cluster with
  Im(lambda) ~ -omega, the quantities (lambda - (-i*omega))/gamma should converge, as gamma->0,
  to the eigenvalues of M_omega (the block projection of D_hat). This is standard degenerate
  first-order PT for a diagonalizable L0 with a (possibly non-normal) perturbation, PROVIDED the
  block M_omega is diagonalizable (verified: cond < 100).

This probe, per omega-group, per pair:
  - count match: does the omega-cluster size at small gamma equal dim(M_omega)? (degeneracy
    preserved => first-order PT applies cleanly)
  - SET match: sort {(lambda+i*omega)/gamma : Im(lambda)~-omega} vs sort(eig(M_omega)); report
    max|diff|. If ~0 for every omega, the central PT step is bit-exact and attack #1 fails.
  - extrapolate to gamma->0 with two gammas (Richardson) to kill O(gamma) curvature.

Run on SOFT + both hard pairs.
"""
from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain, site_op


def Mblocks(H, N=4):
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, "Z") @ V for l in range(N)]
    groups = defaultdict(list)
    for a in range(d):
        for b in range(d):
            groups[round(E[a] - E[b], 6)].append((a, b))
    out = {}
    for omega, modes in groups.items():
        n = len(modes)
        M = np.zeros((n, n), dtype=complex)
        for i, (a, b) in enumerate(modes):
            for j, (ap, bp) in enumerate(modes):
                val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
                if (a, b) == (ap, bp):
                    val -= N
                M[i, j] = val
        out[omega] = M
    return out, E


def cluster_shifts(H, omegas, N=4, g=1e-5):
    """Eigenvalues of L(g), grouped by nearest omega; return per-omega the set
    (lambda + i*omega)/g for the cluster of the right size."""
    L = lindbladian_pauli_dephasing(H, [g] * N, dephase_letter="Z")
    w = np.linalg.eigvals(L)
    # assign each eigenvalue to nearest omega by imaginary part
    om_arr = np.array(sorted(omegas))
    out = defaultdict(list)
    for lam in w:
        j = int(np.argmin(np.abs(lam.imag - (-om_arr))))  # Im(lam) ~ -omega
        om = om_arr[j]
        out[round(om, 6)].append((lam + 1j * om) / g)
    return out


def main():
    N = 4
    cases = [
        ("SOFT XXZ+ZXX", [("X", "X", "Z"), ("Z", "X", "X")]),
        ("FLUX IXY+XIY", [("I", "X", "Y"), ("X", "I", "Y")]),
        ("REAL XXZ+XZX", [("X", "X", "Z"), ("X", "Z", "X")]),
    ]
    for label, pair in cases:
        H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
        blocks, E = Mblocks(H, N)
        omegas = list(blocks.keys())
        # Richardson: two gammas
        g1, g2 = 1e-5, 5e-6
        sh1 = cluster_shifts(H, omegas, N, g1)
        sh2 = cluster_shifts(H, omegas, N, g2)
        worst_set = 0.0
        worst_om = None
        size_mismatch = 0
        n_blocks = 0
        for om, M in blocks.items():
            n_blocks += 1
            bk = round(om, 6)
            s1 = np.array(sorted(np.real(sh1.get(bk, [])))) if bk in sh1 else np.array([])
            s2 = np.array(sorted(np.real(sh2.get(bk, [])))) if bk in sh2 else np.array([])
            dimM = M.shape[0]
            if len(s1) != dimM or len(s2) != dimM:
                size_mismatch += 1
                continue
            # Richardson extrapolate Re shift to g->0: shift(g) = c + O(g); use 2 pts
            # linear in g: c = (g1*s2 - g2*s1)/(g1-g2) per sorted index (assumes same ordering)
            c = (g1 * s2 - g2 * s1) / (g1 - g2)
            bvals = np.sort(np.linalg.eigvals(M).real)
            d = float(np.max(np.abs(np.sort(c) - bvals)))
            if d > worst_set:
                worst_set = d
                worst_om = om
        print(f"\n{label}")
        print(f"  blocks: {n_blocks}   omega-cluster size mismatches: {size_mismatch}")
        print(f"  worst per-omega SET diff (Richardson g->0 shift vs eig(M_omega)) = {worst_set:.3e}  at omega={worst_om}")
        if size_mismatch == 0 and worst_set < 1e-3:
            print("  => degenerate first-order PT VALID: shift-set == M_omega spectrum, every omega.")
        else:
            print("  => MISMATCH: investigate (PT step may not be the M_omega eigenvalues).")


if __name__ == "__main__":
    main()
