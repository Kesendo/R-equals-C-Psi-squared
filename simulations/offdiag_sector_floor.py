#!/usr/bin/env python3
"""Off-diagonal joint-popcount sectors: the 2*gamma*|dp| floor, and which sectors are FLAT.

The claim under test, in two halves.

FLOOR. Every basis coherence |a><b| in the joint-popcount sector (p, q) has
popcount(a) = p and popcount(b) = q, so its Hamming distance is at least |p - q|,
so <n_XY> >= |p - q|, so by the Absorption Theorem Re <= -2*gamma*|p - q|. The
floor is attained, because coherences at distance exactly |p - q| exist in every
such sector. In particular the ADJACENT sectors (k, k+-1) are floored at 2*gamma,
which is why they cannot host the chain dissipation gap: that gap lies below 2*gamma.

FLATNESS. A sector sits ENTIRELY at -2*gamma*|p - q| exactly when it touches the
vacuum or the fully excited state, min(p, q) = 0 or max(p, q) = N. There one index
is pinned to a single basis state, so every coherence is at distance exactly
|p - q|. Everywhere else the sector spreads above its floor. Among the adjacent
sectors this leaves only (0,1) and (N-1,N) flat, which is the statement
CHAIN_GAP_SECTOR_DIAGNOSTIC and SYMMETRY_CENSUS carry.

The gate is two-sided: it fires if a sector leaks, if its minimum |Re| is not the
predicted floor, or if its flatness disagrees with the touching rule in EITHER
direction. Swept over several (J, gamma) so no single lucky parameter passes it.

Docs: experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md (Finding 1),
      experiments/SYMMETRY_CENSUS.md (off-diagonal sector rates),
      docs/proofs/PROOF_ABSORPTION_THEOREM.md (Re = -2*gamma*<n_XY>).
Output: simulations/results/offdiag_sector_floor.txt
"""

import os
import sys

import numpy as np

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "offdiag_sector_floor.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)

TOL_LEAK = 1e-9
TOL_RATE = 1e-9


def op_at(P, site, n):
    """P acting on `site` of an n-qubit register, identity elsewhere."""
    M = np.array([[1]], dtype=complex)
    for i in range(n):
        M = np.kron(M, P if i == site else I2)
    return M


def liouvillian(n, J, gamma):
    """Heisenberg XXX chain, open boundaries, uniform Z-dephasing.

    Column-stacked convention: |a><b| has flat index a*d + b.
    """
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        for P in (X, Y, Z):
            H += J * op_at(P, b, n) @ op_at(P, b + 1, n)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(n):
        Zs = op_at(Z, s, n)
        L += gamma * (np.kron(Zs, Zs.T) - np.kron(Id, Id))
    return L


def popcount(x):
    return bin(x).count("1")


def sectors_of(n):
    """Flat indices grouped by (popcount(row), popcount(col))."""
    d = 2 ** n
    out = {}
    for r in range(d):
        for c in range(d):
            out.setdefault((popcount(r), popcount(c)), []).append(r * d + c)
    return out


def touches_boundary(p, q, n):
    """Does the sector touch the vacuum or the fully excited state?"""
    return min(p, q) == 0 or max(p, q) == n


def check(n, J, gamma):
    """Return (failures, rows) for one (n, J, gamma) point."""
    L = liouvillian(n, J, gamma)
    total = L.shape[0]
    failures, rows = [], []
    for (p, q), idx in sorted(sectors_of(n).items()):
        if p == q:
            continue
        sub = L[np.ix_(idx, idx)]
        outside = np.setdiff1d(np.arange(total), idx, assume_unique=False)
        leak = np.abs(L[np.ix_(outside, idx)]).max() if outside.size else 0.0
        re = np.linalg.eigvals(sub).real
        lo, hi = -re.max(), -re.min()          # smallest and largest |Re|
        floor = 2 * gamma * abs(p - q)
        flat = bool(np.allclose(re, re[0], atol=TOL_RATE))
        expect_flat = touches_boundary(p, q, n)

        if leak > TOL_LEAK:
            failures.append(f"N={n} J={J} g={gamma} ({p},{q}): leaks {leak:.2e}")
        if abs(lo - floor) > TOL_RATE:
            failures.append(f"N={n} J={J} g={gamma} ({p},{q}): min|Re|={lo:.9f} "
                            f"!= floor {floor:.9f}")
        if flat != expect_flat:
            failures.append(f"N={n} J={J} g={gamma} ({p},{q}): flat={flat} but "
                            f"touching rule predicts {expect_flat}")
        rows.append((p, q, len(idx), lo, hi, floor, flat, expect_flat))
    return failures, rows


def main():
    log("=" * 74)
    log("OFF-DIAGONAL JOINT-POPCOUNT SECTORS: FLOOR AND FLATNESS")
    log("=" * 74)
    log()
    log("  floor(p,q)  = 2*gamma*|p-q|, attained (Absorption Theorem + Hamming)")
    log("  flat(p,q)  <=> min(p,q) = 0 or max(p,q) = N (touches vacuum / full)")
    log()

    all_failures = []
    points = [(1.0, 0.5), (1.0, 0.05), (0.3, 0.7), (2.0, 0.11)]
    for n in (3, 4, 5):
        for J, gamma in points:
            failures, rows = check(n, J, gamma)
            all_failures += failures
            log(f"--- N={n}, J={J}, gamma={gamma} "
                f"({len(rows)} off-diagonal sectors) ---")
            for (p, q, dim, lo, hi, floor, flat, exp_flat) in rows:
                if q < p:
                    continue        # print one of each conjugate pair
                mark = "FLAT" if flat else "spread"
                log(f"    ({p},{q}) dim={dim:4d}  |Re| in [{lo:.6f}, {hi:.6f}]  "
                    f"floor={floor:.6f}  {mark}")
            log()

    log("=" * 74)
    if all_failures:
        log(f"FAILURES: {len(all_failures)}")
        for f in all_failures:
            log("  " + f)
    else:
        log("ALL CHECKS PASSED: every off-diagonal sector is exactly closed, every")
        log("one attains 2*gamma*|dp| as its smallest |Re|, and a sector is flat at")
        log("that value exactly when it touches the vacuum or the fully excited")
        log("state. Among the ADJACENT sectors (|dp| = 1) that leaves only (0,1)")
        log("and (N-1,N) flat, at every (N, J, gamma) swept.")
    log("=" * 74)
    return 1 if all_failures else 0


if __name__ == "__main__":
    sys.exit(main())
