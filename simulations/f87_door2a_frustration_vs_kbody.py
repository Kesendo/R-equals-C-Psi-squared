#!/usr/bin/env python3
"""Door 2a thread: why is a frustrated bond hopping still soft? Is hardness a k>=3 phenomenon, not
geometric frustration? (scout, 2026-06-05).

The same GF(2) condition (an odd 𝔽₂-relation among flip-masks) gives HARD on the chain (windowed k>=3
terms, the F87 obstruction) but SOFT on the triangle (bond frustration). So the odd cycle is not the
whole story. Hypothesis: the F87 obstruction is intrinsically a k>=3 / windowed phenomenon; a k=2 bond
hopping is soft or truly regardless of graph topology, so geometric frustration never makes it hard.

This scout sweeps EVERY single bond letter-pair and the symmetric/antisymmetric combos on the frustrated
triangle (and 5-cycle), classifying each. If none is hard, frustration != hardness, and we have learned
that the F87 hard direction is a multi-body (k>=3) structure the chain's windowed placement supplied --
not a property of odd cycles per se.
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def spectral_verdict(N, bonds, bilinear, gamma=0.05, op_tol=1e-10, spec_tol=1e-6):
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [gamma] * N)
    sigma = N * gamma
    M = palindrome_residual(L, sigma, N)
    if np.linalg.norm(M) < op_tol:
        return 'truly'
    ev = np.linalg.eigvals(L)
    used = np.zeros(len(ev), dtype=bool); max_err = 0.0
    for i in range(len(ev)):
        if used[i]: continue
        d = np.abs(ev - (-ev[i] - 2 * sigma)); d[used] = np.inf
        j = int(np.argmin(d)); used[i] = True
        if j != i: used[j] = True
        max_err = max(max_err, float(d[j]))
    return 'soft' if max_err < spec_tol else 'hard'


def main():
    print("=" * 80)
    print("Door 2a thread: is any k=2 bond hopping HARD on a frustrated graph?")
    print("=" * 80)
    letters = "XYZ"
    triangle = (3, [(0, 1), (1, 2), (0, 2)])
    cycle5 = (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])

    for glabel, (N, bonds) in [("triangle K3", triangle), ("5-cycle", cycle5)]:
        print(f"\n  {glabel}:")
        verdicts = {}
        # all 9 single bond terms
        for a, b in product(letters, letters):
            v = spectral_verdict(N, bonds, [(a, b, 1.0)])
            verdicts.setdefault(v, []).append(a + b)
        # symmetric and antisymmetric combos of distinct pairs
        for a, b in [("X", "Y"), ("X", "Z"), ("Y", "Z")]:
            vs = spectral_verdict(N, bonds, [(a, b, 1.0), (b, a, 1.0)])
            va = spectral_verdict(N, bonds, [(a, b, 1.0), (b, a, -1.0)])
            verdicts.setdefault(vs, []).append(f"{a}{b}+{b}{a}")
            verdicts.setdefault(va, []).append(f"{a}{b}-{b}{a}")
        # the full Heisenberg + XY model
        verdicts.setdefault(spectral_verdict(N, bonds, [('X','X',1.),('Y','Y',1.),('Z','Z',1.)]), []).append("Heisenberg")
        for v in ("truly", "soft", "hard"):
            items = verdicts.get(v, [])
            print(f"    {v:<6}: {len(items):>2}  {', '.join(items)}")
        if "hard" not in verdicts:
            print(f"    => NO k=2 bond hopping is hard on {glabel}: frustration does not break the palindrome.")

    print("\n  reading: if no bond hopping is hard on the frustrated graphs, then the F87 obstruction is")
    print("  NOT geometric (bond-cycle) frustration. The hardness lives in the k>=3 windowed structure")
    print("  the chain supplied; the bond odd-cycle is a different cycle that the palindrome survives.")


if __name__ == "__main__":
    main()
