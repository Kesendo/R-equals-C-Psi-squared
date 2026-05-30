#!/usr/bin/env python3
"""Is the windowed-hard break gamma-linear? (the clock's tick as the perturbation parameter)

At gamma=0 the Takt stops: L = -i[H,.] is purely imaginary, hence symmetric about 0 = -sigma,
so SOFT. The break is gamma-driven. If break(gamma)/gamma -> const as gamma->0 (and stays
constant), the break is first-order in gamma, and the windowed converse (non-bipartite => hard)
becomes a first-order perturbation statement: the O(gamma) asymmetry coefficient is nonzero iff
there is an odd cycle. Test break/gamma across gamma for a REAL and a FLUX hard pair, and a soft
pair (should stay ~0).
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain


def break_at(pair, gamma, N=4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    Sigma = N * gamma
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * Sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def main():
    cases = [
        ("REAL  XXZ+XZX", [('X', 'X', 'Z'), ('X', 'Z', 'X')]),
        ("FLUX  IXY+XIY", [('I', 'X', 'Y'), ('X', 'I', 'Y')]),
        ("SOFT  XXZ+ZXX", [('X', 'X', 'Z'), ('Z', 'X', 'X')]),
    ]
    gammas = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]
    for label, pair in cases:
        print(f"{label}:")
        prev = None
        for g in gammas:
            b = break_at(pair, g)
            ratio = b / g
            trend = "" if prev is None else f"  Δ(b/γ)={ratio - prev:+.4f}"
            print(f"    γ={g:.6f}  break={b:.6e}  break/γ={ratio:.5f}{trend}")
            prev = ratio
        print()


if __name__ == "__main__":
    main()
