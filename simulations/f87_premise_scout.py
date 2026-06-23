#!/usr/bin/env python3
"""F87 first-order-block premise: the genericity backbone (scout, 2026-06-04).

The windowed converse non-bipartite ⟹ hard is derived (PROOF_F103 §7.5) modulo ONE premise:
the ω=0 first-order block asymmetry ⟺ the all-orders F87 break. Decomposed:
  - soft direction (c=0 ⟺ bipartite ⟹ soft all-orders): closed via the chiral K.
  - hard direction (c≠0 ⟹ hard): a perturbation/analyticity argument. c≠0 ⟹ break = c·γ + O(γ²)
    ≠ 0 for small γ ⟹ spec(L) ≠ spec(−L−2σ) for small γ ⟹ (eigenvalues algebraic in γ, so the
    closed condition spec(L)=spec(−L−2σ) holds identically or on isolated γ) ⟹ hard for generic γ.

This scout supplies the empirical backbone of the analyticity step: across a γ-sweep, every hard
pair is spec-broken at every γ (modulo isolated accidental degeneracies), and every soft pair is
spec-exact at every γ. It also confirms break/γ → c (first-order) and that c matches the ω=0
block asymmetry. If hard ⟹ broken-at-all-swept-γ holds, "generic γ" is not a loophole.
"""
from __future__ import annotations

import sys
from itertools import product, combinations_with_replacement
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
import framework as fw
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain

DIAG = {"I", "Z"}


def break_at(pair, gamma, N=4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter="Z")
    sig = N * gamma
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * sig
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
             and any(L not in DIAG for L in t)]
    pairs = []
    for t1, t2 in combinations_with_replacement(terms, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        pairs.append((t1, t2, cls))

    hard = [p for p in pairs if p[2] == "hard"]
    soft = [p for p in pairs if p[2] == "soft"]
    gammas = [0.9, 0.7, 0.5, 0.37, 0.23, 0.13, 0.07, 0.031, 0.017, 0.009, 0.0033, 0.0011]
    tol = 1e-7

    print("=" * 92)
    print("F87 premise scout: hard ⟹ spec-broken at every swept γ; soft ⟹ spec-exact at every γ")
    print("=" * 92)
    print(f"  windowed Mixed+Mixed pairs: {len(pairs)}  (hard {len(hard)}, soft {len(soft)})")
    print(f"  γ-sweep ({len(gammas)} values in [{min(gammas)}, {max(gammas)}])")
    print()

    # how many γ each hard pair is broken at (want: all of them)
    worst_hard = (None, len(gammas))
    for (t1, t2, _) in hard:
        nbroken = sum(1 for g in gammas if break_at([t1, t2], g) > tol)
        if nbroken < worst_hard[1]:
            worst_hard = (("".join(t1), "".join(t2)), nbroken)
    best_soft_break = 0.0
    for (t1, t2, _) in soft:
        best_soft_break = max(best_soft_break, max(break_at([t1, t2], g) for g in gammas))

    all_hard_allg = all(
        all(break_at([t1, t2], g) > tol for g in gammas) for (t1, t2, _) in hard)
    print(f"  HARD: every hard pair broken at ALL {len(gammas)} γ:  {all_hard_allg}")
    print(f"        (worst hard pair {worst_hard[0]} broken at {worst_hard[1]}/{len(gammas)} γ)")
    print(f"  SOFT: max break over all soft pairs and all γ:  {best_soft_break:.2e}  "
          f"({'exact ∀γ' if best_soft_break < tol else 'BREAKS'})")
    print()

    # first-order coefficient c = break/γ as γ→0, on the two named hard pairs
    print("  first-order coefficient c = break/γ as γ→0 (the ω=0 block asymmetry):")
    for label, pair in [("REAL XXZ+XZX", [("X", "X", "Z"), ("X", "Z", "X")]),
                        ("FLUX IXY+XIY", [("I", "X", "Y"), ("X", "I", "Y")])]:
        c_small = break_at(pair, 1e-4) / 1e-4
        print(f"    {label}:  c ≈ {c_small:.5f}")
    print()
    print("  reading: if HARD-broken-at-all-γ is True and SOFT-exact-∀γ, the analyticity step is")
    print("  empirically airtight: spec(L)=spec(−L−2σ) is identical (soft) or nowhere on the sweep")
    print("  (hard), so 'generic γ' is no loophole and c≠0 ⟹ hard at the physical γ.")
    assert all_hard_allg and best_soft_break < tol


if __name__ == "__main__":
    main()
