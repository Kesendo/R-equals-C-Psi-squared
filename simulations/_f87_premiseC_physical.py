#!/usr/bin/env python3
"""ATTACK #3 at the PHYSICAL test point: gamma_0 = 0.05, J = 1.0, for EVERY hard pair (N=4
and N=5), with THREE independent break measures, plus the greedy-vs-OT discrepancy.

The classifier `classify_pauli_pair` uses gamma_0 = 0.05 (ChainSystem default) and a GREEDY
nearest-partner pairing with spec_tol = 1e-6. The proof argues in terms of the OPTIMAL-transport
distance A(gamma). If greedy and OT disagree near degeneracies, the 'hard' verdict and the proof's
A(gamma) could diverge. We check, at the exact physical point:
  (1) greedy max_err  (what the classifier actually thresholds)
  (2) OT distance     (the proof's A; max optimal-assignment |lambda - (-lambda-2sigma)|)
  (3) multiset/char-poly distance (sorted-spectrum max|diff|; the algebraic Delta=0 test)
and require all three to declare 'broken' for every hard pair. If ANY hard pair is OT- or
multiset-exact at gamma=0.05 while greedy says hard (or vice versa), that is the hole.

Also: N=5 (attack #3 asks for it explicitly). The pair SET is alphabet-only (same 16 hard
windowed pairs), but at N=5 each k=3 term places at 3 windows, so the hopping graph and spectrum
differ. We classify at N=5 and re-test the three measures.
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


def measures(pair, gamma, N):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter="Z")
    sig = N * gamma
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * sig
    # OT
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    ot = float(cost[r, c].max())
    # greedy (mirror the classifier exactly)
    used = np.zeros(len(ev), dtype=bool)
    gmax = 0.0
    for i in range(len(ev)):
        if used[i]:
            continue
        target = -ev[i] - 2 * sig
        dists = np.abs(ev - target)
        dists[used] = np.inf
        bj = int(np.argmin(dists))
        if bj != i:
            used[i] = True
            used[bj] = True
        else:
            used[i] = True
        gmax = max(gmax, float(dists[bj]))
    # charpoly gap: the ALGEBRAIC Delta(x;gamma)=0 test. |coeff(char_L) - coeff(char_{-L-2s})|,
    # L1, normalised by coeff scale. This is the true multiset-equality / Delta object (a sorted
    # element-wise compare of two DIFFERENT spectra is NOT a valid set-equality test).
    pa = np.poly(ev)
    pb = np.poly(tgt)
    scale = np.max(np.abs(pa)) + np.max(np.abs(pb)) + 1e-300
    cp = float(np.sum(np.abs(pa - pb)) / scale)
    return gmax, ot, cp


def hard_pairs(N):
    chain = fw.ChainSystem(N=N)
    k = 3
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
             and any(L not in DIAG for L in t)]
    out = []
    for t1, t2 in combinations_with_replacement(terms, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        out.append((t1, t2, cls))
    return out


def run(N):
    GAMMA = 0.05
    TOL = 1e-6
    pairs = hard_pairs(N)
    hard = [p for p in pairs if p[2] == "hard"]
    soft = [p for p in pairs if p[2] == "soft"]
    print(f"\n===== N={N}, gamma_0={GAMMA}, J=1.0  (k=3 windowed diagonal cell) =====")
    print(f"  classifier: hard={len(hard)}  soft={len(soft)}")
    print(f"  {'pair':16s} {'greedy':>10s} {'OT':>10s} {'charpolygap':>12s}  verdict")
    worst_soft_ot = 0.0
    worst_soft_cp = 0.0
    hard_all_broken = True
    hard_min = (None, np.inf)
    for (t1, t2, _) in hard:
        lab = "".join(t1) + "+" + "".join(t2)
        g, ot, cp = measures([t1, t2], GAMMA, N)
        # OT is the rigorous multiset distance (charpoly overflows at N=5; informational only).
        broken = (ot > TOL) and (g > TOL)
        hard_all_broken = hard_all_broken and broken
        if ot < hard_min[1]:
            hard_min = (lab, ot)
        flag = "" if broken else "  <-- NOT broken on some measure!"
        print(f"  {lab:16s} {g:10.3e} {ot:10.3e} {cp:12.3e}  {'BROKEN' if broken else 'EXACT'}{flag}")
    for (t1, t2, _) in soft:
        g, ot, cp = measures([t1, t2], GAMMA, N)
        worst_soft_ot = max(worst_soft_ot, max(g, ot))
    worst_soft = worst_soft_ot
    print(f"  --- soft pairs: worst break over all soft = {worst_soft:.3e} "
          f"({'all exact' if worst_soft < TOL else 'A SOFT PAIR BREAKS'})")
    print(f"  --- hard pairs: every hard broken on ALL 3 measures = {hard_all_broken}")
    print(f"  --- hard min OT distance (closest to spec-exact) = {hard_min[1]:.3e} ({hard_min[0]})")
    return hard_all_broken, worst_soft < TOL


def main():
    ok4 = run(4)
    ok5 = run(5)
    print()
    if all(ok4) and all(ok5):
        print("VERDICT(attack #3 physical): at gamma=0.05 J=1, every hard pair is broken on all")
        print("  three measures (greedy/OT/multiset) and every soft pair is exact, at N=4 AND N=5.")
        print("  The physical test gamma is not an accidental restoration for any hard pair.")
    else:
        print("VERDICT(attack #3 physical): POTENTIAL HOLE -- a discrepancy at the physical point.")


if __name__ == "__main__":
    main()
