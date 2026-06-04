#!/usr/bin/env python3
"""ATTACK #3 (genericity) + #4 (higher-order healing): could a HARD pair be
spec-EXACT at some isolated gamma>0?

The proof's analyticity step: Delta(x;gamma) = char_L(x;gamma) - char_{-L-2sigma}(x;gamma)
has coefficients polynomial in gamma (L is affine in gamma), so Delta != 0 (identically)
=> Delta vanishes on only finitely many gamma. The first-order c != 0 shows the spectra
split at small gamma; the claim is then "hard for all but finitely many gamma, and the
physical/test gamma is generic (not an accidental root)".

The hole to realise: a hard pair (c != 0) that is nonetheless spec-EXACT at some isolated
gamma_0 > 0 (an accidental coincidence of the two spectra AS SETS, i.e. a genuine root of
Delta). If such gamma_0 coincides with the physical test gamma for any hard pair, the
converse is false there. Even if not, an accidental root anywhere is a real phenomenon the
proof must (and does) acknowledge; we want to know if it EVER lands on a hard pair, and how
close it can come to a "nice" gamma.

Two independent break measures, to avoid OT-assignment artifacts:
  (A) OT distance A(gamma) = max optimal-assignment |lambda - (-lambda-2sigma)| (the doc's measure)
  (B) char-poly difference: the actual algebraic object. We compute the *sorted-spectrum*
      L2 distance after optimal matching AND the symmetric-function (char-poly coeff) gap.
      The cleanest set-equality test: sort both multisets in the complex plane (by (re,im))
      and take max|diff|. Spectra equal as multisets <=> this is ~0. This is what Delta=0 means.

Sweep: all 16 hard windowed pairs (N=4, Z), fine deterministic grid + 200 random gammas in
(0,2], + specifically the physical test gammas {1.0, 0.5, 0.1, 0.05}. Flag ANY hard pair whose
set-distance dips below 1e-6 at any gamma>1e-3.
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


def L_of(pair, gamma, N=4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    return lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter="Z"), N * gamma


def ot_break(ev, tgt):
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def set_distance(ev, tgt):
    """Multiset distance done RIGHT: the optimal-assignment max. A naive element-wise compare
    of independently-sorted ev and tgt is INVALID (they are two different spectra; dense regions
    misalign). The optimal-assignment distance IS the multiset (Hausdorff/transport) distance and
    is ~0 iff the multisets coincide. We reuse the OT max here; charpoly_gap is the algebraic twin."""
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def charpoly_gap(ev, tgt):
    """Direct |Delta(x;gamma)| proxy: L1 distance between the two characteristic
    polynomials' coefficient vectors (elementary symmetric functions). Equal multisets
    <=> identical coeffs. Normalised by coefficient scale to be gamma-comparable."""
    pa = np.poly(ev)   # monic char poly coeffs of L
    pb = np.poly(tgt)  # of -L-2sigma
    scale = np.max(np.abs(pa)) + np.max(np.abs(pb)) + 1e-300
    return float(np.sum(np.abs(pa - pb)) / scale)


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
        if cls == "hard":
            pairs.append((t1, t2))
    print(f"hard windowed pairs (N=4, Z): {len(pairs)}")

    rng = np.random.default_rng(20260604)
    fine = list(np.linspace(0.005, 2.0, 400))
    rand = list(rng.uniform(1e-3, 2.0, 300))
    physical = [1.0, 0.5, 0.25, 0.1, 0.05, 0.075, 2.0, 1.5]
    gammas = sorted(set(fine + rand + physical))

    SET_TOL = 1e-6
    global_min_setdist = np.inf
    global_min_at = None
    suspicious = []  # (pair, gamma, setdist, otbreak, cpgap)

    for (t1, t2) in pairs:
        lab = "".join(t1) + "+" + "".join(t2)
        local_min = np.inf
        local_at = None
        for g in gammas:
            L, sig = L_of([t1, t2], g)
            ev = np.linalg.eigvals(L)
            tgt = -ev - 2 * sig
            sd = set_distance(ev, tgt)
            if sd < local_min:
                local_min = sd
                local_at = g
            if sd < SET_TOL and g > 1e-3:
                suspicious.append((lab, g, sd, ot_break(ev, tgt), charpoly_gap(ev, tgt)))
        if local_min < global_min_setdist:
            global_min_setdist = local_min
            global_min_at = (lab, local_at)
        # report each pair's closest approach to spec-exact
        print(f"  {lab:14s}  min set-distance over sweep = {local_min:.4e}  at gamma={local_at:.4f}")

    print()
    print(f"GLOBAL min set-distance over all hard pairs & all gamma: {global_min_setdist:.4e}")
    print(f"  achieved by {global_min_at[0]} at gamma={global_min_at[1]:.5f}")
    print(f"Suspicious (set-distance < {SET_TOL} at gamma>1e-3): {len(suspicious)}")
    for s in suspicious[:20]:
        print(f"    {s[0]:14s} gamma={s[1]:.5f}  setdist={s[2]:.2e}  OTbreak={s[3]:.2e}  cpgap={s[4]:.2e}")

    if global_min_setdist > SET_TOL:
        print("\nVERDICT(attack #3/#4): NO hard pair spec-exact anywhere on the dense+random sweep.")
        print("  The closest any hard pair comes to spec-restoration is the GLOBAL min above;")
        print("  it is far from machine zero, so 'generic gamma' is not a loophole here.")
    else:
        print("\nVERDICT(attack #3/#4): POTENTIAL HOLE -- a hard pair approaches spec-exact; inspect above.")


if __name__ == "__main__":
    main()
