#!/usr/bin/env python3
"""ATTACK #2 (cross-block cancellation) + the LOGICAL CORE of section 7.5.

Section 7.5 localizes the break to the omega=0 block specifically:
  soft  <=>  the omega=0 block is centre-symmetric about -N  <=>  -N reflection floor attained
  hard  <=>  omega=0 block asymmetric  <=>  +N population-Perron mode globally unpaired.
The argument: the +N mode (always present, since sum_l Z_l^2 = N*I) pairs ONLY with another
omega=0 mode (partners share omega, 0 = -0), so if -N is absent in the omega=0 block, +N is
unpaired no matter what the other blocks do. This is the crux that makes a PER-BLOCK (indeed
single-block) criterion sufficient.

Potential holes to realise:
  (H2a) A hard pair whose omega=0 block IS symmetric about -N, with the asymmetry living only in
        a nonzero-omega block. That would BREAK the section-7.5 localization (the +N Perron story
        would not be the mechanism), even if the pair is still hard. -> would be HOLE-at-7.5-locus.
  (H2b) A non-bipartite (hard) pair whose TOTAL first-order slope c = 0 by cancellation across
        blocks (so break/gamma -> 0, looking first-order-soft though all-orders-hard). The proof
        does NOT actually need total c; it needs the omega=0 block asymmetric. But if total c can
        be 0 for a hard pair, the f87_break_gamma_scaling 'c -> const != 0' backbone (which the
        scout leans on) would be incomplete. Check both.
  (H2c) A soft (bipartite) pair with SOME block asymmetric about -N (would contradict soft <=>
        every block symmetric, and break the 'c=0 <=> bipartite' equivalence).

For all 42 windowed Mixed+Mixed pairs (16 hard, 26 soft), N=4 Z:
  - omega0_asym  := asymmetry-about-(-N) of the omega=0 block spectrum
  - max_block_asym := max over ALL blocks
  - argmax_omega   := which omega carries the max asymmetry
  - total_c        := size-weighted sum / d^2 (the measured break/gamma)
Report whether: hard <=> omega0_asym>0; soft <=> max_block_asym==0; and whether argmax is ever
a nonzero omega while omega=0 is clean (the H2a hole).
"""
from __future__ import annotations
import sys
from collections import defaultdict
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
from framework.pauli import _build_kbody_chain, site_op

DIAG = {"I", "Z"}


def asym_about(vals, center):
    vals = np.asarray(vals)
    tgt = 2.0 * center - vals
    cost = np.abs(vals[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def block_table(pair, N=4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, "Z") @ V for l in range(N)]
    groups = defaultdict(list)
    for a in range(d):
        for b in range(d):
            groups[round(E[a] - E[b], 6)].append((a, b))
    res = {}
    for omega, modes in groups.items():
        n = len(modes)
        M = np.zeros((n, n), dtype=complex)
        for i, (a, b) in enumerate(modes):
            for j, (ap, bp) in enumerate(modes):
                val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
                if (a, b) == (ap, bp):
                    val -= N
                M[i, j] = val
        s = np.linalg.eigvals(M).real
        res[round(omega, 6)] = (n, asym_about(s, -N), s)
    return res


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

    TOL = 1e-7
    h2a = []  # hard with omega=0 clean but a nonzero-omega block dirty
    h2c = []  # soft with some block dirty
    hard_omega0_all_dirty = True
    soft_all_clean = True
    print(f"{'pair':16s} {'cls':5s} {'om0_asym':>10s} {'max_asym':>10s} {'argmax_om':>10s} {'total_c':>9s}")
    for (t1, t2, cls) in pairs:
        lab = "".join(t1) + "+" + "".join(t2)
        tbl = block_table([t1, t2], N)
        om0 = tbl.get(0.0, (0, 0.0, np.array([])))[1]
        # max over blocks
        argmax_om, max_asym = 0.0, 0.0
        d2 = 0
        total = 0.0
        for om, (n, a, s) in tbl.items():
            d2 += n
            total += a * n
            if a > max_asym:
                max_asym = a
                argmax_om = om
        total_c = total / d2
        print(f"{lab:16s} {cls:5s} {om0:10.4f} {max_asym:10.4f} {argmax_om:10.3f} {total_c:9.4f}")
        if cls == "hard":
            if om0 <= TOL:
                h2a.append((lab, om0, argmax_om, max_asym))
                hard_omega0_all_dirty = False
        if cls == "soft":
            if max_asym > TOL:
                h2c.append((lab, argmax_om, max_asym))
                soft_all_clean = False

    print()
    print(f"H2a (hard with omega=0 block CLEAN, break only at nonzero omega): {len(h2a)}")
    for x in h2a:
        print(f"     {x[0]:16s} om0_asym={x[1]:.2e}  worst at omega={x[2]}  asym={x[3]:.3f}")
    print(f"H2c (soft with SOME block asymmetric about -N): {len(h2c)}")
    for x in h2c:
        print(f"     {x[0]:16s} worst at omega={x[1]}  asym={x[2]:.3f}")
    print()
    print(f"section-7.5 locus claim 'hard => omega=0 block asymmetric': {hard_omega0_all_dirty}")
    print(f"'soft => EVERY block symmetric about -N': {soft_all_clean}")
    if hard_omega0_all_dirty and soft_all_clean:
        print("=> The omega=0 localization holds: every hard pair's omega=0 block is asymmetric")
        print("   (the +N Perron mode is unpaired there), and every soft pair has all blocks")
        print("   symmetric. No cross-block cancellation hole; no nonzero-omega-only break.")
    else:
        print("=> POTENTIAL HOLE in the section-7.5 localization; see H2a/H2c above.")


if __name__ == "__main__":
    main()
