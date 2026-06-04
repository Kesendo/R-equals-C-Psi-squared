#!/usr/bin/env python3
"""F87 spec-B probe 11: is the omega=0 block ALONE the discriminator? (first-order, all 42 pairs)

block_localize showed omega=0 ALWAYS breaks for the 2 hard reference pairs. If the omega=0 block's
(-N)-asymmetry is nonzero  <=>  hard, over ALL 42 pairs, then the first-order converse localizes
entirely to the STATIC block (populations + degenerate-energy coherences), the cleanest possible
target. Test:

  for each of the 42 windowed (N=4,k=3,Z-deph (0,1)) y_par-homog Mixed+Mixed pairs:
    - asym of the FULL omega=0 block about -N
    - asym of the WHOLE first-order spectrum about -N (all blocks)
    - F87 class
  Confirm:  omega=0-block-asymmetric  <=>  hard  (and whether omega!=0 blocks ever break alone).
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


def is_mixed(t):
    return any(L not in DIAG for L in t)


def asym_about(vals, center):
    vals = np.asarray(vals).real
    tgt = 2.0 * center - vals
    cost = np.abs(vals[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def all_blocks(pair, N, tol=1e-6):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
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
    return out


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]

    n_pairs = 0
    o0_eq_hard = 0
    nonzero_omega_breaks_alone = 0
    rows = []
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        n_pairs += 1
        blk = all_blocks([t1, t2], N)
        a0 = asym_about(np.linalg.eigvals(blk[0.0]).real, -N) if 0.0 in blk else 0.0
        # asym of nonzero-omega blocks (max over omega!=0)
        anz = 0.0
        for om, M in blk.items():
            if abs(om) > 1e-9:
                anz = max(anz, asym_about(np.linalg.eigvals(M).real, -N))
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter='Z')
        hard = (cls == 'hard')
        o0_break = a0 > 1e-6
        o0_eq_hard += (o0_break == hard)
        # does a nonzero-omega block break while omega=0 does NOT?
        if anz > 1e-6 and not o0_break:
            nonzero_omega_breaks_alone += 1
        rows.append((t1, t2, hard, a0, anz))

    print("=" * 80)
    print("First-order converse localization: is the omega=0 block alone the discriminator?")
    print(f"  (N=4, k=3, Z-deph (0,1); {n_pairs} y_par-homog Mixed+Mixed pairs)")
    print("=" * 80)
    print(f"  omega=0-block-asymmetric  <=>  hard :  {o0_eq_hard}/{n_pairs}  "
          f"{'ALL' if o0_eq_hard == n_pairs else 'MISMATCH'}")
    print(f"  pairs where a nonzero-omega block breaks but omega=0 does NOT: {nonzero_omega_breaks_alone}")
    if o0_eq_hard == n_pairs:
        print("  => the FIRST-ORDER converse localizes ENTIRELY to the omega=0 (static) block.")
        print("     Target: odd cycle => omega=0 block spectrum asymmetric about -N.")
    # show the omega=0 asym for hard vs soft
    hard_a0 = sorted([r[3] for r in rows if r[2]])
    soft_a0 = max(r[3] for r in rows if not r[2])
    print(f"\n  omega=0 asym:  hard in [{min(hard_a0):.3f}, {max(hard_a0):.3f}],  "
          f"soft max = {soft_a0:.1e}")


if __name__ == "__main__":
    main()
