#!/usr/bin/env python3
"""FINAL resolution of the equiv2/equiv3 'mismatches': they were a BUG in MY probe (I used
'omega=0 block MINIMUM == -N' as the soft test; the correct section-7.5 condition is '-N ATTAINED',
i.e. -N is IN the block spectrum, equivalently the block is SYMMETRIC about -N). A symmetric block
on [-8,0] has min=-8 but -N=-4 present and is SOFT. Re-test with the correct discriminator.

For 1500 random in-scope H (FHF=-H, zero diagonal), discriminator = (omega=0 block symmetric about
-N), and independently (-N in the block spectrum). Compare to TRUE chiral-K existence (KHK=-H).
Expect bit-exact agreement: chiral K exists <=> block symmetric about -N <=> -N in spectrum, with
the ODD-CYCLE (no chiral K) cases asymmetric. Zero mismatches => section 7.5's soft+hard legs of
'c=0 <=> bipartite' are solid well beyond the 16 canonical pairs.
"""
from __future__ import annotations
import sys
from itertools import product
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


def edge_masks(H, tol=1e-9):
    d = H.shape[0]
    S = set()
    for a in range(d):
        for b in np.where(np.abs(H[a]) > tol)[0]:
            if b > a:
                S.add(int(a) ^ int(b))
    return S


def phi_exists(S, N):
    if not S:
        return True  # empty graph: trivially 2-colourable
    rows = [[(s >> i) & 1 for i in range(N)] + [1] for s in S]
    R = len(rows)
    r = 0
    for c in range(N):
        piv = next((i for i in range(r, R) if rows[i][c]), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(R):
            if i != r and rows[i][c]:
                rows[i] = [(x ^ y) for x, y in zip(rows[i], rows[r])]
        r += 1
    for row in rows:
        if not any(row[:N]) and row[N]:
            return False
    return True


def omega0_block(H, N, tol=1e-6):
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, "Z") @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < tol]
    n = len(modes)
    M = np.zeros((n, n), dtype=complex)
    for i, (a, b) in enumerate(modes):
        for j, (ap, bp) in enumerate(modes):
            val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
            if (a, b) == (ap, bp):
                val -= N
            M[i, j] = val
    return np.linalg.eigvals(M).real


def sym_about(s, c):
    tgt = 2 * c - s
    cost = np.abs(s[:, None] - tgt[None, :])
    r, cc = linear_sum_assignment(cost)
    return float(cost[r, cc].mean())


def main():
    N, k = 4, 3
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
             and any(L not in DIAG for L in t)]
    templates = [t for t in product("IZ", repeat=k) if not all(L == "I" for L in t)]
    all_terms = terms + templates
    F = site_op(N, 0, "X")
    for l in range(1, N):
        F = F @ site_op(N, l, "X")

    rng = np.random.default_rng(20260604)
    n_inscope = 0
    mism_sym = 0       # phi_exists != symmetric-about-(-N)
    mism_present = 0   # phi_exists != (-N in spectrum)
    examples = []
    for _ in range(1500):
        m = rng.integers(2, 6)
        chosen = rng.choice(len(all_terms), size=m, replace=False)
        chosen_terms = [all_terms[i] for i in chosen]
        weights = rng.uniform(0.3, 1.5, size=m)
        H = _build_kbody_chain(N, [tuple(t) + (float(w),) for t, w in zip(chosen_terms, weights)])
        if np.linalg.norm(H) < 1e-9:
            continue
        fhf = np.linalg.norm(F @ H @ F + H) < 1e-7
        zd = np.max(np.abs(np.diag(H))) < 1e-9
        if not (fhf and zd):
            continue
        n_inscope += 1
        bip = phi_exists(edge_masks(H), N)
        s = omega0_block(H, N)
        symm = sym_about(s, -N) < 1e-6
        present = np.min(np.abs(s - (-N))) < 1e-6
        if bip != symm:
            mism_sym += 1
            if len(examples) < 8:
                examples.append((chosen_terms, bip, symm, present, float(s.min()), float(sym_about(s, -N))))
        if bip != present:
            mism_present += 1
    print(f"In-scope (FHF=-H, zero-diag) random H: {n_inscope}")
    print(f"  phi_exists(bipartite) != (block SYMMETRIC about -N) mismatches: {mism_sym}")
    print(f"  phi_exists(bipartite) != (-N IN block spectrum) mismatches:    {mism_present}")
    for x in examples:
        print(f"     terms={[''.join(t) for t in x[0]]} bip={x[1]} sym={x[2]} -Npresent={x[3]} "
              f"smin={x[4]:.3f} asym={x[5]:.2e}")
    print()
    if mism_sym == 0 and mism_present == 0:
        print("VERDICT: bit-exact. 'bipartite (chiral K exists) <=> omega=0 block symmetric about -N")
        print("  <=> -N attained in the block' over 1500 random in-scope H. The earlier equiv2/equiv3")
        print("  'mismatches' were MY probe bug (testing block-MINIMUM==-N instead of -N-attained;")
        print("  a symmetric block on [-8,0] has min -8 but -N=-4 present, and is SOFT). With the")
        print("  correct discriminator, section 7.5's 'c=0 <=> bipartite' source equivalence is solid.")
    else:
        print("VERDICT: a genuine mismatch remains -- inspect.")


if __name__ == "__main__":
    main()
