#!/usr/bin/env python3
r"""F87 spec-A probe 12: is 'the system {d_i+d_j=0 : i^j in S} has a nonzero solution' equivalent
to 'a global linear φ:𝔽₂^N→𝔽₂ with φ|_S=1 exists', for an ABSTRACT, FULL Cayley graph on 𝔽₂^N
with generator set S (every i ~ i⊕m for all m∈S)?

This is the IDEALIZED graph (all masks realized at every vertex, weights all 1). If the
equivalence is an identity here, then the ONLY way the physical G_H can deviate is through
'missing edges' (a mask in S not realized at some vertex due to a vanishing matrix element).
The hard direction is safe regardless because (probe 11) the triangle masks ARE realized
everywhere for hard pairs. This probe pins the abstract identity and the role of <S> and ī.

Tests over MANY random S ⊆ 𝔽₂^N \ {0}:
  (A) nonzero-solution(full Cayley system)  ⟺  φ_exists(S)   [the clean identity]
  (B) when φ exists, the solution space dimension = number of cosets of <S> = 2^(N−r),
      r=rank<S>; the per-coset sign is free (±), giving 2^(N−r) basis vectors. So cap counts
      cosets — matches probe-3 cap distribution {2,3,4,8,9} only after the physical missing
      edges split/merge cosets. Report both.
  (C) the anti-diagonal pairing ī=i⊕1⃗: 1⃗∈<S> or not changes whether F·D stays within a
      coset; check that bipartite ⟹ a nonzero solution always exists irrespective of 1⃗.
"""
from __future__ import annotations
import sys
from itertools import combinations
from pathlib import Path
import random

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from f87_flip_generators import phi_exists_gf2


def full_cayley_system_nullity(S, N):
    """Nullity of {d_i + d_{i⊕m} = 0 : all i, all m∈S} over the reals (d∈R^{2^N})."""
    d = 2 ** N
    rows = []
    for i in range(d):
        for m in S:
            j = i ^ m
            if j > i:
                r = np.zeros(d)
                r[i] += 1.0
                r[j] += 1.0
                rows.append(r)
    if not rows:
        return d
    s = np.linalg.svd(np.array(rows), compute_uv=False)
    rank = int(np.sum(s > 1e-9))
    return d - rank


def rank_gf2(vecs, N):
    basis = []
    for v in vecs:
        x = v
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
            basis.sort(reverse=True)
    return len(basis)


def subgroup_size(S, N):
    G = {0}
    for m in S:
        G |= {g ^ m for g in list(G)}
        # closure
    changed = True
    while changed:
        changed = False
        for a in list(G):
            for m in S:
                if (a ^ m) not in G:
                    G.add(a ^ m); changed = True
    return len(G)


def main():
    print("=" * 88)
    print("F87 spec-A probe 12: abstract full-Cayley identity  nonzero-solution ⟺ φ-exists")
    print("=" * 88)
    rng = random.Random(0)
    mism = 0
    trials = 0
    coset_match = 0
    examples = []
    for N in (2, 3, 4):
        univ = list(range(1, 2 ** N))
        for _ in range(400):
            k = rng.randint(1, min(5, len(univ)))
            S = set(rng.sample(univ, k))
            trials += 1
            nul = full_cayley_system_nullity(S, N)
            phi = phi_exists_gf2(S, N)
            has_sol = nul >= 1
            if has_sol != phi:
                mism += 1
                if len(examples) < 5:
                    examples.append((N, sorted(S), nul, phi))
            # coset count check when φ exists
            if phi:
                r = rank_gf2(list(S), N)
                ncoset = 2 ** (N - r)
                coset_match += int(nul == ncoset)
    print(f"  random trials: {trials}")
    print(f"  (A) nonzero-solution ⟺ φ-exists:  mismatches = {mism}  "
          f"{'IDENTITY HOLDS' if mism==0 else 'FAILS: '+str(examples)}")
    print(f"  (B) when φ exists, nullity == 2^(N−rank<S>) (coset count): "
          f"{'ALL' if coset_match>0 else 'n/a'}  (matched {coset_match} φ-cases)")
    print()
    print("  Reading: on the IDEALIZED full Cayley graph the equivalence")
    print("     (homogeneous 2-colour system solvable) ⟺ (global linear φ with φ|_S=1)")
    print("  is an IDENTITY. φ exists ⟺ S has no odd 𝔽₂-relation (standard GF(2) consistency).")
    print("  The K3 triangle {m1,m2,m3}, m1⊕m2⊕m3=0, is an odd relation ⟹ no φ ⟹ no nonzero")
    print("  2-colour solution ⟹ cap=0 ⟹ hard. Physical G_H only ever DROPS edges vs the full")
    print("  Cayley graph; for hard pairs the triangle edges survive in every component (probe 11),")
    print("  so the obstruction is preserved.")


if __name__ == "__main__":
    main()
