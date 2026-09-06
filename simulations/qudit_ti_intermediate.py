#!/usr/bin/env python3
"""
qudit_ti_intermediate.py - finite TI ranks versus the shift-aligned construction
and the combinatorial ceiling of the qudit partial palindrome.

F121 / operator construction (PROOF §6, cap claim retracted):
  shift-aligned Pi_d      -> pairs (2d)^N
  global partial isometry  -> pairs the ceiling Sum_h min(c_h, c_{N-h})
  gap = ceiling - (2d)^N   = gap above that construction, not a locality split

THE QUESTION (the F116 seam, inverted): does a TRANSLATION-INVARIANT but
non-product mirror exceed the shift construction, and if so does it reach the
ceiling? TI is strictly weaker than product (product-with-equal-sites IS
TI, e.g. the restricted map Pi_d P_aligned), so the maximal TI rank lies in
[(2d)^N, ceiling].

METHOD (finite numerical construction, no optimisation search): the palindrome intertwiner W L_D = (-L_D - 2N g) W
forces W to be block-anti-diagonal in the Hamming grading: W maps rung h to
rung N-h. The full (unconstrained) block space has the ceiling
Sum_h min(c_h, c_{N-h}) as its maximum rank. Impose translation invariance [W, T] = 0 (T the
cyclic site shift): allowed entries collapse into T-orbits, one free
coefficient each. A generic point of this linear space attains its maximal
rank, but the finite computation below only exhibits one seeded point. Compare:
   (2d)^N  (shift-aligned rank)  <=  computed TI rank  <=  ceiling.

If the computed TI rank == ceiling: the finite run exhibits a translation-invariant representative (the
"non-locality" is non-product-ness, not non-TI-ness) - the clean answer.
If (2d)^N < computed TI < ceiling: the sampled construction is strictly between.
If computed TI == (2d)^N: this sample buys nothing beyond the shift construction.

The cited finite numerical results are pinned literally and a mutation to the product
construction is required to fail their gate. This is evidence for those constructions,
not an exact or general maximal-rank proof.
"""

import numpy as np
from itertools import product as iprod
from math import comb


def hamming(i, j):
    return sum(1 for a, b in zip(i, j) if a != b)


def cyc_shift(t):
    """cyclic site shift (t0,t1,...,t_{N-1}) -> (t1,...,t_{N-1},t0)."""
    return t[1:] + t[:1]


def ceiling(d, N):
    """Public paired-mode ceiling, including both members of each mirrored rung pair."""
    return ceiling_clean(d, N)[0]


def ceiling_clean(d, N):
    c = [0] * (N + 1)
    for i in iprod(range(d), repeat=N):
        for j in iprod(range(d), repeat=N):
            c[hamming(i, j)] += 1
    paired = 0
    for h in range(N + 1):
        m = N - h
        if h < m:
            paired += 2 * min(c[h], c[m])
        elif h == m:
            paired += c[h]
    return paired, c


def ti_intertwiner_rank(d, N, seed=0, tol=1e-9):
    """Numerical SVD rank of one seeded translation-invariant intertwiner."""
    states = list(iprod(range(d), repeat=N))
    coh = [(i, j) for i in states for j in states]   # d^{2N} coherences
    idx = {c: k for k, c in enumerate(coh)}
    D = len(coh)

    # allowed entries: W[out, in] != 0 only if h(out) = N - h(in)
    allowed = []
    for b, (ib, jb) in enumerate(coh):
        hb = hamming(ib, jb)
        for a, (ia, ja) in enumerate(coh):
            if hamming(ia, ja) == N - hb:
                allowed.append((a, b))

    # T action on coherence index: (i,j) -> (shift i, shift j)
    def Tcoh(c):
        i, j = c
        return (cyc_shift(i), cyc_shift(j))

    # group allowed (a,b) into T-orbits; one free coefficient per orbit
    seen = set()
    orbits = []
    for (a, b) in allowed:
        if (a, b) in seen:
            continue
        orb = []
        ca, cb = coh[a], coh[b]
        for _ in range(N):
            key = (idx[ca], idx[cb])
            if key not in seen:
                seen.add(key)
                orb.append(key)
            ca, cb = Tcoh(ca), Tcoh(cb)
        orbits.append(orb)

    rng = np.random.default_rng(seed)
    W = np.zeros((D, D), dtype=complex)
    for orb in orbits:
        val = rng.standard_normal() + 1j * rng.standard_normal()
        for (a, b) in orb:
            W[a, b] = val
    s = np.linalg.svd(W, compute_uv=False)
    rank = int(np.sum(s > tol * s[0]))
    return rank, len(orbits)


def unconstrained_rank(d, N, seed=0, tol=1e-9):
    """Numerical SVD rank of one seeded unrestricted intertwiner."""
    states = list(iprod(range(d), repeat=N))
    coh = [(i, j) for i in states for j in states]
    D = len(coh)
    h = np.array([hamming(i, j) for (i, j) in coh])
    rng = np.random.default_rng(seed)
    W = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D)))
    mask = (h[:, None] == (N - h[None, :]))
    W *= mask
    s = np.linalg.svd(W, compute_uv=False)
    return int(np.sum(s > tol * s[0]))


CITED_FINITE_RESULTS = {(3, 2): 54, (3, 3): 378, (4, 2): 128}


def validate_cited_finite_rows(rows):
    """Pin the three finite numerical constructions cited by the proof."""
    observed = {(d, N): ti for d, N, _cap, ti, _ceil in rows}
    for key, expected in CITED_FINITE_RESULTS.items():
        actual = observed.get(key)
        assert actual == expected, (
            f"cited finite TI construction {key}: expected rank {expected}, got {actual}"
        )


def main():
    print("=" * 72)
    print("Finite TI rank versus shift-aligned rank (2d)^N and the ceiling")
    print("=" * 72)
    print(f"  {'d':>2}{'N':>2}{'(2d)^N':>9}{'TI rank':>9}{'ceiling':>9}{'d^2N':>8}   verdict")
    rows = []
    for d, N in ((2, 2), (3, 2), (3, 3), (4, 2), (2, 3)):
        cap = (2 * d) ** N
        ti, n_orb = ti_intertwiner_rank(d, N)
        ceil, c = ceiling_clean(d, N)
        unc = unconstrained_rank(d, N)
        assert unc == ceil, f"seeded unrestricted numerical rank {unc} != ceiling {ceil} (d={d},N={N})"
        assert cap <= ti <= ceil, f"TI rank {ti} outside [{cap},{ceil}] (d={d},N={N})"
        if ti == ceil and cap < ceil:
            verdict = "TI REACHES ceiling"
        elif ti == cap and cap < ceil:
            verdict = "TI = product (TI buys nothing)"
        elif cap < ti < ceil:
            verdict = "THIRD LAYER (TI strictly between)"
        else:
            verdict = "d=2: shift=ceiling (full)"
        rows.append((d, N, cap, ti, ceil))
        print(f"  {d:>2}{N:>2}{cap:>9}{ti:>9}{ceil:>9}{d**(2*N):>8}   {verdict}")

    print()
    # sanity: d=2 everything coincides (full mirror)
    for (d, N, cap, ti, ceil) in rows:
        if d == 2:
            assert cap == ti == ceil == (2 * d) ** N, f"d=2 should be full: {(cap,ti,ceil)}"
    validate_cited_finite_rows(rows)
    mutated = list(rows)
    index = next(i for i, row in enumerate(mutated) if row[:2] == (3, 2))
    d, N, cap, _ti, ceil = mutated[index]
    mutated[index] = (d, N, cap, cap, ceil)
    try:
        validate_cited_finite_rows(mutated)
    except AssertionError:
        pass
    else:
        raise AssertionError("mutation control failed: product-rank substitution stayed green")
    print("d=2 columns: (2d)^N = TI = ceiling = d^{2N} (the full mirror). OK")
    print()
    print("READING: the verdict column says whether translation invariance alone")
    print("reaches the ceiling in these seeded finite numerical constructions. It does not prove")
    print("generic maximal rank, and it does not classify the")
    print("optimum over product intertwiners.")
    print("=" * 72)


if __name__ == "__main__":
    main()
