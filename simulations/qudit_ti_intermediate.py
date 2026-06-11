#!/usr/bin/env python3
"""
_qudit_ti_intermediate.py - WIP scout: is there a TI layer between the
product cap (2d)^N and the combinatorial ceiling of the qudit partial palindrome?

F121 / the product-mirror cap (PROOF §6):
  product per-site mirror  -> pairs at most (2d)^N   (local, a theorem)
  global partial isometry  -> pairs the ceiling Sum_h min(c_h, c_{N-h})
  gap = ceiling - (2d)^N   = the non-product part (d=3,N=2: 54 - 36 = 18)

THE QUESTION (the F116 seam, inverted): does a TRANSLATION-INVARIANT but
non-product mirror exceed the product cap, and if so does it reach the
ceiling? TI is strictly weaker than product (product-with-equal-sites IS
TI, e.g. Pi_d = rho^T Shift^N), so TI-rank lies in [(2d)^N, ceiling].

METHOD (exact, no search): the palindrome intertwiner W L_D = (-L_D - 2N g) W
forces W to be block-anti-diagonal in the Hamming grading: W maps rung h to
rung N-h. The full (unconstrained) intertwiner's generic rank is the ceiling
Sum_h min(c_h, c_{N-h}). Impose translation invariance [W, T] = 0 (T the
cyclic site shift): allowed entries collapse into T-orbits, one free
coefficient each. A GENERIC TI intertwiner attains the maximal TI rank
(rank is generic on a linear space). Compare:
   (2d)^N  (product cap)  <=  TI generic rank  <=  ceiling.

If TI rank == ceiling: the non-product part is translation-invariant (the
"non-locality" is non-product-ness, not non-TI-ness) - the clean answer.
If (2d)^N < TI < ceiling: a genuine third layer between local and global.
If TI == (2d)^N: TI buys nothing beyond product.

Self-validating where the targets are known (ceiling, product cap).
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
    c = [0] * (N + 1)
    for i in iprod(range(d), repeat=N):
        for j in iprod(range(d), repeat=N):
            c[hamming(i, j)] += 1
    return sum(min(c[h], c[N - h]) if h != N - h else c[h]
               for h in range(N + 1)) - sum(min(c[h], c[N - h]) for h in range(N + 1) if h < N - h)


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
    """Generic rank of a translation-invariant palindrome intertwiner."""
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
    """Generic rank of a general (non-TI) palindrome intertwiner = the ceiling."""
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


def main():
    print("=" * 72)
    print("Is there a TI layer between the product cap (2d)^N and the ceiling?")
    print("=" * 72)
    print(f"  {'d':>2}{'N':>2}{'(2d)^N':>9}{'TI rank':>9}{'ceiling':>9}{'d^2N':>8}   verdict")
    rows = []
    for d, N in ((2, 2), (3, 2), (3, 3), (4, 2), (2, 3)):
        cap = (2 * d) ** N
        ti, n_orb = ti_intertwiner_rank(d, N)
        ceil, c = ceiling_clean(d, N)
        unc = unconstrained_rank(d, N)
        assert unc == ceil, f"unconstrained generic rank {unc} != ceiling {ceil} (d={d},N={N})"
        assert cap <= ti <= ceil, f"TI rank {ti} outside [{cap},{ceil}] (d={d},N={N})"
        if ti == ceil and cap < ceil:
            verdict = "TI REACHES ceiling (non-product part is TI)"
        elif ti == cap and cap < ceil:
            verdict = "TI = product (TI buys nothing)"
        elif cap < ti < ceil:
            verdict = "THIRD LAYER (TI strictly between)"
        else:
            verdict = "d=2: cap=ceiling (full)"
        rows.append((d, N, cap, ti, ceil, verdict))
        print(f"  {d:>2}{N:>2}{cap:>9}{ti:>9}{ceil:>9}{d**(2*N):>8}   {verdict}")

    print()
    # sanity: d=2 everything coincides (full mirror)
    for (d, N, cap, ti, ceil, _) in rows:
        if d == 2:
            assert cap == ti == ceil == (2 * d) ** N, f"d=2 should be full: {(cap,ti,ceil)}"
    print("d=2 columns: (2d)^N = TI = ceiling = d^{2N} (the full mirror). OK")
    print()
    print("READING: the verdict column says whether translation invariance alone")
    print("recovers the non-product part of the partial palindrome, or whether a")
    print("genuine intermediate (TI-but-below-ceiling) layer exists.")
    print("=" * 72)


if __name__ == "__main__":
    main()
