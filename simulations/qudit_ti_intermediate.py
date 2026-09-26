#!/usr/bin/env python3
"""
qudit_ti_intermediate.py - is there a TI layer between the
product cap and the combinatorial ceiling of the qudit partial palindrome?

F121 / the product-mirror cap (PROOF §6):
  product per-site mirror  -> pairs at most P(d, N) = max_m (2d)^(N-2m)(d^3-d^2)^m
                              (= (2d)^N for d <= 5; local, a theorem)
  global partial isometry  -> pairs the ceiling Sum_h min(c_h, c_{N-h})
  gap = ceiling - P        = the non-product part (d=3,N=2: 54 - 36 = 18)

THE QUESTION (the F116 seam, inverted): does a TRANSLATION-INVARIANT but
non-product mirror exceed the product cap, and if so does it reach the
ceiling?

METHOD (exact, no search): the palindrome intertwiner W L_D = (-L_D - 2N g) W
forces W to be block-anti-diagonal in the Hamming grading: W maps rung h to
rung N-h, so NO intertwiner has rank above the ceiling Sum_h min(c_h, c_{N-h}).
Impose translation invariance [W, T] = 0 (T the cyclic site shift): allowed
entries collapse into T-orbits, one free coefficient each. Give every orbit an
integer coefficient and take the rank over GF(p). Reduction mod p can only
lower a rank, so
   rank_p  <=  rank over Q  <=  maximal TI rank  <=  ceiling,
and rank_p == ceiling proves that translation invariance reaches the ceiling.
The same argument on the unrestricted block space shows the ceiling itself is
attained. Compare with the product cap P(d, N).

If TI rank == ceiling: the non-product part is translation-invariant (the
"non-locality" is non-product-ness, not non-TI-ness) - the clean answer.
If P < TI < ceiling: a genuine third layer between local and global.
If TI == P: TI buys nothing beyond product.

Self-validating: the TI rank must equal the ceiling at every computed case,
d = 2 must be full, and a mutation that replaces a TI rank by the product cap
must fail the gate. Whether TI reaches the ceiling at every (d, N) is open;
this is a rank computation per case.
"""

import random
from itertools import product as iprod
from math import comb

P_MOD = 1_000_000_007


def hamming(i, j):
    return sum(1 for a, b in zip(i, j) if a != b)


def cyc_shift(t):
    """cyclic site shift (t0,t1,...,t_{N-1}) -> (t1,...,t_{N-1},t0)."""
    return t[1:] + t[:1]


def product_cap(d, N):
    """P(d, N) = max_m (2d)^(N-2m) (d^3 - d^2)^m, the product-mirror cap (PROOF §6)."""
    return max((2 * d) ** (N - 2 * m) * (d ** 3 - d ** 2) ** m for m in range(N // 2 + 1))


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


def rank_mod_p(rows, p=P_MOD):
    """Rank over GF(p) of a sparse integer matrix given as a list of {col: value} dicts."""
    rows = [{c: v % p for c, v in r.items() if v % p} for r in rows]
    pivots = {}                     # pivot column -> reduced row (leading coefficient 1)
    rank = 0
    for r in rows:
        r = dict(r)
        while r:
            col = min(r)
            if col not in pivots:
                inv = pow(r[col], p - 2, p)
                pivots[col] = {c: v * inv % p for c, v in r.items()}
                rank += 1
                break
            f = r[col]
            for c, v in pivots[col].items():
                nv = (r.get(c, 0) - f * v) % p
                if nv:
                    r[c] = nv
                else:
                    r.pop(c, None)
    return rank


def _coherences(d, N):
    states = list(iprod(range(d), repeat=N))
    coh = [(i, j) for i in states for j in states]   # d^{2N} coherences
    return coh, {c: k for k, c in enumerate(coh)}


def ti_intertwiner_rank(d, N, seed=0):
    """Exact lower bound (rank mod p) on the maximal rank of a translation-invariant
    palindrome intertwiner: integer coefficient per T-orbit of allowed entries."""
    coh, idx = _coherences(d, N)
    ham = [hamming(i, j) for (i, j) in coh]
    rng = random.Random(seed)
    coeff = {}
    rows = [dict() for _ in coh]
    for a, (ia, ja) in enumerate(coh):
        for b, (ib, jb) in enumerate(coh):
            if ham[a] + ham[b] != N:
                continue
            ca, cb, keys = (ia, ja), (ib, jb), []
            for _ in range(N):
                keys.append((idx[ca], idx[cb]))
                ca = (cyc_shift(ca[0]), cyc_shift(ca[1]))
                cb = (cyc_shift(cb[0]), cyc_shift(cb[1]))
            key = min(keys)
            if key not in coeff:
                coeff[key] = rng.randrange(1, P_MOD)
            rows[a][b] = coeff[key]
    return rank_mod_p(rows), len(coeff)


def unconstrained_rank(d, N, seed=0):
    """Rank mod p of an intertwiner with an independent integer on every allowed entry."""
    coh, _ = _coherences(d, N)
    ham = [hamming(i, j) for (i, j) in coh]
    rng = random.Random(seed)
    rows = [{b: rng.randrange(1, P_MOD) for b in range(len(coh)) if ham[a] + ham[b] == N}
            for a in range(len(coh))]
    return rank_mod_p(rows)


CITED_FINITE_RESULTS = {(3, 2): 54, (3, 3): 378, (4, 2): 128}


def validate_cited_finite_rows(rows):
    """The law at the three cited cases: the TI rank mod p equals the ceiling (so it is
    exact), the ceiling is the closed form, and the product cap sits strictly below."""
    observed = {(d, N): (cap, ti, ceil) for d, N, cap, ti, ceil in rows}
    for key, expected in CITED_FINITE_RESULTS.items():
        cap, ti, ceil = observed.get(key, (None, None, None))
        assert ti == ceil == expected and cap < ceil, (
            f"cited finite TI construction {key}: expected rank {expected} = ceiling, "
            f"got TI {ti}, ceiling {ceil}, product cap {cap}"
        )


def main():
    print("=" * 72)
    print("Translation-invariant rank (exact, mod p) vs product cap and ceiling")
    print("=" * 72)
    print(f"  {'d':>2}{'N':>2}{'P(d,N)':>9}{'TI rank':>9}{'ceiling':>9}{'d^2N':>8}   verdict")
    rows = []
    for d, N in ((2, 2), (3, 2), (3, 3), (4, 2), (2, 3)):
        cap = product_cap(d, N)
        ti, n_orb = ti_intertwiner_rank(d, N)
        ceil, c = ceiling_clean(d, N)
        unc = unconstrained_rank(d, N)
        assert unc == ceil, f"unrestricted rank mod p {unc} != ceiling {ceil} (d={d},N={N})"
        assert ti <= ceil, f"TI rank {ti} above the ceiling {ceil} (d={d},N={N}): impossible"
        if ti == ceil and cap < ceil:
            verdict = "TI REACHES ceiling"
        elif ti == cap and cap < ceil:
            verdict = "TI = product (TI buys nothing)"
        elif cap < ti < ceil:
            verdict = "THIRD LAYER (TI strictly between)"
        else:
            verdict = "d=2: product = ceiling (full)"
        rows.append((d, N, cap, ti, ceil))
        print(f"  {d:>2}{N:>2}{cap:>9}{ti:>9}{ceil:>9}{d**(2*N):>8}   {verdict}")

    print()
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
    print("d=2 columns: P = TI = ceiling = d^{2N} (the full mirror). OK")
    print()
    print("READING: rank mod p <= rank <= ceiling, so each TI rank equal to the ceiling")
    print("is exact: at (3,2), (3,3), (4,2) translation invariance recovers the whole")
    print("non-product part and there is no intermediate layer. Whether it does at")
    print("every (d, N) is open.")
    print("=" * 72)


if __name__ == "__main__":
    main()
