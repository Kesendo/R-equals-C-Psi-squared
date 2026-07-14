"""F129: the level-collision law: the triple level map is injective away from 3|n, 10|n.

STATEMENT (found 2026-07-14 late, the dark-fringe scanner of the F128 arc; proved in
docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md):

  For n >= 5 call a triple tau = {k1<k2<k3} in {1..n-1} CLEAN if no internal pair is
  balanced (k_i + k_j != n).  Two distinct clean triples tau != sigma with EQUAL level
  S = Sum cos(k_i pi/n) exist  ==>  3|n or 10|n  (proved modulo ONE named corner
  family, proof doc SS4: a hypothetical minimal weight-12 sum whose smallest
  ratio-order primes are {2,5,7}, forcing odd 35|n; this census IS its certificate
  through n = 210, next candidate n = 245).  Conversely such pairs exist at every
  3|n >= 9 and every 10|n >= 20 (exhibited per n <= 210 by G3).  Equivalently: away
  from {3|n or 10|n} the level map is INJECTIVE on clean triples: tens of thousands
  of levels, zero accidental collisions.

MECHANISMS: 3|n carries rotated 3-cycle relations (and, at 15|n, the Conway-Jones
sporadics, absorbed); 10|n carries the pentagon + zero-mode family
(a conjugate pair of rotated R5 root-cycles, 10 roots, plus the self-antipodal zero
mode n/2, whose two roots +-i form the one permitted R2).  In the 3-cosine world of
the seed arc the pentagon was chained to 3 (the 15|n sporadic); at six cosines it is
unchained: 10|n, no 3 required.

REDUCTION (SS2 of the proof): S(tau) = S(sigma) <=> the six angles
{k_i} u {n - l_j} have vanishing cosine sum, and "clean + distinct" maps exactly to
"no antipodal pair e + f = n in the 6-set, except at most one self-antipodal zero
mode n/2".  So the law is the classification of vanishing sums of SIX comb cosines,
one weight level above the seed arc's three (Seed.cs, Conway-Jones Thm 6).

GATES (exit 0 iff all pass):
  [G1] the reduction identity: equal-sum <-> 6-cosine zero-sum, exact in Z[x]/Phi_2n,
       verified on the committed example pairs and on random clean pairs.
  [G2] the census certificate, n <= N_HI: at every n with 3∤n and 10∤n, ALL clean-triple
       levels are pairwise distinct modulo a prime p = 1 mod 2n (exact equality over
       Q(zeta) implies equality mod p, so distinctness mod p PROVES injectivity); any
       mod-p candidate collision at such n is refuted exactly in the cyclotomic basis.
  [G3] the converse: at every firing n (3|n >= 9, 10|n >= 20, n <= N_HI) at least one
       exact collision pair is exhibited (confirmed in Z[x]/Phi_2n, no floats).
  [G4] the law equality: {n <= N_HI with confirmed collisions} ==
       {3|n, n >= 9} u {10|n, n >= 20}.
  [G5] the two mechanism anchors, exact: n = 15, (8,12,14) ~ (9,11,13) at a NONZERO
       shared level, its 12 roots partitioning into FOUR rotated 3-cycles; n = 20,
       (1,7,9) ~ (3,5,10) with its 6-set decomposing into an R5 conjugate pair + the
       zero mode (the decompositions are exhibited term by term over Z[x]/Phi_2n).

Run: python simulations/f129_level_collision_law.py            (~40 s, N_HI = 210)
     python simulations/f129_level_collision_law.py --fast     (N_HI = 96, ~3 s)
"""
import itertools
import sys

import numpy as np
import sympy as sp

N_LO = 5
N_HI = 96 if "--fast" in sys.argv else 210
X = sp.symbols("x")
_CYCLO = {}


# ------------------------------------------------------------------ exact layer
def cyclo(m):
    if m not in _CYCLO:
        _CYCLO[m] = sp.Poly(sp.cyclotomic_poly(m, X), X)
    return _CYCLO[m]


def root_sum_vec(n, exps):
    """Sum of zeta^e over exps (zeta of order 2n), exact vector mod Phi_2n."""
    m = 2 * n
    cy = cyclo(m)
    v = np.zeros(cy.degree(), dtype=np.int64)
    for e in exps:
        r = sp.Poly(X ** (e % m), X).rem(cy)
        for mono, c in r.terms():
            v[mono[0]] += int(c)
    return v


def level_vec(n, t):
    """2 * S(t) as an exact cyclotomic vector."""
    m = 2 * n
    return root_sum_vec(n, [e for k in t for e in (k, m - k)])


def is_clean(n, t):
    return t[0] + t[1] != n and t[0] + t[2] != n and t[1] + t[2] != n


# ------------------------------------------------------------------ mod-p layer
def find_prime_and_root(m, seed=0):
    k = (1 << 40) // m + seed * 7919
    while True:
        p = k * m + 1
        if sp.isprime(p):
            fac = sp.factorint(m)
            for h in range(2, 200):
                g = pow(h, (p - 1) // m, p)
                if g != 1 and all(pow(g, m // q, p) != 1 for q in fac):
                    return p, g
        k += 1


def collision_groups_mod_p(n, seed=0):
    """groups of clean triples sharing a level fingerprint mod p (candidates only;
    distinctness mod p is PROOF of exact distinctness)."""
    m = 2 * n
    p, g = find_prime_and_root(m, seed)
    c = np.zeros(n, dtype=np.int64)
    for k in range(1, n):
        c[k] = (pow(g, k, p) + pow(g, m - k, p)) % p
    idx = np.arange(n)
    S = (c[:, None, None] + c[None, :, None] + c[None, None, :]) % p
    i, j, l = np.meshgrid(idx, idx, idx, indexing="ij")
    mask = (i < j) & (j < l) & (i >= 1) & \
           (i + j != n) & (i + l != n) & (j + l != n)
    vals, trips = S[mask], np.stack([i[mask], j[mask], l[mask]], axis=1)
    order = np.argsort(vals, kind="stable")
    vals, trips = vals[order], trips[order]
    groups, a = [], 0
    while a < len(vals):
        b = a + 1
        while b < len(vals) and vals[b] == vals[a]:
            b += 1
        if b - a >= 2:
            groups.append([tuple(int(x) for x in t) for t in trips[a:b]])
        a = b
    return groups


# ------------------------------------------------------------------ gates
def gate1_reduction():
    rng = np.random.default_rng(129)
    checked = 0
    for _ in range(200):
        n = int(rng.integers(7, 60))
        t = tuple(sorted(rng.choice(range(1, n), 3, replace=False).tolist()))
        s = tuple(sorted(rng.choice(range(1, n), 3, replace=False).tolist()))
        if t == s or not (is_clean(n, t) and is_clean(n, s)):
            continue
        m = 2 * n
        equal = not np.any(level_vec(n, t) - level_vec(n, s))
        six = list(t) + [n - l for l in s]
        zero6 = not np.any(root_sum_vec(n, [e for f in six for e in (f, m - f)]))
        assert equal == zero6, (n, t, s)
        checked += 1
    print(f"[G1] reduction: equal-sum <-> 6-cosine zero-sum, exact, {checked} random pairs  PASS")


def gate2_to_4_census():
    firing, silent_pred = [], []
    for n in range(N_LO, N_HI + 1):
        predicted = (n % 3 == 0 and n >= 9) or (n % 10 == 0 and n >= 20)
        groups = collision_groups_mod_p(n)
        if not groups:
            assert not predicted or n in (6, 10) or True, n
            if predicted:
                silent_pred.append(n)
            continue
        if not predicted:
            # law-break candidate: refute every group exactly
            for grp in groups:
                exact = {}
                for t in grp:
                    exact.setdefault(level_vec(n, t).tobytes(), []).append(t)
                for ts in exact.values():
                    assert len(ts) == 1, f"LAW BREAK at n = {n}: {ts}"
            continue
        # firing n: confirm at least one collision exactly (G3).  The [:20] cap is a
        # runtime bound with a false-FAIL risk only: exact collisions always appear
        # inside some mod-p group, so a passing run is sound regardless of the cap.
        confirmed = False
        for grp in groups[:20]:
            exact = {}
            for t in grp:
                exact.setdefault(level_vec(n, t).tobytes(), []).append(t)
            if any(len(ts) >= 2 for ts in exact.values()):
                confirmed = True
                break
        assert confirmed, f"predicted n = {n} has candidates but no exact collision"
        firing.append(n)
    predicted_set = [n for n in range(N_LO, N_HI + 1)
                     if (n % 3 == 0 and n >= 9) or (n % 10 == 0 and n >= 20)]
    assert firing == predicted_set, (firing, predicted_set)
    assert not silent_pred, silent_pred
    print(f"[G2] injectivity PROVED (distinct mod p, candidates refuted exactly) at all "
          f"{N_HI - N_LO + 1 - len(firing)} non-firing n <= {N_HI}  PASS")
    print(f"[G3] one exact collision exhibited at every firing n ({len(firing)} values)  PASS")
    print(f"[G4] law equality: firing set == {{3|n >= 9}} u {{10|n >= 20}}  PASS")


def gate5_mechanism_anchors():
    # n = 15: (8,12,14) ~ (9,11,13), a NONZERO shared level (~ -1.892); the 6-set
    # {8,12,14} u {15-9, 15-11, 15-13} = {2,4,6,8,12,14} has +-root set that
    # partitions into FOUR rotated R3 root-cycles (step m/3 = 10), each exactly zero.
    n, m = 15, 30
    t, s = (8, 12, 14), (9, 11, 13)
    assert not np.any(level_vec(n, t) - level_vec(n, s))
    assert np.any(level_vec(n, t)), "anchor level must be nonzero (else it is F127's world)"
    r3s = [(2, 12, 22), (8, 18, 28), (4, 14, 24), (6, 16, 26)]
    for r in r3s:
        assert (r[1] - r[0]) % m == 10 and (r[2] - r[1]) % m == 10
        assert not np.any(root_sum_vec(n, list(r))), r
    six = [8, 12, 14, 15 - 9, 15 - 11, 15 - 13]
    exps = sorted(e % m for f in six for e in (f, m - f))
    assert sorted(e % m for r in r3s for e in r) == exps
    print("[G5a] n = 15 anchor: (8,12,14) ~ (9,11,13) at nonzero level; 12 roots = four "
          "rotated R3 cycles, each exactly zero, partition exact  PASS")
    # n = 20: (1,7,9) ~ (3,5,10); 6-set = {1,7,9,17,15,10}; R5-pair + zero mode.
    n, m = 20, 40
    t, s = (1, 7, 9), (3, 5, 10)
    assert not np.any(level_vec(n, t) - level_vec(n, s))
    r5a = (1, 9, 17, 25, 33)          # step m/5 = 8: rotated R5
    r5b = tuple((m - e) % m for e in r5a)   # its conjugate
    r2 = (10, 30)                     # the zero mode: zeta^10 = i, zeta^30 = -i
    assert not np.any(root_sum_vec(n, list(r5a)))
    assert not np.any(root_sum_vec(n, list(r5b)))
    assert not np.any(root_sum_vec(n, list(r2)))
    six = [1, 7, 9, 20 - 3, 20 - 5, 20 - 10]
    exps = sorted(e % m for f in six for e in (f, m - f))
    assert sorted(list(r5a) + list(r5b) + list(r2)) == exps, \
        (sorted(list(r5a) + list(r5b) + list(r2)), exps)
    print("[G5b] n = 20 anchor: (1,7,9) ~ (3,5,10), 6-set = R5 conjugate pair + zero mode, "
          "each exactly zero  PASS")


if __name__ == "__main__":
    gate1_reduction()
    gate5_mechanism_anchors()
    gate2_to_4_census()
    print(f"ALL GATES PASS: the level map is injective on clean triples away from "
          f"{{3|n or 10|n}}, n <= {N_HI}; converses exhibited; mechanisms anchored.")
