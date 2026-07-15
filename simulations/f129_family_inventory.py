"""F129 family inventory: the collision counts per firing n, as closed forms per family.

STATEMENT (found 2026-07-15 by exact piece-decomposition of every collision pair;
grade: VERIFIED, not derived -- the ROT3-multiplicity precedent; the owning proof of
the underlying law is docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md):

  Every F129 collision pair (two distinct clean triples with equal level) reduces
  (proof SS2) to a weight-8 or weight-12 vanishing sum W of 2n-th roots of unity,
  and W decomposes into minimal pieces from the CDK weight <= 16 classification
  (arXiv:2008.11268, Table 1).  Classifying every pair by its piece decomposition
  (deterministic FIRST-FIT rule: minimal pieces ordered by size then subset index,
  commit to the first piece whose remainder is still tileable, backtrack only on
  failure; the rule is deterministic but tied to the fixed root enumeration order)
  splits the census into FAMILIES, and each family's count per n is an exact
  closed form.  The family LIST (the 13 signatures) is COMPLETE at every n, a
  one-paragraph consequence of the CDK table + the proof's pairing/complement
  facts (the doc's "Why exactly thirteen" section); the COUNTS are verified:

    label  pieces                    door       count(n)
    A      R3+R3+R3+R3      (d=3)    3|n        (20k+1)(k+1),   k = floor((n-9)/6)
    B      zero+R3+R3       (d=2)    6|n        12(3k+2)(k+1),  k = (n-12)/6
    C      zero+R5-pair     (d=3)    10|n       2(n-10)                [pentagon]
    D      (R5:3R3) w8      (d=2)    15|n       12(n-9)
    E      R3+R3+(R5:R3)    (d=3)    15|n       (20n-264)/3 n odd, (20n-324)/3 n even
    F      (R5:R3)+conj     (d=3)    15|n       10n-149 n odd,  10n-275 n even
    G      zero+(R5:R3)     (d=2)    30|n       6(n-8)
    H      (R7:R3) w8       (d=2)    21|n       6(n-9)
    I      (R7:5R3) w12     (d=3)    21|n       60
    J      zero+(R7:3R3)    (d=3)    42|n       60
    K      (R11:R3) w12     (d=3)    33|n       20
    L      zero+(R7:R5)     (d=3)    70|n       20   [the corner-closure's second
                                                      mechanism, live at n = 70]
    M      w12 order-210    (d=3)    105|n      100  [may bundle up to three CDK
                                                      210-types; not sub-classified]

  d = |A| = |B| is the private-mode count (d=2 <=> overlap-1 <=> weight 8, which is
  why every d=2 family's door is divisible by 3, the F129 sub-law); the polynomial
  DEGREE of each count MATCHES the free rotations of the piece kit (two free
  R3-pair rotations ~ quadratic; one free rotation ~ linear; a rigid single
  weight-12 piece ~ constant) -- an OBSERVATION, not derived.  The formulas encode
  the F129 onsets themselves: A(6) = B(6) = C(10) = 0 (narrative; those n sit
  outside the gated range).

CAVEATS (honest grade): the family LABELS are defined by the deterministic
first-fit decomposition, and the closed forms were FIT to that rule's output, so
gate I2 verifies that the counts follow the forms, not that a label means what its
CDK name says; the piece NAMES are (weight, order) annotations, never re-verified
per pair (family F's one level-0 member per firing n decomposes into two
SELF-conjugate weight-6 pieces, not a conjugate pair).  Counts verified exactly on
every firing n <= 140 plus the capstones n = 150 and n = 210 (all doors at once),
NOT derived.  Thin constants, honestly: L fires at three tested n (70, 140, 210),
M at TWO (105, 210); M may also bundle up to three distinct CDK order-210 types.
Level-0 (F127-world) pairs are counted (F129 asks equal level; equal-to-zero
qualifies): e.g. at n = 30 family A holds 28 level-0 pairs of its 244.

GATES (exit 0 iff all pass):
  [I1] coverage: every exact collision pair at every checked firing n decomposes
       into inventory pieces (first-fit partition succeeds; no UNRESOLVED).
  [I2] the closed forms: per-family counts == the table above, every checked n.
  [I3] completeness: the families partition the census (sum == total pairs, and
       no family key outside the 13 labels appears).
  [I4] the doors: a family is PRESENT at n iff its door divides n and its count
       formula is nonzero (both directions).

Run: python simulations/f129_family_inventory.py            (firing n <= 140, ~15 min)
     python simulations/f129_family_inventory.py --fast     (firing n <= 66, ~1 min)
     python simulations/f129_family_inventory.py --deep     (adds n = 150 and 210;
                                                             the 210 point dominates,
                                                             ~40 min total)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from f129_level_collision_law import (collision_groups_mod_p, find_prime_and_root,
                                      level_vec)

# ------------------------------------------------------------------ the inventory
def _k6(n, off):
    return (n - off) // 6


FAMILIES = {
    # label: (door divisor(s) as predicate, d, piece signature, count closed form)
    "A": (lambda n: n % 3 == 0, 3, ((3, 3), (3, 3), (3, 3), (3, 3)),
          lambda n: (20 * _k6(n, 9) + 1) * (_k6(n, 9) + 1)),
    "B": (lambda n: n % 6 == 0, 2, ((2, 2), (3, 3), (3, 3)),
          lambda n: 12 * (3 * _k6(n, 12) + 2) * (_k6(n, 12) + 1)),
    "C": (lambda n: n % 10 == 0, 3, ((2, 2), (5, 5), (5, 5)),
          lambda n: 2 * (n - 10)),
    "D": (lambda n: n % 15 == 0, 2, ((8, 30),),
          lambda n: 12 * (n - 9)),
    "E": (lambda n: n % 15 == 0, 3, ((3, 3), (3, 3), (6, 30)),
          lambda n: (20 * n - (264 if n % 2 else 324)) // 3),
    "F": (lambda n: n % 15 == 0, 3, ((6, 30), (6, 30)),
          lambda n: 10 * n - (149 if n % 2 else 275)),
    "G": (lambda n: n % 30 == 0, 2, ((2, 2), (6, 30)),
          lambda n: 6 * (n - 8)),
    "H": (lambda n: n % 21 == 0, 2, ((8, 42),),
          lambda n: 6 * (n - 9)),
    "I": (lambda n: n % 21 == 0, 3, ((12, 42),),
          lambda n: 60),
    "J": (lambda n: n % 42 == 0, 3, ((2, 2), (10, 42)),
          lambda n: 60),
    "K": (lambda n: n % 33 == 0, 3, ((12, 66),),
          lambda n: 20),
    "L": (lambda n: n % 70 == 0, 3, ((2, 2), (10, 70)),
          lambda n: 20),
    "M": (lambda n: n % 105 == 0, 3, ((12, 210),),
          lambda n: 100),
}


# ------------------------------------------------------------------ machinery
def exact_collision_pairs(n):
    pairs = []
    for grp in collision_groups_mod_p(n):
        exact = {}
        for t in grp:
            exact.setdefault(level_vec(n, t).tobytes(), []).append(t)
        for ts in exact.values():
            for a in range(len(ts)):
                for b in range(a + 1, len(ts)):
                    pairs.append((ts[a], ts[b]))
    return pairs


def twelve_roots(n, t, s):
    m = 2 * n
    shared = set(t) & set(s)
    A = [k for k in t if k not in shared]
    B = [l for l in s if l not in shared]
    E = A + [n - l for l in B]
    return [x % m for e in E for x in (e, m - e)], len(A)


class _ModP:
    def __init__(self, n, seed):
        self.m = 2 * n
        self.p, self.g = find_prime_and_root(self.m, seed)

    def root(self, e):
        return pow(self.g, e % self.m, self.p)


def piece_decomposition(n, exps, mps):
    """Greedy (smallest-first) partition into minimal vanishing pieces; returns a
    sorted tuple of (size, ratio_order).  Vanishing is tested mod TWO independent
    primes = 1 (mod 2n): a nonzero sum can only read zero mod both by accident with
    probability ~1/p^2; the per-family count assertions against exact totals from
    the committed F129 census make a silent slip loud."""
    m = 2 * n
    k = len(exps)
    full = (1 << k) - 1
    sums = []
    for mp in mps:
        arr = np.zeros(1 << k, dtype=np.int64)
        for i in range(k):
            r = mp.root(exps[i])
            lo = 1 << i
            arr[lo:2 * lo] = (arr[:lo] + r) % mp.p
        sums.append(arr)
    vanish = (sums[0] == 0) & (sums[1] == 0)
    vanish[0] = False
    van_idx = np.nonzero(vanish)[0]
    van_set = set(int(x) for x in van_idx)

    def is_minimal(sub):
        t = (sub - 1) & sub
        while t:
            if t in van_set:
                return False
            t = (t - 1) & sub
        return True

    minimal = sorted((int(s) for s in van_idx if is_minimal(int(s))),
                     key=lambda s: bin(s).count("1"))

    def ratio_order(sub):
        es = [exps[i] for i in range(k) if (sub >> i) & 1]
        g = 0
        for e in es[1:]:
            g = int(np.gcd(g, (e - es[0]) % m))
        g = int(np.gcd(g, m))
        return m // g if g else 1

    def rec(remaining):
        if remaining == 0:
            return []
        for sub in minimal:
            if sub & ~remaining:
                continue
            rest = rec(remaining & ~sub)
            if rest is not None:
                return [(bin(sub).count("1"), ratio_order(sub))] + rest
        return None

    parts = rec(full)
    return tuple(sorted(parts)) if parts is not None else None


def classify(n):
    """Census of one firing n: {family label: count}."""
    mps = (_ModP(n, 0), _ModP(n, 1))
    sig_to_label = {(d, sig): lbl for lbl, (_, d, sig, _) in FAMILIES.items()}
    counts, total = {}, 0
    for t, s in exact_collision_pairs(n):
        exps, d = twelve_roots(n, t, s)
        parts = piece_decomposition(n, exps, mps)
        assert parts is not None, f"[I1] UNRESOLVED pair at n={n}: {t} ~ {s}"
        key = (d, parts)
        lbl = sig_to_label.get(key)
        assert lbl is not None, f"[I3] unknown family at n={n}: d={d} {parts} ({t}~{s})"
        counts[lbl] = counts.get(lbl, 0) + 1
        total += 1
    return counts, total


# ------------------------------------------------------------------ gates
def main():
    hi = 66 if "--fast" in sys.argv else 140
    ns = [n for n in range(9, hi + 1)
          if (n % 3 == 0 and n >= 9) or (n % 10 == 0 and n >= 20)]
    if "--deep" in sys.argv:
        ns += [150, 210]
    for n in ns:
        counts, total = classify(n)
        expected = {lbl: fam[3](n) for lbl, fam in FAMILIES.items()
                    if fam[0](n) and fam[3](n) > 0}
        assert counts == expected, \
            f"[I2/I4] n={n}: census {counts} != closed forms {expected}"
        assert sum(counts.values()) == total
        print(f"  n={n:3d}: {total:5d} pairs, families "
              + " ".join(f"{l}={c}" for l, c in sorted(counts.items())) + "  OK")
        sys.stdout.flush()
    print(f"[I1] every pair decomposes into inventory pieces at all {len(ns)} firing n  PASS")
    print(f"[I2] per-family counts == closed forms at all {len(ns)} firing n  PASS")
    print(f"[I3] the 13 families partition the census everywhere  PASS")
    print(f"[I4] doors: family present iff door divides n and formula nonzero  PASS")
    print("ALL GATES PASS: the F129 family inventory holds on the checked range.")


if __name__ == "__main__":
    main()
