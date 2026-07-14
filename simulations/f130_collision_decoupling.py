"""F130: the collision-decoupling law: equal level implies vanishing cross block.

The statement (the phase-2 physics re-entry of F129).  For the open XX chain with
n = N + 1, mode triples tau = {k1 < k2 < k3} of single-particle modes carry the
hopped Slater-seed cross block

    B(tau, sigma) = (K26 W_tau)^T (K26 W_sigma)
                  = [[U+, U+, 0], [U+, U+ + U-, U-], [0, U-, U-]],

with U+/U- the two half-Gram numbers (experiments/F89_SEED_EXISTENCE_REDUCTION.md,
"The two-number reduction").  The LEVEL of a triple is S(tau) = sum_i cos(k_i pi/n)
(half the energy sum, lam_k = 2 cos(k pi/n)).  F130 says: for ANY two distinct
triples with EQUAL levels S(tau) = S(sigma), the whole cross block vanishes,
B(tau, sigma) = 0.  Resonance (S = 0, the vanishing triples of the seed arc) is the
special case; the decoupling is a property of level COINCIDENCE, not of vanishing.

The proof is a four-cell assembly of committed results (see
docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md):
  disjoint  eps=+1 : Lemma 3 + Lemma 4 (level-free; B = 0 even at unequal levels)
  disjoint  eps=-1 : Lemma 3 + assembly (D) in free angles + F128 on {e1 = f1}
  overlap-1 eps=+1 : the same-two-magnon-energy lemma (equal pair energies of the
                     complements follow from equal triple levels + the shared mode)
  overlap-1 eps=-1 : Lemma 3 + the removable limit of the cross form on {e1 = f1}
  overlap-2        : vacuous (equal levels would force equal third modes)

Gates (all from below, float grade; the exactness lives in the cited proofs):
  [G1] assembly (D) OFF-resonance: at every disjoint equal-level clean pair with
       eps = -1 (n in NS), (U+ - U-)(n/2)^3 == cross_form to 1e-11, and both == 0.
  [G2] the decoupling: at EVERY distinct equal-level pair found numerically
       (disjoint AND overlap-1, clean or not), max(|U+|, |U-|) < 5e-13.
  [G3] non-vacuity controls: unequal-level pairs in each level-sensitive cell
       (disjoint eps=-1, overlap-1 both signs) have |U+-or-U-| > 1e-3; and (D)
       holds there with a value bounded away from 0.
  [G4] the level-free cell: disjoint eps=+1 pairs vanish at UNEQUAL levels too
       (Lemma 3 + 4), max(|U+|, |U-|) < 5e-13.
  [G5] overlap-2 vacuity: no distinct equal-level pair shares two modes.

Run: python simulations/f130_collision_decoupling.py          (~40 s, exit 0)
"""
import itertools
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cross_triple_orthogonality as ct

NS = (12, 13, 14, 15, 20, 21, 24, 30)   # firing n (3|n, 10|n) plus non-firing 13, 14
                                        # (trivial-family collisions exist everywhere)
TOL_ZERO = 5e-13      # |U+|, |U-| at equal levels (observed <= 7e-16)
TOL_D = 1e-11         # assembly (D) residual (values up to ~15)
TOL_LEVEL = 1e-12     # float equal-level grouping (exactness: F129's census)
GENERIC = 1e-3        # controls must exceed this


def all_triples(n):
    return list(itertools.combinations(range(1, n), 3))


def is_clean(n, t):
    return all(t[i] + t[j] != n for i in range(3) for j in range(i + 1, 3))


def level(n, t):
    return sum(math.cos(k * math.pi / n) for k in t)


def equal_level_pairs(n):
    """All distinct triple pairs with equal float level (tol 1e-12), any overlap."""
    trips = sorted(all_triples(n), key=lambda t: level(n, t))
    pairs = []
    for i, t in enumerate(trips):
        for s in trips[i + 1:]:
            if level(n, s) - level(n, t) > 1e-9:
                break
            if abs(level(n, s) - level(n, t)) < TOL_LEVEL:
                pairs.append((t, s))
    return pairs


_GCACHE = {}


def upm(U, n, tau, sig):
    N = n - 1
    for t in (tau, sig):
        if (n, t) not in _GCACHE:
            _GCACHE[(n, t)] = ct.Ggrid(U, t, N)
    return ct.Upm(_GCACHE[(n, tau)], _GCACHE[(n, sig)], N)


def gate1_assembly_off_resonance():
    checked = 0
    worst = 0.0
    for n in NS:
        U = ct.umat(n - 1)
        theta = math.pi / n
        for (t, s) in equal_level_pairs(n):
            if set(t) & set(s) or ct.eps(t, s) != -1:
                continue
            if abs(level(n, t)) < 1e-9:      # off-resonance only in this gate
                continue
            if not (is_clean(n, t) and is_clean(n, s)):
                continue
            up, um = upm(U, n, t, s)
            lhs = (up - um) * (n / 2.0) ** 3
            rhs = ct.cross_form(tuple(k * theta for k in t),
                                tuple(l * theta for l in s))
            worst = max(worst, abs(lhs - rhs), abs(lhs), abs(rhs))
            checked += 1
    assert checked >= 300, f"G1 too few off-resonance pairs: {checked}"
    assert worst < TOL_D, f"G1 assembly (D) off-resonance broken: {worst}"
    print(f"[G1] assembly (D) off-resonance at {checked} disjoint eps=-1 clean pairs, "
          f"worst |residual-or-value| = {worst:.2e}  PASS")


def gate2_decoupling():
    checked = 0
    worst = 0.0
    for n in NS:
        U = ct.umat(n - 1)
        for (t, s) in equal_level_pairs(n):
            if len(set(t) & set(s)) == 2:
                continue                      # G5's cell, vacuous
            up, um = upm(U, n, t, s)
            worst = max(worst, abs(up), abs(um))
            checked += 1
    assert checked >= 2000, f"G2 too few equal-level pairs: {checked}"
    assert worst < TOL_ZERO, f"G2 decoupling broken: {worst}"
    print(f"[G2] B(tau, sigma) = 0 at all {checked} equal-level pairs "
          f"(disjoint + overlap-1, clean or not), worst |U+-or-U-| = {worst:.2e}  PASS")


def gate3_controls():
    hits = {"disjoint_em1": 0, "ov1_ep1": 0, "ov1_em1": 0}
    for n in NS:
        U = ct.umat(n - 1)
        theta = math.pi / n
        trips = all_triples(n)
        for (t, s) in itertools.combinations(trips, 2):
            if all(v > 0 for v in hits.values()):
                break
            if abs(level(n, t) - level(n, s)) < 0.2:
                continue
            ov = len(set(t) & set(s))
            e = ct.eps(t, s)
            if ov == 0 and e == -1 and hits["disjoint_em1"] == 0:
                up, um = upm(U, n, t, s)
                if max(abs(up), abs(um)) > GENERIC:
                    lhs = (up - um) * (n / 2.0) ** 3
                    rhs = ct.cross_form(tuple(k * theta for k in t),
                                        tuple(l * theta for l in s))
                    assert abs(lhs - rhs) < TOL_D and abs(lhs) > GENERIC, \
                        f"G3 (D) control broken at n={n} {t}~{s}: {lhs} vs {rhs}"
                    hits["disjoint_em1"] += 1
            elif ov == 1:
                key = "ov1_ep1" if e == +1 else "ov1_em1"
                if hits[key] == 0:
                    up, um = upm(U, n, t, s)
                    if max(abs(up), abs(um)) > GENERIC:
                        hits[key] += 1
    assert all(v > 0 for v in hits.values()), f"G3 missing controls: {hits}"
    print("[G3] unequal-level controls generic (nonzero) in every level-sensitive "
          "cell; (D) holds there with nonzero value  PASS")


def gate4_level_free_cell():
    checked = 0
    worst = 0.0
    for n in NS[:3]:
        U = ct.umat(n - 1)
        trips = all_triples(n)
        for (t, s) in itertools.combinations(trips, 2):
            if set(t) & set(s) or ct.eps(t, s) != +1:
                continue
            if abs(level(n, t) - level(n, s)) < 0.2:
                continue
            up, um = upm(U, n, t, s)
            worst = max(worst, abs(up), abs(um))
            checked += 1
            if checked >= 200:
                break
        if checked >= 200:
            break
    assert checked >= 200 and worst < TOL_ZERO, \
        f"G4 level-free cell broken: {checked} pairs, worst {worst}"
    print(f"[G4] disjoint eps=+1 pairs vanish at UNEQUAL levels too "
          f"({checked} pairs, worst {worst:.2e}): the level condition is "
          f"load-bearing only in the other three cells  PASS")


def gate5_overlap2_vacuous():
    for n in NS:
        for (t, s) in equal_level_pairs(n):
            assert len(set(t) & set(s)) != 2, f"G5 overlap-2 pair at n={n}: {t}~{s}"
    print("[G5] no distinct equal-level pair shares two modes (vacuous cell)  PASS")


def main():
    gate1_assembly_off_resonance()
    gate2_decoupling()
    gate3_controls()
    gate4_level_free_cell()
    gate5_overlap2_vacuous()
    print("ALL GATES PASS: equal level implies vanishing cross block "
          "(the collision-decoupling law, F130)")


if __name__ == "__main__":
    main()
