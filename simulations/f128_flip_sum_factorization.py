"""F128: the (e1-f1)^2 factorization of the cross form, and the sharper locus.

THEOREM (found 2026-07-14, the handover's lead 1: the 32-term flip-sum's own proof):

    F(x) = -(e1 - f1)^2 * O[ cos(s) * cot(s) * V_a * V_b / P ]        ... (F128)

as an identity of rational functions of the six angles x = (a1,a2,a3,b1,b2,b3), where
F = cross_form (the committed F127 object), e1 = Sum cos a_i, f1 = Sum cos b_j,
s = (Sum a + Sum b)/2, V_a = Prod_{i<j} sin((a_i-a_j)/2), V_b likewise,
P = Prod_{ij} sin((a_i+b_j)/2), and O = the projection onto the part odd in each of
the six angles separately (O[g](x) = 2^-6 Sum_{eps in {+-1}^6} (Prod eps) g(eps o x)).

COROLLARY (the sharper locus): F == 0 on {e1 = f1}, ONE constraint; F127's variety
V = {e1 = 0, f1 = 0} is the codimension-two special case.  Until today every sharper
{e1 = f1} statement in the repo was about the core function T (and needed the sheet);
this one is about F itself and needs nothing else.

PROOF (three committed inputs + one new lemma):
  (1) the SS3A representation (committed, f127_global_representation.py):
        F = -4 * O[ T * cot s ].
  (2) the SS3 closed form (committed, f127_closed_form.py, gate G1):
        T*P = (1/8) [ 2 cos s ((e1-f1)^2 - 2 sin^2 s) + sin s (Sum sin 2x_u) ] V_a V_b,
      which splits T*cot s = h + (e1-f1)^2/4 * cos s cot s V_a V_b / P  with
        h = (1/8) cos s * B * V_a V_b / P,   B = Sum_u sin 2x_u - 2 sin 2s.
      (Term splitting only; no new content beyond G1.)
  (3) e1, f1 are flip-invariant (cos is even), so O passes over (e1-f1)^2.
  (4) NEW LEMMA (the flip lemma, proved exactly over Z below, twice):
        O[ cos s * B * V_a V_b * Ptilde ] == 0   IDENTICALLY,
      Ptilde = Prod_{ij} sin((a_i-b_j)/2).  Since P*Ptilde = 2^-9 Prod_{ij}(cos b_j
      - cos a_i) is flip-invariant, this gives O[h] == 0, and (1)-(3) give F128.

THE FLIP LEMMA'S TWO EXACT PROOFS (both in this gate):
  (A) brute force over Z: the integer trig polynomial cos s * B * V_a V_b * Ptilde
      (8640 monomials, half-angle monomial coordinates) is annihilated by the
      (Z/2)^6 signed character sum -- the projection returns the zero polynomial.
  (B) the trade already knew (Weyl / Murnaghan-Nakayama): V_a V_b Ptilde is the full
      six-angle half-angle Vandermonde Delta = Prod_{u<v} sin((x_u-x_v)/2); the Weyl
      denominator formula folds cos s * Delta into TWO alternants,
        2 (2i)^15 cos s Delta = a_{M1} + a_{M2},  a_M = det[z_u^{m_k}],
        M1 = (3,2,1,0,-1,-2),  M2 = -M1 reversed = (2,1,0,-1,-2,-3);
      the B factor acts by Murnaghan-Nakayama shifts (+-2 on one z-exponent from
      p2(z) - p2(1/z), +-1 on all from the sheet term sin 2s), producing 28
      alternants; every
      shifted exponent set contains a zero, a repeat, or a +- pair; and the odd
      projection of an alternant is det[z_u^{m_k} - z_u^{-m_k}] (the odd Weyl
      numerator), which then has a zero column or equal-up-to-sign columns: zero.

GATES (all exact over Z unless marked numeric):
  [G1] the flip lemma, proof (A): O[integer dict] == the zero polynomial, with a
       self-test that the projector does not annihilate a generic monomial.
  [G2] proof (B) step 1: 2(2i)^15 cos s Delta == a_{M1} + a_{M2} as integer dicts.
  [G3] proof (B) step 2: Btilde * (a_{M1}+a_{M2}) == the signed sum of the 28
       shifted alternants, as integer dicts (the Murnaghan-Nakayama bookkeeping).
  [G4] proof (B) step 3: every shifted exponent set dies structurally (0/repeat/
       +-pair: integer logic), and its exact O-projection is the zero dict.
  [G5] numeric pin: F == -(e1-f1)^2 * O[cos s cot s V_a V_b / P] at generic points
       (brute 64-flip evaluation against the committed cross_form).
  [G6] numeric discriminator: |F| on {e1 = f1 = c}, c != 0, at machine zero while
       generic |F| has median ~10 (single samples span many orders; F has other
       zero sets); plus the flip-sum S = -8F re-pin on V.
  [G7] restriction hygiene: e1 - f1 is not identically zero on any coset of the
       30 pair forms x_u +- x_v == 0 and the 32 sheet forms L.x == 0 -- a
       conservative SUPERSET of the right side's actual polar support (V_a V_b
       sit in the numerator; only P and cot s contribute poles: the 32 sheet
       forms and the 18 cross forms a_i +- b_j); hence no component of
       {e1 = f1} lies inside the polar locus and the factorization restricts to
       every component.  Witness points are floating point with a wide margin
       (> 1.4), certifying the exact fact that the restriction is nonconstant.

Run: python simulations/f128_flip_sum_factorization.py    (~0.3 s, exit 0 iff all pass)
"""
import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "simulations")
from cross_triple_orthogonality import cross_form

RNG = np.random.default_rng(20260714)


# ---------------------------------------------------------------------------
# exact sparse Laurent arithmetic over Z; a polynomial is dict[6-tuple -> int],
# exponents in HALF-ANGLE units (t_u = e^{i x_u / 2}, key k means t_u^k).
# ---------------------------------------------------------------------------
def pmul(p, q):
    r = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2],
                 e1[3] + e2[3], e1[4] + e2[4], e1[5] + e2[5])
            v = r.get(e, 0) + c1 * c2
            if v:
                r[e] = v
            elif e in r:
                del r[e]
    return r


def padd(p, q, sign=1):
    r = dict(p)
    for e, c in q.items():
        v = r.get(e, 0) + sign * c
        if v:
            r[e] = v
        elif e in r:
            del r[e]
    return r


def two_term(u, v):
    """2i * sin((x_u - x_v)/2)  =  t_u/t_v - t_v/t_u."""
    e1, e2 = [0] * 6, [0] * 6
    e1[u], e1[v] = 1, -1
    e2[u], e2[v] = -1, 1
    return {tuple(e1): 1, tuple(e2): -1}


def delta_dict():
    """(2i)^15 * Delta = Prod_{u<v} (2i sin((x_u-x_v)/2)), ordering a1..a3,b1..b3."""
    p = {(0,) * 6: 1}
    for u in range(6):
        for v in range(u + 1, 6):
            p = pmul(p, two_term(u, v))
    return p


COS_S = {(1,) * 6: 1, (-1,) * 6: 1}          # 2 cos s


def b_dict():
    """2i * B = Sum_u (t_u^4 - t_u^-4) - 2 (T^2 - T^-2),  T = t1..t6."""
    B = {}
    for u in range(6):
        e = [0] * 6
        e[u] = 4
        B[tuple(e)] = 1
        e[u] = -4
        B[tuple(e)] = -1
    B[(2,) * 6] = -2
    B[(-2,) * 6] = 2
    return B


def odd_project(p):
    """Prod_u (1 - flip_u) = 64 * O[.] on dicts (flip_u negates exponent u)."""
    for u in range(6):
        r = dict(p)
        for e, c in p.items():
            ef = list(e)
            ef[u] = -ef[u]
            ef = tuple(ef)
            v = r.get(ef, 0) - c
            if v:
                r[ef] = v
            elif ef in r:
                del r[ef]
        p = r
    return p


def perm_sign(perm):
    s, seen = 1, [False] * len(perm)
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def alternant(M):
    """a_M = det[z_u^{m_k}] as a dict in half-angle units (z = t^2)."""
    p = {}
    for perm in itertools.permutations(range(6)):
        sgn = perm_sign(perm)
        e = tuple(2 * M[perm[u]] for u in range(6))
        v = p.get(e, 0) + sgn
        if v:
            p[e] = v
        elif e in p:
            del p[e]
    return p


M1 = (3, 2, 1, 0, -1, -2)
M2 = (2, 1, 0, -1, -2, -3)


def shifted_sets():
    """the 28 signed Murnaghan-Nakayama images of M1, M2 under 2i*B:
    p2(z) shifts one z-exponent by +2, -p2(1/z) by -2, and the sheet term
    -2(Z - 1/Z) with Z = z1..z6 shifts ALL exponents by +-1."""
    out = []
    for M in (M1, M2):
        for k in range(6):
            out.append((+1, tuple(m + (2 if j == k else 0) for j, m in enumerate(M))))
            out.append((-1, tuple(m - (2 if j == k else 0) for j, m in enumerate(M))))
        out.append((-2, tuple(m + 1 for m in M)))
        out.append((+2, tuple(m - 1 for m in M)))
    return out


# ---------------------------------------------------------------------------
def gate1_flip_lemma_bruteforce():
    tst = odd_project({(1, 2, 3, 4, 5, 6): 1})
    assert len(tst) == 64, "projector self-test failed"
    N = pmul(pmul(delta_dict(), COS_S), b_dict())
    assert len(N) == 8640, f"unexpected numerator size {len(N)}"
    G = odd_project(N)
    assert not G, f"flip lemma FAILED: {len(G)} monomials survive"
    print("[G1] flip lemma, proof (A): O[8640-monomial integer polynomial] == 0  PASS")
    return N


def gate2_weyl_denominator():
    lhs = pmul(delta_dict(), COS_S)
    rhs = padd(alternant(M1), alternant(M2))
    assert lhs == rhs, "Weyl-denominator folding failed"
    print("[G2] 2(2i)^15 cos s Delta == a_M1 + a_M2 exactly over Z  PASS")
    return rhs


def gate3_murnaghan_nakayama(cosdelta):
    lhs = pmul(b_dict(), cosdelta)
    rhs = {}
    for sgn, M in shifted_sets():
        rhs = padd(rhs, alternant(M), sgn)
    assert lhs == rhs, "Murnaghan-Nakayama bookkeeping failed"
    print("[G3] 2iB * (a_M1 + a_M2) == the 28 signed shifted alternants exactly  PASS")


def gate4_shifted_alternants_die():
    def dies(M):
        return (0 in M) or (len(set(M)) < 6) or any(-m in M for m in M if m != 0)
    sets = shifted_sets()
    bad = [M for _, M in sets if not dies(M)]
    assert not bad, f"structurally surviving sets: {bad}"
    for _, M in sets:
        assert not odd_project(alternant(M)), f"O[a_{M}] != 0"
    print(f"[G4] all {len(sets)} shifted sets die (0/repeat/+-pair) and O[a_M] == 0  PASS")


def W_brute(a, b):
    """O[cos s cot s V_a V_b / P] by literal 64-term evaluation."""
    x = list(a) + list(b)
    tot = 0.0
    for eps in itertools.product((1, -1), repeat=6):
        xe = [e * v for e, v in zip(eps, x)]
        ae, be = xe[:3], xe[3:]
        se = sum(xe) / 2
        Va = math.prod(math.sin((ae[i] - ae[j]) / 2)
                       for i in range(3) for j in range(i + 1, 3))
        Vb = math.prod(math.sin((be[i] - be[j]) / 2)
                       for i in range(3) for j in range(i + 1, 3))
        P = math.prod(math.sin((ai + bj) / 2) for ai in ae for bj in be)
        tot += math.prod(eps) * math.cos(se) / math.tan(se) * Va * Vb / P
    return tot / 64.0


def rand_angles(n=3):
    return RNG.uniform(0.3, 2 * np.pi - 0.3, n)


def rand_triple(target):
    while True:
        a1, a2 = RNG.uniform(0.2, 2 * np.pi - 0.2, 2)
        c = target - math.cos(a1) - math.cos(a2)
        if abs(c) <= 0.98:
            a3 = math.acos(c) * RNG.choice([1.0, -1.0]) % (2 * np.pi)
            return [a1, a2, a3]


def gate5_factorization_pin():
    worst = 0.0
    for _ in range(10):
        a, b = rand_angles(), rand_angles()
        e1 = sum(math.cos(v) for v in a)
        f1 = sum(math.cos(v) for v in b)
        F = cross_form(tuple(a), tuple(b))
        rhs = -(e1 - f1) ** 2 * W_brute(a, b)
        worst = max(worst, abs(F - rhs) / max(1.0, abs(F)))
    print(f"[G5] numeric pin F == -(e1-f1)^2 W, 10 generic points: rel dev {worst:.2e}  PASS")
    assert worst < 1e-8


def alpha(v):
    out = []
    for i in range(3):
        j, l = [t for t in range(3) if t != i]
        out.append(math.sin(v[l]) - math.sin(v[j]))
    return out


def T_def(a, b):
    aa, bb = alpha(a), alpha(b)
    return sum((-1) ** (i + j) * aa[i] * bb[j] / math.tan((a[i] + b[j]) / 2)
               for i in range(3) for j in range(3))


def gate6_sharper_locus():
    gen = [abs(cross_form(tuple(rand_angles()), tuple(rand_angles()))) for _ in range(10)]
    worst_on = 0.0
    for c in (0.4, -0.7, 1.1, 1.9):
        for _ in range(10):
            a, b = rand_triple(c), rand_triple(c)
            worst_on = max(worst_on, abs(cross_form(tuple(a), tuple(b))))
    # the flip-sum re-pin on V (the F127 special case, via the SS3A representation)
    worstS = 0.0
    for _ in range(5):
        a, b = rand_triple(0.0), rand_triple(0.0)
        x = list(a) + list(b)
        S = 0.0
        for L in [(1,) + t for t in itertools.product((1, -1), repeat=5)]:
            xl = [Li * xi for Li, xi in zip(L, x)]
            S += math.prod(L) * T_def(xl[:3], xl[3:]) / math.tan(sum(xl) / 2)
        worstS = max(worstS, abs(S))
    med = sorted(gen)[len(gen) // 2]
    print(f"[G6] sharper locus: generic |F| median {med:.1e} (range [{min(gen):.1e}, "
          f"{max(gen):.1e}]), max |F| on {{e1=f1=c}} = {worst_on:.1e}; "
          f"max |flip-sum| on V = {worstS:.1e}  PASS")
    assert worst_on < 1e-8 and worstS < 1e-8 and med > 1e-2


def gate7_coset_hygiene():
    """e1 - f1 not identically zero on any denominator-factor coset.  The right
    side's denominator factors (over all 64 flips) are sines of the pair forms
    x_u +- x_v and of the 32 sheet forms L.x; each vanishes on cosets of codim-1
    subtori.  On each coset we exhibit a point with e1 != f1; a nonconstant
    restriction has a dim<=4 zero set inside the dim-5 coset, so no component of
    the hypersurface {e1 = f1} can lie inside it."""
    def ef(x):
        return sum(math.cos(v) for v in x[:3]) - sum(math.cos(v) for v in x[3:])
    worst = math.inf
    # pair forms x_u - x_v == 0 and x_u + x_v == 0 (mod 2pi): parametrize a point
    for u in range(6):
        for v in range(u + 1, 6):
            for sgn in (+1, -1):
                x = [0.3 + 0.4 * k for k in range(6)]        # generic base point
                x[v] = sgn * x[u] % (2 * np.pi)              # sit on the coset
                worst = min(worst, abs(ef(x)))
    # sheet forms L.x == 0 (mod 2pi)
    for L in [(1,) + t for t in itertools.product((1, -1), repeat=5)]:
        x = [0.3 + 0.4 * k for k in range(5)] + [0.0]
        x[5] = (-sum(Li * xi for Li, xi in zip(L, x)) / L[5]) % (2 * np.pi)
        worst = min(worst, abs(ef(x)))
    print(f"[G7] coset hygiene: min |e1 - f1| over the 62 denominator cosets' "
          f"witness points = {worst:.3f}  PASS")
    assert worst > 1e-2


if __name__ == "__main__":
    gate1_flip_lemma_bruteforce()
    cosdelta = gate2_weyl_denominator()
    gate3_murnaghan_nakayama(cosdelta)
    gate4_shifted_alternants_die()
    gate5_factorization_pin()
    gate6_sharper_locus()
    gate7_coset_hygiene()
    print("ALL GATES PASS: F = -(e1-f1)^2 O[cos s cot s Va Vb / P]; F == 0 on {e1 = f1}.")
