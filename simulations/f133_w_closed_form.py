"""F133: the closed form of the F128 cofactor W, via the symplectic character sum.

THEOREM (found 2026-07-17):

    W = -2^9 * (Prod_u sin x_u) * V_c(a) * V_c(b) * K(c) / SP(x)        ... (F133)

as an identity of rational functions of the six angles x = (a1,a2,a3,b1,b2,b3),
where c_u = cos x_u, V_c(a) = Prod_{i<j}(cos a_i - cos a_j) (V_c(b) likewise),
SP = Prod over the 32 canonical sheets L (L in {+-1}^6, L_1 = +1) of sin(L.x/2),
and K is the SYMMETRIC polynomial in (c_1..c_6)

    K = 2^-30 * Sum_lambda n_lambda * chi^{C6}_lambda,

a sum of 143 Weyl characters of Sp(12) (type C6), n_lambda in Z, |n_lambda| <= 8,
committed in simulations/results/f133_w_closed_form/chiC_coeffs.txt (the same K in
the monomial-symmetric basis, 190 dyadic terms, in m_coeffs.txt).  With F128
(F = -(e1-f1)^2 * W) this makes the whole cross form explicit:

    F = +2^9 * (e1-f1)^2 * (Prod sin x_u) * V_c(a) V_c(b) * K(c) / SP(x).

REDUCTION CHAIN (each step an identity; gates in brackets):
  [R1] the P*Ptilde escape (F128's flip-invariant combination, gate G1):
       pairwise sin((a+b)/2) sin((a-b)/2) = (cos b - cos a)/2, so
       W = 2^9 * O[cos s cot s * Dhat] / Prod_{ij}(cos b_j - cos a_i),
       Dhat = the full six-angle half-angle Vandermonde Prod_{u<v} sin((x_u-x_v)/2).
  [R2] cos s cot s = csc s - sin s, and O[sin s * Dhat] == 0 EXACTLY over Z
       (gate G2, the sin-s lemma), so the core is C = O[Dhat / sin s]
       (the trigonometric core; the Lie type is always written C6/Sp(12)).
  [R3] SP is flip-invariant (gate G3, recanonicalization bookkeeping: a sign flip
       permutes the 64 sheet forms; the 32 recanonicalization signs multiply to +1),
       and SP / sin s = Prod over the 31 canonical sheets != 1^6 of sin(s_L).  Hence
       C * SP = O[X],  X = Dhat * Prod_31 sin(s_L),  a trig POLYNOMIAL.
  [R4] X is S6-alternating (gate G4: Dhat alternates as an exact dict; the 31-sheet
       product is S6-invariant, transposition recanonicalization signs even), so
       64*O[X] is W(C6)-antisymmetric and expands in symplectic alternants
       A_mu = det[z_u^{mu_j} - z_u^{-mu_j}] (z = e^{ix}) over strictly dominant mu.
       The coefficients are READ OFF, no fit and no cap:
           n_raw(lambda) = Sum_{eps in {+-1}^6} sgn(eps) [X]_{eps o 2(lambda+rho)},
       rho = (6,5,4,3,2,1), t-exponent units (z = t^2).  Gate G5 recomputes every
       n_raw by exact integer meet-in-the-middle arithmetic and matches the
       committed table (n_raw = 2 n_lambda) on the full even-degree window
       |lambda| <= 18; the flag --full sweeps ALL 18564 dominant lambda with
       parts <= 12 (the window forced by X's exponent support: per-variable
       t-exponents of X lie in [-36, 36], so a contributing alternant needs
       2(lambda_1 + 6) <= 36, i.e. lambda_1 <= 12) and verifies that EXACTLY the
       143 committed terms survive.
  [R5] the C6 Weyl denominator (gate G6, exact dicts):
       A_rho = det[z_u^{7-j} - z_u^{-(7-j)}] = (2i)^6 * 2^15 * Prod sin x_u * Vand_c,
       Vand_c = Prod_{u<v}(c_u - c_v) = V_c(a) V_c(b) * Prod_{ij}(c_{a_i} - c_{b_j}).
       Dividing [R4] by it, with X_dict = (2i)^46 X and 64*O = the dict projector:
       K = C*SP / (Prod sin * Vand_c) = Sum n_lambda chi_lambda * (2i)^6 2^15 /
       (2^5 (2i)^46) = 2^-30 * Sum n_lambda chi^{C6}_lambda.  The committed
       common denominator 2^30 is thereby DERIVED, not fitted.
       Assembling [R1]-[R5] and cancelling the cross factor Prod(c_a - c_b) =
       (-1)^9 Prod(cos b - cos a) gives F133 (the (-1)^9 is the minus sign).

EVIDENCE GRADE: every link above is exact (integer dicts / integer logic) except
the numeric end-to-end pins G8, which corroborate the assembled rational-function
identity against the committed cross_form at generic points.  The 143 coefficients
are derived twice independently (GF(p)+CRT fit, scout grade; and the alternant
read-off of G5); gate G7 pins the two committed tables (chi^C and m-basis) against
each other exactly at random GF(p) points.

GATES:
  [G1] pair identity (2i sin((a+b)/2))(2i sin((a-b)/2)) == 2 cos a - 2 cos b,
       exact dicts, all 9 (i,j) pairs.
  [G2] the sin-s lemma: O[(2i sin s) * (2i)^15 Dhat] == the zero dict, with a
       projector self-test.
  [G3] SP flip-invariance: for each of the 6 generator flips, the canonicalized
       image of the 32-sheet multiset is the same multiset and the sign count is
       even (integer logic); numeric corroboration at a random point.
  [G4] X alternates: Dhat dict is exactly antisymmetric under all 15
       transpositions; the 31-sheet set is stable with even sign count under all
       15 transpositions (integer logic); numeric corroboration for the product.
  [G5] the alternant read-off: meet-in-the-middle X = P1 * P2 (P1 = Dhat * 15
       sheets, P2 = 16 sheets), n_raw(lambda) extracted for all 557 window
       lambda, must equal 2 * (committed n_lambda), zero elsewhere in the window.
       --full: the complete 18564-candidate dominance sweep (support proof).
  [G6] the Weyl denominator identity, exact dicts (46080-term determinant
       expansion vs the product form).
  [G7] chi^C table vs m-basis table: K evaluated both ways at random points
       mod two 30-bit primes, exact equality.
  [G8] end-to-end numeric: W_brute (the literal 64-flip evaluation of the F128
       cofactor) vs the closed form, and cross_form vs -(e1-f1)^2 * W_closed,
       10 generic points each, rel dev < 1e-8.

Run:  python simulations/f133_w_closed_form.py           (~25 s, exit 0 iff all pass)
      python simulations/f133_w_closed_form.py --full    (adds the ~9 min support sweep)
"""
import itertools
import math
import random
import sys
import time
from fractions import Fraction
from functools import lru_cache

import numpy as np

sys.path.insert(0, "simulations")
from f128_flip_sum_factorization import (W_brute, pmul, padd, delta_dict,
                                         odd_project)
from cross_triple_orthogonality import cross_form

T0 = time.time()
RNG = np.random.default_rng(20260717)

CANON = [(1,) + L for L in itertools.product((1, -1), repeat=5)]   # 32 sheets
SHEETS31 = [L for L in CANON if L != (1, 1, 1, 1, 1, 1)]
RHO = (6, 5, 4, 3, 2, 1)
FLIPS = list(itertools.product((1, -1), repeat=6))

DATA = "simulations/results/f133_w_closed_form"


def load_chiC():
    tab = {}
    for line in open(f"{DATA}/chiC_coeffs.txt"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        lam_s, n_s = line.split(";")
        tab[tuple(int(v) for v in lam_s.split(",") if v)] = int(n_s)
    return tab


def load_m():
    tab = {}
    for line in open(f"{DATA}/m_coeffs.txt"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        lam_s, q_s = line.split(";")
        tab[tuple(int(v) for v in lam_s.split(",") if v)] = Fraction(q_s)
    return tab


CHI_TAB = load_chiC()
M_TAB = load_m()


# ---------------------------------------------------------------------------
def gate1_pair_identity():
    """(2i sin((a_i+b_j)/2)) (2i sin((a_i-b_j)/2)) == 2 cos a_i - 2 cos b_j."""
    for i in range(3):
        for j in range(3):
            u, v = i, 3 + j
            plus, minus = [0] * 6, [0] * 6
            plus[u], plus[v] = 1, 1
            minus[u], minus[v] = 1, -1
            f1 = {tuple(plus): 1, tuple(-e for e in plus): -1}
            f2 = {tuple(minus): 1, tuple(-e for e in minus): -1}
            lhs = pmul(f1, f2)
            ea, eb = [0] * 6, [0] * 6
            ea[u] = 2
            eb[v] = 2
            rhs = {tuple(ea): 1, tuple(-e for e in ea): 1,
                   tuple(eb): -1, tuple(-e for e in eb): -1}
            assert lhs == rhs, f"pair identity failed at (i,j)=({i},{j})"
    print(f"[G1] P*Ptilde pair identity, all 9 pairs, exact over Z  PASS"
          f"   [{time.time()-T0:.0f}s]")


def gate2_sin_s_lemma():
    tst = odd_project({(1, 2, 3, 4, 5, 6): 1})
    assert len(tst) == 64, "projector self-test failed"
    SINS = {(1,) * 6: 1, (-1,) * 6: -1}            # 2i sin s
    surv = odd_project(pmul(SINS, delta_dict()))
    assert not surv, f"sin-s lemma FAILED: {len(surv)} monomials survive"
    print(f"[G2] the sin-s lemma: O[sin s * Dhat] == 0 exactly over Z  PASS"
          f"   [{time.time()-T0:.0f}s]")


def canonicalize(L):
    """(sheet form, sign): sin((L.x)/2) = sign * sin((canonical L').x/2)."""
    return (L, 1) if L[0] == 1 else (tuple(-e for e in L), -1)


def gate3_sp_flip_invariance():
    for u in range(6):
        eps = [1] * 6
        eps[u] = -1
        imgs, sgn = [], 1
        for L in CANON:
            Lp, s = canonicalize(tuple(e * l for e, l in zip(eps, L)))
            imgs.append(Lp)
            sgn *= s
        assert sorted(imgs) == sorted(CANON), f"flip {u}: sheet set not stable"
        assert sgn == 1, f"flip {u}: odd recanonicalization sign"
    x = RNG.uniform(0.3, 2 * np.pi - 0.3, 6)
    sp = math.prod(math.sin(sum(l * v for l, v in zip(L, x)) / 2) for L in CANON)
    worst = max(abs(sp - math.prod(
        math.sin(sum(l * e * v for l, e, v in zip(L, eps, x)) / 2) for L in CANON))
        for eps in FLIPS)
    assert worst < 1e-12 * max(1.0, abs(sp))
    print(f"[G3] SP flip-invariant: 6 generator flips exact bookkeeping "
          f"(+ numeric {worst:.1e})  PASS   [{time.time()-T0:.0f}s]")


def gate4_x_alternates():
    DD = delta_dict()
    for u in range(6):
        for v in range(u + 1, 6):
            perm = list(range(6))
            perm[u], perm[v] = v, u
            swapped = {}
            for e, c in DD.items():
                swapped[tuple(e[perm[k]] for k in range(6))] = c
            assert swapped == {e: -c for e, c in DD.items()}, \
                f"Dhat not antisymmetric under ({u},{v})"
            imgs, sgn = [], 1
            for L in SHEETS31:
                Lp, s = canonicalize(tuple(L[perm[k]] for k in range(6)))
                imgs.append(Lp)
                sgn *= s
            assert sorted(imgs) == sorted(SHEETS31), \
                f"31-sheet set not stable under ({u},{v})"
            assert sgn == 1, f"odd sheet sign count under ({u},{v})"
    x = list(RNG.uniform(0.3, 2 * np.pi - 0.3, 6))
    q = math.prod(math.sin(sum(l * v for l, v in zip(L, x)) / 2) for L in SHEETS31)
    xs = [x[1], x[0]] + x[2:]
    qs = math.prod(math.sin(sum(l * v for l, v in zip(L, xs)) / 2) for L in SHEETS31)
    assert abs(q - qs) < 1e-12 * max(1.0, abs(q))
    print(f"[G4] X alternates: Dhat antisymmetric (exact dicts), 31-sheet product "
          f"S6-invariant (bookkeeping + numeric)  PASS   [{time.time()-T0:.0f}s]")


# ---- packed-int meet-in-the-middle (the alternant read-off) -----------------
BIAS = 32768
SHIFTS = [16 * k for k in range(6)]
PACKB = sum(BIAS << s for s in SHIFTS)


def pack(vec):
    return sum((v + BIAS) << s for v, s in zip(vec, SHIFTS))


def mul_factor(poly, vec):
    """poly * (t^vec - t^-vec), poly keyed by packed biased exponents."""
    up = pack(vec) - PACKB
    dn = pack([-v for v in vec]) - PACKB
    out = {}
    for k, c in poly.items():
        for shift, sgn in ((up, 1), (dn, -1)):
            kk = k + shift
            v = out.get(kk, 0) + sgn * c
            if v:
                out[kk] = v
            elif kk in out:
                del out[kk]
    return out


def build_halves():
    P1 = {PACKB: 1}
    for u in range(6):
        for v in range(u + 1, 6):
            e = [0] * 6
            e[u], e[v] = 1, -1
            P1 = mul_factor(P1, e)
    for L in SHEETS31[:15]:
        P1 = mul_factor(P1, list(L))
    P2 = {PACKB: 1}
    for L in SHEETS31[15:]:
        P2 = mul_factor(P2, list(L))
    return P1, list(P2.items())


def n_raw(lam, P1, P2items):
    lam6 = list(lam) + [0] * (6 - len(lam))
    mu = [2 * (lam6[j] + RHO[j]) for j in range(6)]
    g = P1.get
    tot = 0
    for eps in FLIPS:
        sgn = 1
        for e in eps:
            sgn *= e
        tp = pack([e * m for e, m in zip(eps, mu)]) + PACKB
        for k2, c2 in P2items:
            v = g(tp - k2)
            if v:
                tot += sgn * c2 * v
    return tot


def window_lams(full):
    lams = [()]

    def gen(prefix, mx, rem):
        for q in range(min(mx, rem), 0, -1):
            if len(prefix) < 6:
                lam = tuple(prefix + [q])
                if full or sum(lam) % 2 == 0:
                    lams.append(lam)
                gen(prefix + [q], q, rem - q)

    gen([], 12 if full else 18, 100 if full else 18)
    return sorted(set(lams), key=lambda t: (sum(t), t))


def gate5_alternant_readoff(full):
    P1, P2items = build_halves()
    span = max(max(abs((k >> s) % (1 << 16) - BIAS) for s in SHIFTS)
               for k in P1) + 16
    assert span <= 36, f"X exponent span {span} exceeds the stated window"
    lams = window_lams(full)
    label = "full 18564-candidate dominance sweep" if full \
        else "557-candidate even-degree window"
    bad = 0
    for cnt, lam in enumerate(lams):
        want = 2 * CHI_TAB.get(lam, 0)
        got = n_raw(lam, P1, P2items)
        if got != want:
            bad += 1
            if bad <= 10:
                print(f"   MISMATCH {lam}: read-off {got}, table x2 {want}")
        if cnt and cnt % 2000 == 0:
            print(f"   {cnt}/{len(lams)}   [{time.time()-T0:.0f}s]")
    assert bad == 0, f"{bad} coefficient mismatches"
    nz = sum(1 for lam in lams if CHI_TAB.get(lam, 0))
    assert nz == len(CHI_TAB) == 143, "committed table size mismatch"
    print(f"[G5] alternant read-off ({label}): all {len(lams)} n_raw == "
          f"2 n_lambda, exactly the 143 committed terms nonzero  PASS"
          f"   [{time.time()-T0:.0f}s]")


def gate6_weyl_denominator():
    """A_rho = det[z^{7-j} - z^{-(7-j)}] == Prod_u (t_u^2 - t_u^-2) *
    Prod_{u<v} (t_u^2 + t_u^-2 - t_v^2 - t_v^-2), exact dicts (t-units)."""
    det = {}
    for perm in itertools.permutations(range(6)):
        s, seen = 1, [False] * 6
        for i in range(6):
            if seen[i]:
                continue
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                ln += 1
            if ln % 2 == 0:
                s = -s
        for signs in itertools.product((1, -1), repeat=6):
            e = tuple(signs[u] * 2 * (7 - (perm[u] + 1)) for u in range(6))
            v = det.get(e, 0) + s * math.prod(signs)
            if v:
                det[e] = v
            elif e in det:
                del det[e]
    prod = {(0,) * 6: 1}
    for u in range(6):
        e = [0] * 6
        e[u] = 2
        prod = pmul(prod, {tuple(e): 1, tuple(-x for x in e): -1})
    for u in range(6):
        for v in range(u + 1, 6):
            eu, ev = [0] * 6, [0] * 6
            eu[u], ev[v] = 2, 2
            prod = pmul(prod, {tuple(eu): 1, tuple(-x for x in eu): 1,
                               tuple(ev): -1, tuple(-x for x in ev): -1})
    assert det == prod, "Weyl denominator identity failed"
    print(f"[G6] C6 Weyl denominator == product form, exact over Z "
          f"({len(det)} monomials)  PASS   [{time.time()-T0:.0f}s]")


# ---- exact mod-p cross-link of the two committed tables ---------------------
def det_mod(mat, p):
    n = len(mat)
    mat = [r[:] for r in mat]
    det = 1
    for col in range(n):
        piv = next((k for k in range(col, n) if mat[k][col] % p), None)
        if piv is None:
            return 0
        if piv != col:
            mat[col], mat[piv] = mat[piv], mat[col]
            det = -det
        det = det * mat[col][col] % p
        inv = pow(mat[col][col], p - 2, p)
        for k in range(col + 1, n):
            f = mat[k][col] * inv % p
            if f:
                mat[k] = [(v - f * w) % p for v, w in zip(mat[k], mat[col])]
    return det % p


@lru_cache(maxsize=None)
def perms_of(lam6):
    return list(set(itertools.permutations(lam6)))


def gate7_table_crosslink():
    inv2_30 = None
    for p in (2013265921, 1811939329):
        rnd = random.Random(20260717 + p)
        inv2 = pow(2, p - 2, p)
        inv2_30 = pow(pow(2, 30, p), p - 2, p)
        checked = 0
        while checked < 10:
            ts = [rnd.randrange(2, p) for _ in range(6)]
            tinv = [pow(t, p - 2, p) for t in ts]
            c = [(t * t + v * v) % p * inv2 % p for t, v in zip(ts, tinv)]
            if len(set(c)) < 6:
                continue

            def chi(lam6):
                mat = []
                for u in range(6):
                    row = []
                    for j in range(6):
                        e = 2 * (lam6[j] + 7 - (j + 1))
                        row.append((pow(ts[u], e, p) - pow(tinv[u], e, p)) % p)
                    mat.append(row)
                return det_mod(mat, p)

            d0 = chi((0,) * 6)
            if d0 == 0:
                continue
            d0i = pow(d0, p - 2, p)
            k_chi = 0
            for lam, n in CHI_TAB.items():
                lam6 = tuple(list(lam) + [0] * (6 - len(lam)))
                k_chi = (k_chi + n * chi(lam6) * d0i) % p
            k_chi = k_chi * inv2_30 % p
            k_m = 0
            for lam, q in M_TAB.items():
                lam6 = tuple(list(lam) + [0] * (6 - len(lam)))
                mval = sum(math.prod(pow(ci, e, p) for ci, e in zip(c, perm))
                           for perm in perms_of(lam6)) % p
                k_m = (k_m + q.numerator * pow(q.denominator, p - 2, p)
                       * mval) % p
            assert k_chi == k_m, f"table cross-link failed mod {p}"
            checked += 1
    print(f"[G7] chi^C table == m-basis table at 10 exact GF(p) points, "
          f"2 primes  PASS   [{time.time()-T0:.0f}s]")


# ---- end-to-end numeric ------------------------------------------------------
def K_poly(c):
    tot = 0.0
    for lam, q in M_TAB.items():
        lam6 = tuple(list(lam) + [0] * (6 - len(lam)))
        tot += float(q) * sum(math.prod(ci ** e for ci, e in zip(c, perm))
                              for perm in perms_of(lam6))
    return tot


def W_closed(a, b):
    x = list(a) + list(b)
    c = [math.cos(v) for v in x]
    vca = math.prod(c[i] - c[j] for i in range(3) for j in range(i + 1, 3))
    vcb = math.prod(c[3 + i] - c[3 + j] for i in range(3) for j in range(i + 1, 3))
    sp = math.prod(math.sin(sum(l * v for l, v in zip(L, x)) / 2) for L in CANON)
    return (-2 ** 9 * math.prod(math.sin(v) for v in x) * vca * vcb
            * K_poly(c) / sp)


def gate8_end_to_end():
    worst_w = worst_f = 0.0
    for _ in range(10):
        a = list(RNG.uniform(0.3, 2 * np.pi - 0.3, 3))
        b = list(RNG.uniform(0.3, 2 * np.pi - 0.3, 3))
        wb, wc = W_brute(a, b), W_closed(a, b)
        worst_w = max(worst_w, abs(wb - wc) / max(1.0, abs(wb)))
        e1 = sum(math.cos(v) for v in a)
        f1 = sum(math.cos(v) for v in b)
        F = cross_form(tuple(a), tuple(b))
        pred = -(e1 - f1) ** 2 * wc
        worst_f = max(worst_f, abs(F - pred) / max(1.0, abs(F)))
    assert worst_w < 1e-8 and worst_f < 1e-8
    print(f"[G8] end-to-end: W_brute vs closed {worst_w:.1e}, cross_form vs "
          f"-(e1-f1)^2 W {worst_f:.1e}, 10 points each  PASS"
          f"   [{time.time()-T0:.0f}s]")


if __name__ == "__main__":
    full = "--full" in sys.argv
    gate1_pair_identity()
    gate2_sin_s_lemma()
    gate3_sp_flip_invariance()
    gate4_x_alternates()
    gate5_alternant_readoff(full)
    gate6_weyl_denominator()
    gate7_table_crosslink()
    gate8_end_to_end()
    print("ALL GATES PASS: W = -2^9 (Prod sin) V_c(a) V_c(b) K / SP, "
          "K = 2^-30 Sum n_lambda chi^{C6}_lambda (143 terms)."
          + ("" if full else "  (support sweep: run with --full)"))
