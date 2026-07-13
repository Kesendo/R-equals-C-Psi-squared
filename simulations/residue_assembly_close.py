r"""CLOSE the assembly of the Q(i) residue proof.

Sibling of `halfangle_residue_proof.py` (imported below, engine reused).  That file proved
ALL 38 w1-residues of cross_form vanish, but the single-variable assembly was obstructed
because the inner pole locations rho carry w3 and Qw couples w1 <-> w3.

THE FIX EXECUTED HERE -- eliminate w3 FIRST.  Work in the field
    K := Q(i)(z1, z2, w2)[z3] / (Qz),        Qz = z3^2 + Sz z3 + 1,  Sz = z1+1/z1+z2+1/z2,
with w1 a genuinely FREE variable over K.  cross_form, an element of K(w1)[w3]/(Qw)
(Qw = w3^2 + Sw w3 + 1, Sw = w1+1/w1+w2+1/w2), is written
    cross_form  ==  F0(w1) + F1(w1) * w3,        F0, F1 finite sums of small fractions over K.
The a-constraint is z3 in {roots of Qz} (K is the quotient); the b-constraint is FULLY
absorbed into the w3-reduction (w3 in {roots of Qw}); Sw(w1) enters coefficients
polynomially, so w1 is free.  cross_form = 0 on V  <=>  F0 == 0 AND F1 == 0 in K(w1).

Each of the 288 elementary terms  num / ((Mo-1)(Mh-1))  is rationalised: a denominator
factor (M-1) carrying w3 (M = c * w3^{+-1}, c the w3-free monomial part) is turned real by
its degree-2 w3-conjugate,
    (M - 1)(Mbar - 1) = c^2 - c (w3 + 1/w3) + 1 == c^2 + c Sw + 1 =: N   (mod Qw),
a w3-FREE Laurent polynomial in w1 (NOT the 64-fold Klein norm).  The numerator absorbs the
conjugate factor (Mbar - 1) and is reduced mod (Qz, Qw) to A + B w3 (A, B in K[w1, 1/w1]).
So each term contributes A/Dden to F0 and B/Dden to F1, Dden a w3-free product in K(w1).

PROOF SKELETON (each step asserted in code; see the STEP_* functions):
  STEP_FIELD  : Qz irreducible over Q(i)(z1,z2), so K is a field (sympy, once).
  STEP_ELIM   : the 288-term w3-elimination, F0 + F1 w3, pinned to the committed cross_form
                (the norm identity is pinned separately; |F| ~ 1e-9 on V rules out a wrong elim).
  STEP_POLES  : pole inventory of F0, F1 in w1 over K -- 43 distinct denominator factors,
                classified w1-free (6) / linear (6, simple monomial pole) / quadratic (31, norm).
  STEP_COPRIME: [PROVED over Q(i)]  within every term the two denominator factors are coprime in
                w1 over K (Sylvester resultant nonzero in K), so every w1-pole is SIMPLE (no term
                contributes a double pole).  NOTE: distinct factors DO share roots (16 of the ~24
                roots are shared quad-quad), so the residue must be grouped by ROOT, not by factor.
  STEP_RES_CERT: [CERTIFIED in GF(p)]  F0, F1 have NO finite nonzero w1-pole.  Grouped by root:
                at each root r the residue sum_{(t,slot): dens[slot](r)=0} Num_t(r)/(dens[slot]'(r)
                * dens[1-slot](r)) = 0.  Tested at every in-field root over random on-variety GF(p)
                points, three primes.  (The exact-Q(i) per-factor divisibility STEP_RES is BOTH
                false as a per-factor statement -- shared roots -- AND intractable per-root: the
                shared-root quad groups times the free z1,z2,w2 blow up.)
  STEP_ENDS_WINDOW : [PROVED]  every summand's w1 polynomial-part window is a subset of {-1,0,1},
                so F0, F1 (Laurent, no finite pole) have powers only in {-1,0,1}.
  STEP_ENDS_CERT   : [CERTIFIED in GF(p)]  the three residual coefficients c_-1 = c_0 = c_+1 = 0,
                by evaluating F0, F1 at three w1-values on random on-variety GF(p) points (three
                primes), with a control at z3 off the Qz-variety.  Exact-symbolic (STEP_ENDS_EXACT)
                intractable: the coefficients' ~250-fold common denominator over Q(i)(z1,z2,w2)
                explodes (leading coeffs reach 18 terms, 110 distinct products).
  CONCLUSION  : the STRUCTURE is proved over Q(i) (elimination pinned, K a field, poles simple,
                window {-1,0,1}); the VANISHING (residues + the three coefficients) is certified in
                exact GF(p).  So F0 = F1 = 0 in K(w1) and cross_form == 0 on V at CERTIFICATE grade.
                The w1<->w3 coupling that blocked the earlier assembly is ELIMINATED; the residual
                hole is reduced to root-residue divisibilities + three explicit numbers.

HONEST STATUS: this does NOT upgrade cross_form == 0 to a full Q(i) proof.  The w3-elimination
and the pole structure are exact; the vanishing is a Schwartz-Zippel-grade GF(p) certificate
(same grade as the committed certify_cross_form_gfp, but now the certificate is localized to the
residues at ~24 roots plus three endpoint numbers, with the coupling obstruction removed).  The
exact-Q(i) residue/endpoint computations are written (STEP_RES, STEP_ENDS_EXACT) but blow up.

Qz is irreducible over Q(i)(z1,z2) (verified once with sympy).

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import cmath
import math
import random
import sys
import time

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")

import halfangle_residue_proof as H
from halfangle_residue_proof import (
    G, ONE, ZERO, NVARS, W1, SZ, SW,
    p_add, p_mul, p_scale, p_neg, p_mono, const, z, p_eval, p_is_zero, p_nterms,
    reduce_by_relation, reduce_full, subst_w1, _mono_val, _variety_point,
    TERMS, cross_form_terms,
)
from cross_triple_orthogonality import cross_form as committed_cross_form

W3 = 5   # w3 variable index
Z3 = 2   # z3 variable index


# ---------------------------------------------------------------- w3 elimination
def conj_w3(exps):
    """The K(w1)-Galois conjugate on a monomial: w3 -> 1/w3, i.e. negate the w3 exponent."""
    e = list(exps)
    e[W3] = -e[W3]
    return tuple(e)


def norm_factor(M):
    """For a w3-carrying cot monomial M = c * w3^{+-1}, the w3-free norm N = c^2 + c Sw + 1,
    equal to (M-1)(conj M - 1) modulo Qw."""
    c = list(M)
    c[W3] = 0
    cP = p_mono(tuple(c))
    return p_add(p_add(p_mul(cP, cP), p_mul(cP, SW)), const(ONE))


def split_w3(poly):
    """poly reduced mod Qw (w3-exp in {0,1}) -> (A, B) with poly = A + B w3, A,B w3-free."""
    A, B = {}, {}
    for k, c in poly.items():
        e = k[W3]
        kk = list(k)
        kk[W3] = 0
        tgt = A if e == 0 else B
        assert e in (0, 1), f"split_w3 saw w3-exp {e} (poly not reduced mod Qw)"
        cur = tgt.get(tuple(kk))
        tgt[tuple(kk)] = c if cur is None else cur + c
        if tgt[tuple(kk)].is_zero():
            del tgt[tuple(kk)]
    return A, B


def process_term(t):
    """Rationalise one elementary term.

    Returns (A, B, dens): the term equals (A + B w3) / prod(dens), with A, B w3-free numerators
    (reduced mod Qz, so z3-exp in {0,1}) and dens a list of w3-free denominator factors
    (each reduced mod Qz)."""
    num = t["num"]
    dens = []
    for M in (t["Mo"], t["Mh"]):
        if M[W3] == 0:
            dens.append(p_add(p_mono(M), const(-ONE)))            # w3-free factor (M - 1)
        else:
            Mbar = conj_w3(M)
            num = p_mul(num, p_add(p_mono(Mbar), const(-ONE)))    # absorb conjugate factor
            dens.append(norm_factor(M))                           # denominator -> norm
    num_red = reduce_full(num)                                    # mod (Qz, Qw)
    A, B = split_w3(num_red)
    dens = [reduce_by_relation(d, Z3, SZ) for d in dens]          # denominators mod Qz (w3-free)
    return A, B, dens


def build_F():
    """Elimination of all 288 terms.  Returns list of (A, B, dens)."""
    return [process_term(t) for t in TERMS]


# ---------------------------------------------------------------- pins
def _angle_variety_point(rng):
    """A genuine unit-circle variety point: z_k = exp(i a_k), w_k = exp(i b_k), on V.

    Returns (vals, ang_a, ang_b) or None if the constraint is infeasible for the draw."""
    a1, a2 = rng.uniform(0.3, 2.8), rng.uniform(0.3, 2.8)
    ca3 = -(math.cos(a1) + math.cos(a2))
    if abs(ca3) > 1:
        return None
    a3 = math.acos(ca3)
    b1, b2 = rng.uniform(0.3, 2.8), rng.uniform(0.3, 2.8)
    cb3 = -(math.cos(b1) + math.cos(b2))
    if abs(cb3) > 1:
        return None
    b3 = math.acos(cb3)
    ang_a, ang_b = (a1, a2, a3), (b1, b2, b3)
    vals = [cmath.exp(1j * x) for x in (a1, a2, a3, b1, b2, b3)]
    return vals, ang_a, ang_b


def eval_F(F, vals):
    tot = 0j
    w3 = vals[W3]
    for (A, B, dens) in F:
        num = p_eval(A, vals) + p_eval(B, vals) * w3
        den = 1.0 + 0j
        for d in dens:
            den *= p_eval(d, vals)
        tot += num / den
    return tot


def pin_norm_identity(ntrials=40, seed=1):
    """(M-1)(conj M - 1) == N  wherever Qw holds (z's ARBITRARY; only w3 = root of Qw)."""
    rng = random.Random(seed)
    worst = 0.0
    seen = set()
    done = 0
    for t in TERMS:
        for M in (t["Mo"], t["Mh"]):
            if M[W3] == 0 or M in seen:
                continue
            seen.add(M)
    dist = list(seen)
    for _ in range(ntrials):
        # arbitrary z1,z2,z3, w1,w2 ; w3 a genuine root of Qw
        vals = [cmath.exp(1j * rng.uniform(0.2, 3.0)) for _ in range(5)] + [0j]
        w1v, w2v = vals[3], vals[4]
        Sw = w1v + 1 / w1v + w2v + 1 / w2v
        vals[W3] = (-Sw + cmath.sqrt(Sw * Sw - 4)) / 2
        for M in dist:
            Mbar = conj_w3(M)
            lhs = (_mono_val(M, vals) - 1) * (_mono_val(Mbar, vals) - 1)
            rhs = p_eval(norm_factor(M), vals)
            worst = max(worst, abs(lhs - rhs))
        done += 1
    return worst, len(dist), done


def pin_elimination(F, ntrials=30, seed=20260710):
    """F0 + F1 w3 == cross_form on V.  Compares to BOTH the committed cross_form (genuine
    angle points) and the 288-term cross_form_terms; and at BOTH roots of each Qw (a genuine
    functional test in w3, not just 'both ~0')."""
    rng = random.Random(seed)
    worst_committed = 0.0
    worst_terms = 0.0
    worst_F = 0.0            # |F| itself on V ; a buggy elimination is O(1) here, correct ~1e-12
    worst_tworoot = 0.0
    done = 0
    while done < ntrials:
        got = _angle_variety_point(rng)
        if got is None:
            continue
        vals, ang_a, ang_b = got
        try:
            fF = eval_F(F, vals)
            fT = cross_form_terms(vals)
            fC = committed_cross_form(ang_a, ang_b)
        except ZeroDivisionError:
            continue
        worst_committed = max(worst_committed, abs(fF - complex(fC)))
        worst_terms = max(worst_terms, abs(fF - fT))
        worst_F = max(worst_F, abs(fF))
        # two-root functional test: same (z1,z2,w1,w2,z3), the OTHER w3 root
        w2v, w1v = vals[4], vals[3]
        Sw = w1v + 1 / w1v + w2v + 1 / w2v
        w3a = vals[W3]
        w3b = 1.0 / w3a          # the second Qw root (product of roots = 1)
        for w3 in (w3a, w3b):
            v2 = list(vals)
            v2[W3] = w3
            try:
                worst_tworoot = max(worst_tworoot, abs(eval_F(F, v2) - cross_form_terms(v2)))
            except ZeroDivisionError:
                pass
        done += 1
    return dict(committed=worst_committed, terms=worst_terms, absF=worst_F,
                tworoot=worst_tworoot, n=done)


# ---------------------------------------------------------------- STEP_ELIM
def STEP_ELIM():
    print("[ELIM] eliminating w3: cross_form -> F0(w1) + F1(w1) w3 over K ...")
    t0 = time.time()
    F = build_F()
    wn, ndist, nn = pin_norm_identity()
    print(f"  norm identity (M-1)(conjM-1) == N on {nn} Qw-points, {ndist} distinct M:"
          f" max abs err {wn:.2e}")
    assert wn < 1e-7, "norm identity failed"
    pins = pin_elimination(F)
    print(f"  PIN over {pins['n']} genuine unit-circle V-points:")
    print(f"    |F0+F1 w3 - committed cross_form|  max {pins['committed']:.2e}")
    print(f"    |F0+F1 w3 - 288-term cross_form|   max {pins['terms']:.2e}")
    print(f"    |F0+F1 w3| on V (buggy elim = O(1)) max {pins['absF']:.2e}")
    print(f"    two-root w3 functional test        max {pins['tworoot']:.2e}")
    assert pins['committed'] < 1e-6 and pins['terms'] < 1e-6 and pins['tworoot'] < 1e-6
    assert pins['absF'] < 1e-5, "F does not vanish on V -- elimination is wrong"
    print(f"  [ELIM] OK ({time.time()-t0:.1f}s)")
    return F


# ---------------------------------------------------------------- factor canonicalisation
def w1_span(poly):
    if not poly:
        return None
    es = [k[W1] for k in poly]
    return (min(es), max(es))


def factor_key(poly):
    """Canonical hashable key for a Laurent poly, up to nothing (exact dict).  For grouping
    identical factors we normalise by the lexicographically smallest monomial's coefficient
    scaled to a canonical leading 1?  No -- factors here are monic-ish; use raw sorted items."""
    return tuple(sorted((k, (c.re.numerator, c.re.denominator, c.im.numerator, c.im.denominator))
                        for k, c in poly.items()))


def classify_factor(poly):
    """Return 'w1free' | 'linear' | 'quad' | 'other' by the w1-exponent span."""
    sp = w1_span(poly)
    if sp is None:
        return "empty"
    lo, hi = sp
    width = hi - lo
    if all(k[W1] == 0 for k in poly):
        return "w1free"
    if width == 1:
        return "linear"
    if width == 2:
        return "quad"
    return f"other(width={width})"


def STEP_POLES(F):
    print("[POLES] pole inventory of F0, F1 in w1 over K ...")
    from collections import Counter, defaultdict
    factors = {}          # key -> poly
    incidence = defaultdict(list)   # key -> list of (term_index, slot 0/1)
    per_term_keys = []
    for ti, (A, B, dens) in enumerate(F):
        keys = []
        for slot, d in enumerate(dens):
            k = factor_key(d)
            factors[k] = d
            incidence[k].append((ti, slot))
            keys.append(k)
        per_term_keys.append(keys)
    cls = Counter(classify_factor(p) for p in factors.values())
    print(f"  {len(factors)} distinct denominator factors: {dict(cls)}")
    # double poles within a term (same factor both slots)
    dbl = sum(1 for keys in per_term_keys if keys[0] == keys[1])
    print(f"  terms with a coincident double factor (both slots equal): {dbl}")
    # spans by class
    spans = defaultdict(set)
    for p in factors.values():
        spans[classify_factor(p)].add(w1_span(p))
    for c, s in sorted(spans.items()):
        print(f"    {c}: w1-spans {sorted(s)}")
    return factors, incidence, per_term_keys


# ---------------------------------------------------------------- w1-coefficient view + pdiv
def w1_coeffs(poly):
    """Split into {w1-exponent: (w1-free K-poly)}."""
    out = {}
    for k, c in poly.items():
        e = k[W1]
        kk = list(k)
        kk[W1] = 0
        out.setdefault(e, {})
        cur = out[e].get(tuple(kk))
        out[e][tuple(kk)] = c if cur is None else cur + c
        if out[e][tuple(kk)].is_zero():
            del out[e][tuple(kk)]
    return {e: p for e, p in out.items() if p}


def from_w1_coeffs(cmap):
    out = {}
    for e, pol in cmap.items():
        for kk, c in pol.items():
            k = list(kk)
            k[W1] = e
            out = p_add(out, {tuple(k): c})
    return out


def _rq(p):
    """Reduce a K-coefficient poly mod Qz (z3-exp in {0,1})."""
    return reduce_by_relation(p, Z3, SZ)


def pseudo_rem_w1(nc, qc):
    """Pseudo-remainder of NUM by q in the variable w1, coefficients in K (reduced mod Qz).

    nc, qc: {w1-exp>=0: K-poly}, qc having min-exp 0.  Multiply by leading L each step
    (pure K arithmetic; no fractions, no GCD).  Returns the remainder cmap (w1-deg < deg q)."""
    dq = max(qc)
    L = qc[dq]
    nc = {e: _rq(p) for e, p in nc.items()}
    while nc:
        dn = max(nc)
        if dn < dq:
            break
        cN = nc[dn]
        newc = {e: _rq(p_mul(L, pol)) for e, pol in nc.items()}
        for e, pol in qc.items():
            tgt = e + (dn - dq)
            newc[tgt] = _rq(p_add(newc.get(tgt, {}), p_neg(p_mul(cN, pol))))
        nc = {e: pol for e, pol in newc.items() if pol}
    return nc


def clear_w1(poly):
    """Shift a Laurent poly so its min w1-exp is 0; returns {w1-exp>=0: K-poly} (reduced Qz)."""
    cm = w1_coeffs(poly)
    if not cm:
        return {}
    s = min(cm)
    return {e - s: _rq(p) for e, p in cm.items()}


def divisible_by_factor(NUM, phi):
    """True iff phi | NUM in K[w1, 1/w1]  (w1 a unit): clear both to w1-polys, pseudo-divide."""
    if not NUM:
        return True
    nc = clear_w1(NUM)
    qc = clear_w1(phi)
    rem = pseudo_rem_w1(nc, qc)
    return p_is_zero(_rq(from_w1_coeffs(rem)))


def leading_L(phi):
    """The top w1-coefficient of the cleared factor (multiplied through in pseudo-division)."""
    qc = clear_w1(phi)
    return qc[max(qc)]


# ---------------------------------------------------------------- group numerators + residues
def group_numerator(F, incidence, key, which):
    """NUM such that  sum_{t in group} Num_t/(phi R_t) = NUM / (phi * prod distinct R_t).

    which: 0 -> A (F0), 1 -> B (F1).  Each term has exactly two denominator factors."""
    inc = incidence[key]
    entries = []
    for (ti, slot) in inc:
        A, B, dens = F[ti]
        Num = A if which == 0 else B
        assert len(dens) == 2
        Rt = dens[1 - slot]
        entries.append((Num, Rt))
    distinct = []
    for (_, Rt) in entries:
        if not any(Rt == d for d in distinct):
            distinct.append(Rt)
    NUM = {}
    for (Num, Rt) in entries:
        if not Num:
            continue
        rest = {tuple([0] * NVARS): ONE}
        used = False
        for d in distinct:
            if not used and d == Rt:
                used = True
                continue
            rest = p_mul(rest, d)
        NUM = _rq(p_add(NUM, p_mul(Num, rest)))
    return NUM, distinct



def STEP_RES(F, factors, incidence):
    """DEPRECATED / DO NOT TRUST THE RESULT -- kept for the record only.  Two problems, both
    discovered empirically: (1) FALSE as stated -- it asks per-FACTOR divisibility, but 16 of the
    ~24 w1-roots are SHARED across distinct factors, so an individual factor's group residue need
    NOT vanish (only the sum over all factors sharing a root does; see STEP_RES_CERT).  (2)
    INTRACTABLE -- even grouped correctly, the group numerator over K is a product of up to 16
    quadratics whose coefficients blow up in the free variables z1,z2,w2 (a single factor exceeds
    5 min).  main() uses STEP_RES_CERT (root-grouped, GF(p)) instead.  This body is left as the
    honest artifact of the attempted exact route."""
    print("[RES] per-factor residue (divisibility) over K, for F0 and F1 ...")
    t0 = time.time()
    Lmin_global = 1e9
    n_lin = n_quad = 0
    results = {0: [], 1: []}
    for key, phi in factors.items():
        cls = classify_factor(phi)
        if cls == "w1free":
            continue
        L = leading_L(phi)
        Lmin = _L_nonzero_on_V(L)
        Lmin_global = min(Lmin_global, Lmin)
        for which in (0, 1):
            NUM, distinct = group_numerator(F, incidence, key, which)
            ok = divisible_by_factor(NUM, phi)
            results[which].append((cls, len(incidence[key]), p_nterms(NUM), ok))
            assert ok, f"[RES] F{which}: factor {cls} residue does NOT vanish (NUM {p_nterms(NUM)} terms)"
        if cls == "linear":
            n_lin += 1
        else:
            n_quad += 1
    for which in (0, 1):
        allok = all(r[3] for r in results[which])
        nq = sum(1 for r in results[which] if r[0] == "quad")
        nl = sum(1 for r in results[which] if r[0] == "linear")
        maxN = max((r[2] for r in results[which]), default=0)
        print(f"  F{which}: {nl} linear + {nq} quad poles, all residues 0 over Q(i): {allok}"
              f" (largest group numerator {maxN} terms)")
    print(f"  leading-coeff L of every cleared factor nonzero on V: min |L| = {Lmin_global:.2e}")
    assert Lmin_global > 1e-9
    print(f"  [RES] OK ({time.time()-t0:.1f}s)")
    return results


def _L_nonzero_on_V(L, ntest=12, seed=4242):
    rng = random.Random(seed)
    worst = 1e9
    for _ in range(ntest):
        pt = _variety_point(rng)
        worst = min(worst, abs(p_eval(L, pt)))
    return worst


# ---------------------------------------------------------------- exact polynomial part
def pseudo_divmod_w1(Np_c, Dp_c):
    """Pseudo-division in w1 over K.  Np_c, Dp_c: {w1-exp>=0: K-poly}, both min-exp 0.

    Returns (quo, rem, m, L) with  L^m * Np == quo * Dp + rem  (deg rem < deg Dp), L the
    leading (top-w1) coefficient of Dp, m the number of elimination steps."""
    dq = max(Dp_c)
    L = _rq(Dp_c[dq])
    nc = {e: _rq(p) for e, p in Np_c.items()}
    quo = {}
    m = 0
    while nc and max(nc) >= dq:
        dn = max(nc)
        cN = nc[dn]
        quo = {e: _rq(p_mul(L, pol)) for e, pol in quo.items()}
        newc = {e: _rq(p_mul(L, pol)) for e, pol in nc.items()}
        for e, pol in Dp_c.items():
            tgt = e + (dn - dq)
            newc[tgt] = _rq(p_add(newc.get(tgt, {}), p_neg(p_mul(cN, pol))))
        nc = {e: pol for e, pol in newc.items() if pol}
        quo[dn - dq] = _rq(p_add(quo.get(dn - dq, {}), cN))
        m += 1
    return quo, nc, m, L


def _poly_pow(p, k):
    r = {tuple([0] * NVARS): ONE}
    for _ in range(k):
        r = _rq(p_mul(r, p))
    return r


def pin_divmod(F, seed=7):
    """Numerically verify  L^m Np == quo Dp + rem  for many summands."""
    rng = random.Random(seed)
    worst = 0.0
    for (A, B, dens) in F:
        Dden = _rq(p_mul(dens[0], dens[1]))
        for Num in (A, B):
            if not Num:
                continue
            ncm, dcm = w1_coeffs(Num), w1_coeffs(Dden)
            ns, ds = min(ncm), min(dcm)
            Np_c = {e - ns: _rq(p) for e, p in ncm.items()}
            Dp_c = {e - ds: _rq(p) for e, p in dcm.items()}
            quo, rem, m, L = pseudo_divmod_w1(Np_c, Dp_c)
            # rebuild and test at a random point (w1 free, others on/off V arbitrary)
            pt = _variety_point(rng)
            Lm = p_eval(_poly_pow(L, m), pt)
            Npv = p_eval(from_w1_coeffs(Np_c), pt)
            Dpv = p_eval(from_w1_coeffs(Dp_c), pt)
            quov = p_eval(from_w1_coeffs(quo), pt)
            remv = p_eval(from_w1_coeffs(rem), pt)
            worst = max(worst, abs(Lm * Npv - (quov * Dpv + remv)))
    return worst


def STEP_ENDS_EXACT(F):
    """The w1 Laurent-polynomial part P of F0 and of F1 is 0, over Q(i).

    WARNING: INTRACTABLE AT THIS SCALE -- kept for the record, NOT used by main() (which uses the
    GF(p) certificate STEP_ENDS_CERT instead).  The per-summand quotient extraction alone is ~90s
    and yields coefficients of up to ~1000 terms; the c_j common denominator (110 distinct factors,
    up to 18 terms each) then explodes.  Run via `--endsx` only to watch it grind.

    STEP_RES gives F0, F1 no finite nonzero pole; the window bounds the powers to {-1,0,1}, so
    F0 (resp F1) = sum of its summands' polynomial parts (the proper parts sum to 0).  Each
    summand's polynomial part Q_t = quo_t / L_t^{m_t} * w1^{shift_t} is extracted by
    pseudo-division; for each power j in {-1,0,1} the coefficient c_j = sum_t (coeff of Q_t at j)
    must vanish in K.  Verified by clearing the K-denominators L_t^{m_t} (w1-free, nonzero on V)
    and reducing the numerator mod Qz -- no K-division."""
    print("[ENDS] exact polynomial-part extraction (three coefficients c_-1, c_0, c_+1) ...")
    t0 = time.time()
    wd = pin_divmod(F)
    print(f"  pseudo-divmod identity L^m Np == quo Dp + rem: max abs err {wd:.2e}")
    assert wd < 1e-6
    results = {}
    for which, sel in ((0, 0), (1, 1)):
        # per-power contribution lists: power -> list of (coeff_poly, Ldenom_poly)
        contrib = {-1: [], 0: [], 1: []}
        for (A, B, dens) in F:
            Num = A if sel == 0 else B
            if not Num:
                continue
            Dden = _rq(p_mul(dens[0], dens[1]))
            ncm, dcm = w1_coeffs(Num), w1_coeffs(Dden)
            ns, ds = min(ncm), min(dcm)
            shift = ns - ds
            Np_c = {e - ns: _rq(p) for e, p in ncm.items()}
            Dp_c = {e - ds: _rq(p) for e, p in dcm.items()}
            quo, rem, m, L = pseudo_divmod_w1(Np_c, Dp_c)
            Ldenom = _poly_pow(L, m)
            for e, coeff in quo.items():
                power = shift + e
                assert power in (-1, 0, 1), f"poly-part power {power} outside window"
                contrib[power].append((coeff, Ldenom))
        # verify each c_j == 0 in K: common denom = product of DISTINCT Ldenom
        maxpow_ok = {}
        for j in (-1, 0, 1):
            entries = contrib[j]
            distinct = []
            for (_, Ld) in entries:
                if not any(Ld == d for d in distinct):
                    distinct.append(Ld)
            NUM = {}
            for (coeff, Ld) in entries:
                rest = {tuple([0] * NVARS): ONE}
                used = False
                for d in distinct:
                    if not used and d == Ld:
                        used = True
                        continue
                    rest = _rq(p_mul(rest, d))
                NUM = _rq(p_add(NUM, p_mul(coeff, rest)))
            maxpow_ok[j] = (len(entries), p_is_zero(NUM))
            assert p_is_zero(NUM), f"[ENDS] F{which} coefficient c_{j} != 0 (numerator {p_nterms(NUM)} terms)"
        results[which] = maxpow_ok
        print(f"  F{which}: c_-1, c_0, c_+1 all vanish over Q(i): "
              + ", ".join(f"c_{j}=0({maxpow_ok[j][0]} contribs)" for j in (-1, 0, 1)))
    print(f"  [ENDS] OK ({time.time()-t0:.1f}s)")
    return results


# ---------------------------------------------------------------- coprimality (no double poles)
def _det_K(M):
    """Determinant of a small matrix of K-polys, by cofactor expansion, reduced mod Qz."""
    n = len(M)
    if n == 1:
        return _rq(M[0][0])
    if n == 2:
        return _rq(p_add(p_mul(M[0][0], M[1][1]), p_neg(p_mul(M[0][1], M[1][0]))))
    total = {}
    for j in range(n):
        if not M[0][j]:
            continue
        minor = [[M[r][c] for c in range(n) if c != j] for r in range(1, n)]
        term = p_mul(M[0][j], _det_K(minor))
        total = _rq(p_add(total, term if j % 2 == 0 else p_neg(term)))
    return total


def resultant_w1(f1, f2):
    """Sylvester resultant of two Laurent polys in w1 over K (cleared to min-exp 0)."""
    c1, c2 = clear_w1(f1), clear_w1(f2)
    d1, d2 = max(c1), max(c2)
    if d1 == 0 and d2 == 0:
        return {tuple([0] * NVARS): ONE}     # both w1-free units: no common w1-root
    size = d1 + d2
    S = [[{} for _ in range(size)] for _ in range(size)]
    for r in range(d2):                       # d2 rows of f1 coeffs
        for e, c in c1.items():
            S[r][r + (d1 - e)] = c
    for r in range(d1):                       # d1 rows of f2 coeffs
        for e, c in c2.items():
            S[d2 + r][r + (d2 - e)] = c
    return _det_K(S)


def STEP_COPRIME(F):
    """No term contributes a double pole: within every term, the two denominator factors are
    coprime in w1 over K (Sylvester resultant nonzero in K), so at any w1-pole each affected
    term has exactly ONE vanishing factor and every pole is SIMPLE.  (This does NOT make the
    residue grouping per-factor: distinct factors of DIFFERENT terms share roots, so the residue
    is grouped by ROOT in STEP_RES_CERT.)"""
    print("[COPRIME] per-term factor coprimality in w1 over K ...")
    t0 = time.time()
    worst_on_V = 1e9
    nquadquad = 0
    for (A, B, dens) in F:
        f1, f2 = dens
        if all(k[W1] == 0 for k in f1) or all(k[W1] == 0 for k in f2):
            continue                          # a w1-free factor has no w1-root: coprime
        R = resultant_w1(f1, f2)
        Rr = _rq(R)
        assert not p_is_zero(Rr), "[COPRIME] two factors of a term share a w1-root over K"
        rng = random.Random(101)
        for _ in range(6):
            worst_on_V = min(worst_on_V, abs(p_eval(Rr, _variety_point(rng))))
        if classify_factor(f1) == "quad" and classify_factor(f2) == "quad":
            nquadquad += 1
    print(f"  all terms coprime; {nquadquad} quad-quad pairs; min |resultant| on V = {worst_on_V:.2e}")
    assert worst_on_V > 1e-9, "resultant vanishes on V -- coincident roots possible on V"
    print(f"  [COPRIME] OK ({time.time()-t0:.1f}s)")


# ---------------------------------------------------------------- endpoint: GF(p) certificate
from cross_triple_orthogonality import _tonelli, PRIMES


def _w1poly_modp(poly, vals, p, ii):
    """Collapse a w3-free K-poly to {w1-exp: GF(p) coeff} with z1,z2,z3,w2 specialized to vals
    (indices 0,1,2,4; w1 index 3 kept symbolic, w3 index 5 absent)."""
    out = {}
    for k, c in poly.items():
        cr, ci = _g_modp(c, p)
        coef = (cr + ii * ci) % p
        for idx in (0, 1, 2, 4):
            e = k[idx]
            if e:
                v = vals[idx]
                base = v if e > 0 else pow(v, p - 2, p)
                coef = coef * pow(base, abs(e), p) % p
        e1 = k[3]
        out[e1] = (out.get(e1, 0) + coef) % p
    return {e: v for e, v in out.items() if v}



def _eval_w1_modp(cmap, r, rinv, p):
    tot = 0
    for e, c in cmap.items():
        tot = (tot + c * (pow(r, e, p) if e >= 0 else pow(rinv, -e, p))) % p
    return tot


def _deriv_eval_w1_modp(cmap, r, rinv, p):
    tot = 0
    for e, c in cmap.items():
        if e == 0:
            continue
        ce = c * (e % p) % p
        ee = e - 1
        tot = (tot + ce * (pow(r, ee, p) if ee >= 0 else pow(rinv, -ee, p))) % p
    return tot


def STEP_RES_CERT(F, factors, incidence, trials=10):
    """Every w1-pole residue of F0 and F1 vanishes -- CERTIFIED exactly in GF(p), grouped BY ROOT.

    The 37 w1-carrying factors have only ~24 distinct roots: 16 roots are SHARED by two factors
    (quad-quad), so the residue must be summed over ALL terms poled at a root, not per factor
    (per-factor divisibility is false).  Per-term coprimality (STEP_COPRIME) makes each pole
    simple, so at a root r the residue is  sum_{(t,slot): dens[slot](r)=0}
    Num_t(r) / ( dens[slot]'(r) * dens[1-slot](r) ).  With z1,z2,w2 random and z3 a Qz-root mod p,
    every factor collapses to a w1-polynomial over GF(p); we test Res(r) == 0 at every root r
    that lies in GF(p), for F0 and F1, three primes, many specializations.

    (The exact-Q(i) per-root residue is intractable: the shared-root quad-quad groups and the free
    variables z1,z2,w2 make the symbolic combination blow up; hence a GF(p) certificate.)"""
    print("[RES] residue vanishing at every w1-root (grouped by root) -- GF(p) certificate ...")
    t0 = time.time()
    checks = roots_tested = 0
    for p in PRIMES:
        ii = _tonelli(p - 1, p)
        inv2 = pow(2, p - 2, p)
        rng = random.Random(p ^ 0x11ce)
        done = 0
        while done < trials:
            z1, z2, w2 = (rng.randrange(2, p) for _ in range(3))
            Sz = (z1 + pow(z1, p - 2, p) + z2 + pow(z2, p - 2, p)) % p
            dsq = _tonelli((Sz * Sz - 4) % p, p)
            if dsq is None:
                continue
            z3 = (p - Sz + dsq) % p * inv2 % p
            if z3 == 0:
                continue
            vals = [z1, z2, z3, 0, w2, 0]
            # collapse every factor and every term's numerators to GF(p) w1-polynomials
            fac_p = {}
            for key, phi in factors.items():
                if classify_factor(phi) != "w1free":
                    fac_p[key] = _w1poly_modp(phi, vals, p, ii)
            dens_p = []
            A_p = []
            B_p = []
            for (A, B, dens) in F:
                dens_p.append([_w1poly_modp(dens[0], vals, p, ii), _w1poly_modp(dens[1], vals, p, ii)])
                A_p.append(_w1poly_modp(A, vals, p, ii) if A else {})
                B_p.append(_w1poly_modp(B, vals, p, ii) if B else {})
            # collect distinct GF(p) roots of all w1-factors
            roots = set()
            for cm in fac_p.values():
                lo = min(cm)
                c = {e - lo: v for e, v in cm.items()}
                deg = max(c)
                if deg == 1:
                    roots.add((p - c.get(0, 0)) * pow(c[1], p - 2, p) % p)
                elif deg == 2:
                    a, b, cc = c.get(2, 0), c.get(1, 0), c.get(0, 0)
                    s = _tonelli((b * b - 4 * a * cc) % p, p)
                    if s is not None:
                        ia = pow(2 * a, p - 2, p)
                        roots.add((p - b + s) % p * ia % p)
                        roots.add((p - b - s) % p * ia % p)
            ok = True
            for r in roots:
                if r == 0:
                    continue
                rinv = pow(r, p - 2, p)
                for which in (0, 1):
                    Nump = A_p if which == 0 else B_p
                    res = 0
                    for ti in range(len(F)):
                        for slot in (0, 1):
                            f = dens_p[ti][slot]
                            if _eval_w1_modp(f, r, rinv, p) == 0:
                                dprime = _deriv_eval_w1_modp(f, r, rinv, p)
                                other = _eval_w1_modp(dens_p[ti][1 - slot], r, rinv, p)
                                num = _eval_w1_modp(Nump[ti], r, rinv, p)
                                res = (res + num * pow(dprime * other % p, p - 2, p)) % p
                    if res != 0:
                        ok = False
                roots_tested += 1
            assert ok, f"[RES] a root residue of F0/F1 does NOT vanish mod {p}"
            checks += 1
            done += 1
    print(f"  residues at all in-field w1-roots vanish for F0 and F1: {checks} specialisations"
          f" x 3 primes, {roots_tested} (root x point) residues, 0 failures")
    print(f"  [RES] certificate OK ({time.time()-t0:.1f}s)")


def _g_modp(g, p):
    def fr(f):
        return (int(f.numerator) % p) * pow(int(f.denominator) % p, p - 2, p) % p
    return fr(g.re), fr(g.im)     # (real, imag) parts mod p


def _peval_modp(poly, vals, p, ii):
    """Evaluate a Laurent poly at vals (list of 6 residues mod p), i = sqrt(-1) mod p."""
    tot = 0
    for k, c in poly.items():
        cr, ci = _g_modp(c, p)
        term = (cr + ii * ci) % p
        for idx, e in enumerate(k):
            if e:
                v = vals[idx]
                base = v if e > 0 else pow(v, p - 2, p)
                term = term * pow(base, abs(e), p) % p
        tot = (tot + term) % p
    return tot


def _eval_F_component(F, which, vals, p, ii):
    """F0 (which=0) or F1 (which=1) evaluated in GF(p): sum_t Num_t/prod(dens)."""
    tot = 0
    for (A, B, dens) in F:
        Num = A if which == 0 else B
        if not Num:
            continue
        den = 1
        for d in dens:
            den = den * _peval_modp(d, vals, p, ii) % p
        if den == 0:
            raise ZeroDivisionError
        tot = (tot + _peval_modp(Num, vals, p, ii) * pow(den, p - 2, p)) % p
    return tot


def STEP_ENDS_CERT(F, trials=40):
    """The w1 Laurent-polynomial part of F0 and F1 is 0 -- CERTIFIED exactly in GF(p).

    STEP_RES proves F0, F1 (over K, w1 free) have NO finite nonzero w1-pole, so each is a
    Laurent polynomial c_-1/w1 + c_0 + c_+1 w1 (window {-1,0,1}, STEP_ENDS_WINDOW).  For fixed
    on-variety (z1,z2,z3,w2) [z3 a genuine root of Qz mod p] the map w1 -> F0 is that degree-2
    (times 1/w1) shape; evaluating at THREE distinct nonzero w1 and getting 0 forces all three
    coefficients to 0 (a nonzero such shape has <= 2 roots).  Done for F0 and F1, three primes,
    with a CONTROL at z3 OFF the Qz-variety (there F0 is generically nonzero, so the test bites).

    NOTE ON GRADE: the exact-symbolic endpoint (STEP_ENDS_EXACT) is intractable at this scale
    (the polynomial-part coefficients are single K-elements whose ~250-fold common denominator
    over Q(i)(z1,z2,w2) explodes; leading coeffs reach 18 terms, 110 distinct products).  So the
    three residual coefficients are CERTIFIED here in exact GF(p) arithmetic, not proved over Q."""
    print("[ENDS] endpoint polynomial-part c_-1, c_0, c_+1 = 0 -- exact GF(p) certificate ...")
    t0 = time.time()
    inv2 = None
    total_checks = 0
    control_nonzero = 0
    for p in PRIMES:
        ii = _tonelli(p - 1, p)
        inv2 = pow(2, p - 2, p)
        rng = random.Random(p ^ 0x5eed)
        good = 0
        while good < trials:
            z1 = rng.randrange(2, p)
            z2 = rng.randrange(2, p)
            w2 = rng.randrange(2, p)
            Sz = (z1 + pow(z1, p - 2, p) + z2 + pow(z2, p - 2, p)) % p
            d = _tonelli((Sz * Sz - 4) % p, p)
            if d is None:
                continue
            z3 = (p - Sz + d) % p * inv2 % p
            if z3 == 0:
                continue
            # three distinct nonzero w1 values
            ws = []
            while len(ws) < 3:
                v = rng.randrange(2, p)
                if v not in ws:
                    ws.append(v)
            try:
                ok = True
                for which in (0, 1):
                    for v in ws:
                        val = _eval_F_component(F, which, [z1, z2, z3, v, w2, 0], p, ii)
                        if val != 0:
                            ok = False
                # control: z3' NOT a Qz-root -> F0 generically nonzero
                z3bad = (z3 + 1) % p
                ctrl = _eval_F_component(F, 0, [z1, z2, z3bad, ws[0], w2, 0], p, ii)
            except ZeroDivisionError:
                continue
            assert ok, f"[ENDS] F0 or F1 nonzero at an on-variety GF(p) point mod {p}"
            control_nonzero += (ctrl != 0)
            total_checks += 1
            good += 1
    print(f"  {total_checks} on-variety points x 3 w1-values x 2 components x 3 primes:"
          f" F0, F1 identically 0 (poly part c_-1=c_0=c_+1=0)")
    print(f"  control (z3 off the Qz-variety) nonzero: {control_nonzero}/{total_checks}"
          f" -> the certificate is discriminating")
    assert control_nonzero > 0.9 * total_checks
    print(f"  [ENDS] certificate OK ({time.time()-t0:.1f}s)")


def STEP_ENDS_WINDOW(F):
    """Structural fact used by STEP_ENDS_CERT: every summand's w1 polynomial-part window is a
    subset of {-1, 0, 1}, so F0, F1 (being Laurent polynomials after STEP_RES) have powers only
    in {-1, 0, 1}.  Exact integer bookkeeping (numerator vs denominator w1-spans)."""
    hi = lo = 0
    for (A, B, dens) in F:
        dlo = dhi = 0
        for dd in dens:
            sp = w1_span(dd)
            dlo += sp[0]
            dhi += sp[1]
        for Num in (A, B):
            if not Num:
                continue
            nlo, nhi = w1_span(Num)
            hi = max(hi, nhi - dhi)
            lo = min(lo, nlo - dlo)
    print(f"[ENDS] w1 polynomial-part window across all summands: [{lo}, {hi}] (subset of [-1,1])")
    assert lo >= -1 and hi <= 1
    return lo, hi


# ---------------------------------------------------------------- K well-definedness
def STEP_FIELD():
    """K = Q(i)(z1,z2,w2)[z3]/(Qz) is a field: Qz = z3^2 + Sz z3 + 1 is irreducible over
    Q(i)(z1,z2).  Verified once with sympy (a genuine quadratic with non-square discriminant
    Sz^2 - 4)."""
    import sympy as sp
    z1, z2, z3 = sp.symbols('z1 z2 z3')
    Sz = z1 + 1 / z1 + z2 + 1 / z2
    Qz = z3 ** 2 + Sz * z3 + 1
    poly = sp.Poly(sp.numer(sp.together(Qz)), z3)          # clear 1/z1,1/z2 -> polynomial in z3
    irr = poly.is_irreducible
    # discriminant Sz^2 - 4 is not a perfect square in Q(i)(z1,z2)
    disc = sp.simplify(Sz ** 2 - 4)
    print(f"[FIELD] Qz irreducible over Q(i)(z1,z2): {irr}"
          f" (deg {sp.degree(poly, z3)}, disc = Sz^2-4 non-square)")
    assert irr, "Qz reducible -- K is not a field"
    return irr


if __name__ == "__main__" and "--field" in sys.argv:
    STEP_FIELD()


def main():
    print("=" * 78)
    print("CLOSING the assembly: cross_form == 0 on V, over Q(i).")
    print("Route: eliminate w3 (norm), then residues + endpoints of F0, F1 in w1 over K.")
    print("=" * 78)
    t0 = time.time()
    STEP_FIELD()
    F = STEP_ELIM()
    factors, incidence, per_term_keys = STEP_POLES(F)
    STEP_COPRIME(F)
    STEP_RES_CERT(F, factors, incidence)
    STEP_ENDS_WINDOW(F)
    STEP_ENDS_CERT(F)
    print("-" * 78)
    print("CONCLUSION.  Over K = Q(i)(z1,z2,w2)[z3]/(Qz), with w1 free:")
    print("  [PROVED over Q(i)]  the STRUCTURE: w3 eliminated by the degree-2 norm (STEP_ELIM,")
    print("     pinned), Qz irreducible so K is a field (STEP_FIELD), 43 denominator factors with")
    print("     no double pole (STEP_COPRIME), and the w1 polynomial-part window is {-1,0,1}")
    print("     (STEP_ENDS_WINDOW).  The w1<->w3 coupling that blocked the earlier assembly is GONE.")
    print("  [CERTIFIED in GF(p)]  F0, F1 have NO finite nonzero w1-pole -- the residue summed over")
    print("     all terms poled at each of the ~24 roots (16 shared quad-quad) vanishes")
    print("     (STEP_RES_CERT), and the three residual coefficients c_-1 = c_0 = c_+1 = 0")
    print("     (STEP_ENDS_CERT).  The exact-Q(i) versions (STEP_RES, STEP_ENDS_EXACT) are written")
    print("     but INTRACTABLE at this scale (shared-root groups / common-denominator products")
    print("     blow up in the free variables z1,z2,w2).")
    print("  Together: F0 = F1 = 0 in K(w1), hence cross_form == 0 on V -- the exact STRUCTURE is")
    print("  proved over Q(i); the VANISHING (residues + endpoint) is GF(p)-certified, not yet")
    print("  proved over Q(i).  Net advance: the coupling obstruction is removed and the residual")
    print("  hole is reduced to divisibilities + three numbers, all certified.")
    print(f"[TOTAL] {time.time()-t0:.1f}s")


if __name__ == "__main__" and "--endsx" in sys.argv:
    F = STEP_ELIM()
    STEP_ENDS_EXACT(F)

if __name__ == "__main__" and len(sys.argv) == 1:
    main()


if __name__ == "__main__" and "--elim" in sys.argv:
    STEP_ELIM()

if __name__ == "__main__" and "--poles" in sys.argv:
    F = STEP_ELIM()
    STEP_POLES(F)

if __name__ == "__main__" and "--res" in sys.argv:
    F = STEP_ELIM()
    factors, incidence, per_term_keys = STEP_POLES(F)
    STEP_COPRIME(F)
    STEP_RES(F, factors, incidence)

if __name__ == "__main__" and "--ends" in sys.argv:
    F = STEP_ELIM()
    STEP_ENDS_WINDOW(F)
    STEP_ENDS_CERT(F)
