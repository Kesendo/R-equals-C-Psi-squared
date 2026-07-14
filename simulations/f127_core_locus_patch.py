r"""A9 -- the core-identity degenerate-locus patch (hostile-referee round).

BACKGROUND.  `f127_core_identity.py` proves  T == 0  on the variety
    V = { Qz(z3) = 0 } cap { P(z3) = 0 },
    Qz = z3^2 + Sz z3 + 1,   Sz = z1+1/z1+z2+1/z2,
    P  = A z3^2 + B z3 + C,   A = m0^2,  B = Sw m0,  C = 1,   m0 = z1 z2 w1 w2,
over K = Q(i)(z1,z2,w1,w2), by the route
    N mod Qz = r0 + r1 z3,   E := r0 (B - A Sz) - r1 (C - A),   Res := Res_z3(Qz, P),
and the exact divisibility  Res | E  (re-confirmed: baseline remainder is zero).

THE REFEREE'S HOLE.  The divisibility only settles T=0 where the common root can be SOLVED,
    z3 = (C - A) / (A Sz - B),   which needs  A Sz - B != 0.
Indeed  N|_{z3=common root} = E / (B - A Sz)  (derived below), so N=0 <=> E=0 ONLY off
{A Sz = B}.  On the sublocus  {A Sz = B} cap V  the argument is vacuous: E == 0 identically
there and proves nothing.

THIS PATCH delivers three separate, independently-asserted results.
  (a) THE FORCING CHAIN, derived symbolically (not assumed):
        A Sz = B  AND common root  ==>  A = C  ==>  m0^2 = 1  ==>  Sz = m0 Sw.
        (And then B = A Sz = Sz and C = A = 1, so P == Qz: the two quadratics COINCIDE.)
  (b) B - A Sz = m0 (Sw - m0 Sz), and {Sw = m0 Sz} cap V is a PROPER closed sublocus of every
        irreducible component (branch) of {Res = 0}: EXACT certificate -- no irreducible factor
        f of Res divides g := (Sw - m0 Sz), so {f=0} is NOT contained in {g=0}, the good locus
        {B - A Sz != 0} is dense on every branch, and T (rational, denominator not identically
        zero) = 0 on the whole branch by Zariski closure.  This already covers the degenerate
        points (they lie in the closure of the good locus, being codim >= 2 in each branch).
  (c) T = 0 DIRECTLY on the degenerate sublocus  {m0^2 = 1, Sz = m0 Sw, Qz = 0}, as a genuinely
        different (smaller) identity: there P == Qz, so BOTH roots of Qz lie on the variety and
        T=0 forces  r0 == r1 == 0.  Proved exactly per sign eps = m0 in {+1,-1}: eliminate w2 by
        w2 = eps/(z1 z2 w1) (which bakes in m0 = eps, hence m0^2 = 1), leaving the single relation
        R_eps := Sz - eps Sw; then every irreducible factor of R_eps divides the numerators of
        r0 and r1 (exact ideal membership), i.e. r0 = r1 = 0 on the sublocus.  A numeric gate on
        genuine sublocus points confirms |T| ~ 0 first.

If (c) had FAILED (r0, r1 not vanishing on a sublocus component), that would be reported LOUDLY;
in that event (b)'s closure argument alone is the fallback (the sublocus is codim >= 2 in each
Res-branch, so still inside the closure of the good locus).  Both routes in fact succeed.

Runs the same product-form assembly as `f127_core_identity.build_numerator`; NO
sympy cancel/together is called on the raw unassembled 9-product sum (documented 2h hang).

Authors: Thomas Wicht and Claude, 2026-07-14.
"""
import cmath
import importlib
import random
import sys
import time

import sympy as sp

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
core = importlib.import_module("f127_core_identity")

z1, z2, z3, w1, w2 = core.z1, core.z2, core.z3, core.w1, core.w2
I = sp.I

m0 = z1 * z2 * w1 * w2
Sz = z1 + 1 / z1 + z2 + 1 / z2
Sw = w1 + 1 / w1 + w2 + 1 / w2
A = m0 ** 2
B = sp.together(Sw * m0)
C = sp.Integer(1)
Qz = z3 ** 2 + Sz * z3 + 1
P = A * z3 ** 2 + B * z3 + C
g = Sw - m0 * Sz            # B - A Sz = m0 * g ; the degenerate factor


def numer_cleared(expr):
    """Clear all Laurent denominators (monomial units on the torus) -> a genuine polynomial."""
    return sp.expand(sp.numer(sp.cancel(sp.together(expr))))


# ============================================================ (a) the forcing chain
def part_a():
    print("=" * 78)
    print("(a) THE FORCING CHAIN  [A Sz = B & common root ==> A=C ==> m0^2=1 ==> Sz=m0 Sw]")
    print("=" * 78)

    # The z3^2-eliminating combination P - A Qz.
    lin = sp.expand(P - A * Qz)
    target = (B - A * Sz) * z3 + (C - A)
    assert sp.expand(lin - target) == 0, "P - A Qz != (B - A Sz) z3 + (C - A)"
    print("  [1] P - A*Qz = (B - A*Sz) z3 + (C - A)                         [identity, verified]")

    # coefficient of z3 in that combination is exactly B - A Sz
    assert sp.expand(sp.Poly(lin, z3).coeff_monomial(z3) - (B - A * Sz)) == 0
    # constant term is exactly C - A
    assert sp.expand(sp.Poly(lin, z3).coeff_monomial(1) - (C - A)) == 0
    print("  [2] its z3-coeff = (B - A*Sz), its constant = (C - A)          [verified]")

    # At a common root z3*, Qz=0 and P=0, so (P - A Qz)|_{z3*} = 0.  The ONLY z3-dependence of
    # (P - A Qz) is its leading term (B - A Sz) z3; on {B - A Sz = 0} it drops and the whole
    # combination is the CONSTANT (C - A), so a common root forces C - A = 0, i.e. A = C.
    assert sp.expand(lin - (B - A * Sz) * z3 - (C - A)) == 0   # P - A Qz minus its z3-term = C - A
    print("  [3] on {B = A*Sz}: P - A*Qz = (C - A), a CONSTANT; a common root forces C - A = 0")
    print("      ==> A = C.                                                  [forced, verified]")

    # A = C  with A = m0^2, C = 1  <=>  m0^2 = 1.
    assert sp.expand((C - A) - (1 - m0 ** 2)) == 0
    print("  [4] C - A = 1 - m0^2, so A = C  <=>  m0^2 = 1.                  [verified]")

    # B - A Sz = m0 (Sw - m0 Sz)  [same identity used in (b)].  With m0^2 = 1 (so m0 = 1/m0):
    #   B - A Sz = 0  <=>  Sw - m0 Sz = 0  <=>  Sz - m0 Sw = 0   (multiply by m0).
    assert numer_cleared((B - A * Sz) - m0 * g) == 0
    # the m0^2 = 1 pivot:  m0 (Sw - m0 Sz) = m0 Sw - m0^2 Sz  (identity), = m0 Sw - Sz mod (m0^2-1)
    assert numer_cleared(m0 * g - (m0 * Sw - m0 ** 2 * Sz)) == 0
    print("  [5] B - A*Sz = m0 (Sw - m0 Sz) = m0 Sw - m0^2 Sz; under m0^2 = 1 this is")
    print("      m0 Sw - Sz = -(Sz - m0 Sw).  So B = A*Sz forces  Sz = m0 Sw.   [forced, verified]")

    # consequence: on the sublocus B = A Sz = m0^2 Sz and C = A = m0^2, so
    #   P - Qz = (A - 1) z3^2 + (A Sz - Sz) z3 + (C - 1) = (m0^2 - 1)(z3^2 + Sz z3),
    # which is 0 exactly when m0^2 = 1: the two quadratics COINCIDE (P == Qz).
    diff = (A * z3 ** 2 + (A * Sz) * z3 + C) - Qz
    assert numer_cleared(diff - (m0 ** 2 - 1) * (z3 ** 2 + Sz * z3)) == 0
    print("  [6] P - Qz = (m0^2 - 1)(z3^2 + Sz z3) on {B = A*Sz}: with m0^2 = 1, P == Qz")
    print("      (quadratics COINCIDE).                                       [verified]")
    print("  (a) OK.\n")


# ============================================================ shared: r0, r1 of N mod Qz
def build_r0_r1():
    print("[assemble] building N (product-form) and reducing N mod Qz over K ...")
    t0 = time.time()
    N, poles = core.build_numerator()
    N = sp.expand(N)
    K = sp.FractionField(sp.QQ_I, [z1, z2, w1, w2])
    Sz_k = z1 + 1 / z1 + z2 + 1 / z2
    Np = sp.Poly(N, z3, domain=K)
    Qzp = sp.Poly(z3 ** 2 + Sz_k * z3 + 1, z3, domain=K)
    _q, r = sp.div(Np, Qzp)
    r0 = r.as_expr().coeff(z3, 0)
    r1 = r.as_expr().coeff(z3, 1)
    print(f"  N: {len(sp.Add.make_args(N))} terms; r0, r1 obtained ({time.time()-t0:.1f}s)")
    return N, r0, r1, poles


# ============================================================ (b) properness of the good locus
def part_b(r0, r1):
    print("=" * 78)
    print("(b) B - A Sz = m0 (Sw - m0 Sz);  {Sw = m0 Sz} PROPER on every Res-branch")
    print("=" * 78)

    # the identity behind the whole reduction: the common root of Qz and P (when B - A Sz != 0)
    # is z3* = (C - A)/(A Sz - B) = -(C - A)/(B - A Sz), the root of the linear reduction
    # (B - A Sz) z3 + (C - A) = 0; and there N|_{z3*} = (r0 + r1 z3*) = E/(B - A Sz).
    z3_star = (C - A) / (A * Sz - B)
    assert sp.simplify((B - A * Sz) * z3_star + (C - A)) == 0
    Ediff = (B - A * Sz)
    E = r0 * Ediff - r1 * (C - A)
    assert sp.cancel((r0 + r1 * z3_star) * Ediff - E) == 0
    print("  [1] common root z3* = (C-A)/(A Sz - B); N|_{z3*} (B - A Sz) = E")
    print("      (so N=0 <=> E=0 exactly where B - A Sz != 0).                 [verified]")

    # the factorization B - A Sz = m0 (Sw - m0 Sz)
    assert numer_cleared((B - A * Sz) - m0 * g) == 0
    print("  [2] B - A Sz = m0 (Sw - m0 Sz);  m0 a unit, so {B - A Sz = 0} = {Sw = m0 Sz}.  [ok]")

    # Res and g cleared to polynomials over Q(i)
    Res = sp.resultant(Qz, P, z3)
    Res_poly = numer_cleared(Res)
    g_poly = numer_cleared(g)
    print(f"  [3] Res cleared: {len(sp.Add.make_args(sp.expand(Res_poly)))} terms;"
          f" g = Sw - m0 Sz cleared: {len(sp.Add.make_args(sp.expand(g_poly)))} terms")

    # PROPERNESS CERTIFICATE by a single exact gcd (no full factorization needed).
    #   A branch f of {Res=0} lies in {g=0}  <=>  f is an irreducible factor of BOTH Res and g
    #   <=>  f | gcd(Res, g).  So SOME branch lies in {g=0}  <=>  gcd(Res, g) has a non-monomial
    #   factor.  Hence: NO branch lies in {g=0}  <=>  gcd(Res, g) is a MONOMIAL (its only zeros
    #   are off the torus, where z1..w2 = 0, i.e. not genuine components of our variety).
    Rp = sp.Poly(Res_poly, z1, z2, w1, w2, domain=sp.QQ_I)
    Gp = sp.Poly(g_poly, z1, z2, w1, w2, domain=sp.QQ_I)
    d = sp.gcd(Rp, Gp)
    d_is_monomial = d.is_monomial            # constant (unit) or pure monomial -> off-torus only
    print(f"  [4] gcd(Res, g) = {d.as_expr()}  ->  {'monomial/unit' if d_is_monomial else 'NON-monomial'}")
    assert d_is_monomial, ("gcd(Res, g) has a non-monomial factor -- a Res-branch lies in"
                           " {Sw = m0 Sz}; the closure argument would fail on it")
    print("  [5] no non-monomial branch of {Res=0} lies in {Sw = m0 Sz}: the good locus"
          " {B - A Sz != 0} is")
    print("      DENSE (open, nonempty) on every branch.  E = 0 on {Res=0} (baseline Res|E) gives")
    print("      N = 0 on the good locus; N polynomial => N = 0 on the closure = whole branch;")
    print("      T = N/denom (denom not identically 0) => T = 0 on every branch of {Res=0}.")
    print("  (b) OK.\n")


# ============================================================ (c) direct on the degenerate sublocus
def _alpha_z(zs, i):
    j, l = [t for t in range(3) if t != i]
    return (zs[l] - zs[j]) * (zs[j] * zs[l] + 1) / (2j * zs[j] * zs[l])


def _T_numeric(zv, wv):
    """T(a;b) in z/w coordinates: sum (-1)^{i+j} alpha^z_i alpha^w_j * i (z_i w_j+1)/(z_i w_j-1)."""
    tot = 0j
    for i in range(3):
        ai = _alpha_z(zv, i)
        for j in range(3):
            aj = _alpha_z(wv, j)
            u = zv[i] * wv[j]
            tot += ((-1) ** (i + j)) * ai * aj * 1j * (u + 1) / (u - 1)
    return tot


def _sublocus_point(rng, eps):
    """A genuine point on {m0=eps, Sz=eps Sw, Qz=0}: z1,z2 free on the circle; z3 a Qz-root;
    w2 = eps/(z1 z2 w1); w1 a root of (1+eps c) w1^2 - (eps Sz) w1 + (1+eps/c) = 0, c = z1 z2;
    w3 = 1/(m0 z3) = 1/(eps z3).  Returns (zv, wv) or None."""
    z1v = cmath.exp(1j * rng.uniform(0.3, 2.8))
    z2v = cmath.exp(1j * rng.uniform(0.3, 2.8))
    Szv = z1v + 1 / z1v + z2v + 1 / z2v
    z3v = (-Szv + cmath.sqrt(Szv * Szv - 4)) / 2
    c = z1v * z2v
    a2 = 1 + eps * c
    a0 = 1 + eps / c
    if abs(a2) < 1e-9:
        return None
    T0 = eps * Szv                       # target Sw = eps Sz  (eps = 1/eps)
    disc = cmath.sqrt(T0 * T0 - 4 * a2 * a0)
    for sgn in (+1, -1):
        w1v = (T0 + sgn * disc) / (2 * a2)
        if abs(w1v) < 1e-9:
            continue
        w2v = eps / (c * w1v)
        # verify the relations numerically
        m0v = z1v * z2v * w1v * w2v
        Swv = w1v + 1 / w1v + w2v + 1 / w2v
        if abs(m0v - eps) > 1e-8 or abs(Szv - eps * Swv) > 1e-7:
            continue
        w3v = 1 / (eps * z3v)
        # verify the sheet w3 = 1/(m0 z3)
        if abs(w3v - 1 / (m0v * z3v)) > 1e-8:
            continue
        return [z1v, z2v, z3v], [w1v, w2v, w3v]
    return None


def part_c(r0, r1):
    print("=" * 78)
    print("(c) T = 0 DIRECTLY on the degenerate sublocus {m0^2=1, Sz=m0 Sw, Qz=0}")
    print("=" * 78)

    # ---- numeric gate: the POLE-FREE witnesses r0, r1 ~ 0 on the sublocus ----
    # SUBTLETY: T ITSELF is singular on the sublocus.  There z3 w3 = 1/eps (w3 = 1/(m0 z3),
    # m0 = eps), so for eps = +1 the (i=3,j=3) cotangent's denominator z3 w3 - 1 vanishes:
    # evaluating T directly blows up (~1e31).  The pole-free statement is that T, as a rational
    # function in K[z3]/(Qz), i.e. its reduction coefficients r0, r1 in K, vanishes on the
    # sublocus.  T = 0 on the variety means exactly r0 = r1 = 0 there.  We gate on r0, r1.
    f0 = sp.lambdify((z1, z2, w1, w2), r0, "numpy")
    f1 = sp.lambdify((z1, z2, w1, w2), r1, "numpy")
    rng = random.Random(1234)
    ctrl = 0.0
    for _ in range(30):
        pp = [cmath.exp(1j * rng.uniform(0.3, 2.8)) for _ in range(4)]
        ctrl = max(ctrl, abs(complex(f0(*pp))), abs(complex(f1(*pp))))
    worst = 0.0
    tpole = 0.0
    counts = {1: 0, -1: 0}
    tries = 0
    while min(counts.values()) < 12 and tries < 6000:
        tries += 1
        eps = 1 if tries % 2 else -1
        got = _sublocus_point(rng, eps)
        if got is None:
            continue
        zv, wv = got
        pt = (zv[0], zv[1], wv[0], wv[1])
        worst = max(worst, abs(complex(f0(*pt))), abs(complex(f1(*pt))))
        counts[eps] += 1
        if eps == 1:                       # record the (3,3) cot blow-up for the record
            try:
                tpole = max(tpole, abs(_T_numeric(zv, wv)))
            except ZeroDivisionError:
                tpole = float("inf")
    print(f"  [gate] control |r0|,|r1| off-sublocus ~ {ctrl:.2e};"
          f" on {counts} sublocus points worst |r0|,|r1| = {worst:.2e}")
    print(f"         (for the record, |T| at eps=+1 sublocus points ~ {tpole:.2e}:"
          f" the (3,3) cotangent pole -- T is singular, r0,r1 are the pole-free witnesses)")
    assert ctrl > 1.0 and worst < 1e-6, "r0, r1 do NOT vanish on the sublocus -- STOP, report loudly"

    # ---- exact: r0 == r1 == 0 on the sublocus, per sign eps ----
    # On the sublocus P == Qz, so BOTH roots of Qz lie on the variety; T=0 for both forces
    # r0 + r1 z3 = 0 at two distinct z3, i.e. r0 = r1 = 0 (as functions on the sublocus).
    n0 = numer_cleared(r0)
    n1 = numer_cleared(r1)
    print(f"  [exact] numerators: r0 -> {len(sp.Add.make_args(sp.expand(n0)))} terms,"
          f" r1 -> {len(sp.Add.make_args(sp.expand(n1)))} terms")

    for eps in (1, -1):
        # eliminate w2 by w2 = eps/(z1 z2 w1): bakes in m0 = eps (=> m0^2 = 1).
        subs = {w2: sp.Integer(eps) / (z1 * z2 * w1)}
        # remaining relation R_eps := Sz - eps Sw, cleared to a polynomial in z1,z2,w1
        R_eps = numer_cleared((Sz - eps * Sw).subs(subs))
        n0e = numer_cleared(n0.subs(subs))
        n1e = numer_cleared(n1.subs(subs))
        # factor R_eps; r_i vanishes on {R_eps=0} iff every irreducible factor divides r_i.
        Rfl = sp.factor_list(R_eps, gaussian=True)
        rfacs = []
        for fac, mult in Rfl[1]:
            pf = sp.Poly(fac, z1, z2, w1, domain=sp.QQ_I)
            if pf.is_monomial:
                continue
            rfacs.append(fac)
        assert rfacs, f"R_eps (eps={eps}) has no non-monomial factor -- unexpected"
        print(f"  [eps={eps:+d}] R_eps = Sz - eps Sw has {len(rfacs)} non-monomial factor(s);"
              f" testing r0, r1 membership ...")
        for name, ne in (("r0", n0e), ("r1", n1e)):
            nep = sp.Poly(ne, z1, z2, w1, domain=sp.QQ_I)
            for k, fac in enumerate(rfacs):
                fp = sp.Poly(fac, z1, z2, w1, domain=sp.QQ_I)
                _, rem = sp.reduced(nep, [fp], z1, z2, w1, order="lex", domain=sp.QQ_I)
                ok = (rem == sp.Poly(0, z1, z2, w1, domain=sp.QQ_I))
                assert ok, (f"[STOP-LOUD] {name} NOT divisible by R_eps factor[{k}] at eps={eps}"
                            f" -- {name} does NOT vanish on a sublocus component!")
            print(f"      {name}: every R_eps factor divides it -> {name} == 0 on the sublocus. [ok]")
    print("  [conclude] r0 = r1 = 0 on {m0^2=1, Sz=m0 Sw}; N mod Qz = 0 there, so T = 0")
    print("  on the degenerate sublocus for BOTH roots of Qz.  (c) OK.\n")


def main():
    t0 = time.time()
    part_a()
    _N, r0, r1, _poles = build_r0_r1()
    part_b(r0, r1)
    part_c(r0, r1)
    print("=" * 78)
    print("A9 PATCH COMPLETE.  The degenerate locus {A Sz = B} is closed three ways:")
    print("  (a) the forcing chain A=C ==> m0^2=1 ==> Sz=m0 Sw (and P == Qz) is DERIVED;")
    print("  (b) no Res-branch lies in {Sw=m0 Sz}, so the good-locus closure argument covers all")
    print("      of {Res=0} (T rational, vanishing on a dense open subset of each branch);")
    print("  (c) INDEPENDENTLY, r0 = r1 = 0 on {m0^2=1, Sz=m0 Sw} exactly, so T = 0 there directly.")
    print(f"[TOTAL] {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
