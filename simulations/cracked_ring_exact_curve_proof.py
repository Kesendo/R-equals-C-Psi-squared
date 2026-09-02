"""Symbolic gates for docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md (F160), exact in sympy, no floats.

Every trigonometric identity is decided by ONE exact route: expand the multiple angles, write sin k -> s and
cos k -> c, clear the (monomial) denominators, and reduce the numerator modulo s^2 + c^2 - 1; the identity holds
iff the remainder is the zero polynomial in (s, c, u). A first version of P3 and P4 fell back to a 40-digit
evaluation at one rational point where sympy's simplifier did not close the form; that was a tolerance where an
exact route exists, and a reviewer said so.

The cracked ring: N sites, bonds 1 (units of J), the wrap bond N-1 <-> 0 carrying u. Everything below is an
identity in the symbols x, u (and k, delta), decided by sympy's exact simplification to 0.

  [P1] Theorem A: det(x I - H) == U_N(x/2) - u^2 U_{N-2}(x/2) - 2u, the symbolic determinant of the N x N pencil
       against the Chebyshev recursion, N = 3..9, as polynomials in (x, u).
  [P2] the Cassini identity U_{n}^2 - U_{n+1} U_{n-1} = 1 in the 2cos normalization (the step of the second route,
       the matrix determinant lemma), n = 1..12.
  [P3] Corollary B: P(2 cos k) * sin k == G(k) in both printed forms, N = 3..9, exact (reduction modulo s^2 + c^2 - 1).
  [P4] Theorem D, the join: 1/t == (1+u^2)/(2u) + i(1-u^2) cot k/(2u) with e^{+-ik} = c +- i s, and
       Re[e^{-iNk}/t(k)] - 1 == G(k)/(2u sin k) with Re[(cos Nk - i sin Nk)(A + iB)] = A cos Nk + B sin Nk for real A, B;
       N = 3..9, exact by the same reduction.
  [P5] Theorem E, the split's next order: the series of the pair's two roots to O(delta^2) and
       Delta E/(4 delta/N) == 1 + delta (1/2 - 1/(N s^2)) + O(delta^2), symbolic in (N, s, c) with c^2 = 1 - s^2.
  [P6] Theorem F, the band-edge factor forms: P(2) and P(-2) against the two linear factors, N = 3..12 (both parities),
       and the odd-N threshold u = (N+1)/(N-1) as the exact root of P(-2).
  [P7] Theorem G, the chain-end velocity, exact: dP/du|0 = -2; dP/dx|0 at x = 2 cos theta is the displayed
       trigonometric form (reduction); substituting sin((N+1)theta) = 0, cos((N+1)theta) = (-1)^k gives
       (-1)^(k+1) (4/(N+1)) sin^2 theta from both 2/P' and 2 psi_k(0) psi_k(N-1); N = 3..9.
  [P8] Corollary B's simplicity clause: G == 2AB symbolic in N, the two reflection sectors through the crack, and
       the Bezout identity (cos a + u cos b)A + (sin a - u sin b)B == 1 - u^2 symbolic in N, so a common zero of A
       and B forces u^2 = 1.
  [P9] the prefactor relations behind the clause: the ring folded along the reflection through the crack gives two
       Jacobi blocks with characteristic polynomials chi_e, chi_o, and A = cos(k/2) chi_e(2 cos k), B = sin(k/2)
       chi_o(2 cos k) at even N, A = chi_e/2, B = sin k chi_o at odd N, with chi_e chi_o = P; N = 3..9, u symbolic.

Run: python simulations/cracked_ring_exact_curve_proof.py     (a few seconds)
"""
import sys

import sympy as sp

x, u, k, d = sp.symbols('x u k delta', real=True)
FAIL = []


def gate(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def cheb_u_half(n):
    """U_n(x/2) by the recursion p_0 = 1, p_1 = x, p_n = x p_{n-1} - p_{n-2}."""
    p0, p1 = sp.Integer(1), x
    if n == 0:
        return p0
    for _ in range(n - 1):
        p0, p1 = p1, sp.expand(x * p1 - p0)
    return p1


def road_poly(n):
    return sp.expand(cheb_u_half(n) - u**2 * cheb_u_half(n - 2) - 2 * u)


def pencil(n):
    m = sp.zeros(n, n)
    for i in range(n):
        m[i, i] = x
    for i in range(n - 1):
        m[i, i + 1] = -1
        m[i + 1, i] = -1
    m[0, n - 1] = -u
    m[n - 1, 0] = -u
    return m


def curve_first(n):
    return (1 - u**2) * sp.sin(n * k) * sp.cos(k) + ((1 + u**2) * sp.cos(n * k) - 2 * u) * sp.sin(k)


def curve_second(n):
    return sp.sin((n + 1) * k) - u**2 * sp.sin((n - 1) * k) - 2 * u * sp.sin(k)


s_, c_ = sp.symbols('s c')


def trig_zero(expr):
    """Exact decision of a trigonometric identity in k with polynomial coefficients in u: expand the multiple
    angles, map sin k -> s and cos k -> c (tan and cot as their quotients), clear the monomial denominators, and
    reduce the numerator modulo s^2 + c^2 - 1. True iff the remainder is the zero polynomial."""
    e = sp.expand_trig(expr)
    e = e.subs({sp.sin(k): s_, sp.cos(k): c_, sp.tan(k): s_ / c_, sp.cot(k): c_ / s_})
    e = sp.together(sp.expand(e))
    num, den = sp.fraction(e)
    assert set(den.free_symbols) <= {s_, c_, u}, den      # denominators are powers of sin k, cos k, u only
    _, r = sp.reduced(sp.expand(num), [s_**2 + c_**2 - 1], s_, c_, u)
    return sp.expand(r) == 0


# ---------------------------------------------------------------- P1, P2
ok = True
for n in range(3, 10):
    det = sp.expand(pencil(n).det(method='berkowitz'))
    if sp.expand(det - road_poly(n)) != 0:
        ok = False
gate("P1 Theorem A: det(x I - H) == U_N(x/2) - u^2 U_{N-2}(x/2) - 2u as polynomials in (x, u), N = 3..9", ok)

ok = all(sp.expand(cheb_u_half(n)**2 - cheb_u_half(n + 1) * cheb_u_half(n - 1) - 1) == 0 for n in range(1, 13))
gate("P2 Cassini: U_n^2 - U_{n+1} U_{n-1} == 1 in the 2cos normalization, n = 1..12", ok)

# ---------------------------------------------------------------- P3
ok = True
for n in range(3, 10):
    lhs = road_poly(n).subs(x, 2 * sp.cos(k)) * sp.sin(k)
    for rhs in (curve_first(n), curve_second(n)):
        if not trig_zero(lhs - rhs):
            ok = False
gate("P3 Corollary B: P(2 cos k) sin k == G(k), both printed forms, N = 3..9, exact (reduction modulo s^2 + c^2 - 1)", ok)

# ---------------------------------------------------------------- P4
# 1/t with e^{+-ik} = cos k +- i sin k: 1/t = (e^{-ik} - u^2 e^{ik}) / (-2iu sin k) against the displayed form.
denom_t = (c_ - sp.I * s_) - u**2 * (c_ + sp.I * s_)
inv_t_target = (1 + u**2) / (2 * u) + sp.I * (1 - u**2) * c_ / (2 * u * s_)
e = sp.together(sp.expand(denom_t / (-2 * sp.I * u * s_) - inv_t_target))
num, den = sp.fraction(e)
_, r = sp.reduced(sp.expand(num), [s_**2 + c_**2 - 1], s_, c_, u)
ok = sp.expand(r) == 0
# Re[e^{-iNk}/t]: with 1/t = A + iB (A, B real), Re[(cos Nk - i sin Nk)(A + iB)] = A cos Nk + B sin Nk, so the
# left side is cos(Nk)(1+u^2)/(2u) + sin(Nk)(1-u^2) cot k/(2u); minus 1 minus G/(2u sin k) must vanish identically.
A_, B_ = sp.symbols('A B', real=True)
ok = ok and sp.simplify(sp.re((sp.cos(k) - sp.I * sp.sin(k)) * (A_ + sp.I * B_)) - (A_ * sp.cos(k) + B_ * sp.sin(k))) == 0
for n in range(3, 10):
    re_f = sp.cos(n * k) * (1 + u**2) / (2 * u) + sp.sin(n * k) * (1 - u**2) * sp.cot(k) / (2 * u)
    if not trig_zero(re_f - 1 - curve_first(n) / (2 * u * sp.sin(k))):
        ok = False
gate("P4 Theorem D: 1/t == (1+u^2)/(2u) + i(1-u^2) cot k/(2u), and Re[e^{-iNk}/t] - 1 == G/(2u sin k), N = 3..9, exact (reduction modulo s^2 + c^2 - 1)", ok)

# ---------------------------------------------------------------- P5
s, Nn = sp.symbols('s N', positive=True)
c = sp.symbols('c', real=True)          # cos k_m takes both signs on the pairs; only s = sin k_m > 0 is assumed
q = x / Nn
cosk = c * sp.cos(q) - s * sp.sin(q)
sink = s * sp.cos(q) + c * sp.sin(q)
uu = 1 - d
G = (1 - uu**2) * sp.sin(x) * cosk + ((1 + uu**2) * sp.cos(x) - 2 * uu) * sink
a1, a2 = sp.symbols('a1 a2')
Gs = sp.expand(sp.series(G.subs(x, a1 * d + a2 * d**2), d, 0, 4).removeO())
e2 = sp.factor(Gs.coeff(d, 2))
e3 = sp.expand(Gs.coeff(d, 3))
ok = sp.expand(e2 - (-a1**2 * s + 2 * a1 * c + s)) == 0
roots = sp.solve(e2, a1)
sols = []
for r1 in roots:
    r2 = sp.solve(e3.subs(a1, r1), a2)
    sols.append((r1, r2[0]))
Es = []
for r1, r2 in sols:
    xx = r1 * d + r2 * d**2
    Es.append(sp.expand(sp.series(2 * (c * sp.cos(xx / Nn) - s * sp.sin(xx / Nn)), d, 0, 3).removeO()))
dE = sp.expand((Es[0] - Es[1]).subs(c**2, 1 - s**2))
if sp.simplify(dE.coeff(d, 1)) < 0:
    dE = -dE
ratio = sp.expand(sp.simplify(dE / (4 * d / Nn)).subs(c**2, 1 - s**2))
ok = ok and sp.simplify(ratio - (1 + d * (sp.Rational(1, 2) - 1 / (Nn * s**2)))) == 0
gate("P5 Theorem E: delta^2 coefficient is -a1^2 s + 2 a1 c + s (roots (c +- 1)/s), and Delta E/(4 delta/N) == 1 + delta(1/2 - 1/(N s^2)) + O(delta^2)",
     ok, f"a1 = {[sp.simplify(r.subs(sp.sqrt(c**2 + s**2), 1)) for r in roots]} (with c^2 + s^2 = 1)")

# ---------------------------------------------------------------- P6
ok = True
for n in range(3, 13):
    p = road_poly(n)
    top = sp.expand(p.subs(x, 2))
    bot = sp.expand(p.subs(x, -2))
    top_f = sp.expand(-((n - 1) * u + (n + 1)) * (u - 1))
    bot_f = sp.expand(((n - 1) * u - (n + 1)) * (u + 1)) if n % 2 == 1 else sp.expand(-((n - 1) * u + (n + 1)) * (u - 1))
    if top - top_f != 0 or bot - bot_f != 0:
        ok = False
    if n % 2 == 1:
        rts = sp.solve(bot, u)
        if sp.Rational(n + 1, n - 1) not in rts:
            ok = False
gate("P6 Theorem F: P(+2) and P(-2) equal the two linear factors, N = 3..12; the odd-N threshold (N+1)/(N-1) is the exact positive root of P(-2)", ok)

# ---------------------------------------------------------------- P7, exact
# (1) dP/du at u = 0 is -2. (2) dP/dx at u = 0 is (d/dx) U_N(x/2), which at x = 2 cos theta equals
#     [(N+1) cos((N+1)theta) sin theta - sin((N+1)theta) cos theta] / (-2 sin^3 theta): a trigonometric identity
#     in theta, decided by the reduction. (3) At theta_k = k pi/(N+1), sin((N+1)theta_k) = 0 and
#     cos((N+1)theta_k) = (-1)^k, so 2/P' = (-1)^(k+1) 4 sin^2(theta_k)/(N+1): pure substitution, done symbolically
#     with S = sin((N+1)theta) -> 0 and C = cos((N+1)theta) -> (-1)^k, k an integer symbol.
# (4) The eigenvector route: 2 psi_k(0) psi_k(N-1) = (4/(N+1)) sin theta sin(N theta), and
#     sin(N theta) = sin((N+1)theta) cos theta - cos((N+1)theta) sin theta (an identity) -> -(-1)^k sin theta_k.
kk = sp.symbols('k_label', integer=True)
S_, C_ = sp.symbols('S C', real=True)
ok = True
for n in range(3, 10):
    p = road_poly(n)
    ok = ok and sp.expand(sp.diff(p, u).subs(u, 0) + 2) == 0
    dpdx = sp.diff(p, x).subs(u, 0).subs(x, 2 * sp.cos(k))
    target = ((n + 1) * sp.cos((n + 1) * k) * sp.sin(k) - sp.sin((n + 1) * k) * sp.cos(k)) / (-2 * sp.sin(k)**3)
    ok = ok and trig_zero(dpdx - target)
    # (3) the substitution at theta_k, symbolic in k_label
    dpdx_k = ((n + 1) * C_ * sp.sin(k) - S_ * sp.cos(k)) / (-2 * sp.sin(k)**3)
    vel = (2 / dpdx_k).subs({S_: 0, C_: (-1)**kk})
    closed = (-1)**(kk + 1) * sp.Rational(4, n + 1) * sp.sin(k)**2
    ok = ok and sp.simplify(vel - closed) == 0
    # (4) the eigenvector route
    sinN = sp.sin((n + 1) * k) * sp.cos(k) - sp.cos((n + 1) * k) * sp.sin(k)
    ok = ok and trig_zero(sp.sin(n * k) - sinN)
    eig = sp.Rational(4, n + 1) * sp.sin(k) * (S_ * sp.cos(k) - C_ * sp.sin(k))
    ok = ok and sp.simplify(eig.subs({S_: 0, C_: (-1)**kk}) - closed) == 0
gate("P7 Theorem G, exact: dP/du|0 = -2; dP/dx|0 at x = 2 cos theta is the displayed trig form (reduction); at sin((N+1)theta) = 0, cos((N+1)theta) = (-1)^k both 2/P' and 2 psi_k(0) psi_k(N-1) equal (-1)^(k+1) (4/(N+1)) sin^2 theta; N = 3..9", ok)

# ---------------------------------------------------------------- P8 the simplicity clause (found by a review lens, 2026-09-02)
Ns = sp.symbols('N', positive=True)
aa, bb = k * (Ns + 1) / 2, k * (Ns - 1) / 2
Afac = sp.cos(aa) - u * sp.cos(bb)
Bfac = sp.sin(aa) + u * sp.sin(bb)
Gsym = sp.sin((Ns + 1) * k) - u**2 * sp.sin((Ns - 1) * k) - 2 * u * sp.sin(k)
ok = sp.simplify(sp.expand_trig(sp.expand(2 * Afac * Bfac - Gsym))) == 0
# a common zero of A and B is impossible unless u^2 = 1, by the Bezout identity
#   (cos a + u cos b) A + (sin a - u sin b) B = 1 - u^2,  symbolic in N.
# (A first version of this half compared two Pythagorean identities typed by hand, f(x) == f(x); a reviewer caught it.)
bezout = (sp.cos(aa) + u * sp.cos(bb)) * Afac + (sp.sin(aa) - u * sp.sin(bb)) * Bfac - (1 - u**2)
ok = ok and sp.simplify(sp.expand(bezout)) == 0
gate("P8 Corollary B's simplicity clause: G == 2AB symbolic in N (A = cos(k(N+1)/2) - u cos(k(N-1)/2), B = sin(k(N+1)/2) + u sin(k(N-1)/2)), and (cos a + u cos b)A + (sin a - u sin b)B == 1 - u^2 symbolic in N, so A = B = 0 forces u^2 = 1", ok)

# ---------------------------------------------------------------- P9 the prefactor relations (the step a review lens found missing)
# Fold the ring along the reflection through the crack bond. Even N = 2M: R-even basis (e_j + e_{N-1-j})/sqrt2,
# j = 0..M-1, giving a tridiagonal block with off-diagonals 1, diagonal u at j = 0 (the wrap bond joins a site to
# its own partner) and +1 at j = M-1 (the middle bond does the same); the R-odd block has -u and -1 there. Odd
# N = 2M+1: the middle site M is R-fixed and joins the even block through a sqrt2 (no fold-end diagonal); the odd
# block is M x M with -u at j = 0 and nothing at the fold end. With chi_e, chi_o their characteristic polynomials:
#   even N:  A(k) = cos(k/2) chi_e(2 cos k),   B(k) = sin(k/2) chi_o(2 cos k)
#   odd N:   A(k) = (1/2) chi_e(2 cos k),      B(k) = sin(k)  chi_o(2 cos k)
# decided exactly in the half angle phi = k/2 (so that cos(k/2), sin(k/2), cos k, sin(Nk/2) are all integer multiples).
phi = sp.symbols('phi', real=True)


def trig_zero_in(expr, var):
    e = sp.expand_trig(expr)
    e = e.subs({sp.sin(var): s_, sp.cos(var): c_})
    e = sp.together(sp.expand(e))
    num, den = sp.fraction(e)
    assert set(den.free_symbols) <= {s_, c_, u}, den
    _, r = sp.reduced(sp.expand(num), [s_**2 + c_**2 - 1], s_, c_, u)
    return sp.expand(r) == 0


def folded_blocks(n):
    M = n // 2
    if n % 2 == 0:
        e = sp.zeros(M, M); o = sp.zeros(M, M)
        for j in range(M - 1):
            e[j, j + 1] = e[j + 1, j] = 1
            o[j, j + 1] = o[j + 1, j] = 1
        e[0, 0] = u; o[0, 0] = -u
        e[M - 1, M - 1] = 1; o[M - 1, M - 1] = -1
    else:
        e = sp.zeros(M + 1, M + 1); o = sp.zeros(M, M)
        for j in range(M - 1):
            e[j, j + 1] = e[j + 1, j] = 1
            o[j, j + 1] = o[j + 1, j] = 1
        e[M - 1, M] = e[M, M - 1] = sp.sqrt(2)
        e[0, 0] = u; o[0, 0] = -u
    return e, o


ok = True
for n in range(3, 10):
    e_blk, o_blk = folded_blocks(n)
    # charpoly() returns a PurePoly over its own generator; carry it back onto THIS script's x (the
    # generator prints as x but is a different Symbol, and a difference of the two never cancels)
    cp_e, cp_o = e_blk.charpoly(), o_blk.charpoly()
    chi_e = cp_e.as_expr().subs(cp_e.gens[0], x)
    chi_o = cp_o.as_expr().subs(cp_o.gens[0], x)
    kk_ = 2 * phi
    A_ = sp.cos(kk_ * (n + 1) / 2) - u * sp.cos(kk_ * (n - 1) / 2)
    B_ = sp.sin(kk_ * (n + 1) / 2) + u * sp.sin(kk_ * (n - 1) / 2)
    xk = 2 * sp.cos(kk_)
    if n % 2 == 0:
        okA = trig_zero_in(A_ - sp.cos(phi) * chi_e.subs(x, xk), phi)
        okB = trig_zero_in(B_ - sp.sin(phi) * chi_o.subs(x, xk), phi)
    else:
        okA = trig_zero_in(A_ - sp.Rational(1, 2) * chi_e.subs(x, xk), phi)
        okB = trig_zero_in(B_ - sp.sin(kk_) * chi_o.subs(x, xk), phi)
    # and the fold is consistent with the whole: chi_e * chi_o == P (the blocks are a similarity of H)
    okP = sp.expand(chi_e * chi_o - road_poly(n)) == 0
    if not (okA and okB and okP):
        ok = False
        print("   P9 fails at N =", n, okA, okB, okP)
gate("P9 the prefactor relations of the simplicity clause: A == cos(k/2) chi_e(2 cos k) and B == sin(k/2) chi_o(2 cos k) at even N, A == chi_e/2 and B == sin k chi_o at odd N, with chi_e chi_o == P, N = 3..9, u symbolic, exact in the half angle", ok)

print("ALL GATES PASS" if not FAIL else f"FAILED: {FAIL}")
sys.exit(0 if not FAIL else 1)
