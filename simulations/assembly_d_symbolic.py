"""The assembly (D), made symbolic: (U+ - U-) * (n/2)^3 = cross_form(a; b), given D1-D4.

THE STATEMENT. The object-identification bridge of the twinning proof
(experiments/F89_SEED_EXISTENCE_REDUCTION.md; registered as F127 in docs/ANALYTICAL_FORMULAS.md):
the discrete Gram difference of the -6-rung Slater lifts equals the six-angle cross form,

    (U+ - U-) * (n/2)^3  =  cross_form(a; b),   a_i = k_i*theta, b_j = l_j*theta, theta = pi/n,

for mode-disjoint triples tau = {k_1<k_2<k_3}, sigma = {l_1<l_2<l_3} with eps(tau,sigma) = -1
(at eps = +1 both sides are 0: D4's even branch kills every term, re-deriving Lemma 4). Until
2026-07-14 this equation was checked numerically (16200 pairs, 3e-13, the gate_angle_forms
endpoint in cross_triple_orthogonality.py); this script makes the assembly SYMBOLIC, so the
identification is derived, not measured.

THE CHAIN (each numbered step is D1-D4 of cross_triple_orthogonality.py, already proved and
gated there, or one of the elementary steps S1-S3 verified here):

    (U+ - U-) * (n/2)^3
      = sum_{b,c} sgn(c-b) G_tau(b,c) G_sigma(b,c)
            [definitional: Ggrid = raw/||D||, ||D_tau||^2 = (n/2)^3 proved by sine
             orthogonality + Cauchy-Binet; the (n/2)^3 clears the two norms]
      = sum_{i,j} (-1)^{i+j} sum_b M_i(b) M_j(b) * sum_c sgn(c-b) u_{k_i}(c) u_{l_j}(c)
            [D1, Laplace along the c column; the sign is (-1)^{i+1}(-1)^{j+1} = (-1)^{i+j}]
      = sum_{i,j} (-1)^{i+j} sum_b M_i M_j * (1/2)[Theta_{P-}(b) - Theta_{P+}(b)]
            [S2 product-to-sum sin*sin = (cos(dif) - cos(sum))/2, then D3 per cosine;
             P+ = k_i + l_j, P- = k_i - l_j; preconditions for the cotangent branch:
             0 < P+ <= 2n-2 and 0 < |P-| < n for mode-disjoint pairs, so neither is
             0 mod 2n and the Theta_0 branch never fires]
      = sum_{i,j} (-1)^{i+j} (1/2) sum_b M_i M_j [cot(mu+/2) s_{P+}(b) - cot(mu-/2) s_{P-}(b)]
            [S2 parity: P+ - P- = 2*l_j is even, so the two Theta constants
             [1-(-1)^P]/2 are equal and cancel in the difference; mu+- = P+- * theta]
      = sum_{i,j} (-1)^{i+j} (1/2) [cot(mu+/2) Xs(mu+) - cot(mu-/2) Xs(mu-)]
            [expand M_i M_j into the four sine products (D2's closed form) and apply D4
             per product. Parity uniformity: every x + y + P == k(tau) + k(sigma) mod 2,
             so at eps = -1 ALL sixteen D4 calls take the cotangent branch (gated below).
             Negative P: D4's both sides are ODD in P (cot odd), so the identity extends
             from the gated positive range (gated symbolically below).]
      = cross_form(a; b)
            [S1: the D2 half-angle coefficients ARE _pieces' alpha/beta
             (sin(q t) - sin(p t) = 2 cos((p+q)t/2) sin((q-p)t/2), and the beta twin);
             with them the assembled expression is term-for-term cross_form's formula.
             S3 verifies the equality of the two six-angle expressions symbolically.]

WHAT THIS CLOSES. With D1-D4 proved in cross_triple_orthogonality.py and S1-S3 here, the
assembly (D) is a theorem; combined with F127 (cross_form == 0 on V, the grid+CRT wall) the
full-spectrum twinning is proof grade over Q(i) WITHOUT the modulo-(D) qualifier. The remaining
caveat on the pair is code trust (single implementations), no longer a missing mathematical step.

GATES (default run, ~40 s). The four summation identities D1-D4 are each proved
SYMBOLICALLY here (no sampling; the committed numeric gates in cross_triple_orthogonality.py
become corroboration), the bookkeeping lemmas are symbolic, and the endpoint is pinned to the
committed cross_form. The composition of the steps is linearity of finite sums (reordering a
finite double sum), which is not a gateable step; the exact N=9 anchor (G12) closes end-to-end.
  G1   S1 coefficient identities, symbolic.
  G2   S2 product-to-sum + the Theta-constant parity cancellation, symbolic.
  G3   D4 oddness extension RHS(-P) = -RHS(P), symbolic.
  G4   parity uniformity: every (x, y, P) of the expansion sits in the eps branch. SYMBOLIC:
       the exponent cancellation (-1)^{(k(tau)-k_i)+(k(sigma)-l_j)+(k_i +- l_j)} = (-1)^{k(tau)+k(sigma)}
       with integer symbols, plus the exhaustive integer sweep (all mode-disjoint pairs,
       N = 9, 11, 12) as corroboration.
  G5   D3 preconditions: P+ and P- never 0 mod 2n for mode-disjoint pairs. The one-line bound:
       distinct modes in 1..N = n-1 give 3 <= P+ <= 2n-3 and 0 < |P-| <= n-2, so neither is
       0 mod 2n; gated by the same exhaustive sweep.
  G6   S3: the half-angle and sin-difference COEFFICIENT CONVENTIONS agree across the whole
       assembled six-angle sum (free angles, exact simplify). This is S1 lifted through the
       sum; the discrete side enters via G7/G12, the committed cross_form via G11.
  G7   the definitional normalization reading pinned against the committed engine:
       (Upm difference)*(n/2)^3 == sum sgn * G_raw * G_raw at N = 9 (float, 1e-12).
  G8   D1 SYMBOLIC: the cofactor expansion along the c column, generic 3x3.
  G9   D2 SYMBOLIC: the telescoping closed form of M_pq, free variables.
  G10  D3 + D4 SYMBOLIC: the signed cosine sum (both parities of P) and the odd-Q full-range
       sine sum, both from the geometric series with z^n = (-1)^Q as the ONLY integrality;
       plus the triple-product expansion and the match to the committed D4 cotangent form.
  G11  transcription pin: this file's _form(halfangle=False) == the committed
       cross_triple_orthogonality.cross_form on a deterministic angle grid (float, 1e-8).
--slow adds:
  G12 the end-to-end exact-arithmetic anchor: N = 9, pair ({1,2,3},{4,5,6}), LHS from exact
      Slater sums (algebraic numbers), RHS = cross_form exact, difference == 0 by sympy
      simplify. ~90 minutes of exact algebra; the anchor that no float enters anywhere.

Run:  python simulations/assembly_d_symbolic.py [--slow]
"""
from __future__ import annotations

import itertools
import sys
import time

import sympy as sp

sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

T0 = time.time()


def _pass(name):
    print(f"PASS {name}   [{time.time() - T0:.1f}s]", flush=True)


COT = lambda z: sp.cos(z) / sp.sin(z)


def gate1_coefficients():
    p, q, t = sp.symbols('p q t', positive=True)
    s1a = sp.sin(q*t) - sp.sin(p*t) - 2*sp.cos((p+q)*t/2)*sp.sin((q-p)*t/2)
    s1b = -(sp.sin(q*t) + sp.sin(p*t)) + 2*sp.sin((p+q)*t/2)*sp.cos((q-p)*t/2)
    assert sp.simplify(sp.expand_trig(s1a)) == 0
    assert sp.simplify(sp.expand_trig(s1b)) == 0
    _pass("G1 S1: the D2 half-angle coefficients == _pieces alpha/beta (symbolic)")


def gate2_product_to_sum_and_parity():
    x, y = sp.symbols('X Y', real=True)
    assert sp.simplify(sp.sin(x)*sp.sin(y) - (sp.cos(x-y) - sp.cos(x+y))/2) == 0
    k, l = sp.symbols('k l', integer=True)
    assert sp.simplify((-1)**(k+l) - (-1)**(k-l)) == 0   # Theta constants equal => cancel
    _pass("G2 S2: product-to-sum + Theta-constant parity cancellation (symbolic)")


def gate3_d4_oddness():
    p, x, y, th = sp.symbols('P x y theta', real=True)
    rhs = sp.Rational(1, 4)*(COT((p+x-y)*th/2) + COT((p-x+y)*th/2)
                             - COT((p+x+y)*th/2) - COT((p-x-y)*th/2))
    assert sp.simplify(rhs.subs(p, -p) + rhs) == 0
    _pass("G3 D4 extends to negative P: RHS(-P) = -RHS(P) (cot odd, symbolic)")


def _disjoint_pairs(n_sites):
    for tau in itertools.combinations(range(1, n_sites + 1), 3):
        for sig in itertools.combinations(range(1, n_sites + 1), 3):
            if tau < sig and not (set(tau) & set(sig)):
                yield tau, sig


def gate4_parity_uniformity(ns=(9, 11, 12)):
    # the symbolic half: the exponent cancellation with integer symbols, both P branches
    kt, ks, ki, lj = sp.symbols('ktau ksigma ki lj', integer=True)
    for p_exp in (ki + lj, ki - lj):
        exponent = (kt - ki) + (ks - lj) + p_exp
        assert sp.simplify((-1)**exponent - (-1)**(kt + ks)) == 0
    # the exhaustive integer sweep as corroboration
    for n_sites in ns:
        for tau, sig in _disjoint_pairs(n_sites):
            eps_par = (sum(tau) + sum(sig)) % 2
            for i in range(3):
                rest_t = [k for t, k in enumerate(tau) if t != i]
                xs = (rest_t[0] + rest_t[1], rest_t[1] - rest_t[0])
                for j in range(3):
                    rest_s = [l for t2, l in enumerate(sig) if t2 != j]
                    ys = (rest_s[0] + rest_s[1], rest_s[1] - rest_s[0])
                    for x in xs:
                        for y in ys:
                            for p in (tau[i] + sig[j], tau[i] - sig[j]):
                                assert (x + y + p) % 2 == eps_par
    _pass(f"G4 parity uniformity: SYMBOLIC exponent cancellation + exhaustive sweep (N={ns})")


def gate5_d3_preconditions(ns=(9, 11, 12)):
    for n_sites in ns:
        n = n_sites + 1
        for tau, sig in _disjoint_pairs(n_sites):
            for i in range(3):
                for j in range(3):
                    p_plus, p_minus = tau[i] + sig[j], tau[i] - sig[j]
                    assert 0 < p_plus <= 2*n - 2 and p_plus % (2*n) != 0
                    assert p_minus != 0 and abs(p_minus) < n
    _pass(f"G5 D3 preconditions: P+ and P- never 0 mod 2n on mode-disjoint pairs (N={ns})")


def _pieces_sym(ang, halfangle):
    out = []
    for i in range(3):
        j2, l2 = [x for x in range(3) if x != i]
        aj, al = ang[j2], ang[l2]
        if halfangle:
            alpha = 2*sp.cos((aj+al)/2)*sp.sin((al-aj)/2)
            beta = -2*sp.sin((aj+al)/2)*sp.cos((al-aj)/2)
        else:
            alpha = sp.sin(al) - sp.sin(aj)
            beta = -(sp.sin(al) + sp.sin(aj))
        out.append(dict(psum=aj+al, pdif=al-aj, alpha=alpha, beta=beta))
    return out


def _xh(mu, xi, up):
    return sp.Rational(1, 4)*(COT((mu+xi-up)/2) + COT((mu-xi+up)/2)
                              - COT((mu+xi+up)/2) - COT((mu-xi-up)/2))


def _form(ang_tau, ang_sigma, halfangle):
    pa, pb = _pieces_sym(ang_tau, halfangle), _pieces_sym(ang_sigma, halfangle)
    tot = 0
    for i in range(3):
        pi = pa[i]
        for j in range(3):
            pj = pb[j]
            def xs(mu):
                return (pi['alpha']*pj['alpha']*_xh(mu, pi['psum'], pj['psum'])
                        + pi['alpha']*pj['beta']*_xh(mu, pi['psum'], pj['pdif'])
                        + pi['beta']*pj['alpha']*_xh(mu, pi['pdif'], pj['psum'])
                        + pi['beta']*pj['beta']*_xh(mu, pi['pdif'], pj['pdif']))
            mup, mum = ang_tau[i] + ang_sigma[j], ang_tau[i] - ang_sigma[j]
            tot += (-1)**(i+j)*(COT(mup/2)*xs(mup) - COT(mum/2)*xs(mum))
    return tot/2


def gate6_structural_identity():
    a = sp.symbols('a1 a2 a3', real=True)
    b = sp.symbols('b1 b2 b3', real=True)
    diff = _form(a, b, halfangle=True) - _form(a, b, halfangle=False)
    d = sp.cancel(sp.together(diff.rewrite(sp.exp).expand()))
    assert sp.simplify(d) == 0, "coefficient conventions disagree across the assembled sum"
    _pass("G6 S3: half-angle == sin-difference coefficient convention across the assembled sum (symbolic)")


def gate7_normalization_reading(n_sites=9, max_pairs=40):
    import cross_triple_orthogonality as ct
    n = n_sites + 1
    u = ct.umat(n_sites)
    checked = 0
    worst = 0.0
    for tau, sig in _disjoint_pairs(n_sites):
        if (sum(tau) + sum(sig)) % 2 == 0:
            continue
        g_t, g_s = ct.Ggrid(u, tau, n_sites), ct.Ggrid(u, sig, n_sites)
        up, um = ct.Upm(g_t, g_s, n_sites)
        lhs = (up - um) * (n / 2) ** 3
        rhs = sum((1 if c > b else -1)
                  * (ct.slater(u, tau, b-1, b, c) + ct.slater(u, tau, b, b+1, c))
                  * (ct.slater(u, sig, b-1, b, c) + ct.slater(u, sig, b, b+1, c))
                  for b in range(1, n_sites+1) for c in range(1, n_sites+1) if c != b)
        worst = max(worst, abs(lhs - rhs))
        checked += 1
        if checked >= max_pairs:
            break
    assert worst < 1e-12, f"normalization reading off: {worst:.2e}"
    _pass(f"G7 definitional normalization vs the committed engine: {worst:.1e} over {checked} pairs (N={n_sites})")


def gate8_d1_symbolic():
    m = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'm{i}{j}'))
    cof = sum((-1)**i * m[i, 2] * m.minor_submatrix(i, 2).det() for i in range(3))
    assert sp.simplify(m.det() - cof) == 0
    _pass("G8 D1 symbolic: cofactor expansion along the c column (generic 3x3)")


def gate9_d2_symbolic():
    p, q, b, t = sp.symbols('p q b t', real=True)
    u = lambda k, x: sp.sin(k*x*t)
    nab = lambda y: u(p, y)*u(q, y+1) - u(q, y)*u(p, y+1)
    lhs = nab(b-1) + nab(b)
    rhs = (2*sp.cos((p+q)*t/2)*sp.sin((q-p)*t/2)*sp.sin((p+q)*b*t)
           - 2*sp.sin((p+q)*t/2)*sp.cos((q-p)*t/2)*sp.sin((q-p)*b*t))
    assert sp.simplify(sp.expand_trig(lhs - rhs)) == 0
    _pass("G9 D2 symbolic: the telescoping closed form of M_pq (free variables)")


def gate10_d3_d4_symbolic():
    # D3: sgn-sum of cosines. Variables phi = P*theta, psi = P*b*theta free; the ONLY
    # integrality is z^n = (-1)^P =: s. The geometric series does the c-sums.
    phi, psi = sp.symbols('phi psi', real=True)
    z, zb = sp.exp(sp.I*phi), sp.exp(sp.I*psi)
    for s in (1, -1):
        full = (s - z)/(z - 1)                    # sum_{c=1}^{N} z^c   (z^n = s)
        upto_b = z*(zb - 1)/(z - 1)               # sum_{c=1}^{b} z^c   (z^b = zb)
        sgn_sum = sp.re(sp.expand_complex(full - 2*upto_b)) + sp.cos(psi)
        target = sp.Rational(1, 2)*(1 - s) - (sp.cos(phi/2)/sp.sin(phi/2))*sp.sin(psi)
        assert sp.simplify(sp.expand_complex(sgn_sum - target)) == 0
    # D4: for ODD Q the full-range sine sum collapses to a cotangent
    al = sp.symbols('alpha', real=True)
    zq = sp.exp(sp.I*al)
    assert sp.simplify(sp.expand_complex((-1 - zq)/(zq - 1) - sp.I*sp.cos(al/2)/sp.sin(al/2))) == 0
    # the triple-product expansion feeding it
    aa, bb, cc = sp.symbols('A B C', real=True)
    expd = sp.Rational(1, 4)*(sp.sin(aa+bb-cc) + sp.sin(cc+aa-bb) + sp.sin(bb+cc-aa) - sp.sin(aa+bb+cc))
    assert sp.simplify(sp.expand_trig(sp.sin(aa)*sp.sin(bb)*sp.sin(cc) - expd)) == 0
    # and the assembled four-cotangent form equals the committed D4 statement
    x, y, p, th = sp.symbols('x y P theta', real=True)
    mine = sp.Rational(1, 4)*(COT((x+y-p)*th/2) + COT((p+x-y)*th/2)
                              + COT((p-x+y)*th/2) - COT((x+y+p)*th/2))
    committed = sp.Rational(1, 4)*(COT((p+x-y)*th/2) + COT((p-x+y)*th/2)
                                   - COT((p+x+y)*th/2) - COT((p-x-y)*th/2))
    assert sp.simplify(mine - committed) == 0
    _pass("G10 D3+D4 symbolic: the geometric series with z^n = (-1)^Q as the only integrality")


def gate11_transcription_pin():
    import cross_triple_orthogonality as ct
    a_sym = sp.symbols('a1 a2 a3', real=True)
    b_sym = sp.symbols('b1 b2 b3', real=True)
    mine = _form(a_sym, b_sym, halfangle=False)
    worst = 0.0
    for da in (0.31, 0.57, 0.83):
        ang_a = [0.4, 0.4 + da, 0.4 + 2.1*da]
        for db in (0.29, 0.61):
            ang_b = [0.7, 0.7 + db, 0.7 + 1.7*db]
            committed = ct.cross_form(ang_a, ang_b)
            val = float(mine.subs(dict(zip(list(a_sym) + list(b_sym), ang_a + ang_b))).evalf())
            worst = max(worst, abs(val - committed))
    assert worst < 1e-8, f"transcription drift {worst:.2e}"
    _pass(f"G11 transcription pin vs the committed cross_form: {worst:.1e} on a deterministic angle grid")


def gate12_exact_anchor():
    n_sites = 9
    n = n_sites + 1
    th = sp.pi/n
    tau, sig = (1, 2, 3), (4, 5, 6)      # mode-disjoint, k(tau)+k(sigma) = 21 odd => eps = -1
    u = [[sp.sin(k*x*th) for x in range(n+1)] for k in range(1, n_sites+1)]

    def slater(t3, x, y, z):
        return sp.Matrix(3, 3, lambda r, c: u[t3[r]-1][(x, y, z)[c]]).det()

    # robustness note: sympy's simplify settles the algebraic difference at n = 10 (odd N = 9);
    # at even n the nested dyadic radicals can defeat simplify (a surd that IS zero but does not
    # reduce), so an even-n anchor would need nsimplify/high-precision instead. N = 9 is chosen
    # deliberately.
    def g(t3, b, c):
        return slater(t3, b-1, b, c) + slater(t3, b, b+1, c)

    lhs = 0
    for b in range(1, n_sites+1):
        for c in range(1, n_sites+1):
            if c != b:
                lhs += (1 if c > b else -1)*g(tau, b, c)*g(sig, b, c)
    rhs = _form([k*th for k in tau], [l*th for l in sig], halfangle=False)
    assert sp.simplify(lhs - rhs) == 0, "exact anchor failed"
    _pass("G12 end-to-end EXACT anchor at N=9, pair ({1,2,3},{4,5,6}): difference == 0, no float")


if __name__ == "__main__":
    gate1_coefficients()
    gate2_product_to_sum_and_parity()
    gate3_d4_oddness()
    gate4_parity_uniformity()
    gate5_d3_preconditions()
    gate6_structural_identity()
    gate7_normalization_reading()
    gate8_d1_symbolic()
    gate9_d2_symbolic()
    gate10_d3_d4_symbolic()
    gate11_transcription_pin()
    if "--slow" in sys.argv:
        gate12_exact_anchor()
    else:
        print("(G12, the ~90 min exact end-to-end anchor, runs with --slow; it PASSED 2026-07-14)")
    print("THE ASSEMBLY (D) IS SYMBOLIC: given D1-D4 (proved in cross_triple_orthogonality.py),")
    print("(U+ - U-)*(n/2)^3 = cross_form(a; b) on mode-disjoint pairs; at eps=+1 both sides vanish.")
