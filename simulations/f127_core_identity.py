"""The fragile-thing hunt, scout 4: exact core-identity proof, NO sympy cancel.

Object:  T(a;b) = Sum_{i,j} (-1)^(i+j) alpha_i(a) alpha_j(b) cot((a_i+b_j)/2)
Claim:   T = 0 on  Ca (Sum cos a = 0),  Cb (Sum cos b = 0),  S (Sum a + Sum b = 0 mod 2pi).

z-coords: z_k = e^{i a_k}, w_k = e^{i b_k}; the sheet S is MONOMIAL: w3 = 1/(m0 z3),
m0 = z1 z2 w1 w2.  Factored building blocks (all small):
    alpha_i(z)  = (z_l - z_j)(z_j z_l + 1) / (2 I z_j z_l)
    cot((x+y)/2) with e^{ix}=u, e^{iy}=v:  I (u v + 1)/(u v - 1)
Numerator N assembled as a sum of 9 PRODUCTS over the common denominator
(no cancel; monomial denominators are units and dropped).

Then, in K[z3] with K = QQ_I(z1,z2,w1,w2):
    Qz = z3^2 + Sz z3 + 1,       Sz = z1 + 1/z1 + z2 + 1/z2      (Ca)
    P  = A z3^2 + B z3 + C,      A = m0^2, B = Sw m0, C = 1      (Cb via sheet)
    N mod Qz = r0 + r1 z3;   E := r0 (B - A Sz) - r1 (C - A);
    Res := Res_z3(Qz, P).
Claim reduces to: E vanishes on {Res = 0}   (both quadratic roots are variety pts).

Steps:
  [A] assemble N; sanity gate vs float T at random angle points on the sheet.
  [B] NUMERIC DECISION on the complex variety: random complex z1,z2,w1;
      solve Res(w2)=0; test E at every root.  (Also test the torus-slice
      branches separately.)
  [C] exact divisibility E mod Res over QQ_I (only if [B] says yes).
"""
import math
import random
import sys

import sympy as sp

sys.path.insert(0, "simulations")

z1, z2, z3, w1, w2 = sp.symbols("z1 z2 z3 w1 w2")
I = sp.I
V5 = (z1, z2, z3, w1, w2)


def build_numerator():
    """N and the common denominator's pole factors, exactly, product-form only."""
    m0 = z1 * z2 * w1 * w2
    # b-side monomials with w3 = 1/(m0 z3)
    wmon = [w1, w2, None]

    # alpha numerator/denominator per slot, a-side
    def alpha_z(i):
        j, l = [t for t in range(3) if t != i]
        zs = [z1, z2, z3]
        return (zs[l] - zs[j]) * (zs[j] * zs[l] + 1), zs[j] * zs[l]   # num, den(monomial)

    # a-side: plain.  b-side: slots may involve w3 = 1/(m0 z3).
    # sin b_k = (w_k - 1/w_k)/(2I).  For k=3: (1/(m0 z3) - m0 z3)/(2I)
    #   = (1 - (m0 z3)^2) / (2I m0 z3).
    def sin_w(k):
        if k < 2:
            v = [w1, w2][k]
            return v ** 2 - 1, v                     # num/(2I den)
        return 1 - (m0 * z3) ** 2, m0 * z3

    def alpha_w(j):
        p, q = [t for t in range(3) if t != j]
        np_, dp = sin_w(q)
        nq_, dq = sin_w(p)
        # alpha_j = sin b_q' ... committed convention: alpha = sin(ang[l]) - sin(ang[j])
        # with (j_, l_) = complement ascending; here complement = (p, q), p < q:
        # alpha_j(b) = sin b_q - sin b_p
        return np_ * dq - nq_ * dp, dp * dq          # num/(2I den)

    # cot((a_i + b_j)/2) = I (z_i w_j + 1)/(z_i w_j - 1); for j=3:
    #   z_i w3 = z_i/(m0 z3):  I (z_i + m0 z3)/(z_i - m0 z3)
    def cot_ij(i, j):
        zs = [z1, z2, z3]
        if j < 2:
            u = zs[i] * [w1, w2][j]
            return (u + 1), (u - 1)                  # times I; num, den(POLE)
        return (zs[i] + m0 * z3), (zs[i] - m0 * z3)  # times I

    poles = {}
    for i in range(3):
        for j in range(3):
            poles[(i, j)] = cot_ij(i, j)[1]

    N = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            an, ad = alpha_z(i)
            bn, bd = alpha_w(j)
            cn, cd = cot_ij(i, j)
            # term = (-1)^(i+j) * an/(2I ad) * bn/(2I bd) * I cn/cd
            #      = (-1)^(i+j) * I/(4 I^2) * an bn cn / (ad bd cd)
            #      = (-1)^(i+j) * (-I/4) ... constants global, drop 1/4, keep sign+I
            other = sp.Integer(1)
            for kl, pf in poles.items():
                if kl != (i, j):
                    other = sp.expand(other * pf)
            # monomial denominators ad, bd are units: multiply the WHOLE sum by
            # the global monomial LCM implicitly -- equivalently, multiply this
            # term by (AD/ad)(BD/bd) where AD, BD are fixed monomial caps.
            term = (-1) ** (i + j) * an * bn * cn * other
            # clear the term's own monomial dens by multiplying the others in:
            N = N + term * MON_CAP_A / ad * MON_CAP_B(j) / bd
    return N, poles


# monomial caps: lcm of all ad = (z1 z2 z3)^2? ad = z_j z_l over complement:
# products over i: ad_i = z_j z_l; lcm over i = z1 z2 z3 (each pair misses one).
# Use cap = z1 z2 z3 (divides evenly: cap/ad = z_i).
MON_CAP_A = z1 * z2 * z3


def MON_CAP_B(_j):
    # bd_j = dp dq with d in {w1, w2, m0 z3}; lcm over j = w1 w2 m0 z3.
    return w1 * w2 * (z1 * z2 * w1 * w2) * z3


def naked_falsification(n=25):
    """The sheet relation is LOAD-BEARING: naked T on Ca+Cb alone is O(1), not 0."""
    def alphaf(ang, i):
        j, l = [t for t in range(3) if t != i]
        return math.sin(ang[l]) - math.sin(ang[j])

    def Tf(a, b):
        return sum(((-1) ** (i + j)) * alphaf(a, i) * alphaf(b, j)
                   / math.tan((a[i] + b[j]) / 2)
                   for i in range(3) for j in range(3))

    def rand_on_ca():
        while True:
            a1 = random.uniform(-math.pi, math.pi)
            a2 = random.uniform(-math.pi, math.pi)
            c = -math.cos(a1) - math.cos(a2)
            if abs(c) <= 0.98:
                return [a1, a2, random.choice([1, -1]) * math.acos(c)]

    random.seed(4)
    biggest = max(abs(Tf(rand_on_ca(), rand_on_ca())) for _ in range(n))
    print(f"[0] naked T on Ca+Cb only: max|T| = {biggest:.2e} (must be O(1): "
          f"the sheet constraint is not removable)")
    assert biggest > 1.0, "naked T unexpectedly small -- falsification gate broken"


def main():
    naked_falsification()
    print("[A] assembling N ...")
    N, poles = build_numerator()
    N = sp.expand(N)
    print(f"    N: {len(sp.Add.make_args(N))} terms, deg_z3 = {sp.degree(N, z3)}")

    # ---- sanity gate: N == 0 exactly where float T == 0? compare N against
    # T * (denominators) at random SHEET angle points (no Ca/Cb):
    def alphaf(ang, i):
        j, l = [t for t in range(3) if t != i]
        return math.sin(ang[l]) - math.sin(ang[j])

    def Tf(a, b):
        return sum(((-1) ** (i + j)) * alphaf(a, i) * alphaf(b, j)
                   / math.tan((a[i] + b[j]) / 2)
                   for i in range(3) for j in range(3))

    random.seed(11)
    worst = 0.0
    for _ in range(6):
        a = [random.uniform(-3, 3) for _ in range(3)]
        b = [random.uniform(-3, 3) for _ in range(2)]
        b3 = -(sum(a) + sum(b))                     # the sheet
        b = b + [b3]
        zv = [complex(math.cos(x), math.sin(x)) for x in a]
        wv = [complex(math.cos(x), math.sin(x)) for x in b]
        sub = {z1: zv[0], z2: zv[1], z3: zv[2], w1: wv[0], w2: wv[1]}
        Nf = complex(sp.N(N.subs(sub), 20))
        # reconstruct the full prefactor: T = (-I/4) * N / (monomials * poles)
        polef = 1
        for pf in poles.values():
            polef *= complex(sp.N(pf.subs(sub), 20))
        monf = complex(sp.N((MON_CAP_A * MON_CAP_B(0) / 1).subs(sub), 20))
        # global constant: each term carried I/(2I)(2I) = -I/4
        Tpred = (-1j / 4) * Nf / (monf * polef)
        err = abs(Tpred - Tf(a, b)) / max(1.0, abs(Tf(a, b)))
        worst = max(worst, err)
    print(f"    gate N vs float T: worst rel err {worst:.2e}")
    assert worst < 1e-9, "N assembly does NOT match T"

    # ---- K[z3] reduction
    K = sp.FractionField(sp.QQ_I, [z1, z2, w1, w2])
    Sz = z1 + 1 / z1 + z2 + 1 / z2
    m0 = z1 * z2 * w1 * w2
    Sw = w1 + 1 / w1 + w2 + 1 / w2
    A_, B_, C_ = m0 ** 2, sp.together(Sw * m0), sp.Integer(1)

    Np = sp.Poly(N, z3, domain=K)
    Qzp = sp.Poly(z3 ** 2 + Sz * z3 + 1, z3, domain=K)
    Pp = sp.Poly(A_ * z3 ** 2 + B_ * z3 + C_, z3, domain=K)
    q, r = sp.div(Np, Qzp)
    r0 = r.as_expr().coeff(z3, 0)
    r1 = r.as_expr().coeff(z3, 1)
    print(f"    reduced mod Qz: r0/r1 obtained")

    E = sp.together(r0 * (B_ - A_ * Sz) - r1 * (C_ - A_))
    En = sp.expand(sp.fraction(sp.cancel(E))[0])
    Res = sp.resultant(Qzp.as_expr(), Pp.as_expr(), z3)
    Resn = sp.expand(sp.fraction(sp.cancel(sp.together(Res)))[0])
    print(f"    E numerator: {len(sp.Add.make_args(En))} terms; "
          f"Res numerator: {len(sp.Add.make_args(Resn))} terms")

    # ---- [B] numeric decision on the COMPLEX variety
    print("[B] testing E on random branches of {Res = 0} ...")
    import numpy as np
    rng = random.Random(23)
    bad = good = 0
    for trial in range(8):
        sub = {z1: complex(rng.uniform(-1, 1), rng.uniform(-1, 1)),
               z2: complex(rng.uniform(-1, 1), rng.uniform(-1, 1)),
               w1: complex(rng.uniform(-1, 1), rng.uniform(-1, 1))}
        Rw = sp.Poly(Resn.subs(sub), w2)
        coeffs = [complex(c) for c in Rw.all_coeffs()]
        roots = np.roots(coeffs)
        Ecoeffs = [complex(c) for c in sp.Poly(En.subs(sub), w2).all_coeffs()]
        scale = max(abs(c) for c in Ecoeffs)

        def Ew(x, cs=Ecoeffs):
            acc = 0j
            for c in cs:
                acc = acc * x + c
            return acc

        for rt in roots:
            if abs(rt) < 1e-8:
                continue
            val = abs(Ew(rt)) / max(scale, 1e-30)
            if val < 1e-7:
                good += 1
            else:
                bad += 1
                if bad <= 5:
                    print(f"    NONZERO branch: |E|/scale = {val:.2e} at w2 = {rt:.4f}")
    print(f"    branches: {good} vanish, {bad} do not")
    if bad == 0:
        print("[B] E vanishes on ALL branches of Res=0  ->  divisibility expected.")
        print("[C] exact division E mod Res over QQ_I ...")
        quo, rem = sp.div(sp.Poly(En, z1, z2, w1, w2, domain=sp.QQ_I),
                          sp.Poly(Resn, z1, z2, w1, w2, domain=sp.QQ_I))
        print(f"    remainder zero: {rem.is_zero}")
        if rem.is_zero:
            print("*** CORE IDENTITY PROVED over Q(i): Res | E, i.e. T = 0 on the variety. ***")
    else:
        print("[B] E does NOT vanish on some complex branches -> the theorem lives on")
        print("    the torus-relevant component only; identify the factor of Res:")
        print("    factoring Res ...")
        for fac, mult in sp.factor_list(Resn)[1]:
            print(f"      factor (mult {mult}, {len(sp.Add.make_args(sp.expand(fac)))} terms)")


if __name__ == "__main__":
    main()
