"""Gate for the triangle's palindromic locus with planar fields: Heisenberg bonds on the triangle (d, u, v), the jump Z
on d alone, no field on d, and fields f_u, f_v of any direction perpendicular to z.

PROOF_PALINDROME_COMPLEMENT_CONNECTION, section "The triangle with planar fields", proves: such a connected row pairs
exactly when it lies in A (f_u parallel to f_v, one of them zero included), B (J_du = J_dv, |f_u| = |f_v|), C (e2 = 0,
|f_u| = |f_v|) or on the three relations R1: J_uv (J_du - J_dv)^2 = J_du J_dv (J_du + J_dv), R2: f_u.f_v = -2 e2,
R3: (J_du + J_uv)|f_u|^2 + (J_dv + J_uv)|f_v|^2 + 2(J_du + J_dv + J_uv) f_u.f_v = 0. The computational steps are checked here
in exact arithmetic (sympy), the main ones with a control that sends a wrong input through the same door. Notation: H = sum_x,y |x><y| (x) H_xy
over the bit of d, H00 the dimer u-v with fields a = f_u + J_du z, b = f_v + J_dv z, H11 with a' = f_u - J_du z,
b' = f_v - J_dv z, V = H01 = 2(J_du s-_u + J_dv s-_v) the one hop, K = V U with U a twin map, H11 = U H00 U^-1 (up to scale).

G1  (symbolic, every parameter) The twin map U = u (x) u, u = (n.sigma) Z, n = a' x b', carries H00 onto H11; Theta =
    (time reversal) o (pi turn about a x b) fixes H00, carries K to its adjoint and squares to +1; U^dag U = |n|^4.
    Control: a z part on f_u, the twin map rebuilt for it, breaks the twin.
G2  (symbolic) The dimer's characteristic polynomial in secular form, the basis of the simple-spectrum lemma: with
    g = a + b, k = a - b, mu = x - r, |g|^2 det(x - H) = |g|^2 mu (mu + 4r)(mu^2 - |g|^2) - (k.g)^2 (mu^2 - |g|^2)
    - (|k|^2 |g|^2 - (k.g)^2) mu^2.
G3  (symbolic, chart J_du = 1, f_u = (A, 0), f_v = (B1, B2)) Four alternating words Im tr(H^a K H^b K^dag H^c K H^d K^dag):
    (0,1,1,2) = 128 (t-1) A B2 |n|^8 e2 R1 exactly; (1,3,0,1) = const A B2 |n|^8 R1 Q with
    Q = s t (|f_v|^2 - |f_u|^2) + (t - s) e2 (s + t + r); on R1, (2,1,0,2) carries R2 and on R1 and R2, (0,3,1,2) carries R3,
    every other factor from {t, t-1, t^2+1, A, B2, |n|^2}. Control: without the ZZ part of the bond u-v, R1 does not divide
    the first word. The charts: every R1 triple with J_du != 0, J_du != J_dv is a multiple of (w, w+1, w(w+1)(2w+1)), every
    e2 = 0 triple with J_du != 0, J_du + J_dv != 0 one of (1+t, t(1+t), -t); on E, e2 = 2 J_du J_dv (J_du^2 + J_dv^2)/(J_du - J_dv)^2
    and J_uv = 0 exactly when J_dv = -J_du.
G4  (symbolic) Sufficiency on R1, R2, R3: in the chart J = (w, w+1, w(w+1)(2w+1)), f_u = (A, 0), f_v = (-2 e2/A, B2) with B2^2
    from R3, the system "F K and K F Hermitian" for F = sum (x_k + i B2 y_k) H00^k has the explicit kernel vector v(w, A)
    (M v = 0 identically on every row), F F^dag = sigma 1 identically, and the factors of sigma that could vanish (Q1, Q2)
    meet no real point of the locus (the discriminant of Q2 in A^2 has no real zero, so its roots are simple at every
    breakpoint; signs read exactly by counting roots). Control: R3 shifted, the system has no kernel at a random point.
G5  (symbolic) The degenerate stratum J_uv = 0, |a| = |b|: the levels are +-2 alpha and 0 twice; the six commutant
    elements commute with H00 and are independent; the E+/E- invariant Npm factors as stated; on the branch
    t = -s (symmetric gauge f_u = (p, q), f_v = (p, -q)) maximal minors share only p, q, p^2+q^2, p^2+q^2+1, p^2-q^2-2;
    on the branch Q5 the minors share only factors that vanish at excluded rows (t = 0, +-1, mu = 0, +-t) or nowhere
    real, and two of them no common real zero; the corner
    p^2 - q^2 = 2 has a one-dimensional kernel and det F with only factors that vanish at q = 0. Control: p^2 - q^2 = 3
    has no kernel; the common-zero test reports a planted pair of factors meeting at a real point, on both branches.
G6  (symbolic) C at every angle: chart J = (1+t, t(1+t), -t), symmetric gauge; kernel one-dimensional over Q(t, p, q),
    F F^dag = sigma 1 with sigma positive off pq = 0. Control: |f_u| != |f_v| on e2 = 0 has no kernel.
G7  (exact rows, a second route sharing only the block builders) Near and far counted straight from the block equations
    over Q(i), with no K, no U and no words: the row
    J = (1,2,6), f_u = 8X, f_v = -5X + 3Y and the corner J = (2,-2,0), f_u = (1,4), f_v = (4,1) pair; their perturbations
    break; and the carrier G4 builds at w = 1, A = 8 is minus the eighteen-string element of n3_offletter_locus_gate.py G5
    (its strings copied here as literals).
Not checked here, argued in the proof: the reduction to a Hermitian far element (step 2), the reality of the walks
(step 3) and of the block product on the stratum (step 6), the monotonicity of the secular function (step 4), the case
logic of steps 5 and 6 around the computed factors (with the conjugation eta2 -> -eta2 and the reach of the charts), the
sign of sigma's denominator in G4 (powers of A, w+1, 2w^2+w+1), and the sufficiency of A and B (the section above).
Output: simulations/results/n3_planar_field_locus_gate.txt. Run: python simulations/n3_planar_field_locus_gate.py
(about twenty minutes)."""
import itertools, random, sys, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
OUT = []
FAILS = 0
T0 = time.time()


def check(name, ok, detail=""):
    global FAILS
    line = ("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else "")
    if not ok:
        FAILS += 1
    OUT.append(line)
    print(line, f"({time.time() - T0:.0f} s)", flush=True)


I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Z = sp.Matrix([[1, 0], [0, -1]])
SM = sp.Matrix([[0, 0], [1, 0]])          # |1><0|, the lowering of bit 0 -> 1
kp = sp.kronecker_product
S3 = (X, Y, Z)


def dsig(v):
    return sum((v[k] * S3[k] for k in range(3)), sp.zeros(2))


def dimer(r, a, b, zz=True):
    bond = r * (kp(X, X) + kp(Y, Y) + (kp(Z, Z) if zz else sp.zeros(4)))
    return bond + kp(dsig(a), I2) + kp(I2, dsig(b))


def twin_U(s, t, fu, fv):
    n = sp.Matrix([fu[0], fu[1], -s]).cross(sp.Matrix([fv[0], fv[1], -t]))
    u = dsig(n) * Z
    return kp(u, u), sp.expand((n.T * n)[0])


def hop(s, t):
    return 2 * (s * kp(SM, I2) + t * kp(I2, SM))


def conjI(e):
    return e.subs(sp.I, -sp.I)


def factor_set(e):
    return {f for f, _ in sp.factor_list(sp.expand(e))[1]}


# ---------------------------------------------------------------- G1
s, t, r, A1, A2, B1, B2, c3 = sp.symbols("s t r A1 A2 B1 B2 c3", real=True)
fu = (A1, A2, 0); fv = (B1, B2, 0)
H00 = dimer(r, (A1, A2, s), (B1, B2, t)); H11 = dimer(r, (A1, A2, -s), (B1, B2, -t))
U, n2 = twin_U(s, t, fu, fv)
K = sp.expand(hop(s, t) * U)
m = sp.Matrix([A1, A2, s]).cross(sp.Matrix([B1, B2, t])); mm = sp.expand((m.T * m)[0])
Th = kp(dsig(m) * (sp.I * Y), dsig(m) * (sp.I * Y))      # Theta = Th o complex conjugation
check("G1 the twin map carries H00 onto H11: U H00 U^dag = |n|^4 H11, identically in s, t, r and planar f_u, f_v",
      sp.expand(U * H00 * U.H - n2 ** 2 * H11) == sp.zeros(4))
check("G1 Theta fixes H00 and carries K = V U to its adjoint, identically",
      sp.expand(Th * H00.conjugate() * Th.H - mm ** 2 * H00) == sp.zeros(4)
      and sp.expand(Th * K.conjugate() * Th.H - mm ** 2 * K.H) == sp.zeros(4))
H00c = dimer(r, (A1, A2, s + c3), (B1, B2, t)); H11c = dimer(r, (A1, A2, c3 - s), (B1, B2, -t))
check("G1 Theta squares to +1 and U^dag U = |n|^4 1, identically",
      sp.expand(Th * Th.conjugate() - mm ** 2 * sp.eye(4)) == sp.zeros(4) and sp.expand(U.H * U - n2 ** 2 * sp.eye(4)) == sp.zeros(4))
Uc, n2c3 = twin_U(s - c3, t, fu, fv)          # the same construction for the row with f_u's z part: a' = (A1, A2, c3 - s)
check("G1 control: a z part c3 on f_u breaks the twin identity, with the twin map rebuilt for that row",
      sp.expand(Uc * H00c * Uc.H - n2c3 ** 2 * H11c) != sp.zeros(4))

# ---------------------------------------------------------------- G2
x = sp.symbols("x")
a = sp.symbols("a1:4", real=True); b = sp.symbols("b1:4", real=True)
Hd = dimer(r, a, b)
g = [a[i] + b[i] for i in range(3)]; kv = [a[i] - b[i] for i in range(3)]
g2 = sum(e * e for e in g); k2 = sum(e * e for e in kv); kg = sum(kv[i] * g[i] for i in range(3))
mu = x - r
lhs = sp.expand(g2 * (x * sp.eye(4) - Hd).det())
rhs = sp.expand(g2 * mu * (mu + 4 * r) * (mu ** 2 - g2) - kg ** 2 * (mu ** 2 - g2) - (k2 * g2 - kg ** 2) * mu ** 2)
check("G2 secular form of the dimer's characteristic polynomial (singlet at mu = -4r, triplet in the field g, coupled by k)",
      sp.expand(lhs - rhs) == 0)
check("G2 control: the same form with the singlet at mu = -2r fails",
      sp.expand(lhs - (rhs + g2 * mu * (-2 * r) * (mu ** 2 - g2))) != 0)

# ---------------------------------------------------------------- G3
t, r, A, B1, B2 = sp.symbols("t r A B1 B2", real=True)


def words(wds, zz=True):
    H = dimer(r, (A, 0, 1), (B1, B2, t), zz)
    U, _ = twin_U(1, t, (A, 0, 0), (B1, B2, 0))
    Kw = sp.expand(hop(1, t) * U); Kd = Kw.H
    P = [sp.eye(4), H, sp.expand(H * H)]; P.append(sp.expand(P[2] * H))
    out = {}
    for wd in wds:
        w_ = sp.eye(4)
        for i, e in enumerate(wd):
            w_ = sp.expand(w_ * P[e] * (Kw if i % 2 == 0 else Kd))
        out[wd] = sp.expand(sp.im(sp.expand(w_.trace())))
    return out


n2c = A ** 2 * B2 ** 2 + A ** 2 * t ** 2 - 2 * A * B1 * t + B1 ** 2 + B2 ** 2     # |n|^2 in the chart
e2c = t + r + t * r
R1c = r * (1 - t) ** 2 - t * (1 + t)
Qc = t * (B1 ** 2 + B2 ** 2 - A ** 2) + (t - 1) * e2c * (1 + t + r)
W = words([(0, 1, 1, 2), (1, 3, 0, 1), (2, 1, 0, 2), (0, 3, 1, 2)])
check("G3 word (0,1,1,2) = 128 (t-1) A B2 |n|^8 e2 R1, exactly",
      sp.expand(W[(0, 1, 1, 2)] - 128 * (t - 1) * A * B2 * n2c ** 4 * e2c * R1c) == 0)
q13 = sp.cancel(W[(1, 3, 0, 1)] / (A * B2 * n2c ** 4 * R1c * Qc))
check("G3 word (1,3,0,1) = const A B2 |n|^8 R1 Q, Q = s t (|f_v|^2 - |f_u|^2) + (t - s) e2 (s + t + r) at s = 1",
      q13 == -128, f"const {q13}")
allowed = {t, t - 1, t ** 2 + 1, A, B2}
rr = t * (1 + t) / (1 - t) ** 2
num = sp.numer(sp.together(W[(2, 1, 0, 2)].subs(r, rr)))
R2b = sp.numer(sp.together((A * B1 + 2 * e2c).subs(r, rr)))
quo, rem = sp.div(sp.Poly(num, t, A, B1, B2), sp.Poly(R2b, t, A, B1, B2))
n2fs = factor_set(n2c)
check("G3 on R1, word (2,1,0,2) = R2 times factors from {t, t-1, t^2+1, A, B2, |n|^2}",
      rem.is_zero and not quo.is_zero and factor_set(quo.as_expr()) <= allowed | n2fs, f"cofactor factors {sorted(map(str, factor_set(quo.as_expr())))}")
b1 = -4 * t * (t ** 2 + 1) / (A * (t - 1) ** 2)                       # R2 solved on R1
num = sp.numer(sp.together(W[(0, 3, 1, 2)].subs({r: rr, B1: b1})))
R3b = sp.numer(sp.together(((1 + r) * A ** 2 + (t + r) * (B1 ** 2 + B2 ** 2) + 2 * (1 + t + r) * A * B1).subs({r: rr, B1: b1})))
quo, rem = sp.div(sp.Poly(num, t, A, B2), sp.Poly(R3b, t, A, B2))
n2b = factor_set(sp.numer(sp.together(n2c.subs({r: rr, B1: b1}))))
check("G3 on R1 and R2, word (0,3,1,2) = R3 times factors from {t, t-1, t^2+1, A, B2, |n|^2 on the branch}",
      rem.is_zero and not quo.is_zero and factor_set(quo.as_expr()) <= allowed | n2b, f"cofactor factors {len(factor_set(quo.as_expr()))}")
Wc = words([(0, 1, 1, 2)], zz=False)[(0, 1, 1, 2)]
check("G3 control: without the ZZ part of the bond u-v, R1 does not divide word (0,1,1,2)",
      not sp.div(sp.Poly(Wc, t, r, A, B1, B2), sp.Poly(R1c, t, r, A, B1, B2))[1].is_zero)
jd, jv, ju, wc_, tc_ = sp.symbols("jd jv ju wc tc")
R1g = ju * (jd - jv) ** 2 - jd * jv * (jd + jv)
juR1 = jd * jv * (jd + jv) / (jd - jv) ** 2
check("G3 the charts: R1 at (w, w+1, w(w+1)(2w+1)) vanishes, and R1 solved for J_uv at (1, (w+1)/w) is (w+1)(2w+1), the chart's "
      "ratio; e2 = 0 solved for J_uv at (1+t, t(1+t)) is -t",
      sp.expand(R1g.subs({jd: wc_, jv: wc_ + 1, ju: wc_ * (wc_ + 1) * (2 * wc_ + 1)})) == 0
      and sp.cancel(juR1.subs({jd: 1, jv: (wc_ + 1) / wc_}) - (wc_ + 1) * (2 * wc_ + 1)) == 0
      and sp.cancel((-jd * jv / (jd + jv)).subs({jd: 1 + tc_, jv: tc_ * (1 + tc_)}) + tc_) == 0)
check("G3 on E: e2 = 2 J_du J_dv (J_du^2 + J_dv^2)/(J_du - J_dv)^2, and J_uv = 0 exactly when J_du J_dv (J_du + J_dv) = 0",
      sp.cancel(jd * jv + juR1 * (jd + jv) - 2 * jd * jv * (jd ** 2 + jv ** 2) / (jd - jv) ** 2) == 0
      and sp.expand(sp.numer(sp.together(juR1)) - jd * jv * (jd + jv)) == 0 and sp.expand(sp.denom(sp.together(juR1)) - (jd - jv) ** 2) == 0)

# ---------------------------------------------------------------- G4
w, A, B2, B1s, rs, ts, ss = sp.symbols("w A B2 B1s rs ts ss", real=True)
xs = sp.symbols("x0:4", real=True); ys = sp.symbols("y0:4", real=True); unk = list(xs) + list(ys)
Sw, Tw, Rw = w, w + 1, w * (w + 1) * (2 * w + 1)
e2w = sp.expand(Sw * Tw + Rw * (Sw + Tw)); cw = -2 * e2w; kw = -2 * (Sw + Tw + Rw) * cw; Dw = Tw + Rw


def e_system(shift=0):
    """rows of 'F K and K F Hermitian' for F = sum (x_k + i B2 y_k) H00^k, reduced by parity in B2 and put on the chart"""
    H = dimer(rs, (A, 0, ss), (B1s, B2, ts))
    U, _ = twin_U(ss, ts, (A, 0, 0), (B1s, B2, 0))
    Kq = sp.expand(hop(ss, ts) * U); Kd = Kq.H
    P = [sp.eye(4), H, sp.expand(H * H)]; P.append(sp.expand(P[2] * H))
    F = sum(((xs[k] + sp.I * B2 * ys[k]) * P[k] for k in range(4)), sp.zeros(4))
    beta = (kw + shift - (Sw + Rw) * A ** 2) / Dw - cw ** 2 / A ** 2
    rows = []
    for E in (sp.expand(F * Kq - Kd * F.H), sp.expand(Kq * F - F.H * Kd)):
        for e in E:
            for part in ((e + conjI(e)) / 2, (e - conjI(e)) / (2 * sp.I)):
                part = sp.expand(part)
                if part == 0:
                    continue
                ev = sp.expand(part.subs(B2, -B2))
                if sp.expand(ev - part) == 0:
                    pass
                elif sp.expand(ev + part) == 0:
                    part = sp.expand(sp.cancel(part / B2))
                else:
                    raise SystemExit("mixed parity in B2")
                out = sum(cf * beta ** (kk // 2) for (kk,), cf in sp.Poly(part, B2).terms())
                out = sp.numer(sp.together(out.subs({ss: Sw, ts: Tw, rs: Rw, B1s: cw / A})))
                rows.append([sp.expand(sp.diff(out, v)) for v in unk])
    return sp.Matrix(rows), beta


Mp, betaE = e_system()
fu2, fv2, fuv = A ** 2, (cw / A) ** 2 + betaE, cw
check("G4 the chart satisfies R2 (f_u.f_v = -2 e2) and R3 with B2^2 = beta, identically",
      sp.expand(fuv + 2 * e2w) == 0 and sp.cancel((Sw + Rw) * fu2 + (Tw + Rw) * fv2 + 2 * (Sw + Tw + Rw) * fuv) == 0)
pt = {w: sp.Rational(3, 2), A: sp.Rational(7, 3)}
Mn = Mp.subs(pt)
check("G4 the parity-reduced system on R1, R2, R3 has rank 7 at a random point (a one-dimensional kernel)",
      Mn.rank() == 7, f"{Mp.shape}")
deg = lambda i: max((sp.Poly(e, w, A).total_degree() for e in Mp.row(i) if e != 0), default=-1)
sel = []
for i in sorted(range(Mp.rows), key=deg):
    if deg(i) >= 0 and Mn[sel + [i], :].rank() == len(sel) + 1:
        sel.append(i)
    if len(sel) == 7:
        break
R7 = Mp[sel, :]
v = []
for j in range(8):
    d = DomainMatrix.from_Matrix(R7[:, [c for c in range(8) if c != j]]).convert_to(sp.ZZ[w, A]).det()
    v.append((-1) ** j * sp.ZZ[w, A].to_sympy(d))
gg = v[0]
for e in v[1:]:
    gg = sp.gcd(gg, e)
v = [sp.cancel(e / gg) for e in v]
prod = (DomainMatrix.from_Matrix(Mp).convert_to(sp.ZZ[w, A]) * DomainMatrix.from_Matrix(sp.Matrix(v)).convert_to(sp.ZZ[w, A])).to_Matrix()
vdeg = [(sp.Poly(e, w, A).degree(w), sp.Poly(e, w, A).degree(A)) for e in v]
check("G4 the cofactor vector v(w, A) solves every row identically (M v = 0 over Z[w, A]), degrees at most 25 in w and 6 in A",
      all(sp.expand(e) == 0 for e in prod) and any(e != 0 for e in v) and max(d_[0] for d_ in vdeg) <= 25 and max(d_[1] for d_ in vdeg) <= 6,
      f"degrees {[(sp.Poly(e, w, A).degree(w), sp.Poly(e, w, A).degree(A)) for e in v]}")
# F F^dag through the characteristic polynomial: F = p(H), F F^dag = p(H) pbar(H), reduced mod charpoly(H)
lam = sp.symbols("lam")
HE = dimer(Rw, (A, 0, Sw), (cw / A, B2, Tw))
cpE = sp.Poly(sp.expand((lam * sp.eye(4) - HE).det()), lam)
cks = [v[k] + sp.I * B2 * v[4 + k] for k in range(4)]
pp = sp.expand(sum(cks[k] * lam ** k for k in range(4)) * sum(conjI(cks[k]) * lam ** k for k in range(4)))
rem = sp.Poly(pp, lam).rem(cpE)


def red_b2(e):
    e = sp.expand(e)
    return sp.numer(sp.together(sum(cf * betaE ** (kk // 2) * (B2 if kk % 2 else 1) for (kk,), cf in sp.Poly(e, B2).terms())))


coefs = [red_b2(rem.coeff_monomial(lam ** k)) for k in range(4)]
check("G4 F F^dag = sigma 1 identically on the locus (F F^dag reduced modulo the characteristic polynomial is a constant)",
      all(sp.expand(c_) == 0 for c_ in coefs[1:]) and sp.expand(coefs[0]) != 0)
sig = sp.factor(coefs[0])
aa = sp.symbols("aa", real=True)
Q1 = -aa + 16 * w ** 5 + 48 * w ** 4 + 64 * w ** 3 + 52 * w ** 2 + 24 * w + 4
Q2 = (2 * aa ** 2 * w ** 3 + 3 * aa ** 2 * w ** 2 + 2 * aa ** 2 * w - 32 * aa * w ** 7 - 112 * aa * w ** 6 - 192 * aa * w ** 5
      - 204 * aa * w ** 4 - 136 * aa * w ** 3 - 56 * aa * w ** 2 - 12 * aa * w - aa + 128 * w ** 11 + 704 * w ** 10
      + 1728 * w ** 9 + 2560 * w ** 8 + 2560 * w ** 7 + 1808 * w ** 6 + 912 * w ** 5 + 320 * w ** 4 + 72 * w ** 3 + 8 * w ** 2)
safe = {A, w + 1, 2 * w + 1, 2 * w ** 2 + w + 1, 2 * w ** 2 + 2 * w + 1}
sfs = factor_set(sig)
Q1A = sp.expand(Q1.subs(aa, A ** 2)); Q2A = sp.expand(Q2.subs(aa, A ** 2))
ok_f = all(f in safe or sp.expand(f - Q1A) == 0 or sp.expand(f + Q1A) == 0 or sp.expand(f - Q2A) == 0 or sp.expand(f + Q2A) == 0 for f in sfs)
check("G4 the factors of sigma: A, w+1, 2w+1, 2w^2+w+1, 2w^2+2w+1 and the two factors Q1, Q2 (in a = A^2), nothing else",
      ok_f and any(sp.expand(f - Q1A) == 0 or sp.expand(f + Q1A) == 0 for f in sfs) and any(sp.expand(f - Q2A) == 0 or sp.expand(f + Q2A) == 0 for f in sfs),
      f"{len(sfs)} factors")
beta_a = sp.together(betaE.subs(A ** 2, aa).subs(A, sp.sqrt(aa)))
a_on_Q1 = sp.solve(Q1, aa)[0]
check("G4 Q1 meets no real point: on Q1, A^2 = 4(w+1)^2 (4w^3+4w^2+4w+1) and B2^2 (4w^3+4w^2+4w+1) = -16 w^2 (2w+1)^2 (w^2+w+1)^2",
      sp.expand(a_on_Q1 - 4 * (w + 1) ** 2 * (4 * w ** 3 + 4 * w ** 2 + 4 * w + 1)) == 0
      and sp.simplify(beta_a.subs(aa, a_on_Q1) * (4 * w ** 3 + 4 * w ** 2 + 4 * w + 1) + 16 * w ** 2 * (2 * w + 1) ** 2 * (w ** 2 + w + 1) ** 2) == 0)
numb = sp.numer(beta_a); denb = sp.denom(beta_a)


def sign_at(P, root):
    """exact sign of the polynomial P (in aa, rational coefficients) at a real algebraic root: LC sign times (-1) to the
    number of P's real roots above it; 0 when root is one of them"""
    Pp = sp.Poly(P, aa)
    if Pp.degree() < 1:
        return sp.sign(Pp.as_expr())
    rts = Pp.real_roots()
    if any(r_ == root for r_ in rts):
        return 0
    return sp.sign(Pp.LC()) * (-1) ** sum(1 for r_ in rts if r_ > root)


def beta_positive(xv, root):
    return sign_at(numb.subs(w, xv), root) * sign_at(denb.subs(w, xv), root) > 0


sext = sp.expand(192 * w ** 6 + 576 * w ** 5 + 688 * w ** 4 + 416 * w ** 3 + 130 * w ** 2 + 18 * w + 1)
check("G4 the discriminant of Q2 in a is (2w^2+2w+1)^3 times a sextic with no real zero, so Q2's roots in a are real and simple "
      "wherever Q2 is quadratic (leading coefficient w(2w^2+3w+2), zero only at w = 0); its resultant with B2^2's numerator "
      "has real zeros only at w in {0, -1, -1/2}",
      sp.expand(sp.discriminant(Q2, aa) - (2 * w ** 2 + 2 * w + 1) ** 3 * sext) == 0 and sp.Poly(sext, w).count_roots() == 0
      and sp.expand(sp.Poly(Q2, aa).LC() - w * (2 * w ** 2 + 3 * w + 2)) == 0
      and set(sp.Poly(sp.resultant(Q2, numb, aa), w).real_roots()) <= {0, -1, sp.Rational(-1, 2)})
brk = {sp.Integer(-1), sp.Rational(-1, 2), sp.Integer(0)}
nonzero_pols = True
for pol in [sp.discriminant(Q2, aa), Q2.subs(aa, 0), sp.Poly(Q2, aa).LC(), sp.resultant(Q2, numb, aa), sp.resultant(Q2, denb, aa)]:
    pe = sp.expand(pol)
    nonzero_pols = nonzero_pols and pe != 0
    if pe.free_symbols:
        brk |= set(sp.Poly(pe, w).real_roots())
brk = sorted(brk, key=float); fl = [float(z) for z in brk]
samples = [sp.Rational(int(fl[0]) - 2)] + [sp.nsimplify((fl[i] + fl[i + 1]) / 2, rational=True, tolerance=1e-9) for i in range(len(fl) - 1)] + [sp.Rational(int(fl[-1]) + 2)]
ordered = bool(samples[0] < brk[0]) and bool(samples[-1] > brk[-1]) and all(bool(brk[i] < samples[i + 1]) and bool(samples[i + 1] < brk[i + 1]) for i in range(len(brk) - 1))
hits = 0; seen = 0
for xv in samples:
    for root in sp.Poly(Q2.subs(w, xv), aa).real_roots():
        if root > 0:
            seen += 1
            hits += 1 if beta_positive(xv, root) else 0
check("G4 Q2 meets no real point: on every interval between the breakpoints (discriminant, leading coefficient, a = 0, "
      "the resultants with B2^2's numerator and denominator) no positive root a has B2^2 > 0 (signs exact), and at a breakpoint "
      "every root is simple, so a limit of roots of the neighbouring intervals",
      hits == 0 and seen > 0 and ordered and nonzero_pols, f"{len(brk)} breakpoints, {seen} positive roots read")
hits_c = 0
for xv in samples:
    for root in sp.Poly(sp.expand((aa - 64) * (w ** 2 + 1)).subs(w, xv), aa).real_roots():
        if root > 0 and beta_positive(xv, root):
            hits_c += 1
check("G4 control: the same interval read with the curve A^2 = 64 in place of Q2 (it meets the locus, B2^2 = 9 at w = 1) "
      "reports real points", hits_c > 0, f"{hits_c}")
Mc, _ = e_system(shift=1)
check("G4 control: with R3 shifted by one the system has full rank 8 at the same point (no carrier)",
      Mc.subs(pt).rank() == 8)

# ---------------------------------------------------------------- G5
def r0_system(s_, t_, fu_, fv_):
    a_ = sp.Matrix([fu_[0], fu_[1], s_]); b_ = sp.Matrix([fv_[0], fv_[1], t_])
    Au = kp(dsig(a_), I2); Bv = kp(I2, dsig(b_)); al2 = sp.expand((a_.T * a_)[0])
    P0 = sp.expand(al2 * sp.eye(4) - Au * Bv)
    zh = sp.Matrix([0, 0, 1]); c1 = a_.cross(zh); c2 = a_.cross(c1); d1 = b_.cross(zh)
    basis = [sp.eye(4), Au, Bv, sp.expand(Au * Bv), sp.expand(kp(dsig(c1), dsig(d1)) * P0), sp.expand(kp(dsig(c2), dsig(d1)) * P0)]
    Ur = kp(dsig(sp.Matrix([fu_[0], fu_[1], 0])), dsig(sp.Matrix([fv_[0], fv_[1], 0])))      # separate pi turns
    Kr = sp.expand(hop(s_, t_) * Ur); Kd = Kr.H
    cols = []
    for Bm in basis:
        for ph in (1, sp.I):
            Fb = ph * Bm; col = []
            for E in (sp.expand(Fb * Kr - Kd * Fb.H), sp.expand(Kr * Fb - Fb.H * Kd)):
                for e in E:
                    col += [sp.expand(sp.re(e)), sp.expand(sp.im(e))]
            cols.append(col)
    M = sp.Matrix(cols).T
    rows = []; seenr = set()
    for i in range(M.rows):
        rw = [sp.together(e) for e in M.row(i)]
        if all(e == 0 for e in rw):
            continue
        dn = sp.lcm([sp.denom(e) for e in rw]); rw = tuple(sp.expand(sp.cancel(e * dn)) for e in rw)
        if rw not in seenr:
            seenr.add(rw); rows.append(list(rw))
    return sp.Matrix(rows), basis, Kr, Ur


# the E+/E- invariant Npm on the stratum
t, A, B1, B2, al, lm = sp.symbols("t A B1 B2 alpha lm", real=True)
Hs = dimer(0, (A, 0, 1), (B1, B2, t))
Us = kp(A * X, B1 * X + B2 * Y)
H11s = dimer(0, (A, 0, -1), (B1, B2, -t))
check("G5 at J_uv = 0 the separate pi turns about f_u and f_v are a twin map: U' H00 U'^dag = |f_u|^2 |f_v|^2 H11",
      sp.expand(Us * Hs * Us.H - A ** 2 * (B1 ** 2 + B2 ** 2) * H11s) == sp.zeros(4))
Ks = sp.expand(hop(1, t) * Us); Kds = Ks.H; Hs2 = sp.expand(Hs * Hs)
Pp = sp.expand(Hs2 + 2 * al * Hs); Pm = sp.expand(Hs2 - 2 * al * Hs)
wv = sp.eye(4)
for Mx in (Pp, Ks, Pp, Kds, Pm, Ks, Pm, Kds):
    wv = sp.expand(wv * Mx)
N1 = sp.expand(sp.im(wv.trace()))
PN = sp.Poly(N1, al)
check("G5 Npm = Im tr(P+ K P+ K^dag P- K P- K^dag) is even in alpha (P+ and P- trade places under alpha -> -alpha)",
      all(k_ % 2 == 0 for (k_,) in PN.monoms()))
N1e = sum(cf * (A ** 2 + 1) ** (k_ // 2) for (k_,), cf in PN.terms())
mq = B1 ** 2 + t ** 2 - 1
N1s = sp.numer(sp.together(N1e.subs({A: (lm + mq / lm) / 2, B2: (mq / lm - lm) / 2})))
Q5lit = (B1 ** 4 + 2 * B1 ** 3 * lm + 2 * B1 ** 2 * lm ** 2 + 2 * B1 ** 2 * t ** 2 - 2 * B1 ** 2 + 2 * B1 * lm ** 3
         + 2 * B1 * lm * t ** 2 - 2 * B1 * lm + lm ** 4 - 2 * lm ** 2 * t ** 2 + 4 * lm ** 2 * t - 2 * lm ** 2 + t ** 4 - 2 * t ** 2 + 1)
fsN = factor_set(N1s)
expect = {t, t - 1, t + 1, sp.expand(mq - lm ** 2), sp.expand(mq + lm ** 2), sp.expand((mq + lm ** 2) ** 2 + 4 * lm ** 2),
          sp.expand((mq + lm ** 2) ** 2 + 4 * lm ** 2 * (1 - t ** 2)), sp.expand(Q5lit)}
check("G5 Npm on the stratum factors into t, t-1, t+1, mq -+ lm^2 (B2 = 0, A = 0), (mq+lm^2)^2 + 4lm^2 = 4lm^2(A^2+1), "
      "(mq+lm^2)^2 + 4lm^2(1-t^2) = 4lm^2|f_v|^2, and Q5", fsN == expect, f"{len(fsN)} factors")
tr5 = sp.expand(sp.expand(Q5lit.subs(lm, A - B2)).subs(t ** 2, A ** 2 - B2 ** 2 - B1 ** 2 + 1))
cpS = sp.Poly(sp.expand((lm * sp.eye(4) - Hs).det()), B2)
cpR = sp.expand(cpS.rem(sp.Poly(B2 ** 2 - (A ** 2 + 1 - B1 ** 2 - t ** 2), B2)).as_expr())
check("G5 on the stratum the levels are +-2 alpha and 0 twice: det(lm - H00) = lm^2 (lm^2 - 4 alpha^2), alpha^2 = |f_u|^2 + s^2",
      sp.expand(cpR - lm ** 2 * (lm ** 2 - 4 * (A ** 2 + 1))) == 0)
check("G5 Q5 read back on the stratum: 4 (A - B2)^2 (f_u.f_v + |f_v|^2 + s t - s^2) at s = 1",
      sp.expand(tr5 - 4 * (A - B2) ** 2 * (A * B1 + B1 ** 2 + B2 ** 2 + t - 1)) == 0)


def minors(M, sym, k_, seed, pt_):
    Mn_ = M.subs(pt_); random.seed(seed); found = []
    while len(found) < k_:
        rw = sorted(random.sample(range(M.rows), M.cols))
        if Mn_[rw, :].rank() == M.cols and rw not in found:
            found.append(rw)
    R = sp.ZZ[sym]
    return [sp.expand(R.to_sympy(DomainMatrix.from_Matrix(M[rw, :]).convert_to(R).det())) for rw in found]


def common_real_zeros(f, g_, xv, yv, excluded_x):
    """exact: for each pair of distinct irreducible factors of f and g, every real root of their resultant in y must be
    one of the excluded x values (a list of rational numbers), so no common real zero has another x. Returns the pairs
    whose resultant has a real root elsewhere (empty: no common real zero outside the excluded x values)."""
    out_ = []
    ff = [h for h, _ in sp.factor_list(sp.expand(f))[1] if h.has(yv)]
    gf = [h for h, _ in sp.factor_list(sp.expand(g_))[1] if h.has(yv)]
    fx = [h for h, _ in sp.factor_list(sp.expand(f))[1] if h.free_symbols and not h.has(yv)]
    gx = [h for h, _ in sp.factor_list(sp.expand(g_))[1] if h.free_symbols and not h.has(yv)]
    for h in fx + gx:                       # factors in x alone vanish on whole lines
        for x0 in sp.Poly(h, xv).real_roots():
            if x0 not in excluded_x:
                out_.append(("factor in x alone", str(h)))
    for h1 in ff:
        for h2 in gf:
            if sp.expand(h1 - h2) == 0 or sp.expand(h1 + h2) == 0:
                out_.append(("shared factor", str(h1))); continue
            res = sp.Poly(sp.resultant(h1, h2, yv), xv)
            if res.is_zero:
                out_.append(("zero resultant", str(h1))); continue
            if res.degree() < 1:
                continue
            bad = [x0 for x0 in res.real_roots() if x0 not in excluded_x]
            if bad:
                out_.append((str(h1)[:40], str(h2)[:40], [float(x0) for x0 in bad[:4]]))
    return out_


p, q = sp.symbols("p q", real=True)
Mm, basis_m, _, _ = r0_system(1, -1, (p, q, 0), (p, -q, 0))
H00m = dimer(0, (p, q, 1), (p, -q, -1))
flat = sp.Matrix([[e for e in Bm] for Bm in basis_m]).subs({p: sp.Rational(7, 3), q: sp.Rational(2, 5)})
check("G5 on t = -s the six commutant elements commute with H00 and are independent (rank 6 at a point)",
      all(sp.expand(Bm * H00m - H00m * Bm) == sp.zeros(4) for Bm in basis_m) and flat.rank() == 6)
mins = minors(Mm, (p, q), 5, 4, {p: sp.Rational(7, 3), q: sp.Rational(2, 5)})
common = set.intersection(*[factor_set(e) for e in mins])
check("G5 branch t = -s: five maximal minors share exactly p, q, p^2+q^2, p^2+q^2+1, p^2-q^2-2",
      common == {p, q, p ** 2 + q ** 2, p ** 2 + q ** 2 + 1, p ** 2 - q ** 2 - 2})
strip = lambda e: sp.cancel(e / sp.prod([f ** mlt for f, mlt in sp.factor_list(e)[1] if f in common]))
hz = common_real_zeros(sp.expand(strip(mins[0])), sp.expand(strip(mins[1])), p, q, [sp.Integer(0)])
check("G5 branch t = -s: every real root in p of the pairwise resultants of the rest of two minors is p = 0, so they have no "
      "common real zero with p q != 0 and a palindrome needs f_u.f_v = 2 s^2 (exact)", hz == [], f"{hz}")
hzc = common_real_zeros(sp.expand(strip(mins[0]) * (q - p - 1)), sp.expand(strip(mins[1]) * (p ** 2 + q ** 2 - 5)), p, q, [sp.Integer(0)])
check("G5 control: with q - p - 1 planted on one minor and p^2 + q^2 - 5 on the other (they meet at p = 1 and p = -2) the same "
      "test reports their resultant's real roots", any(len(e_) == 3 for e_ in hzc), f"{hzc}")
ka = sp.symbols("kappa", real=True)
pk = (ka + 2 / ka) / 2; qk = (ka - 2 / ka) / 2
Mk, basis_k, _, _ = r0_system(1, -1, (pk, qk, 0), (pk, -qk, 0))
ns = DomainMatrix.from_Matrix(Mk).convert_to(sp.QQ.frac_field(ka)).nullspace().to_Matrix()
zk = list(ns.row(0)) if ns.rows == 1 else None
if zk:
    Fk = sum(((zk[2 * j] + sp.I * zk[2 * j + 1]) * basis_k[j] for j in range(6)), sp.zeros(4)).applyfunc(sp.cancel)
    dF = sp.factor(sp.cancel(Fk.det()))
    poles = [x0 for e in zk for f in factor_set(sp.denom(sp.together(e))) if f.free_symbols
             for x0 in sp.Poly(f, ka).real_roots() if x0 != 0]          # kappa = 0 lies outside the chart
    bad = [f for f in factor_set(sp.numer(dF)) if f not in {ka, ka ** 2 - 2, ka ** 2 + 2, ka ** 2 - 2 * ka + 2, ka ** 2 + 2 * ka + 2, ka ** 4 + 2 * ka ** 2 + 4}]
else:
    bad = ["no one-dimensional kernel"]
check("G5 the corner p^2 - q^2 = 2 (p, q = (kappa +- 2/kappa)/2): kernel one-dimensional over Q(kappa), det F vanishing only at "
      "kappa^2 = 2 (q = 0, collinear) and at no real kappa otherwise, the kernel vector with no pole at a real kappa != 0",
      bad == [] and zk is not None and dF != 0 and poles == [], f"{bad}")
Mx3, _, _, _ = r0_system(1, -1, (p, q, 0), (p, -q, 0))
pt3 = {p: sp.Rational(2), q: sp.Integer(1)}       # p^2 - q^2 = 3
check("G5 control: on t = -s with f_u.f_v = 3 s^2 the system has full rank 12 (no carrier)", Mx3.subs(pt3).rank() == 12)
# the Q5 branch
tq, mq_ = sp.symbols("tq mq", real=True)
Aq = (mq_ + tq ** 2 / mq_) / 2; wq = (mq_ - tq ** 2 / mq_) / 2
B1q = sp.together((tq ** 2 - tq - Aq ** 2) / Aq); B2q = sp.together((tq - 1) * wq / Aq)
check("G5 the chart (t, mu) lies on the stratum and on Q5",
      sp.simplify(Aq ** 2 + 1 - (B1q ** 2 + B2q ** 2 + tq ** 2)) == 0 and sp.simplify(Aq * B1q + B1q ** 2 + B2q ** 2 + tq - 1) == 0)
Mq, _, _, _ = r0_system(1, tq, (Aq, 0, 0), (B1q, B2q, 0))
minsq = minors(Mq, (tq, mq_), 2, 9, {tq: sp.Rational(5, 3), mq_: sp.Rational(7, 2)})
cq = factor_set(minsq[0]) & factor_set(minsq[1])
positive = {mq_ ** 2 + tq ** 2, sp.expand((mq_ ** 2 + tq ** 2) ** 2 + 4 * mq_ ** 2), sp.expand((mq_ ** 2 - tq ** 2) ** 2 + 4 * mq_ ** 2)}
G6f = sp.expand(mq_ ** 4 * (tq + 1) ** 2 - 2 * mq_ ** 2 * tq ** 2 * (tq ** 2 + 2 * tq - 7) + tq ** 4 * (tq + 1) ** 2)
excluded_f = {tq, tq - 1, tq + 1, mq_, mq_ - tq, mq_ + tq}
check("G5 branch Q5: the two minors share only excluded factors (t, t +- 1, mu, mu +- t), positive ones, and P6 (quadratic in mu^2, "
      "both roots negative where real)", cq <= excluded_f | positive | {G6f, -mq_ + tq}, f"{sorted(map(str, cq - excluded_f - positive))}")
xx = sp.symbols("xx")
G6x = sp.Poly(sp.expand(G6f.subs(mq_ ** 2, xx)), xx)
check("G5 factor P6 in x = mu^2: discriminant -64 t^4 (t+3)(t-1), so real roots only on -3 <= t <= 1, where the roots' sum "
      "2 t^2 (t^2+2t-7)/(t+1)^2 is negative (t^2+2t-7 = (t+1)^2 - 8 <= -4 there) and their product t^4 positive: no real mu",
      sp.expand(tq ** 2 + 2 * tq - 7 - ((tq + 1) ** 2 - 8)) == 0 and
      sp.expand(sp.discriminant(G6x) + 64 * tq ** 4 * (tq + 3) * (tq - 1)) == 0
      and sp.expand(G6x.all_coeffs()[1] + 2 * tq ** 2 * (tq ** 2 + 2 * tq - 7)) == 0
      and sp.expand(G6x.all_coeffs()[0] - (tq + 1) ** 2) == 0 and sp.expand(G6x.all_coeffs()[2] - tq ** 4 * (tq + 1) ** 2) == 0)
fq = [sp.expand(sp.cancel(e / sp.prod([f ** mlt for f, mlt in sp.factor_list(e)[1] if f in cq]))) for e in minsq]
hq = common_real_zeros(fq[0], fq[1], tq, mq_, [sp.Integer(0), sp.Integer(1), sp.Integer(-1)])
check("G5 branch Q5: every real root in t of the pairwise resultants of the rest of the two minors is t = 0 or t = +-1 "
      "(excluded), so they have no common real zero on a connected row and Q5 never pairs (exact)", hq == [], f"{hq}")
hqc = common_real_zeros(sp.expand(fq[0] * (mq_ - tq - 3)), sp.expand(fq[1] * (mq_ ** 2 + tq ** 2 - 25)), tq, mq_,
                        [sp.Integer(0), sp.Integer(1), sp.Integer(-1)])
check("G5 control: with mu - t - 3 planted on one minor and mu^2 + t^2 - 25 on the other (they meet at real t = (-3 +- sqrt 41)/2) "
      "the same test reports their resultant's real roots", any(len(e_) == 3 for e_ in hqc), f"{len(hqc)}")

# ---------------------------------------------------------------- G6
t, p, q = sp.symbols("t p q", real=True)


def c_system(fv_):
    s_, t_, r_ = 1 + t, t * (1 + t), -t
    H = dimer(r_, (p, q, s_), (fv_[0], fv_[1], t_))
    U, _ = twin_U(s_, t_, (p, q, 0), fv_)
    Kc = sp.expand(hop(s_, t_) * U); Kd = Kc.H
    P = [sp.eye(4), H, sp.expand(H * H)]; P.append(sp.expand(P[2] * H))
    cols = []
    for k_ in range(4):
        for ph in (1, sp.I):
            Fb = ph * P[k_]; col = []
            for E in (sp.expand(Fb * Kc - Kd * Fb.H), sp.expand(Kc * Fb - Fb.H * Kd)):
                for e in E:
                    col += [sp.expand(sp.re(e)), sp.expand(sp.im(e))]
            cols.append(col)
    M = sp.Matrix(cols).T
    return M[[i for i in range(M.rows) if any(e != 0 for e in M.row(i))], :], P


Mcl, Pc = c_system((p, -q, 0))
check("G6 the chart J = (1+t, t(1+t), -t) has e2 = 0 identically", sp.expand((1 + t) * t * (1 + t) - t * (1 + t + t * (1 + t))) == 0)
nsC = DomainMatrix.from_Matrix(Mcl).convert_to(sp.QQ.frac_field(t, p, q)).nullspace().to_Matrix()
okC = nsC.rows == 1
if okC:
    zc = [sp.together(e) for e in nsC.row(0)]
    dn = sp.lcm([sp.denom(e) for e in zc]); zc = [sp.cancel(e * dn) for e in zc]
    gC = zc[0]
    for e in zc[1:]:
        gC = sp.gcd(gC, e)
    zc = [sp.cancel(e / gC) for e in zc]          # the primitive kernel vector; without this F F^dag is unwieldy
    Fc =sum(((zc[2 * k_] + sp.I * zc[2 * k_ + 1]) * Pc[k_] for k_ in range(4)), sp.zeros(4))
    Gc = sp.expand(Fc * Fc.H)
    okC = all(sp.expand(Gc[i, j]) == 0 for i in range(4) for j in range(4) if i != j) and all(sp.expand(Gc[i, i] - Gc[0, 0]) == 0 for i in range(4))
    sgC = factor_set(Gc[0, 0])
    expC = {sp.expand(p ** 2 * (t - 1) ** 2 + q ** 2 * (t + 1) ** 2), sp.expand(p ** 2 * (t + 1) ** 2 + q ** 2 * (t - 1) ** 2),
            sp.expand(4 * p ** 2 * q ** 2 + p ** 2 * (t ** 2 - 1) ** 2 + q ** 2 * (t + 1) ** 4)}
    extra = {f for f in sgC if not f.is_number} - expC
else:
    extra = {"kernel"}
check("G6 C at every angle: kernel one-dimensional over Q(t, p, q), F F^dag = sigma 1, sigma a product of "
      "p^2(t-1)^2 + q^2(t+1)^2, p^2(t+1)^2 + q^2(t-1)^2, 4p^2q^2 + p^2(t^2-1)^2 + q^2(t+1)^4 (positive off pq = 0) and constants",
      okC and sp.expand(Gc[0, 0]) != 0 and extra <= {p, q}, f"extra {extra}")
Mcc, _ = c_system((sp.Rational(6, 5) * p, -sp.Rational(6, 5) * q, 0))
check("G6 control: on e2 = 0 with |f_v| = 6/5 |f_u| the system has full rank 8 at a random point (no carrier)",
      Mcc.subs({t: sp.Rational(3, 2), p: sp.Rational(5, 3), q: sp.Rational(2, 7)}).rank() == 8)

# ---------------------------------------------------------------- G7
def ends(J, fu_, fv_):
    s_, t_, r_ = [sp.nsimplify(e) for e in J]
    H0 = dimer(r_, (fu_[0], fu_[1], s_), (fv_[0], fv_[1], t_)); H1 = dimer(r_, (fu_[0], fu_[1], -s_), (fv_[0], fv_[1], -t_))
    Vh = hop(s_, t_)
    Hb = {(0, 0): H0, (1, 1): H1, (0, 1): Vh, (1, 0): Vh.H}
    Hbar = {k_: Hb[(1 - k_[0], 1 - k_[1])] for k_ in Hb}

    def nul(L):
        rows = []
        for (x_, y_), Hxy in Hb.items():
            Bk = sp.zeros(16, 32)
            Bk[:, 16 * y_:16 * y_ + 16] += kp(sp.eye(4), L[(x_, y_)]); Bk[:, 16 * x_:16 * x_ + 16] -= kp(Hxy.T, sp.eye(4))
            rows.append(Bk)
        return 32 - sp.Matrix.vstack(*rows).rank(simplify=True)
    return nul(Hb), nul(Hbar)


rowsG7 = [((1, 2, 6), (8, 0), (-5, 3), (1, 1)), ((1, 2, 6), (8, 0), (-5, 4), (1, 0)), ((1, 2, 7), (8, 0), (-5, 3), (1, 0)),
          ((2, -2, 0), (1, 4), (4, 1), (1, 1)), ((2, -2, 0), (1, 4), (4, -1), (1, 0)), ((2, -2, 0), (1, 4), (-1, 4), (1, 0))]
got = [ends(J, f1, f2) for J, f1, f2, _ in rowsG7]
check("G7 the block route (no K, no U): the E row and the corner pair with near = far = 1, their perturbations break",
      got == [e[3] for e in rowsG7], f"{got}")
w0, A0 = sp.Integer(1), sp.Integer(8)
b2v = sp.sqrt(betaE.subs({w: w0, A: A0})); b1v = (cw / A).subs({w: w0, A: A0})
S0, T0_, R0 = Sw.subs(w, w0), Tw.subs(w, w0), Rw.subs(w, w0)
H0 = dimer(R0, (A0, 0, S0), (b1v, b2v, T0_))
U0, _ = twin_U(S0, T0_, (A0, 0, 0), (b1v, b2v, 0))
vv = [e.subs({w: w0, A: A0}) for e in v]
F0 = sum(((vv[j] + sp.I * b2v * vv[4 + j]) * H0 ** j for j in range(4)), sp.zeros(4))
d0 = sp.expand(U0 * F0)
W8 = sp.zeros(8); W8[4:8, 0:4] = d0; W8[0:4, 4:8] = d0.H      # X_d (x) diag(d0, d0^dag): row = complement of column
W8 = W8.applyfunc(sp.expand)
lam0 = sp.expand((W8 * W8)[0, 0])
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}
WE = {"XII": 12, "XIX": 30, "XIY": -10, "XXI": 15, "XXX": -30, "XXY": -2, "XYI": 5, "XYX": -8, "XZZ": -12, "YII": -4,
      "YIX": -10, "YIY": 30, "YXI": 5, "YXY": 24, "YYI": 15, "YYX": 6, "YYY": 10, "YZZ": 4}
str3 = lambda L: kp(kp(PAULI[L[0]], PAULI[L[1]]), PAULI[L[2]])
coef = {}
for L in itertools.product("IXYZ", repeat=3):
    cf = sp.expand((str3(L) * W8).trace() / 8)
    if cf != 0:
        coef["".join(L)] = cf
lamW = sp.cancel(coef.get("XII", 0) / -12)
WEm = sum((v_ * str3(k_) for k_, v_ in WE.items()), sp.zeros(8))
check("G7 at w = 1, A = 8 the locus gives B2^2 = 9, and the carrier built from v is a real multiple of minus the eighteen-string "
      "element of n3_offletter_locus_gate.py G5 (sites d, u, v), the multiple positive, the element squaring to 4500 I; "
      "W^2 = lam0 I with lam0 = 4500 times the multiple squared (exact, no square root of lam0 taken)",
      betaE.subs({w: w0, A: A0}) == 9 and b2v == 3 and lamW.is_real and lamW > 0
      and set(coef) == set(WE) and all(sp.expand(coef[k_] + lamW * WE[k_]) == 0 for k_ in WE)
      and sp.expand(WEm * WEm - 4500 * sp.eye(8)) == sp.zeros(8) and sp.expand(W8 * W8 - lam0 * sp.eye(8)) == sp.zeros(8)
      and sp.expand(lam0 - 4500 * lamW ** 2) == 0, f"{len(coef)} strings")

OUT.append(f"{'ALL PASS' if FAILS == 0 else str(FAILS) + ' FAIL'}  ({time.time() - T0:.0f} s)")
print(OUT[-1])
open("simulations/results/n3_planar_field_locus_gate.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
