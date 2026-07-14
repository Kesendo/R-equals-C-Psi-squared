"""The fragile-thing hunt, scout 5: EXACT proofs of the four finite gates over Q(i).

Everything below is proved as an identity of RATIONAL FUNCTIONS in the unit-circle
monomials z_k = e^{i a_k}, w_k = e^{i b_k} over Q(i).  No floats enter a proof; floats
are only 3-point sanity pins fired before each exact argument.  cot, cos, sin, tan are
the Weierstrass transcriptions:

    cos x = (X + 1/X)/2,   sin x = (X - 1/X)/(2 I),
    cot(theta/2) = I (W + 1)/(W - 1),   tan(theta/2) = -I (W - 1)/(W + 1),   W = e^{i theta}.

The four statements (coordinator brief + addendum):

  S1  cross_form = the 72-atom S5 decomposition, per (i,j) block.
  S2  per canonical sheet L: sum of the 9 event residues = -(prod L) * 1/4 * T(L o x),
      T the 9-term core function, on the sheet {L . x = 0 mod 2pi}.
  S3  oddness: cross_form(a;-b1,b2,b3) = -cross_form(a;b) and cross_form(a;b2,b1,b3) = -cross_form(a;b).
  S4  mirror specialisation: cross_form(a; pi-a3, pi-a2, pi-a1) = 1/2 mirror_form(a).
  S5  oddness: cross_form(a; b1, b2, -b3) = -cross_form(a; b)  (anchor's second Qw-branch).
  S6  oddness: cross_form(a; b3, b2, b1) = -cross_form(a; b)  (transport lemma's relabel).

THE HANG TRAP (documented) AND HOW IT IS AVOIDED.  together()+expand() on the numerator
of a difference of two rational functions with DIFFERENT pole structures explodes even in
3 variables (cot-poles vs cos-poles vs 1/sin-poles never share a common denominator, so
the common numerator has astronomically many terms).  The cure is to make the two sides
share pole structure BEFORE combining, then the check is per-term / per-event / per-block
and each numerator is tiny.  Two half-angle lemmas do the folding (both proved here, 2
variables, instant):

    LEM      cot(mu/2)[cot((mu+rho)/2)+cot((mu-rho)/2)] = 2(1+cos mu)/(cos rho - cos mu)
    tanLEM   tan(mu/2)[tan((mu+rho)/2)+tan((mu-rho)/2)] = 2(1-cos mu)/(cos mu + cos rho)

LEM rewrites cross_form's cot-block into a cos-pole "reduced" block (same poles as the S5
atoms); tanLEM does the same for mirror_form's tan-block.  Every exact check below is then
a small identity: is_zero_rat = "numerator(together(.)) == 0" (gcd-free cross-multiply),
sound because P/dP = Q/dQ  <=>  numer(together(P/dP - Q/dQ)) == 0 for nonzero denominators.
"""
import math
import random
import sys
import time

import sympy as sp

sys.path.insert(0, "simulations")
from cross_triple_orthogonality import cross_form as cf_committed, mirror_form as mf_committed

I = sp.I
z1, z2, z3, w1, w2, w3 = sp.symbols("z1 z2 z3 w1 w2 w3")
Msym, Rsym = sp.symbols("Msym Rsym")
Z = [z1, z2, z3]
W = [w1, w2, w3]
QUARTER = sp.Rational(1, 4)
HALF = sp.Rational(1, 2)


# =====================================================================  monomial algebra
def cosm(X):
    return (X + 1 / X) / 2


def sinm(X):
    return (X - 1 / X) / (2 * I)


def cot_half(E):
    return I * (E + 1) / (E - 1)


def tan_half(E):
    return -I * (E - 1) / (E + 1)


def compl(i):
    return [t for t in range(3) if t != i]


def Esum_a(i):
    j, l = compl(i)
    return Z[j] * Z[l]


def Edif_a(i):
    j, l = compl(i)
    return Z[l] / Z[j]


def alpha_a(i):
    j, l = compl(i)
    return sinm(Z[l]) - sinm(Z[j])


def beta_a(i):
    j, l = compl(i)
    return -(sinm(Z[l]) + sinm(Z[j]))


def Esum_b(j):
    p, q = compl(j)
    return W[p] * W[q]


def Edif_b(j):
    p, q = compl(j)
    return W[q] / W[p]


def alpha_b(j):
    p, q = compl(j)
    return sinm(W[q]) - sinm(W[p])


def beta_b(j):
    p, q = compl(j)
    return -(sinm(W[q]) + sinm(W[p]))


def _a_side(i, xi):
    return (Esum_a(i), alpha_a(i)) if xi == "psum" else (Edif_a(i), beta_a(i))


def _b_side(j, up):
    return (Esum_b(j), alpha_b(j)) if up == "psum" else (Edif_b(j), beta_b(j))


def is_zero_rat(expr):
    """Gcd-free exact zero test for a rational function: numerator(together(expr)) == 0."""
    num, _den = sp.fraction(sp.together(expr))
    return sp.expand(num) == 0


# =====================================================================  the two half-angle lemmas
def prove_LEM():
    """cot(mu/2)[cot((mu+rho)/2)+cot((mu-rho)/2)] = 2(1+cos mu)/(cos rho - cos mu)."""
    lhs = cot_half(Msym) * (cot_half(Msym * Rsym) + cot_half(Msym / Rsym))
    rhs = 2 * (1 + cosm(Msym)) / (cosm(Rsym) - cosm(Msym))
    assert is_zero_rat(lhs - rhs), "LEM failed"


def prove_tanLEM():
    """tan(mu/2)[tan((mu+rho)/2)+tan((mu-rho)/2)] = 2(1-cos mu)/(cos mu + cos rho)."""
    lhs = tan_half(Msym) * (tan_half(Msym * Rsym) + tan_half(Msym / Rsym))
    rhs = 2 * (1 - cosm(Msym)) / (cosm(Msym) + cosm(Rsym))
    assert is_zero_rat(lhs - rhs), "tanLEM failed"


# =====================================================================  symbolic objects
def committed_block_sym(i, j):
    """The (i,j) contribution to the committed cross_form (cot form; carries the outer 1/2).
    Used only for numeric transcription pins; the exact proof uses its LEM-reduced twin."""
    Emp, Emm = Z[i] * W[j], Z[i] / W[j]
    tot = sp.Integer(0)
    for muE, sign in ((Emp, 1), (Emm, -1)):
        Xs = sp.Integer(0)
        for xi in ("psum", "pdif"):
            Exi, cxi = _a_side(i, xi)
            for up in ("psum", "pdif"):
                Eup, cup = _b_side(j, up)
                Xh = QUARTER * (cot_half(muE * Exi / Eup) + cot_half(muE / Exi * Eup)
                                - cot_half(muE * Exi * Eup) - cot_half(muE / Exi / Eup))
                Xs += cxi * cup * Xh
        tot += sign * cot_half(muE) * Xs
    return HALF * sp.Integer((-1) ** (i + j)) * tot


def committed_reduced_block(i, j):
    """committed_block(i,j) folded by LEM: cot(mu/2)*Xh(mu,xi,up) -> the cos-pole form.
    Each fold is the substitution instance of LEM at M = e^{i mu}, R = e^{i(xi -+ up)}."""
    Emp, Emm = Z[i] * W[j], Z[i] / W[j]
    tot = sp.Integer(0)
    for muE, sign in ((Emp, 1), (Emm, -1)):
        cmu = cosm(muE)
        inner = sp.Integer(0)
        for xi in ("psum", "pdif"):
            Exi, cxi = _a_side(i, xi)
            for up in ("psum", "pdif"):
                Eup, cup = _b_side(j, up)
                # cot(mu/2)*Xh = 1/4 [ 2(1+cos mu)/(cos(xi-up)-cos mu) - 2(1+cos mu)/(cos(xi+up)-cos mu) ]
                inner += cxi * cup * QUARTER * (2 * (1 + cmu) / (cosm(Exi / Eup) - cmu)
                                                - 2 * (1 + cmu) / (cosm(Exi * Eup) - cmu))
        tot += sign * inner
    return HALF * sp.Integer((-1) ** (i + j)) * tot


def atom_value_sym(atom):
    """One S5 atom (i, j, xi, up, e) as a cos-pole rational function in the six monomials."""
    i, j, xi, up, e = atom
    Exi, cxi = _a_side(i, xi)
    Eup, cup = _b_side(j, up)
    Ephi = Exi * Eup ** e
    cphi = cosm(Ephi)
    cmp_, cmm = cosm(Z[i] * W[j]), cosm(Z[i] / W[j])
    return (QUARTER * sp.Integer((-1) ** (i + j)) * sp.Integer(-e) * cxi * cup
            * (1 + cphi) * (1 / (cphi - cmp_) - 1 / (cphi - cmm)))


def mf_reduced_block(i, j):
    """mirror_form's (i,j) tan-block folded by tanLEM into cos-pole form."""
    ja, la = compl(i)
    jb, lb = compl(j)
    Ai = dict(Es=Z[ja] * Z[la], Ed=Z[la] / Z[ja],
              al=sinm(Z[la]) - sinm(Z[ja]), be=-(sinm(Z[la]) + sinm(Z[ja])), lap=(-1) ** (i + 1))
    Aj = dict(Es=Z[jb] * Z[lb], Ed=Z[lb] / Z[jb],
              al=sinm(Z[lb]) - sinm(Z[jb]), be=-(sinm(Z[lb]) + sinm(Z[jb])), lap=(-1) ** (j + 1))

    def cval(muE):
        cmu = cosm(muE)
        s = sp.Integer(0)
        for Exi, cxi in ((Ai["Es"], Ai["al"]), (Ai["Ed"], Ai["be"])):
            for Eup, cup in ((Aj["Es"], Aj["al"]), (Aj["Ed"], Aj["be"])):
                # tan(mu/2)*Xi = 1/4 [ 2(1-cos mu)/(cos mu + cos(xi-up)) - 2(1-cos mu)/(cos mu + cos(xi+up)) ]
                s += cxi * cup * QUARTER * (2 * (1 - cmu) / (cmu + cosm(Exi / Eup))
                                            - 2 * (1 - cmu) / (cmu + cosm(Exi * Eup)))
        return s
    return Ai["lap"] * Aj["lap"] * (-cval(Z[i] * Z[j]) + cval(Z[i] / Z[j]))


def mirror_form_sym():
    """Committed mirror_form transcribed to z-monomials (tan form), for the numeric pin."""
    P = []
    for i in range(3):
        j, l = compl(i)
        P.append(dict(Es=Z[j] * Z[l], Ed=Z[l] / Z[j],
                      al=sinm(Z[l]) - sinm(Z[j]), be=-(sinm(Z[l]) + sinm(Z[j])), lap=(-1) ** (i + 1)))

    def Xi(muE, xiE, upE):
        return QUARTER * (tan_half(muE / xiE * upE) + tan_half(muE * xiE / upE)
                          - tan_half(muE / xiE / upE) - tan_half(muE * xiE * upE))

    expr = sp.Integer(0)
    for i in range(3):
        pi = P[i]
        for j in range(3):
            pj = P[j]

            def Xs(muE, pi=pi, pj=pj):
                return (pi["al"] * pj["al"] * Xi(muE, pi["Es"], pj["Es"])
                        + pi["al"] * pj["be"] * Xi(muE, pi["Es"], pj["Ed"])
                        + pi["be"] * pj["al"] * Xi(muE, pi["Ed"], pj["Es"])
                        + pi["be"] * pj["be"] * Xi(muE, pi["Ed"], pj["Ed"]))
            expr += pi["lap"] * pj["lap"] * (-tan_half(Z[i] * Z[j]) * Xs(Z[i] * Z[j])
                                             + tan_half(Z[i] / Z[j]) * Xs(Z[i] / Z[j]))
    return expr


def T_ij_flip_sym(i, j, L):
    """The (i,j) term of T(L o x): (-1)^(i+j) alpha_i(L o a) alpha_j(L o b) cot((L_i a_i + L_j b_j)/2)."""
    sub = {Z[0]: Z[0] ** L[0], Z[1]: Z[1] ** L[1], Z[2]: Z[2] ** L[2],
           W[0]: W[0] ** L[3], W[1]: W[1] ** L[4], W[2]: W[2] ** L[5]}
    return (sp.Integer((-1) ** (i + j)) * alpha_a(i) * alpha_b(j) * cot_half(Z[i] * W[j])).subs(
        sub, simultaneous=True)


def event_weight_sym(atom, s, eps):
    i, j, xi, up, e = atom
    Exi, cxi = _a_side(i, xi)
    Eup, cup = _b_side(j, up)
    Ephi = Exi * Eup ** e
    A = QUARTER * sp.Integer((-1) ** (i + j)) * sp.Integer(-e) * cxi * cup * (1 + cosm(Ephi))
    c_s = sp.Integer(1) if s == 1 else sp.Integer(-1)
    return -A * c_s * sp.Integer(eps) / sinm(Ephi)


# =====================================================================  float mirrors (pins only)
def alpha_f(ang, i):
    j, l = compl(i)
    return math.sin(ang[l]) - math.sin(ang[j])


def T_f(a, b):
    return sum(((-1) ** (i + j)) * alpha_f(a, i) * alpha_f(b, j) / math.tan((a[i] + b[j]) / 2)
               for i in range(3) for j in range(3))


def coeff_f(ang, i, which):
    j, l = compl(i)
    if which == "psum":
        return math.sin(ang[l]) - math.sin(ang[j])
    return -(math.sin(ang[l]) + math.sin(ang[j]))


def xival_f(ang, i, which):
    j, l = compl(i)
    return ang[j] + ang[l] if which == "psum" else ang[l] - ang[j]


def atom_value_f(atom, a, b):
    i, j, xi, up, e = atom
    phi = xival_f(a, i, xi) + e * xival_f(b, j, up)
    A = 0.25 * ((-1) ** (i + j)) * (-e) * coeff_f(a, i, xi) * coeff_f(b, j, up) * (1 + math.cos(phi))
    mp, mm = a[i] + b[j], a[i] - b[j]
    return A * (1.0 / (math.cos(phi) - math.cos(mp)) - 1.0 / (math.cos(phi) - math.cos(mm)))


def event_weight_f(atom, s, tau, eps, a, b):
    i, j, xi, up, e = atom
    phi = xival_f(a, i, xi) + e * xival_f(b, j, up)
    A = 0.25 * ((-1) ** (i + j)) * (-e) * coeff_f(a, i, xi) * coeff_f(b, j, up) * (1 + math.cos(phi))
    c_s = 1 if s == 1 else -1
    return -A * c_s * eps / math.sin(phi)


# =====================================================================  event / sheet combinatorics
def enum_atoms():
    return [(i, j, xi, up, e)
            for i in range(3) for j in range(3)
            for xi in ("psum", "pdif") for up in ("psum", "pdif") for e in (1, -1)]


def phi_vec(atom):
    i, j, xi, up, e = atom
    v = [0] * 6
    ja, la = compl(i)
    if xi == "psum":
        v[ja] += 1
        v[la] += 1
    else:
        v[ja] -= 1
        v[la] += 1
    jb, lb = compl(j)
    if up == "psum":
        v[3 + jb] += e
        v[3 + lb] += e
    else:
        v[3 + jb] -= e
        v[3 + lb] += e
    return v


def enum_events():
    out = []
    for atom in enum_atoms():
        i, j = atom[0], atom[1]
        base = phi_vec(atom)
        for s in (1, -1):
            for tau in (1, -1):
                L = list(base)
                L[i] -= tau
                L[3 + j] -= tau * s
                assert all(x in (-1, 1) for x in L), (atom, s, tau, L)
                out.append((atom, s, tau, tuple(L)))
    return out


def canon(L):
    eps = 1 if L[0] > 0 else -1
    return tuple(eps * x for x in L), eps


def sheets():
    out = {}
    for atom, s, tau, L in enum_events():
        key, eps = canon(L)
        out.setdefault(key, []).append((atom, s, tau, eps))
    return out


def sheet_relation_subst(L):
    """Eliminate w3 via the sheet relation prod_k X_k^{L_k} = 1  (L[0] = +1)."""
    rest = Z[0] ** L[0] * Z[1] ** L[1] * Z[2] ** L[2] * W[0] ** L[3] * W[1] ** L[4]
    return rest ** (-L[5])


FAM_ALL = enum_events()


# =====================================================================  STATEMENT 1
def statement1():
    """S1  cross_form = 1/4 sum_{i,j,xi,up,e} (-1)^(i+j)(-e) c_xi c_up (1+cos phi)
             [ 1/(cos phi - cos(a_i+b_j)) - 1/(cos phi - cos(a_i-b_j)) ],  phi = xi + e up.

    Route: LEM folds the committed cot-block into committed_reduced_block (cos-pole);
    then reduced == atoms is proved PER TERM (each (i,j,xi,up,e) is a 3-monomial identity).
    """
    print("=" * 78)
    print("STATEMENT 1  (S5 atom decomposition, exact)")

    rng = random.Random(101)
    worst = 0.0
    for _ in range(3):
        a = [rng.uniform(-3, 3) for _ in range(3)]
        b = [rng.uniform(-3, 3) for _ in range(3)]
        cf = cf_committed(a, b)
        worst = max(worst, abs(sum(atom_value_f(t, a, b) for t in enum_atoms()) - cf) / max(1.0, abs(cf)))
    print(f"  float pin: |atom_sum - cross_form| worst rel {worst:.2e}")
    assert worst < 1e-9, "S1 float pin FAILED"

    # numeric wiring pin: committed_block_sym (cot transcription) sums to cross_form,
    # and committed_reduced_block sums to the same -- confirms the LEM fold is wired right.
    ang_a = [0.7, 1.9, 2.6]
    ang_b = [1.1, 0.4, 2.2]
    sub = {z1: sp.exp(I * sp.Float(ang_a[0])), z2: sp.exp(I * sp.Float(ang_a[1])),
           z3: sp.exp(I * sp.Float(ang_a[2])), w1: sp.exp(I * sp.Float(ang_b[0])),
           w2: sp.exp(I * sp.Float(ang_b[1])), w3: sp.exp(I * sp.Float(ang_b[2]))}
    cot_sum = complex(sp.N(sum(committed_block_sym(i, j) for i in range(3) for j in range(3)).subs(sub), 25))
    red_sum = complex(sp.N(sum(committed_reduced_block(i, j) for i in range(3) for j in range(3)).subs(sub), 25))
    cf_ref = cf_committed(ang_a, ang_b)
    print(f"  wiring pin: committed {cot_sum.real:.6f}, reduced {red_sum.real:.6f}, cross_form {cf_ref:.6f}")
    assert abs(cot_sum.real - cf_ref) < 1e-8 and abs(red_sum.real - cf_ref) < 1e-8, "S1 wiring pin FAILED"

    t0 = time.time()
    prove_LEM()
    print(f"  EXACT: LEM proved  =>  committed_block == committed_reduced_block (per fold)")

    def reduced_term(i, j, xi, up, e):
        Exi, cxi = _a_side(i, xi)
        Eup, cup = _b_side(j, up)
        cphi = cosm(Exi * Eup ** e)
        cmp_, cmm = cosm(Z[i] * W[j]), cosm(Z[i] / W[j])
        K = QUARTER * sp.Integer((-1) ** (i + j)) * sp.Integer(-e) * cxi * cup
        return K * ((1 + cmp_) / (cphi - cmp_) - (1 + cmm) / (cphi - cmm))

    for atom in enum_atoms():
        i, j, xi, up, e = atom
        assert is_zero_rat(reduced_term(i, j, xi, up, e) - atom_value_sym(atom)), \
            f"S1 term {atom} reduced != atom"
    dt = time.time() - t0
    print(f"  EXACT: all 72 terms  reduced_term == atom  over Q(i)  ({dt:.1f}s)")
    print("  [S1 PASS] cross_form equals the 72-atom S5 decomposition, exactly.")


# =====================================================================  STATEMENT 2
def statement2():
    """S2  per canonical sheet L: sum of the 9 event residues = -(prod L)/4 * T(L o x).

    REFINEMENT proved (stronger than the sum): the 9 events biject to the 9 (i,j) terms of
    T, and PER EVENT  weight = -(prod L)/4 * T_ij(L o x)  on the sheet.  Hence the sum too.
    """
    print("=" * 78)
    print("STATEMENT 2  (per-sheet residue sum = flip-T, exact, all 32 sheets)")
    S = sheets()
    mult = sorted(set(len(v) for v in S.values()))
    print(f"  {len(S)} canonical sheets, events per sheet: {mult}")
    assert len(S) == 32 and mult == [9], "sheet combinatorics unexpected"
    for key, fam in S.items():
        ijs = sorted((a[0], a[1]) for a, s, tau, eps in fam)
        assert ijs == [(i, j) for i in range(3) for j in range(3)], f"sheet {key} not one event per (i,j)"
    print("  every sheet carries exactly one event per (i,j) -> per-event correspondence")

    # ---- float pins at 3 random sheets/points; print the empirical constant once
    rng = random.Random(202)
    pinned = 0
    for key in random.Random(7).sample(list(S), 3):
        pt = _sample_sheet_point(key, rng)
        assert pt is not None, "could not sample a clean sheet point"
        a, b = pt
        lhs = sum(event_weight_f(atom, s, tau, eps, a, b) for atom, s, tau, eps in S[key])
        aL = [key[0] * a[0], key[1] * a[1], key[2] * a[2]]
        bL = [key[3] * b[0], key[4] * b[1], key[5] * b[2]]
        prodL = 1
        for c in key:
            prodL *= c
        rhs = -prodL * 0.25 * T_f(aL, bL)
        if pinned == 0:
            tcore = 0.25 * T_f(aL, bL)
            ratio = lhs / tcore if abs(tcore) > 1e-9 else float("nan")
            print(f"  float pin: sum_w / (1/4 T) = {ratio:+.4f}   (expected -prod L = {-prodL:+d})")
        assert abs(lhs - rhs) < 1e-6 * max(1.0, abs(rhs)), f"S2 float pin FAILED at sheet {key}"
        pinned += 1

    # ---- exact, per event, per sheet
    t0 = time.time()
    for key, fam in sorted(S.items()):
        prodL = 1
        for c in key:
            prodL *= c
        sub_w3 = sheet_relation_subst(key)
        for atom, s, tau, eps in fam:
            i, j = atom[0], atom[1]
            lhs = event_weight_sym(atom, s, eps)
            rhs = -sp.Integer(prodL) * QUARTER * T_ij_flip_sym(i, j, key)
            assert is_zero_rat((lhs - rhs).subs(W[2], sub_w3)), f"S2 sheet {key} event {atom} fails"
    dt = time.time() - t0
    print(f"  EXACT: all 32 sheets x 9 events  weight == -(prod L)/4 * T_ij(L o x)  ({dt:.1f}s)")
    print("  [S2 PASS] every sheet's residue sum is the flipped core T, exactly (term by term).")


def _sample_sheet_point(L, rng, tries=200):
    for _ in range(tries):
        a = [rng.uniform(-3, 3) for _ in range(3)]
        b = [rng.uniform(-3, 3) for _ in range(2)]
        rest = (L[0] * a[0] + L[1] * a[1] + L[2] * a[2] + L[3] * b[0] + L[4] * b[1])
        b = b + [-rest / L[5]]
        ok = True
        for atom, s, tau, eps in FAM_ALL:
            i, j, xi, up, e = atom
            phi = xival_f(a, i, xi) + e * xival_f(b, j, up)
            if abs(math.sin(phi)) < 5e-2:
                ok = False
                break
        if ok:
            return a, b
    return None


# =====================================================================  STATEMENT 3
def _fingerprint(expr, points):
    out = []
    for pt in points:
        val = complex(sp.N(expr.subs(pt, simultaneous=True), 30))
        out.append((round(val.real, 8), round(val.imag, 8)))
    return tuple(out)


def _prove_oddness(label, subst, atoms, V, points):
    """sum_atoms V|subst = - sum_atoms V, via a sign-flipping bijection sigma of the 72 atoms."""
    Vt = {k: V[k].subs(subst, simultaneous=True) for k in atoms}
    table = {}
    for m in atoms:
        table.setdefault(_fingerprint(V[m], points), []).append(m)
    used = set()
    bijection = {}
    for k in atoms:
        target = tuple((-r, -im) for r, im in _fingerprint(Vt[k], points))
        cands = [m for m in table.get(target, []) if m not in used]
        assert cands, f"{label}: no sign-flip partner for atom {k}"
        m = next((c for c in cands if is_zero_rat(Vt[k] + V[c])), None)
        assert m is not None, f"{label}: fingerprint partner(s) {cands} fail the exact check for {k}"
        used.add(m)
        bijection[k] = m
    assert len(used) == len(atoms), f"{label}: not a bijection ({len(used)}/{len(atoms)})"
    return bijection


def statement3():
    """S3  cross_form(a;-b1,b2,b3) = -cross_form(a;b) and cross_form(a;b2,b1,b3) = -cross_form(a;b).

    cross_form = sum of the 72 atoms (S1); for each b-operation we exhibit an exact
    sign-flipping bijection sigma on the atoms with atom|subst = -atom(sigma).  Then
    sum atom|subst = -sum atom(sigma) = -sum atom = -cross_form.
    """
    print("=" * 78)
    print("STATEMENT 3  (oddness under b1->-b1 and b1<->b2, exact bijection)")
    atoms = enum_atoms()
    V = {k: atom_value_sym(k) for k in atoms}

    rng = random.Random(303)
    w_flip = w_swap = 0.0
    for _ in range(3):
        a = [rng.uniform(-3, 3) for _ in range(3)]
        b = [rng.uniform(-3, 3) for _ in range(3)]
        d = max(1.0, abs(cf_committed(a, b)))
        w_flip = max(w_flip, abs(cf_committed(a, [-b[0], b[1], b[2]]) + cf_committed(a, b)) / d)
        w_swap = max(w_swap, abs(cf_committed(a, [b[1], b[0], b[2]]) + cf_committed(a, b)) / d)
    print(f"  float pin: b1->-b1 residual {w_flip:.2e}, b1<->b2 residual {w_swap:.2e}")
    assert w_flip < 1e-9 and w_swap < 1e-9, "S3 float pin FAILED"

    pts = [{z1: 2, z2: 3, z3: 5, w1: 7, w2: 11, w3: 13},
           {z1: 3, z2: 7, z3: 2, w1: 5, w2: 13, w3: 11},
           {z1: 5, z2: 2, z3: 11, w1: 3, w2: 7, w3: 13}]

    t0 = time.time()
    bij_flip = _prove_oddness("b1->-b1", {w1: 1 / w1}, atoms, V, pts)
    bij_swap = _prove_oddness("b1<->b2", {w1: w2, w2: w1}, atoms, V, pts)
    dt = time.time() - t0
    print(f"  EXACT: both sign-flipping bijections verified atom-by-atom ({dt:.1f}s)")

    def show(bij, name):
        fixed = sum(1 for k in bij if bij[k] == k)
        print(f"    {name}: 72 atoms, {fixed} fixed points; sample (atom -> partner):")
        for k in list(bij)[:4]:
            print(f"      {k}  ->  {bij[k]}")
    show(bij_flip, "b1->-b1")
    show(bij_swap, "b1<->b2")
    print("  [S3 PASS] cross_form is odd under both b-operations, exactly.")
    return bij_flip, bij_swap


# =====================================================================  STATEMENT 5
def statement5():
    """S5  cross_form(a; b1, b2, -b3) = -cross_form(a; b), as an identity of functions.

    Load-bearing anchor: the second Qw-branch at the mirror point is reached by b3 -> -b3,
    so this oddness must be proved EXACTLY, not inferred by composing S3's flips.  Same
    method as S3: cross_form = sum of the 72 atoms (S1); exhibit an exact sign-flipping
    bijection sigma with atom|_{w3->1/w3} = -atom(sigma).
    """
    print("=" * 78)
    print("STATEMENT 5  (b3->-b3 oddness, exact bijection)")
    atoms = enum_atoms()
    V = {k: atom_value_sym(k) for k in atoms}

    rng = random.Random(505)
    worst = 0.0
    for _ in range(3):
        a = [rng.uniform(-3, 3) for _ in range(3)]
        b = [rng.uniform(-3, 3) for _ in range(3)]
        worst = max(worst, abs(cf_committed(a, [b[0], b[1], -b[2]]) + cf_committed(a, b))
                    / max(1.0, abs(cf_committed(a, b))))
    print(f"  float pin: b3->-b3 residual {worst:.2e}")
    assert worst < 1e-9, "S5 float pin FAILED"

    pts = [{z1: 2, z2: 3, z3: 5, w1: 7, w2: 11, w3: 13},
           {z1: 3, z2: 7, z3: 2, w1: 5, w2: 13, w3: 11},
           {z1: 5, z2: 2, z3: 11, w1: 3, w2: 7, w3: 13}]

    t0 = time.time()
    bij = _prove_oddness("b3->-b3", {w3: 1 / w3}, atoms, V, pts)
    dt = time.time() - t0
    print(f"  EXACT: sign-flipping bijection verified atom-by-atom ({dt:.1f}s)")
    fixed = sum(1 for k in bij if bij[k] == k)
    print(f"    b3->-b3: 72 atoms, {fixed} fixed points; sample (atom -> partner):")
    for k in list(bij)[:4]:
        print(f"      {k}  ->  {bij[k]}")
    print("  [S5 PASS] cross_form is odd under b3->-b3, exactly.")
    return bij


# =====================================================================  STATEMENT 6
def statement6():
    """S6  cross_form(a; b3, b2, b1) = -cross_form(a; b), as an identity of functions.

    Load-bearing: the transport lemma relabels the 3 sparse a_i + b3 outer poles onto the
    identically-empty a_i + b1 outer poles via exactly this transposition; it is currently
    the only numerically-pinned link, so it must be made exact.  Same method as S3/S5:
    cross_form = sum of the 72 atoms (S1); exhibit an exact sign-flipping bijection sigma
    with atom|_{w1<->w3} = -atom(sigma).
    """
    print("=" * 78)
    print("STATEMENT 6  (b1<->b3 transposition oddness, exact bijection)")
    atoms = enum_atoms()
    V = {k: atom_value_sym(k) for k in atoms}

    rng = random.Random(606)
    worst = 0.0
    for _ in range(3):
        a = [rng.uniform(-3, 3) for _ in range(3)]
        b = [rng.uniform(-3, 3) for _ in range(3)]
        worst = max(worst, abs(cf_committed(a, [b[2], b[1], b[0]]) + cf_committed(a, b))
                    / max(1.0, abs(cf_committed(a, b))))
    print(f"  float pin: b1<->b3 residual {worst:.2e}")
    assert worst < 1e-9, "S6 float pin FAILED"

    pts = [{z1: 2, z2: 3, z3: 5, w1: 7, w2: 11, w3: 13},
           {z1: 3, z2: 7, z3: 2, w1: 5, w2: 13, w3: 11},
           {z1: 5, z2: 2, z3: 11, w1: 3, w2: 7, w3: 13}]

    t0 = time.time()
    bij = _prove_oddness("b1<->b3", {w1: w3, w3: w1}, atoms, V, pts)
    dt = time.time() - t0
    print(f"  EXACT: sign-flipping bijection verified atom-by-atom ({dt:.1f}s)")
    fixed = sum(1 for k in bij if bij[k] == k)
    print(f"    b1<->b3: 72 atoms, {fixed} fixed points; sample (atom -> partner):")
    for k in list(bij)[:4]:
        print(f"      {k}  ->  {bij[k]}")
    print("  [S6 PASS] cross_form is odd under the b1<->b3 transposition, exactly.")
    return bij


# =====================================================================  STATEMENT 4
def statement4():
    """S4  cross_form(a; pi-a3, pi-a2, pi-a1) = 1/2 mirror_form(a).

    In monomials the mirror is w1=-1/z3, w2=-1/z2, w3=-1/z1.  LEM folds cross_form to
    cos-pole blocks, tanLEM folds mirror_form to cos-pole blocks; then the identity holds
    PER BLOCK with the reindexing sigma(i,j)=(i,2-j):
        committed_reduced_block(i,j)|mirror  ==  1/2 mf_reduced_block(i, 2-j).
    Since sigma is a bijection of the 9 blocks, summing gives the whole gate.
    """
    print("=" * 78)
    print("STATEMENT 4  (mirror specialisation, exact, per block)")

    rng = random.Random(404)
    worst = 0.0
    pinned = 0
    while pinned < 3:
        a = [rng.uniform(0.4, 2.7) for _ in range(3)]
        b = [math.pi - a[2], math.pi - a[1], math.pi - a[0]]
        try:
            cf, mf = cf_committed(tuple(a), tuple(b)), mf_committed(tuple(a))
        except ZeroDivisionError:
            continue
        if abs(mf) < 1e-9:
            continue
        worst = max(worst, abs(cf - 0.5 * mf) / abs(mf))
        pinned += 1
    print(f"  float pin: |cross_form(mirror) - 1/2 mirror_form| worst rel {worst:.2e}")
    assert worst < 1e-9, "S4 float pin FAILED"

    # numeric wiring pin: mf_reduced (tanLEM fold) reproduces mirror_form
    ang = (0.7, 1.3, 2.1)
    sub = {z1: sp.exp(I * sp.Float(ang[0])), z2: sp.exp(I * sp.Float(ang[1])), z3: sp.exp(I * sp.Float(ang[2]))}
    red = complex(sp.N(sum(mf_reduced_block(i, j) for i in range(3) for j in range(3)).subs(sub), 25))
    full = complex(sp.N(mirror_form_sym().subs(sub), 25))
    print(f"  wiring pin: mf_reduced {red.real:.5f}, mirror_form_sym {full.real:.5f}, mirror_form {mf_committed(ang):.5f}")
    assert abs(red.real - mf_committed(ang)) < 1e-7 and abs(full.real - mf_committed(ang)) < 1e-7, "S4 wiring pin FAILED"

    mirror = {w1: -1 / z3, w2: -1 / z2, w3: -1 / z1}
    t0 = time.time()
    prove_LEM()
    prove_tanLEM()
    for i in range(3):
        for j in range(3):
            lhs = committed_reduced_block(i, j).subs(mirror, simultaneous=True)
            rhs = HALF * mf_reduced_block(i, 2 - j)
            assert is_zero_rat(lhs - rhs), f"S4 block ({i},{j})->(({i},{2 - j})) fails"
    dt = time.time() - t0
    print(f"  EXACT: LEM + tanLEM + all 9 blocks under sigma(i,j)=(i,2-j)  ({dt:.1f}s)")
    print("  [S4 PASS] cross_form(a; mirror) == 1/2 mirror_form(a), exactly.")


# =====================================================================  driver
def main():
    t0 = time.time()
    statement1()
    statement2()
    statement3()
    statement4()
    statement5()
    statement6()
    print("=" * 78)
    print(f"ALL SIX STATEMENTS PROVED EXACTLY OVER Q(i)   (total {time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()
