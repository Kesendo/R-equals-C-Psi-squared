"""F127: the closed form of the core identity T (the Cauchy-kernel lens find).

Found 2026-07-14 by the borrowing-a-discipline scout (idea 1 of the F127
handover): classical determinant theory owns T's shape.

THE TWO THEOREMS (both proved exactly here, over Q(i)):

[G1] THE CLOSED FORM.  With P = prod_{ij} sin((a_i+b_j)/2),
     Va = prod_{i<j} sin((a_i-a_j)/2), Vb likewise in b,
     s = (sum a + sum b)/2, e1 = sum_i cos a_i, f1 = sum_j cos b_j:

       T(a;b) * P = (1/8) [ 2 cos(s) ( (e1 - f1)^2 - 2 sin(s)^2 )
                             + sin(s) ( sum_i sin 2a_i + sum_j sin 2b_j ) ] * Va * Vb

     COROLLARY 1 (= the F127 core identity, section 3 of
     PROOF_F127_RESIDUE_COLLAPSE.md): on the variety
     {e1 = 0, f1 = 0, s = 0 mod pi} every closed-form term vanishes, so
     T*P = 0 identically on the variety (the pole-free assembly T*P is
     what gate G1 proves; pointwise T = 0 wherever P != 0, and on the
     P = 0 sublocus the pole-free witness of section 3.2 stands, now read
     as the degenerate-denominator case of this same identity).
     SCOPE: this replaces ONLY the section-3 core-identity argument
     (the 540-term assembly + Res|E divisibility).  Sections 1/2/4/5/6
     of the residue-collapse proof (atoms, sheet lattice, transport,
     window, anchor) and the grid+CRT wall are untouched and remain
     load-bearing for F127 itself.
     COROLLARY 2 (sharper, NEW): T = 0 already on {e1 = f1, s = 0 mod pi};
     on {s = 0 mod pi} alone, T = +-(Va Vb / 4P) (e1 - f1)^2, a signed square.
     COROLLARY 3 (the naked falsification explained): on {e1 = f1 = 0} alone
     8 T P/(Va Vb) = -4 cos s sin^2 s + sin s (s2a + s2b) != 0 generically.

[G2] THE STRUCTURAL ANCHOR.  T is a bordered Frobenius determinant:

       T = det [[ C, U ], [ W^T, 0 ]],   C_ij = cot((a_i+b_j)/2),
       U = [1 | sin a_i] (3x2),  W = [1 | sin b_j] (3x2).

     (Laplace expansion along the border columns gives T's (-1)^{i+j}
     alpha_i alpha_j sum; the closed form is the confluent evaluation of the
     classical Frobenius product formula det[cot(x_i+y_j)] =
     cos(sum x + sum y) prod sin(x_i-x_j) sin(y_i-y_j) / prod sin(x_i+y_j),
     which gate [G0] pins numerically.)

Conventions match simulations/f127_core_identity.py exactly:
alpha_i = sin(ang_l) - sin(ang_j), (j,l) = complement of i, ascending.

Run: python simulations/f127_closed_form.py   (~75 s, exits 0 iff all gates pass)
"""
import numpy as np
import sympy as sp

I = sp.I


# ---------- shared numeric pieces (committed convention) ----------
def alpha_num(ang):
    out = np.empty(3)
    for i in range(3):
        j, l = [t for t in range(3) if t != i]
        out[i] = np.sin(ang[l]) - np.sin(ang[j])
    return out


def T_direct(a, b):
    al, be = alpha_num(a), alpha_num(b)
    t = 0.0
    for i in range(3):
        for j in range(3):
            t += (-1) ** (i + j) * al[i] * be[j] / np.tan((a[i] + b[j]) / 2)
    return t


def closed_rhs_over_P(a, b):
    """The closed form's RHS divided by P: G * Va * Vb / P with 8G explicit."""
    x, y = a / 2, b / 2
    P = np.prod(np.sin(x[:, None] + y[None, :]))
    Va = np.prod([np.sin(x[i] - x[j]) for i in range(3) for j in range(i + 1, 3)])
    Vb = np.prod([np.sin(y[i] - y[j]) for i in range(3) for j in range(i + 1, 3)])
    s = (a.sum() + b.sum()) / 2
    e1, f1 = np.cos(a).sum(), np.cos(b).sum()
    G8 = (2 * np.cos(s) * ((e1 - f1) ** 2 - 2 * np.sin(s) ** 2)
          + np.sin(s) * (np.sin(2 * a).sum() + np.sin(2 * b).sum()))
    return (G8 / 8) * Va * Vb / P


# ---------- [G0] Frobenius baseline (numeric pin) ----------
def gate0():
    rng = np.random.default_rng(20260714)
    worst = 0.0
    for _ in range(300):
        a = rng.uniform(0.2, 2 * np.pi - 0.2, 3)
        b = rng.uniform(0.2, 2 * np.pi - 0.2, 3)
        x, y = a / 2, b / 2
        C = 1.0 / np.tan(x[:, None] + y[None, :])
        V = np.prod([np.sin(x[i] - x[j]) for i in range(3) for j in range(i + 1, 3)])
        W = np.prod([np.sin(y[i] - y[j]) for i in range(3) for j in range(i + 1, 3)])
        P = np.prod(np.sin(x[:, None] + y[None, :]))
        rhs = np.cos(x.sum() + y.sum()) * V * W / P
        worst = max(worst, abs(np.linalg.det(C) - rhs) / max(1.0, abs(rhs)))
    assert worst < 1e-11, f"G0 Frobenius baseline: {worst}"
    print(f"[G0] Frobenius product formula, 300 pts: worst rel dev {worst:.2e}  PASS")


# ---------- [G2] bordered determinant == T (exact, Q(i)) ----------
def gate2():
    Z = sp.symbols("zeta1 zeta2 zeta3")
    W_ = sp.symbols("omega1 omega2 omega3")

    def sin_of(m):
        return (m - 1 / m) / (2 * I)

    def cot_of(m):  # cot(theta) with e^{2 i theta} = m^2, i.e. e^{i theta} = m
        return I * (m ** 2 + 1) / (m ** 2 - 1)

    C = sp.Matrix(3, 3, lambda i, j: cot_of(Z[i] * W_[j]))
    U = sp.Matrix(3, 2, lambda i, k: 1 if k == 0 else sin_of(Z[i] ** 2))
    V = sp.Matrix(3, 2, lambda j, k: 1 if k == 0 else sin_of(W_[j] ** 2))
    M = sp.zeros(5, 5)
    M[:3, :3] = C
    M[:3, 3:] = U
    M[3:, :3] = V.T
    detM = M.det(method="berkowitz")

    def alpha_sym(vars_):
        out = []
        for i in range(3):
            j, l = [t for t in range(3) if t != i]
            out.append(sin_of(vars_[l] ** 2) - sin_of(vars_[j] ** 2))
        return out

    alp, bet = alpha_sym(Z), alpha_sym(W_)
    T = sum((-1) ** (i + j) * alp[i] * bet[j] * cot_of(Z[i] * W_[j])
            for i in range(3) for j in range(3))
    nT = len(sp.Add.make_args(sp.expand(sp.together(T).as_numer_denom()[0])))
    assert nT > 50, "vacuous G2 expansion"
    diff = sp.simplify(sp.together(detM - T))
    assert diff == 0, "G2 bordered determinant failed"
    print("[G2] T == det[[C,U],[W^T,0]] exactly over Q(i)  PASS")


# ---------- [G1] the closed form (exact, Q(i)) ----------
def gate1():
    Z = sp.symbols("zeta1 zeta2 zeta3")
    W_ = sp.symbols("omega1 omega2 omega3")

    def sin_of(m):
        return (m - 1 / m) / (2 * I)

    def cos_of(m):
        return (m + 1 / m) / 2

    def alpha_sym(vars_):
        out = []
        for i in range(3):
            j, l = [t for t in range(3) if t != i]
            out.append(sin_of(vars_[l] ** 2) - sin_of(vars_[j] ** 2))
        return out

    alp, bet = alpha_sym(Z), alpha_sym(W_)
    sin_ij = {(i, j): sin_of(Z[i] * W_[j]) for i in range(3) for j in range(3)}
    cos_ij = {(i, j): cos_of(Z[i] * W_[j]) for i in range(3) for j in range(3)}

    LHS = sp.Integer(0)  # T * P, assembled product-form (cot*sin = cos clears poles)
    for i in range(3):
        for j in range(3):
            term = (-1) ** (i + j) * alp[i] * bet[j] * cos_ij[(i, j)]
            for kl, sv in sin_ij.items():
                if kl != (i, j):
                    term *= sv
            LHS += sp.expand(term)

    mu = Z[0] * Z[1] * Z[2] * W_[0] * W_[1] * W_[2]  # e^{i s}
    cs, ss = cos_of(mu), sin_of(mu)
    e1 = sum(cos_of(z ** 2) for z in Z)
    f1 = sum(cos_of(w ** 2) for w in W_)
    s2a = sum(sin_of(z ** 4) for z in Z)
    s2b = sum(sin_of(w ** 4) for w in W_)
    Va = sin_of(Z[0] / Z[1]) * sin_of(Z[0] / Z[2]) * sin_of(Z[1] / Z[2])
    Vb = sin_of(W_[0] / W_[1]) * sin_of(W_[0] / W_[2]) * sin_of(W_[1] / W_[2])
    G8 = 2 * cs * ((e1 - f1) ** 2 - 2 * ss ** 2) + ss * (s2a + s2b)
    RHS = sp.expand(G8 * Va * Vb / 8)

    nL = len(sp.Add.make_args(LHS))
    nR = len(sp.Add.make_args(RHS))
    assert nL > 100 and nR > 100, "vacuous expansion"
    diff = sp.expand(LHS - RHS)
    assert diff == 0, "G1 closed form failed"
    print(f"[G1] T*P == closed form exactly over Q(i) ({nL} = {nR} monomials)  PASS")


# ---------- [G3] numeric gates: agreement + all three corollaries ----------
def gate3():
    rng = np.random.default_rng(2026)
    worst = 0.0
    for _ in range(1000):
        a = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        b = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        worst = max(worst,
                    abs(T_direct(a, b) - closed_rhs_over_P(a, b))
                    / max(1.0, abs(T_direct(a, b))))
    assert worst < 1e-10, f"G3 agreement: {worst}"
    print(f"[G3a] closed form == direct T, 1000 pts: worst rel dev {worst:.2e}  PASS")

    # Corollary 1: exact variety points -> T = 0 (float tolerance)
    hits = worst_v = 0
    while hits < 50:
        a = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        c3 = -(np.cos(a[0]) + np.cos(a[1]))
        if abs(c3) > 0.999:
            continue
        a[2] = np.arccos(c3)  # e1 = 0 exactly (up to float)
        b0, b1 = rng.uniform(0.15, 2 * np.pi - 0.15, 2)
        b2 = (-(a.sum() + b0 + b1)) % (2 * np.pi)  # sheet: sum a + sum b = 0 mod 2pi
        b = np.array([b0, b1, b2])
        if abs(np.cos(b).sum()) > 1e-12:  # need f1 = 0 too: solve instead
            # solve f1 = 0 with the sheet: one equation, adjust b1 by bisection
            def g(t):
                bb2 = (-(a.sum() + b0 + t)) % (2 * np.pi)
                return np.cos(b0) + np.cos(t) + np.cos(bb2)
            lo, hi = 0.15, 2 * np.pi - 0.15
            tv = np.linspace(lo, hi, 400)
            gv = [g(t) for t in tv]
            k = next((k for k in range(399) if gv[k] * gv[k + 1] < 0), None)
            if k is None:
                continue
            lo, hi = tv[k], tv[k + 1]
            for _ in range(80):
                mid = (lo + hi) / 2
                if g(lo) * g(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            b1 = (lo + hi) / 2
            b = np.array([b0, b1, (-(a.sum() + b0 + b1)) % (2 * np.pi)])
        if abs(np.cos(b).sum()) > 1e-10:
            continue
        worst_v = max(worst_v, abs(T_direct(a, b)))
        hits += 1
    assert worst_v < 1e-8, f"G3b variety vanishing: {worst_v}"
    print(f"[G3b] Corollary 1 (F127 core identity), 50 variety pts: "
          f"max |T| = {worst_v:.2e}  PASS")

    # Corollary 2: sharper locus e1 = f1 (nonzero), s = 0 mod 2pi -> T = 0
    hits = worst_s = 0
    while hits < 50:
        a = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        b0 = rng.uniform(0.15, 2 * np.pi - 0.15)

        def h(t):
            bb2 = (-(a.sum() + b0 + t)) % (2 * np.pi)
            return (np.cos(b0) + np.cos(t) + np.cos(bb2)) - np.cos(a).sum()

        tv = np.linspace(0.15, 2 * np.pi - 0.15, 400)
        hv = [h(t) for t in tv]
        k = next((k for k in range(399) if hv[k] * hv[k + 1] < 0), None)
        if k is None:
            continue
        lo, hi = tv[k], tv[k + 1]
        for _ in range(80):
            mid = (lo + hi) / 2
            if h(lo) * h(mid) <= 0:
                hi = mid
            else:
                lo = mid
        b1 = (lo + hi) / 2
        b = np.array([b0, b1, (-(a.sum() + b0 + b1)) % (2 * np.pi)])
        if abs(np.cos(a).sum()) < 0.05:  # want e1 = f1 NONZERO (sharper than F127)
            continue
        worst_s = max(worst_s, abs(T_direct(a, b)))
        hits += 1
    assert worst_s < 1e-8, f"G3c sharper locus: {worst_s}"
    print(f"[G3c] Corollary 2 (SHARPER: e1 = f1 != 0, sheet), 50 pts: "
          f"max |T| = {worst_s:.2e}  PASS")

    # Corollary 3: naked falsification preserved: e1 = f1 = 0 WITHOUT sheet -> O(1)
    big = 0.0
    for _ in range(200):
        a = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        c3 = -(np.cos(a[0]) + np.cos(a[1]))
        if abs(c3) > 0.999:
            continue
        a[2] = np.arccos(c3)
        b = rng.uniform(0.15, 2 * np.pi - 0.15, 3)
        c3 = -(np.cos(b[0]) + np.cos(b[1]))
        if abs(c3) > 0.999:
            continue
        b[2] = np.arccos(c3)
        big = max(big, abs(T_direct(a, b)))
    assert big > 1.0, f"G3d naked non-vanishing lost: {big}"
    print(f"[G3d] Corollary 3 (naked T on e1 = f1 = 0 alone is O(1)): "
          f"max |T| = {big:.2f}  PASS")


if __name__ == "__main__":
    gate0()
    gate2()
    gate1()
    gate3()
    print("ALL GATES PASS: the closed form stands.")
