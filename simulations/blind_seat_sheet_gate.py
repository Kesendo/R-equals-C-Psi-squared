"""Gate: the blind seat's sheet at N = 3 and N = 4 (PROOF_ONE_DEPHASED_SEAT, section on blind dephased
seats, paragraph "The sheet and blindness, at N = 3 and N = 4").

Exact throughout (sympy: Gaussian rationals, with square roots where a half-angle needs one), no threshold. Four parts.

1. N = 3, planar fields of one length, neither parallel nor antiparallel, no field on d. The hand reduction: W = |1><0|_d (x) M + h.c. with
   M = beta·X'X' + 2m·P'_s in the frame where each field points along x' (Sigma = theta_u + theta_v, Delta = theta_u -
   theta_v), and the three equations
      (A')  2i·J_uv·sin(Delta)·p = (J_du - J_dv)·m,                 p = beta - m,
      (C1)  J_dv·(m + p)·e^{-i Sigma} = J_du·conj(m) + J_dv·conj(p)·e^{-i Delta},
      (C2)  J_du·(m + p)·e^{-i Sigma} = J_du·conj(p)·e^{i Delta} + J_dv·conj(m).
   At every row: the solution space of (A'), (C1), (C2) is compared with the far elements of the full commutator
   [H, W] = 0 on the 8-dim space (M unknown on all sixteen entries, and on the pair-parity pattern alone); every
   solution of the hand system is rebuilt as W and checked: [H, W] = 0, {Z_d, W} = 0, W^2 a positive multiple of 1.
   Predicted: one solution exactly when (J_du - J_dv)·e2 = 0, none otherwise; p = 0 on the swap branch J_du = J_dv,
   and, read back from (A'), m = 2i·J_uv·sin(Delta)/(J_du - J_dv)·p on the branch e2 = 0 (a readout, not a test); no far
   element off the pattern. Blindness is read
   at every row from the Krylov rank of seat d as well as from the factor (J_du - J_dv)·e2, and symbolically an
   eigenvector of the Laplacian vanishing at d is (0, J_dv, -J_du) up to a factor, an eigenvector exactly where that
   factor vanishes. The field step is checked in the turned frame, symbolic in the lengths f, g > 0: f·X_u + g·X_v has
   eigenvalues +-f +-g, so off f = g its commutant keeping the pair's parity is span{1, X_uX_v}, and at f = g that
   commutant is three-dimensional with Y_uY_v + Z_uZ_v. The reduction is checked symbolically in the couplings and the half-angles (e^{i theta_u/2} = wu,
   e^{i theta_v/2} = wv as symbols with conj(w) = 1/w): on M = c·1 + beta·X'X' + 2m·P'_s every entry of [H_b, W] is a
   multiple of c plus 0 or +-2 times (A'), (C1), (C2) or a conjugate; the hop equation H_01·d0 = d1·H_10 applied to |00>
   is 2c·(J_dv|01> + J_du|10>), applied to |11> it is (C1) and (C2) plus multiples of c, and the far equation
   H_11·d0 = d0·H_00 is (A') plus a multiple of c on the pair's one-excitation sector and a multiple of c on |00>, |11>;
   and four identities carry (A'), (C1), (C2) into the
   displayed form in U, V (four equations, the middle two the Laplacian of the pair with d struck) and back, dividing by
   sin(Delta/2) alone (reading zeta_prod back from psi2 divides by cos(Delta/2), as the document says). Symbolically on
   e2 = 0, the closed form solves all four displayed equations, and e^{i Sigma/2}·xi_prod is (J_du + J_uv)·e^{i theta_u}
   + J_uv·e^{i theta_v} up to the real sin(Delta/2)·cos(Delta/2); symbolically too, W^2 = (|zeta_swap|^2 + |zeta_prod|^2)·1 + 2(n_swap·n_prod)·SWAP'. The displayed form also holds at every solution, and at rows whose field directions are squares of
   Gaussian rationals (rational half-angles) its solution count is the full commutator's. At those rows the element is
   also built from its closed form at two rows of C (zeta_prod from the struck block's null vector (psi1, psi2) =
   (J_du + 2J_uv, J_du), zeta_swap from the first equation) and checked against the full commutator, its rank across
   d | uv two there, and B's element from the hand system has rank one; the census row's element has twelve strings. Controls on the hypotheses: disconnected rows (two couplings
   zero), one with a field on d and one without it at unequal lengths, have a far element of odd height (u is isolated
   there), so both parts need at most one coupling zero; at parallel fields of one length the pair's singlet gives the
   far element (f_d·sigma_d) (x) P_s at J_du = J_dv with a field on d off the bisector, so part 2 needs the fields
   not parallel.
2. A planar field f_d on d: one axis at d, along f_d. Predicted: the swap class has a far element of odd height exactly
   with f_d on the bisector f_u + f_v, its axis at d along f_d; e2 = 0 couplings outside B with f_d != 0 have none.
   A table of far elements by height (odd, even, all) at the rows TrianglePalindromeLocusWitnessTests pins, rows
   outside the setting among them (parallel and antiparallel fields, the singlet, a disconnected row), is the second
   route for those pins.
3. N = 4: the singlet block. With the pair (1, 2) seen by sites 0 and 3 through its difference only,
   K = b0·s0 + b3·s3 - h·x, K·K = 3(b0^2 + b3^2) + h^2 + 2b0·b3·s0·s3 - 2h(b0·X0 + b3·X3), and the condition
   [S, K·K] = 0 on S in span{X0, X0Y3, Y0, Y0Y3} (read string by string, so at every nonzero b0, b3, h) has no solution for general
   b0, b3, and exactly the line h·X0 + b0·Y0Y3 at b3 = -b0 (symbolic; at h = 0 Y0Y3 solves at every b0, b3, a control
   on the hypothesis h != 0). The pair's block of the Laplacian with 0 and 3 struck has determinant -(b0 + b3)^2, and a
   vector vanishing at 0 and 3 is an eigenvector only as (0, 1, 1, 0) and only at b3 = -b0. For a
   swap-even W = diag(S, T) the singlet entry of [H, W] is g·[Y3, S] (symbolic, on the family; the code embeds
   diag(S, 0), T entering the singlet entry of a block-diagonal W not at all), and on the widened family H's blocks in the
   pair's singlet-triplet basis are -3e + g·Y3, (e + g·Y3) on each triplet state and 2K_alpha between (symbolic), which
   gives the singlet entry g·[Y3, S] and the triplet entries g·[Y3, T_alpha beta]. Rows: the family (one
   far element, swap-even); difference-only coupling with b3 != -b0 (seat 0 blind by another vector, no far element);
   a coupling of site 0 to the pair's sum (seat 0 not blind there, no far element); at g = 0 the difference-only rows
   carry nine dimensions of swap-even far elements, a control on the hypothesis g != 0; with the pair's opposite fields
   along y instead of x, a swap-even far element exists at b3 = -2·b0 (Y^(x)4 is one), a control on the fields' direction.
4. The join: at the letter configuration f_u = h·X_u, f_v = h·Y_v the far element of a blind row is the common carrier
   W(s, t) of the section "A common carrier for the triangle's whole field line", up to a factor, it commutes with
   H_b + h·F at h = 5, and Z^(x)3 splits H into the bonds (kept) and the fields (reversed).

Controls: a sign flipped in (A') makes the hand system lose every C solution, so it would disagree with the full
commutator; a coupling detuned off e2 = 0 gives no solution of either route; unequal field lengths give no far element.

Notation: the code writes m for the document's zeta_swap, p for zeta_prod, Sigma and Delta for theta_plus and theta_minus,
U and V for psi1 and psi2, and (A') for the first displayed equation.

Run: python simulations/blind_seat_sheet_gate.py   (about twenty seconds; writes results/blind_seat_sheet_gate.txt)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")   # the console's default code page cannot print every character below

OUT: list[str] = []
FAILS: list[str] = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    OUT.append(line)


def check(cond, what):
    say(("   PASS  " if cond else "   FAIL  ") + what)
    if not cond:
        FAILS.append(what)


I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
LOW = sp.Matrix([[0, 0], [1, 0]])     # |1><0| on d
UP = LOW.T
SP = sp.Matrix([[0, 1], [0, 0]])      # sigma^+ = |0><1|
SM = SP.T                             # sigma^- = |1><0|


def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.kronecker_product(out, m)
    return out


def op(n, at):
    return kron(*[PAULI[at.get(l, 'I')] for l in range(n)])


def heis(n, bonds, fields):
    h = sp.zeros(2 ** n)
    for (i, j), w in bonds.items():
        if w:
            for c in 'XYZ':
                h += w * op(n, {i: c, j: c})
    for s, f in fields.items():
        for c, x in zip('XYZ', f):
            if x:
                h += x * op(n, {s: c})
    return h


def far_count(H, n, pattern=None, extra=None):
    """Real dimension of the Hermitian far elements W = |1><0|_0 (x) M + h.c. (M on sites 1..n-1) with [H, W] = 0;
    pattern restricts M's entries; extra(E) returns further linear constraints on M (a matrix that must vanish)."""
    m = 2 ** (n - 1)
    pos = [(i, j) for i in range(m) for j in range(m) if pattern is None or pattern(i, j)]
    cols = []
    for (i, j) in pos:
        for unit in (1, sp.I):
            E = sp.zeros(m); E[i, j] = unit
            W = kron(LOW, E) + kron(UP, E.H)
            C = H * W - W * H
            col = [sp.re(x) for x in C] + [sp.im(x) for x in C]
            if extra is not None:
                D = extra(E)
                col += [sp.re(x) for x in D] + [sp.im(x) for x in D]
            cols.append(col)
    A = sp.Matrix(cols).T
    return A.shape[1] - DomainMatrix.from_Matrix(A).convert_to(sp.QQ).rank()


def even_pair(i, j):
    return (i in (0, 3)) == (j in (0, 3))


def blind_seat0(n, bonds):
    L = sp.zeros(n)
    for (i, j), w in bonds.items():
        L[i, i] += w; L[j, j] += w; L[i, j] -= w; L[j, i] -= w
    e = sp.zeros(n, 1); e[0] = 1
    return sp.Matrix.hstack(*[L ** k * e for k in range(n)]).rank() < n


def unit(f):
    z = sp.sympify(f[0]) + sp.I * sp.sympify(f[1])          # exact: the inputs are integers or rationals
    return sp.simplify(z / sp.sqrt(sp.sympify(f[0]) ** 2 + sp.sympify(f[1]) ** 2))


def frame_ops(eu, ev):
    """X'_u X'_v and P'_s on (u, v) for the field directions e^{i theta_u}, e^{i theta_v}."""
    Xu = sp.conjugate(eu) * SP + eu * SM
    Xv = sp.conjugate(ev) * SP + ev * SM
    Yu = -sp.I * (sp.conjugate(eu) * SP - eu * SM)
    Yv = -sp.I * (sp.conjugate(ev) * SP - ev * SM)
    XX = kron(Xu, Xv)
    Ps = (sp.eye(4) - XX - kron(Yu, Yv) - kron(Z, Z)) / 4
    return sp.simplify(XX), sp.simplify(Ps)


def hand_solutions(J, fu, fv, sign=1):
    """Solve (A'), (C1), (C2) over the reals (m, p complex); returns a basis of (m, p) solutions."""
    a, b, r = J
    eu, ev = unit(fu), unit(fv)
    eS = sp.simplify(eu * ev)
    eD = sp.simplify(eu * sp.conjugate(ev))
    sD = sp.im(eD)
    m1, m2, p1, p2 = sp.symbols('m1 m2 p1 p2', real=True)
    m, p = m1 + sp.I * m2, p1 + sp.I * p2
    eqs = [sign * 2 * sp.I * r * sD * p - (a - b) * m,
           b * (m + p) * sp.conjugate(eS) - a * sp.conjugate(m) - b * sp.conjugate(p) * sp.conjugate(eD),
           a * (m + p) * sp.conjugate(eS) - a * sp.conjugate(p) * eD - b * sp.conjugate(m)]
    real_eqs = []
    for e in eqs:
        e = sp.expand(e)
        real_eqs += [sp.re(e), sp.im(e)]
    A, _ = sp.linear_eq_to_matrix(real_eqs, [m1, m2, p1, p2])
    sols = []
    for v in A.nullspace():
        sols.append((v[0] + sp.I * v[1], v[2] + sp.I * v[3]))
    return sols, eu, ev, eS, eD, sD


def displayed_system(J, m, p, eu, ev):
    """The displayed form of the hand system, in x = e^{-i Sigma/2}·m, y = e^{-i Sigma/2}·p = u1 + i·u2,
    U = sin(Delta/2)·u1, V = cos(Delta/2)·u2, Delta = theta_u - theta_v. Returns the four left-hand sides:
    (J_du - J_dv)·m - 2i·J_uv·sin(Delta)·p, J_dv·U + (J_dv + 2J_uv)·V, J_du·U - (J_du + 2J_uv)·V,
    (J_du + J_dv)·Im(x) + 2J_dv·cos(Delta/2)·(U + V). The branch of e^{i Sigma/2} flips x, y, and both half-angle
    factors together, which leaves every line homogeneous and the zero set unchanged."""
    a, b, r = J
    hS = sp.sqrt(sp.simplify(eu * ev))                 # e^{i Sigma/2}, one branch
    hD = sp.simplify(eu / hS)                          # e^{i Delta/2} = e^{i theta_u} / e^{i Sigma/2}
    sD = sp.im(sp.simplify(eu * sp.conjugate(ev)))
    x, y = sp.simplify(m / hS), sp.simplify(p / hS)
    u1, u2 = sp.re(y), sp.im(y)
    U, V = sp.im(hD) * u1, sp.re(hD) * u2
    rows = [(a - b) * m - 2 * sp.I * r * sD * p, b * U + (b + 2 * r) * V, a * U - (a + 2 * r) * V,
            (a + b) * sp.im(x) + 2 * b * sp.re(hD) * (U + V)]
    return [sp.simplify(sp.radsimp(sp.expand(z))) for z in rows]      # exact simplification, no numerical identification


def w_from(m, p, eu, ev):
    XX, Ps = frame_ops(eu, ev)
    M = (m + p) * XX + 2 * m * Ps
    return kron(LOW, M) + kron(UP, M.H)


tri = lambda J: {(0, 1): J[0], (0, 2): J[1], (1, 2): J[2]}
e2 = lambda J: J[0] * J[1] + J[0] * J[2] + J[1] * J[2]

# ---- 1. N = 3, planar fields of one length, no field on d -----------------------------------------------------------
say("1. N = 3: the hand system against the full commutator, planar fields of one length, no field on d")
S_ = sp.Integer
Q_ = sp.Rational
rows1 = [
    # name, J, f_u, f_v (one length), class
    ("C, t = 2", (S_(3), S_(6), S_(-2)), (S_(39), S_(52)), (S_(25), S_(60)), "C"),
    ("C, t = -3", (S_(-2), S_(6), S_(3)), (S_(104), S_(-195)), (S_(204), S_(85)), "C"),
    ("C, t = 1/2", (Q_(3, 2), Q_(3, 4), Q_(-1, 2)), (S_(87), S_(-116)), (S_(105), S_(100)), "C"),
    ("C, t = -5/3", (Q_(-2, 3), Q_(10, 9), Q_(5, 3)), (S_(260), S_(273)), (S_(348), S_(145)), "C"),
    ("C, census row", (S_(3), S_(6), S_(-2)), (S_(3), S_(4)), (S_(5), S_(0)), "C"),
    ("C, letter configuration", (S_(3), S_(6), S_(-2)), (S_(1), S_(0)), (S_(0), S_(1)), "C"),
    ("B, (2, 2, 3)", (S_(2), S_(2), S_(3)), (S_(75), S_(100)), (S_(120), S_(35)), "B"),
    ("B, (1, 1, -5)", (S_(1), S_(1), S_(-5)), (S_(25), S_(-60)), (S_(52), S_(39)), "B"),
    ("B, (3, 3, 1)", (S_(3), S_(3), S_(1)), (S_(25), S_(60)), (S_(-52), S_(39)), "B"),
    ("B and C at once, (2, 2, -1)", (S_(2), S_(2), S_(-1)), (S_(3), S_(4)), (S_(5), S_(0)), "BC"),
    ("not blind, (1, 2, 3)", (S_(1), S_(2), S_(3)), (S_(91), S_(312)), (S_(-300), S_(125)), "-"),
    ("not blind, (2, 5, 1)", (S_(2), S_(5), S_(1)), (S_(260), S_(273)), (S_(348), S_(145)), "-"),
    ("not blind, (1, -3, 4)", (S_(1), S_(-3), S_(4)), (S_(40), S_(75)), (S_(68), S_(51)), "-"),
    ("J_du + J_dv = 0, (2, -2, 1)", (S_(2), S_(-2), S_(1)), (S_(39), S_(52)), (S_(25), S_(60)), "-"),
    ("B at J_uv = 0, (3, 3, 0)", (S_(3), S_(3), S_(0)), (S_(39), S_(52)), (S_(25), S_(60)), "B"),
    ("d at an end, (0, 2, 3)", (S_(0), S_(2), S_(3)), (S_(39), S_(52)), (S_(25), S_(60)), "-"),
    ("J_du + J_dv = 0, (1, -1, 3)", (S_(1), S_(-1), S_(3)), (S_(3), S_(4)), (S_(5), S_(0)), "-"),
    ("control: C detuned, (3, 6, -3)", (S_(3), S_(6), S_(-3)), (S_(39), S_(52)), (S_(25), S_(60)), "-"),
]
for name, J, fu, fv, cls in rows1:
    H = heis(3, tri(J), {1: (fu[0], fu[1], 0), 2: (fv[0], fv[1], 0)})
    sols, eu, ev, eS, eD, sD = hand_solutions(J, fu, fv)
    full = far_count(H, 3)
    patt = far_count(H, 3, even_pair)
    blind = sp.simplify((J[0] - J[1]) * e2(J)) == 0
    say(f"   {name:32s} blind={blind!s:5}  hand system {len(sols)}  far elements: pattern {patt}, all {full}")
    check(blind_seat0(3, tri(J)) == blind, f"{name}: the Krylov rank of seat d reads blind = {blind}, as the factor does")
    pred = 1 if blind else 0            # at J_du = J_dv, (A') forces p = 0, so the two branches meet in one line
    check(len(sols) == patt == full == pred, f"{name}: hand = pattern = all = {pred}")
    for (m, p) in sols:
        W = w_from(m, p, eu, ev)
        Zd = op(3, {0: 'Z'})
        W2 = sp.simplify(W * W)
        ok = sp.simplify(H * W - W * H) == sp.zeros(8) and sp.simplify(Zd * W + W * Zd) == sp.zeros(8) \
            and W2 == W2[0, 0] * sp.eye(8) and sp.simplify(W2[0, 0]) > 0
        check(ok, f"{name}: the rebuilt W commutes with H, anticommutes with Z_d, W^2 = {sp.factor(W2[0, 0])}·1")
        check(all(z == 0 for z in displayed_system(J, m, p, eu, ev)),
              f"{name}: the four displayed equations in U, V hold at the solution")
        if cls in ("B", "BC"):
            check(sp.simplify(p) == 0, f"{name}: swap branch, p = 0")
            check(sp.simplify(sp.im(m * sp.conjugate(sp.sqrt(sp.simplify(eu * ev))))) == 0,
                  f"{name}: m lies on e^(i Sigma/2), d's axis on the bisector of the fields")
        if cls == "C":
            ratio = sp.simplify(m / p)
            check(sp.simplify(ratio - 2 * sp.I * J[2] * sD / (J[0] - J[1])) == 0,
                  f"{name}: readout of (A') at the solution, m/p = {ratio}")
        if name == "C, census row":
            nstr = sum(1 for lab in itertools.product('IXYZ', repeat=3)
                       if sp.simplify((op(3, dict(enumerate(lab))) * W).trace()) != 0)
            check(nstr == 12, f"{name}: the element has {nstr} Pauli strings")

say("   the two middle equations are the Laplacian of the pair with d struck, applied to (U - V, U + V)")
sa, sb, sr, sU, sV = sp.symbols('a b r U V')
Lstruck = sp.Matrix([[sa + sr, -sr], [-sr, sb + sr]])
lhs = Lstruck * sp.Matrix([sU - sV, sU + sV])
check(sp.expand(lhs[0] - (sa * sU - (sa + 2 * sr) * sV)) == 0 and sp.expand(lhs[1] - (sb * sU + (sb + 2 * sr) * sV)) == 0
      and sp.expand(Lstruck.det() - (sa * sb + sa * sr + sb * sr)) == 0,
      "[[J_du + J_uv, -J_uv], [-J_uv, J_dv + J_uv]]·(U - V, U + V) is the pair (J_du·U - (J_du + 2J_uv)·V, "
      "J_dv·U + (J_dv + 2J_uv)·V), of determinant e2")
Ja, Jb, Jr = sp.symbols('J_du J_dv J_uv', real=True)
Lsym = sp.Matrix([[Ja + Jb, -Ja, -Jb], [-Ja, Ja + Jr, -Jr], [-Jb, -Jr, Jb + Jr]])
wd = sp.Matrix([0, Jb, -Ja])            # d's row of L forces (y_u, y_v) proportional to (J_dv, -J_du) on a vector vanishing at d
Lwd = (Lsym * wd).applyfunc(sp.expand)
check(Lwd[0] == 0 and sp.expand(Lwd[1] * wd[2] - Lwd[2] * wd[1] + (Ja - Jb) * (Ja * Jb + Ja * Jr + Jb * Jr)) == 0,
      "(0, J_dv, -J_du) spans the vectors vanishing at d that d's row of L sends to zero (J_du, J_dv not both zero), "
      "and it is an eigenvector of L exactly where (J_du - J_dv)·e2 = 0")

say("   the field step in the turned frame, an exact relabelling of the letters: f·X_u + g·X_v with lengths f, g > 0")
fl, gl = sp.symbols('f g', positive=True)
lam_ = sp.Symbol('lambda')
Fp = fl * kron(X, I2) + gl * kron(I2, X)
PARITY = [('I', 'I'), ('Z', 'I'), ('I', 'Z'), ('Z', 'Z'), ('X', 'X'), ('X', 'Y'), ('Y', 'X'), ('Y', 'Y')]


def parity_commutant(Fm):
    """The operators in the span of the pair-parity-keeping strings that commute with Fm."""
    basis = [kron(PAULI[p_], PAULI[q_]) for p_, q_ in PARITY]
    A_ = sp.Matrix([list((Fm * B_ - B_ * Fm).applyfunc(sp.expand)) for B_ in basis]).T
    return [sum((v[k] * basis[k] for k in range(8)), sp.zeros(4)) for v in A_.nullspace()]


four = [sp.eye(4), kron(X, I2), kron(I2, X), kron(X, X)]
gen_comm = parity_commutant(Fp)
check(sp.expand((Fp - lam_ * sp.eye(4)).det() - (lam_ ** 2 - (fl + gl) ** 2) * (lam_ ** 2 - (fl - gl) ** 2)) == 0
      and all((Fp * B_ - B_ * Fp) == sp.zeros(4) for B_ in four)
      and len(gen_comm) == 2 and sp.Matrix([list(B_) for B_ in gen_comm + [sp.eye(4), kron(X, X)]]).rank() == 2,
      "eigenvalues +-f +-g, simple off f = g, so there the commutant is span{1, X_u, X_v, X_uX_v}; its parity-keeping "
      "part, computed over Q(f, g), is span{1, X_uX_v}")
eq_comm = parity_commutant(Fp.subs(gl, fl))
YYZZ = kron(Y, Y) + kron(Z, Z)
check(len(eq_comm) == 3 and sp.Matrix([list(B_) for B_ in eq_comm + [YYZZ]]).rank() == 3,
      "at f = g the parity-keeping commutant is three-dimensional and holds Y_uY_v + Z_uZ_v")

say("   the reduction, symbolic in the couplings and the half-angles e^(i theta_u/2) = wu, e^(i theta_v/2) = wv")
wu, wv = sp.symbols('wu wv')          # no assumptions: they stand for e^{i theta/2}, conjugated by BAR alone
mS, mSb, mP, mPb, cS, cSb = sp.symbols('zeta_swap zeta_swap_bar zeta_prod zeta_prod_bar c c_bar')
BAR = {sp.I: -sp.I, wu: 1 / wu, wv: 1 / wv, mS: mSb, mSb: mS, mP: mPb, mPb: mP, cS: cSb, cSb: cS}


def bar(z):
    """Complex conjugation on these symbols: |wu| = |wv| = 1, the couplings real, each unknown paired with its bar."""
    return sp.expand(sp.sympify(z).xreplace(BAR))


Hsym = sp.zeros(8)
for (i, j), w in {(0, 1): Ja, (0, 2): Jb, (1, 2): Jr}.items():
    for ch in 'XYZ':
        Hsym += w * op(3, {i: ch, j: ch})
eu_s, ev_s = wu ** 2, wv ** 2
XXs = kron(SP / eu_s + eu_s * SM, SP / ev_s + ev_s * SM)
Ps_s = (sp.eye(4) - XXs - kron(-sp.I * (SP / eu_s - eu_s * SM), -sp.I * (SP / ev_s - ev_s * SM)) - kron(Z, Z)) / 4
M_s = (cS * sp.eye(4) + (mS + mP) * XXs + 2 * mS * Ps_s).applyfunc(sp.expand)
W_s = kron(LOW, M_s) + kron(UP, M_s.T.applyfunc(bar))
C_s = (Hsym * W_s - W_s * Hsym).applyfunc(sp.expand)
eS_s, eD_s = eu_s * ev_s, eu_s / ev_s
sinD_s = sp.expand((eD_s - 1 / eD_s) / (2 * sp.I))
forms = {"(A')": sp.expand(2 * sp.I * Jr * sinD_s * mP - (Ja - Jb) * mS),
         "(C1)": sp.expand(Jb * (mS + mP) / eS_s - Ja * mSb - Jb * mPb / eD_s),
         "(C2)": sp.expand(Ja * (mS + mP) / eS_s - Ja * mPb * eD_s - Jb * mSb)}
forms.update({k + "*": bar(v) for k, v in list(forms.items())})
ZETA0 = {mS: 0, mSb: 0, mP: 0, mPb: 0}


def form_of(zx):
    """(name, factor) with zx = factor·form, factor = +-2, or None."""
    for k, f in forms.items():
        for sgn in (2, -2):
            if sp.expand(zx - sgn * f) == 0:
                return k, sgn
    return None


ent = {}
for i in range(8):
    for j in range(8):
        zc = sp.expand(C_s[i, j].subs(ZETA0))
        zx = sp.expand(C_s[i, j] - zc)
        ent[(i, j)] = (zc, zx, form_of(zx) if zx != 0 else None)
check(all(zx == 0 or f is not None for (_, zx, f) in ent.values()),
      "every entry of [H_b, W] on M = c·1 + beta·X'X' + 2m·P'_s is a multiple of c plus 0 or +-2 times (A'), (C1), (C2) or a conjugate")
check(sp.expand(ent[(1, 0)][0] - 2 * Jb * cS) == 0 and sp.expand(ent[(2, 0)][0] - 2 * Ja * cS) == 0
      and all(ent[(k, 0)][1] == 0 for k in range(4)) and ent[(0, 0)][0] == 0 and ent[(3, 0)][0] == 0,
      "the hop equation H_01·d0 = d1·H_10 applied to |00> is 2c·(J_dv|01> + J_du|10>): c = 0 unless J_du = J_dv = 0")
check(ent[(1, 3)][2] == ("(C1)", 2) and ent[(2, 3)][2] == ("(C2)", 2),
      "applied to |11> it is (C1) on <01| and (C2) on <10|")
hop_blocks = [(i, j) for i in range(4) for j in range(4)] + [(i, j) for i in range(4, 8) for j in range(4, 8)]
check(all(ent[k][2] is None or ent[k][2][0] in ("(C1)", "(C2)", "(C1)*", "(C2)*") for k in hop_blocks),
      "every other entry of the two hop equations is, at c = 0, +-2 times (C1), (C2) or a conjugate, or zero")
far_odd = [(4 + p, q) for p in (1, 2) for q in (1, 2)]
far_even = [(4 + p, q) for p in (0, 3) for q in (0, 3)]
far_mixed = [(4 + p, q) for p in range(4) for q in range(4) if (p in (0, 3)) != (q in (0, 3))]
check(ent[(5, 1)][2] == ("(A')", 2) and ent[(6, 2)][2] == ("(A')", -2) and ent[(5, 2)][1] == 0 == ent[(6, 1)][1]
      and all(ent[k][1] == 0 for k in far_even) and all(ent[k][0] == 0 == ent[k][1] for k in far_mixed),
      "the far equation H_11·d0 = d0·H_00 is (A') on the pair's one-excitation sector and a multiple of c on |00>, |11>")
x1, x2, u1, u2 = sp.symbols('x1 x2 u1 u2', real=True)             # xi_swap = x1 + i·x2, xi_prod = u1 + i·u2
hS_s, hD_s = wu * wv, wu / wv                                      # e^{i Sigma/2}, e^{i Delta/2}
XI = {mS: hS_s * (x1 + sp.I * x2), mSb: (x1 - sp.I * x2) / hS_s, mP: hS_s * (u1 + sp.I * u2), mPb: (u1 - sp.I * u2) / hS_s}
re_ = lambda z: sp.expand((z + bar(z)) / 2)
im_ = lambda z: sp.expand((z - bar(z)) / (2 * sp.I))
xs, xsb, xp, xpb = x1 + sp.I * x2, x1 - sp.I * x2, u1 + sp.I * u2, u1 - sp.I * u2
h1 = sp.expand(hS_s * forms["(C1)"].subs(XI))
h2 = sp.expand(hS_s * forms["(C2)"].subs(XI))
first = sp.expand((Ja - Jb) * xs - 2 * sp.I * Jr * sinD_s * xp)    # e^{-i Sigma/2} times the first equation, left minus right
s2, c2 = sp.expand((hD_s - 1 / hD_s) / (2 * sp.I)), sp.expand((hD_s + 1 / hD_s) / 2)
psi1, psi2 = s2 * u1, c2 * u2
eq2 = Ja * psi1 - (Ja + 2 * Jr) * psi2
eq3 = Jb * psi1 + (Jb + 2 * Jr) * psi2
eq4 = (Ja + Jb) * x2 + 2 * Jb * c2 * (psi1 + psi2)
zero = lambda z: sp.cancel(sp.expand(z)) == 0
check(zero(h1 - (Jb * (xs + xp) - Ja * xsb - Jb * xpb / eD_s)) and zero(h2 - (Ja * (xs + xp) - Jb * xsb - Ja * xpb * eD_s))
      and zero(forms["(A')"].subs(XI) + hS_s * first) and zero(sinD_s - 2 * s2 * c2),
      "times e^{i Sigma/2}, (C1) and (C2) are J_dv·(xi_s + xi_p) = J_du·conj(xi_s) + J_dv·e^{-i Delta}·conj(xi_p) and "
      "J_du·(xi_s + xi_p) = J_dv·conj(xi_s) + J_du·e^{i Delta}·conj(xi_p); (A') is -e^{i Sigma/2} times the first equation")
check(zero(im_(h1) - eq4) and zero(im_(h2) - (eq4 - 2 * c2 * (eq2 + eq3)))
      and zero(re_(h1) + re_(first) - 2 * s2 * eq3) and zero(re_(h2) - re_(first) - 2 * s2 * eq2),
      "Im of the first is the last equation, Im of the second the last minus 2cos(Delta/2)·(second + third); Re of the "
      "first plus, and Re of the second minus, Re of the first equation are 2sin(Delta/2) times the third and the second")
xp_gen = (Ja + 2 * Jr) / s2 + sp.I * Ja / c2                       # xi_prod from (psi1, psi2) = (J_du + 2J_uv, J_du)
check(zero(hS_s * xp_gen * s2 * c2 - ((Ja + Jr) * wu ** 2 + Jr * wv ** 2)),
      "the null vector read back: e^{i Sigma/2}·xi_prod is (J_du + J_uv)·e^{i theta_u} + J_uv·e^{i theta_v} times the real "
      "sin(Delta/2)·cos(Delta/2)")
r_C = -Ja * Jb / (Ja + Jb)                                          # e2 = 0
xp_C = xp_gen.subs(Jr, r_C)
xs_C = 2 * sp.I * r_C * sinD_s * xp_C / (Ja - Jb)                   # the first equation
on_C = {x1: re_(xs_C), x2: im_(xs_C), u1: re_(xp_C), u2: im_(xp_C), Jr: r_C}
check(all(zero(sp.together(z_.subs(on_C))) for z_ in (re_(first), im_(first), eq2, eq3, eq4)),
      "on e2 = 0, symbolically in J_du, J_dv and the half-angles: the closed form solves all four displayed equations")
W0 = W_s.subs({cS: 0, cSb: 0})
SWAPs = sp.eye(4) - 2 * Ps_s
check(((W0 * W0).applyfunc(sp.expand) - (mS * mSb + mP * mPb) * sp.eye(8)
       - (mS * mPb + mSb * mP) * kron(sp.eye(2), SWAPs)).applyfunc(lambda z_: sp.cancel(sp.expand(z_))) == sp.zeros(8),
      "W^2 = (|zeta_swap|^2 + |zeta_prod|^2)·1 + 2(n_swap·n_prod)·SWAP'_uv: invertible wherever n_swap is perpendicular to n_prod")
check(not zero(im_(h2) - eq4) and not zero(re_(h1) - 2 * s2 * eq3)
      and sp.expand(ent[(1, 3)][1] - 2 * (Ja * (mS + mP) / eS_s - Jb * mSb - Ja * mPb / eD_s)) != 0,
      "controls: without the corrections Im of the second is not the last equation and Re of the first not 2sin(Delta/2) "
      "times the third; (C1) with J_du and J_dv exchanged is not the hop equation's entry on <01|, |11>")
say("   the displayed equations against the full commutator at rows: field directions zu^2, zv^2 with rational half-angles")


def displayed_count(J, zu, zv, mutate=False):
    """Solution count of the displayed four equations alone, in x = e^{-i Sigma/2}·m and u1 + i·u2 = e^{-i Sigma/2}·p,
    for field directions e^{i theta_u} = zu^2, e^{i theta_v} = zv^2 (e^{i Sigma/2} = zu·zv, e^{i Delta/2} = zu·conj(zv))."""
    a, b, r = J
    hS, hD = zu * zv, sp.expand(zu * sp.conjugate(zv))
    sD = sp.im(sp.expand(hD ** 2))
    x1, x2, u1, u2 = sp.symbols('x1 x2 u1 u2', real=True)
    m, p = hS * (x1 + sp.I * x2), hS * (u1 + sp.I * u2)
    U, V = sp.im(hD) * u1, sp.re(hD) * u2
    first = sp.expand((a - b) * m - 2 * sp.I * r * sD * p)
    last_coef = a if mutate else b                    # the control writes J_du for J_dv in the last equation
    eqs = [sp.re(first), sp.im(first), b * U + (b + 2 * r) * V, a * U - (a + 2 * r) * V,
           (a + b) * x2 + 2 * last_coef * sp.re(hD) * (U + V)]
    A_, _ = sp.linear_eq_to_matrix(eqs, [x1, x2, u1, u2])
    return 4 - A_.rank()


zu, zv = (3 + 4 * sp.I) / 5, (5 + 12 * sp.I) / 13
fu_sq = tuple(4225 * c for c in (sp.re(sp.expand(zu ** 2)), sp.im(sp.expand(zu ** 2))))     # one length, 4225
fv_sq = tuple(4225 * c for c in (sp.re(sp.expand(zv ** 2)), sp.im(sp.expand(zv ** 2))))
for name, J in (("C, (3, 6, -2)", (S_(3), S_(6), S_(-2))), ("C, (-2, 6, 3)", (S_(-2), S_(6), S_(3))),
                ("B, (2, 2, 3)", (S_(2), S_(2), S_(3))), ("not blind, (1, 2, 3)", (S_(1), S_(2), S_(3)))):
    H = heis(3, tri(J), {1: (fu_sq[0], fu_sq[1], 0), 2: (fv_sq[0], fv_sq[1], 0)})
    dc = displayed_count(J, zu, zv)
    hc = len(hand_solutions(J, fu_sq, fv_sq)[0])
    pc, ac = far_count(H, 3, even_pair), far_count(H, 3)
    say(f"   {name:22s} f_u = {fu_sq}, f_v = {fv_sq}: displayed {dc}, hand system {hc}, far elements of odd height {pc}, all {ac}")
    check(dc == hc == pc == ac, f"{name}: displayed system, hand system and full commutator agree, the far end of dimension {ac}")
zu2, zv2 = (5 + 12 * sp.I) / 13, (8 + 15 * sp.I) / 17
fu2 = tuple(48841 * c for c in (sp.re(sp.expand(zu2 ** 2)), sp.im(sp.expand(zu2 ** 2))))   # one length, 13^2·17^2
fv2 = tuple(48841 * c for c in (sp.re(sp.expand(zv2 ** 2)), sp.im(sp.expand(zv2 ** 2))))
for name, J in (("C, (3, 6, -2)", (S_(3), S_(6), S_(-2))), ("B, (1, 1, -5)", (S_(1), S_(1), S_(-5))),
                ("not blind, (2, 5, 1)", (S_(2), S_(5), S_(1)))):
    H = heis(3, tri(J), {1: (fu2[0], fu2[1], 0), 2: (fv2[0], fv2[1], 0)})
    dc, pc, ac = displayed_count(J, zu2, zv2), far_count(H, 3, even_pair), far_count(H, 3)
    say(f"   {name:22s} second angle pair: displayed {dc}, far elements of odd height {pc}, all {ac}")
    check(dc == pc == ac, f"{name}, second angle pair: displayed system and full commutator agree, the far end of dimension {ac}")
check(displayed_count((S_(3), S_(6), S_(-2)), zu, zv, mutate=True) != displayed_count((S_(3), S_(6), S_(-2)), zu, zv),
      "control: J_du written for J_dv in the last displayed equation changes the count at C, (3, 6, -2)")
say("   the closed form at those rows, and the rank across d | uv")


def schmidt_d_uv(W):
    """Operator Schmidt rank of an 8 x 8 operator across d | uv (d the first factor)."""
    R = sp.zeros(4, 16)
    for i in range(8):
        for j in range(8):
            R[(i >> 2) * 2 + (j >> 2), (i & 3) * 4 + (j & 3)] = W[i, j]
    return R.rank()


hS_, hD_ = zu * zv, sp.expand(zu * sp.conjugate(zv))
for name, J in (("C, (3, 6, -2)", (S_(3), S_(6), S_(-2))), ("C, (-2, 6, 3)", (S_(-2), S_(6), S_(3)))):
    a_, b_, r_ = J
    sD_ = sp.im(sp.expand(hD_ ** 2))
    p_cf = hS_ * ((a_ + 2 * r_) / sp.im(hD_) + sp.I * a_ / sp.re(hD_))
    m_cf = 2 * sp.I * r_ * sD_ * p_cf / (a_ - b_)
    eu_, ev_ = sp.expand(zu ** 2), sp.expand(zv ** 2)
    W = w_from(m_cf, p_cf, eu_, ev_)
    H = heis(3, tri(J), {1: (fu_sq[0], fu_sq[1], 0), 2: (fv_sq[0], fv_sq[1], 0)})
    Zd = op(3, {0: 'Z'})
    ok = sp.simplify(H * W - W * H) == sp.zeros(8) and sp.simplify(Zd * W + W * Zd) == sp.zeros(8)
    check(ok, f"{name}: the closed form (zeta_prod from (psi1, psi2) = (J_du + 2J_uv, J_du), zeta_swap from the first equation) is a far element")
    check(schmidt_d_uv(W) == 2, f"{name}: rank two across d | uv")
    check(far_count(H, 3) == 1, f"{name}: the far end there is one-dimensional")
sols_b, eu_b, ev_b, *_ = hand_solutions((S_(2), S_(2), S_(3)), fu_sq, fv_sq)
check(schmidt_d_uv(w_from(*sols_b[0], eu_b, ev_b)) == 1, "B, (2, 2, 3): rank one across d | uv (the swap times a product)")
say("   controls on the hypotheses")
H = heis(3, tri((0, 1, 0)), {0: (1, 0, 0), 1: (3, 4, 0), 2: (5, 0, 0)})
check(far_count(H, 3, even_pair) == 1, "disconnected (0, 1, 0) with f_d = (1, 0): a far element of odd height, so part 2 needs at most one coupling zero")
H = heis(3, tri((0, 1, 0)), {1: (3, 4, 0), 2: (2, 0, 0)})
check(far_count(H, 3, even_pair) == 1,
      "disconnected (0, 1, 0), no field on d, unequal lengths |f_u| = 5, |f_v| = 2: a far element of odd height, so part 1 needs it too")
Hs = heis(3, tri((2, 2, 3)), {0: (1, 2, 0), 1: (3, 1, 0), 2: (3, 1, 0)})
Ps_ = (sp.eye(4) - kron(X, X) - kron(Y, Y) - kron(Z, Z)) / 4
Wsing = kron(1 * X + 2 * Y, Ps_)
check(sp.simplify(Hs * Wsing - Wsing * Hs) == sp.zeros(8) and sp.simplify(op(3, {0: 'Z'}) * Wsing + Wsing * op(3, {0: 'Z'})) == sp.zeros(8),
      "parallel fields of one length, (2, 2, 3), f_u = f_v = (3, 1), f_d = (1, 2): (f_d·sigma_d) (x) P_s is a far element of odd height")
say("   control: unequal lengths at e2 = 0, (3, 6, -2), f_u = (39, 52), f_v = (50, 120)")
H = heis(3, tri((3, 6, -2)), {1: (39, 52, 0), 2: (50, 120, 0)})
check(far_count(H, 3) == 0 and far_count(H, 3, even_pair) == 0, "unequal lengths: no far element")
H = heis(3, tri((2, 2, 3)), {1: (39, 52, 0), 2: (50, 120, 0)})
check(far_count(H, 3) == 0, "unequal lengths at B, (2, 2, 3): no far element")
say("   control: the sign of (A') flipped")
bad = 0
for name, J, fu, fv, cls in rows1[:4]:
    sols, *_ = hand_solutions(J, fu, fv, sign=-1)
    bad += len(sols)
check(bad == 0, "with (A') sign-flipped the hand system loses every C solution (it would disagree with the full commutator)")

# ---- 2. a planar field on d ------------------------------------------------------------------------------------------
say("2. N = 3 with a planar field on d: one axis at d, along f_d")
rows2 = [
    ("B, f_d on the bisector", (2, 2, 3), (-9, 33), (25, 60), (-52, 39), 1),
    ("B, f_d on the bisector, (1, 1, -5)", (1, 1, -5), (77, -21), (25, -60), (52, 39), 1),
    ("B, f_d across the bisector", (2, 2, 3), (77, 21), (25, 60), (-52, 39), 0),
    ("B, f_d generic", (2, 2, 3), (3, 1), (25, 60), (-52, 39), 0),
    ("C, e2 = 0, f_d generic", (3, 6, -2), (1, 2), (39, 52), (25, 60), 0),
    ("C, e2 = 0, f_d on the bisector", (3, 6, -2), (64, 112), (39, 52), (25, 60), 0),
]
for name, J, fd, fu, fv, pred in rows2:
    H = heis(3, tri(J), {0: (fd[0], fd[1], 0), 1: (fu[0], fu[1], 0), 2: (fv[0], fv[1], 0)})
    patt, allc = far_count(H, 3, even_pair), far_count(H, 3)
    om = fd[0] + sp.I * fd[1]
    one_letter = far_count(H, 3, even_pair, extra=lambda E: sp.conjugate(om) * E - om * E.H) if patt else 0
    say(f"   {name:36s} pattern {patt}, all {allc}, of them with the axis at d along f_d: {one_letter}")
    check(patt == pred == allc and one_letter == patt, f"{name}: predicted {pred}, the far end of that dimension")

say("   far elements by height (odd, even, all) at the rows TrianglePalindromeLocusWitnessTests pins")
rows_h = [
    ("C", (3, 6, -2), (0, 0), (3, 4), (5, 0), (1, 0, 1)), ("C, unequal lengths", (3, 6, -2), (0, 0), (3, 4), (4, 4), (0, 0, 0)),
    ("B", (2, 2, 3), (0, 0), (3, 4), (0, 5), (1, 0, 1)), ("not blind", (1, 2, 3), (0, 0), (3, 4), (5, 0), (0, 0, 0)),
    ("E", (1, 2, 6), (0, 0), (8, 0), (-5, 3), (0, 0, 1)), ("d at an end", (0, 2, 3), (0, 0), (3, 4), (5, 0), (0, 0, 0)),
    ("B at J_uv = 0", (3, 3, 0), (0, 0), (3, 4), (5, 0), (1, 0, 1)), ("B, f_d on the line", (2, 2, 3), (1, 1), (3, 4), (4, 3), (1, 0, 1)),
    ("B, f_d off it", (2, 2, 3), (1, 0), (3, 4), (4, 3), (0, 0, 0)), ("C, f_d on the bisector", (3, 6, -2), (8, 4), (3, 4), (5, 0), (0, 0, 0)),
    ("D", (1, 0, 3), (1, 0), (-1, -2), (0, 2), (0, 0, 1)), ("B (1, 1, -5), f_d on the line", (1, 1, -5), (1, 1), (3, 4), (4, 3), (1, 0, 1)),
    ("B and C at once, (2, 2, -1)", (2, 2, -1), (0, 0), (3, 4), (5, 0), (1, 0, 1)),
    ("C (-2, 6, 3) at another angle", (-2, 6, 3), (0, 0), (1, 2), (2, 1), (1, 0, 1)),
    ("outside: parallel fields", (1, 2, 3), (0, 0), (1, 2), (2, 4), (1, 0, 1)),
    ("outside: antiparallel fields", (1, 2, 3), (0, 0), (1, 0), (-1, 0), (1, 0, 1)),
    ("outside: the singlet", (2, 2, 3), (1, 2), (3, 1), (3, 1), (1, 0, 1)),
    ("outside: disconnected", (0, 1, 0), (0, 0), (3, 4), (2, 0), (1, 1, 2)),
]
for name, J, fd, fu, fv, want in rows_h:
    H = heis(3, tri(J), {0: (fd[0], fd[1], 0), 1: (fu[0], fu[1], 0), 2: (fv[0], fv[1], 0)})
    got = (far_count(H, 3, even_pair), far_count(H, 3, lambda i, j: not even_pair(i, j)), far_count(H, 3))
    check(got == want, f"{name}: odd {got[0]}, even {got[1]}, all {got[2]}")

# ---- 3. N = 4: the singlet block ------------------------------------------------------------------------------------
say("3. N = 4: the singlet block and the blind vector (0, 1, 1, 0)")
b0, b3, h, g, a, e = sp.symbols('b0 b3 h g a e', nonzero=True)
k2 = sp.kronecker_product
K = [b0 * k2(X, I2) + b3 * k2(I2, X) - h * sp.eye(4), b0 * k2(Y, I2) + b3 * k2(I2, Y), b0 * k2(Z, I2) + b3 * k2(I2, Z)]
KK = sp.expand(sum((Ka * Ka for Ka in K), sp.zeros(4)))
check(sp.expand(KK - (3 * (b0 ** 2 + b3 ** 2) + h ** 2) * sp.eye(4)
                - 2 * b0 * b3 * (k2(X, X) + k2(Y, Y) + k2(Z, Z)) + 2 * h * (b0 * k2(X, I2) + b3 * k2(I2, X))) == sp.zeros(4),
      "K·K = 3(b0^2 + b3^2) + h^2 + 2b0·b3·sigma0·sigma3 - 2h(b0·X0 + b3·X3)")
cS1, cS2, cS3, cS4 = sp.symbols('c1:5')
S_gen = cS1 * k2(X, I2) + cS2 * k2(X, Y) + cS3 * k2(Y, I2) + cS4 * k2(Y, Y)
SKK = (S_gen * KK - KK * S_gen).applyfunc(sp.expand)
coef = {a_ + b_: sp.expand((k2(PAULI[a_], PAULI[b_]) * SKK).trace() / 4) for a_ in 'IXYZ' for b_ in 'IXYZ'}
want = {'IZ': -4 * sp.I * b0 * b3 * cS2, 'ZX': -4 * sp.I * b0 * b3 * cS3, 'YZ': -4 * sp.I * b3 * (b0 * cS1 - h * cS4),
        'ZY': 4 * sp.I * b0 * (b3 * cS1 + h * cS4), 'XZ': 4 * sp.I * b3 * (b0 * cS3 + h * cS2), 'ZI': 4 * sp.I * b0 * (b3 * cS2 + h * cS3)}
check(all(sp.expand(coef[k_] - want.get(k_, 0)) == 0 for k_ in coef),
      "[S, K·K] has the strings IZ, ZX, YZ, ZY with -4i·b0·b3·c2, -4i·b0·b3·c3, -4i·b3·(b0·c1 - h·c4), 4i·b0·(b3·c1 + h·c4) "
      "and XZ, ZI that follow: at every nonzero b0, b3, h it vanishes exactly when c2 = c3 = 0, h·c4 = b0·c1, (b0 + b3)·c1 = 0")
span = [k2(X, I2), k2(X, Y), k2(Y, I2), k2(Y, Y)]
cols = []
for B in span:
    C = sp.expand(B * KK - KK * B)
    cols.append([sp.expand(sp.re(x)) for x in C] + [sp.expand(sp.im(x)) for x in C])
A = sp.Matrix(cols).T
check(A.nullspace() == [], "general b0, b3: [S, K·K] = 0 has no solution in span{X0, X0Y3, Y0, Y0Y3}")
ns = A.subs(b3, -b0).nullspace()
check(len(ns) == 1 and sp.simplify(ns[0][1]) == 0 and sp.simplify(ns[0][2]) == 0
      and sp.simplify(ns[0][0] / ns[0][3] - h / b0) == 0, "b3 = -b0: exactly the line h·X0 + b0·Y0Y3")
ns0 = A.subs(h, 0).nullspace()
check(any(sp.simplify(v[0]) == 0 and sp.simplify(v[1]) == 0 and sp.simplify(v[2]) == 0 and v[3] != 0 for v in ns0),
      "control, h = 0 (the Y^(x)4-coloured row): Y0Y3 solves for every b0, b3, so the step needs h != 0")
# the pair's Laplacian block with 0 and 3 struck, difference-only coupling J01 = b0, J02 = -b0, J13 = b3, J23 = -b3
L4 = sp.zeros(4)
for (i, j), w in {(0, 1): b0, (0, 2): -b0, (1, 3): b3, (2, 3): -b3, (1, 2): e}.items():
    L4[i, i] += w; L4[j, j] += w; L4[i, j] -= w; L4[j, i] -= w
check(sp.factor(L4.extract([1, 2], [1, 2]).det()) == -(b0 + b3) ** 2
      and list(L4 * sp.Matrix([0, 1, 1, 0])) == [0, b0 + b3, -(b0 + b3), 0]
      and list(sp.expand(L4 * sp.Matrix([0, b3, b3, b0 + b3]))) == [0, 0, 0, 0],
      "the pair's block with 0 and 3 struck has determinant -(b0 + b3)^2; (0, 1, 1, 0) is in ker L exactly at "
      "b3 = -b0, and (0, b3, b3, b0 + b3), which seat 3 sees unless b3 = -b0, keeps seat 0 blind at every b0, b3")
y1, y2 = sp.symbols('y1 y2')
check(sp.expand((L4.row(0) * sp.Matrix([0, y1, y2, 0]))[0] + b0 * (y1 - y2)) == 0,
      "a vector (0, y1, y2, 0) vanishing at 0 and 3 is an eigenvector only with y1 = y2 (row 0 is -b0·(y1 - y2)), so only "
      "as (0, 1, 1, 0), and L·(0, 1, 1, 0) = (b0 + b3)·(0, 1, -1, 0) makes it one exactly at b3 = -b0")
# the singlet entry of [H, W] for a block-diagonal W, on the family (e and g symbolic): g·[Y3, S]
s_pair = (kron(sp.Matrix([1, 0]), sp.Matrix([0, 1])) - kron(sp.Matrix([0, 1]), sp.Matrix([1, 0]))) / sp.sqrt(2)
Hfam = sp.zeros(16)
for c in 'XYZ':
    Hfam += e * op(4, {1: c, 2: c}) + a * (op(4, {0: c, 1: c}) - op(4, {0: c, 2: c}) - op(4, {1: c, 3: c}) + op(4, {2: c, 3: c}))
Hfam += -h * op(4, {1: 'X'}) + h * op(4, {2: 'X'}) + g * op(4, {3: 'Y'})
Sx = sp.symbols('Sx0:16')
Sop = sp.Matrix(4, 4, lambda i, j: Sx[4 * i + j])      # an arbitrary operator on (0, 3)
# embed S on the singlet: (0, 1, 2, 3) with the pair's singlet; reorder (0,3) x pair -> (0, 1, 2, 3)
perm = []
for q0 in range(2):
    for q1_ in range(2):
        for q2_ in range(2):
            for q3 in range(2):
                perm.append(((q0 * 2 + q3), (q1_ * 2 + q2_)))
Psing = s_pair * s_pair.T
Wemb = sp.zeros(16)
for I_, (o1, q1) in enumerate(perm):
    for J_, (o2, q2) in enumerate(perm):
        Wemb[I_, J_] = Sop[o1, o2] * Psing[q1, q2]
comm = sp.expand(Hfam * Wemb - Wemb * Hfam)
proj = sp.zeros(4)                                      # <s| comm |s> as an operator on (0, 3)
for I_, (o1, q1) in enumerate(perm):
    for J_, (o2, q2) in enumerate(perm):
        proj[o1, o2] += s_pair[q1] * comm[I_, J_] * s_pair[q2]
gY3S = g * (k2(I2, Y) * Sop - Sop * k2(I2, Y))
check(sp.simplify(proj - gY3S) == sp.zeros(4), "on the family the singlet entry of [H, diag(S, 0)] is g·[Y3, S]: e drops out")
# the widened family's blocks in the pair's singlet-triplet basis, symbolic in b0, b3, e, h, g
Hwid = sp.zeros(16)
for c in 'XYZ':
    Hwid += e * op(4, {1: c, 2: c}) + b0 * (op(4, {0: c, 1: c}) - op(4, {0: c, 2: c})) \
        + b3 * (op(4, {1: c, 3: c}) - op(4, {2: c, 3: c}))
Hwid += -h * op(4, {1: 'X'}) + h * op(4, {2: 'X'}) + g * op(4, {3: 'Y'})
where = {pq: I_ for I_, pq in enumerate(perm)}
pstate = {'s': s_pair, 'x': kron(X, I2) * s_pair, 'y': kron(Y, I2) * s_pair, 'z': kron(Z, I2) * s_pair}


def blk(a_, b_):
    """<a| Hwid |b> on the pair, as an operator on (0, 3)."""
    out = sp.zeros(4)
    for o1 in range(4):
        for o2 in range(4):
            out[o1, o2] = sp.expand(sum(sp.conjugate(pstate[a_][q1]) * Hwid[where[(o1, q1)], where[(o2, q2)]] * pstate[b_][q2]
                                        for q1 in range(4) for q2 in range(4)))
    return out


Kw = dict(zip('xyz', K))
Y3 = k2(I2, Y)
check(blk('s', 's') == -3 * e * sp.eye(4) + g * Y3
      and all(blk(a_, b_) == ((e * sp.eye(4) + g * Y3) if a_ == b_ else sp.zeros(4)) for a_ in 'xyz' for b_ in 'xyz')
      and all((blk(a_, 's') - 2 * Kw[a_]).applyfunc(sp.expand) == sp.zeros(4) for a_ in 'xyz'),
      "widened family: H is -3e + g·Y3 on the singlet, e + g·Y3 on each triplet state and 2K_alpha from the singlet to |alpha>, "
      "so for W = diag(S, T) the singlet entry of [H, W] is g·[Y3, S] and the triplet entries are g·[Y3, T_alpha beta]")


def n4(a_, b_, e_, h_, g_, c0=0, c3=0, j03=0):
    bonds = {(1, 2): e_, (0, 1): a_ + c0, (0, 2): -a_ + c0, (1, 3): -b_ + c3, (2, 3): b_ + c3, (0, 3): j03}
    return bonds, {1: (-h_, 0, 0), 2: (h_, 0, 0), 3: (0, g_, 0)}


SW12 = sp.zeros(8)
for c in 'IXYZ':
    SW12 += kron(PAULI[c], PAULI[c], I2) / 2
rows3 = [
    ("the family a = 1, e = 2, h = 1, g = 3", n4(1, 1, 2, 1, 3), True, 1, 1),
    ("the family a = 2, e = -1, h = 3, g = 1", n4(2, 2, -1, 3, 1), True, 1, 1),
    ("difference only, b3 = -2·b0", n4(1, 2, 2, 1, 3), True, 0, 0),
    ("difference only, b3 = b0", n4(1, -1, 2, 1, 3), True, 0, 0),
    ("difference only, (b0, b3) = (2, -3)", n4(2, 3, 2, 1, 3), True, 0, 0),
    ("site 0 sees the pair's sum, c0 = 1", n4(1, 1, 2, 1, 3, c0=1), False, 0, 0),
]
for name, (bonds, fields), pb, pall, peven in rows3:
    H = heis(4, bonds, fields)
    b = blind_seat0(4, bonds)
    full = far_count(H, 4)
    even = far_count(H, 4, extra=lambda E: SW12 * E - E * SW12)
    say(f"   {name:40s} blind(0)={b!s:5}  far elements: all {full}, swap-even {even}")
    check(b == pb and full == pall and even == peven, f"{name}: blind {pb}, all {pall}, swap-even {peven}")
for name, (bonds, fields) in (("g = 0, b3 = -2·b0", n4(1, 2, 2, 1, 0)), ("g = 0, b3 = b0", n4(1, -1, 2, 1, 0))):
    H = heis(4, bonds, fields)
    even = far_count(H, 4, extra=lambda E: SW12 * E - E * SW12)
    check(even == 9, f"control, {name}: {even} dimensions of swap-even far elements, so the step needs g != 0")
bonds_y = {(1, 2): 2, (0, 1): 1, (0, 2): -1, (1, 3): -2, (2, 3): 2}
H = heis(4, bonds_y, {1: (0, -1, 0), 2: (0, 1, 0), 3: (0, 3, 0)})
check(far_count(H, 4, extra=lambda E: SW12 * E - E * SW12) >= 1,
      "control, the pair's opposite fields along y at b3 = -2·b0: a swap-even far element (Y^(x)4), so the step needs them along x")

# ---- 4. the join with the common carrier ------------------------------------------------------------------------------
say("4. the join: the common carrier W(s, t) at the letter configuration")
s, t, r = sp.symbols('s t r', nonzero=True)
Hb = sp.zeros(8)
for c in 'XYZ':
    Hb += s * op(3, {0: c, 1: c}) + t * op(3, {0: c, 2: c}) + r * op(3, {1: c, 2: c})
F = op(3, {1: 'X'}) + op(3, {2: 'Y'})
Z3 = op(3, {0: 'Z', 1: 'Z', 2: 'Z'})
check(sp.simplify(Z3 * Hb * Z3 - Hb) == sp.zeros(8) and sp.simplify(Z3 * F * Z3 + F) == sp.zeros(8),
      "Z^(x)3 keeps the bonds and reverses the fields: its pencil H_b + lambda·F is the field line")
for (sv, tv, rv) in ((3, 6, -2), (2, 2, 3), (1, -2, -2)):
    Wst = sum((coef * op(3, at) for coef, at in (
        (sv * tv ** 2, {0: 'X'}), (sv ** 3, {0: 'X', 1: 'X', 2: 'Y'}), (sv * tv ** 2, {0: 'X', 1: 'Y', 2: 'X'}),
        (-sv * tv ** 2, {0: 'X', 1: 'Z', 2: 'Z'}), (sv ** 2 * tv, {0: 'Y'}), (tv ** 3, {0: 'Y', 1: 'X', 2: 'Y'}),
        (sv ** 2 * tv, {0: 'Y', 1: 'Y', 2: 'X'}), (-sv ** 2 * tv, {0: 'Y', 1: 'Z', 2: 'Z'}))), sp.zeros(8))
    H = Hb.subs({s: sv, t: tv, r: rv}) + 5 * F
    sols, eu, ev, *_ = hand_solutions((S_(sv), S_(tv), S_(rv)), (1, 0), (0, 1))
    W = w_from(*sols[0], eu, ev)
    i0, j0 = next((i, j) for i in range(8) for j in range(8) if Wst[i, j] != 0)
    k0 = sp.simplify(W[i0, j0] / Wst[i0, j0])
    prop = k0 != 0 and (W - k0 * Wst).applyfunc(lambda x: sp.simplify(sp.expand(x))) == sp.zeros(8)
    check(len(sols) == 1 and prop and sp.simplify(Z3 * W * Z3 + W) == sp.zeros(8)
          and sp.simplify(H * W - W * H) == sp.zeros(8),
          f"(s, t, r) = {(sv, tv, rv)}: the hand element is a multiple of W(s, t), odd under Z^(x)3, a far element at h = 5")

say("ALL CHECKS PASS" if not FAILS else f"{len(FAILS)} CHECK(S) FAILED")
Path(__file__).parent.joinpath("results", "blind_seat_sheet_gate.txt").write_bytes(("\n".join(OUT) + "\n").encode("utf-8"))
raise SystemExit(1 if FAILS else 0)
