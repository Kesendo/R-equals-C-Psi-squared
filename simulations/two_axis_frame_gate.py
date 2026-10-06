"""Gate for the frame of a palindrome with two jumps on one site (THE_PALINDROME_AS_A_COLOURING, "Two letters on a site").

A frame of a palindrome is a unitary V commuting with every jump such that a lit string F (one anticommuting with every
jump) commutes with V H V^dagger; equivalently V^dagger F V lies in W, the anticommutant {W : [H, W] = 0, A W A = -W for
every jump A} of F158. With single-letter jumps, the page's theorem builds a frame whenever one jump commutes with all the
others, so the open configuration is: every dephased site carries two letters, and some site is undephased. There, with
D the dephased sites, O the undephased ones (m = |O|, not F158's jump count) and G_D the product of the third letters on D:
  (a) an operator anticommuting with both letters on a site is a multiple of the third, and one commuting with both is a
      multiple of 1, so the unitaries commuting with the jumps are 1_D (x) v, the lit strings are G_D (x) f, and
      W = G_D (x) W', W' = {M on U : [H, G_D (x) M] = 0};
  (b) the page's polar step (U |U|^-1, then K = e^{it} W0 + e^{-it} W0^dagger, then K |K|^-1) uses only that the
      commutant is a *-algebra, that it multiplies W into itself, that W W lies in it and that W is closed under the
      adjoint, all true for any Hermitian involutive jumps, so W' holds a Hermitian unitary s whenever the palindrome holds;
  (c) a frame image is G_D (x) v^dagger f v, with v^dagger f v a Hermitian unitary on U of trace 0 or +-2^m, so a frame exists exactly when W'
      holds a Hermitian unitary of trace 0 or +-2^m;
  (d) at m = 1 every Hermitian unitary on the qubit is +-1 or n.sigma, so every such palindrome has a frame.
G1 (exact) the one-site facts of (a) and the qubit fact of (d), symbolically.
G2 (exact) m = 1: for rational unit vectors n (the frame built up to scale; A = 0 at n = -z, not sampled), H = H0 + Q H0 Q with Q = Y (x) n.sigma (jumps X, Z on site 0) pairs
    (dim of the commutant = dim W = 1, W spanned by Q), and A = Z + n.sigma carries the lit string Y (x) Z onto Q:
    A Z A = 2(1 + n_z) n.sigma, so (1 (x) A) H (1 (x) A) commutes with Y (x) Z. Controls: H0 alone does not pair.
G3 (exact) m = 2, the page's row: commutant = span(1), W = span(Y (x) diag(1, 1, 1, -1)) (exact elements for the lower
    bounds, ranks modulo two primes for the upper); the only Hermitian unitaries
    in W' are +-diag(1, 1, 1, -1), of trace +-2, so there is no frame. Rows symmetrised on s = v0^dagger (Z (x) 1) v0 and
    on v0^dagger diag(1, 1, 1, -1) v0, v0 a rational unitary (Cayley transform): the first has the frame 1 (x) v0, exactly,
    the second has dim W = 1 with s of trace 2, so no frame.
G4 (exact) the extension: jumps X, Z on site 0 and Z on site 1, site 2 undephased; H = H0 + S H0 S with S = Y (x) s,
    s = v^dagger (X (x) 1) v, v commuting with Z_1 (rational, Cayley blocks): the page's V = P + (1 - P) F S with
    P = (1 + Z_1)/2 is unitary, commutes with all three jumps and makes every lit string F commute with V H V^dagger.
    Control: P built from Z_0 (which X_0 anticommutes with) breaks [V, X_0] = 0.
(The 'carried' and 'framed' checks of G2 and G3 and the frame of G4 hold by construction on rows built for them; what the
data decides there are the dimensions.)
Output: simulations/results/two_axis_frame_gate.txt. Run: python simulations/two_axis_frame_gate.py"""
import itertools, os, random, sys
import sympy as sp

if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass
OUT = []
_print = print
def print(*a, **k):
    _print(*a, **k); OUT.append(" ".join(str(x) for x in a))
FAILS = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok: FAILS.append(name)

I2 = sp.eye(2); X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Z = sp.Matrix([[1, 0], [0, -1]])
def kron(*ms):
    out = ms[0]
    for m in ms[1:]: out = sp.kronecker_product(out, m)
    return out
def solve_space(dim, conds):
    """basis of {M (dim x dim) : every cond(M) = 0}, conds linear"""
    syms = sp.symbols(f"m0:{dim * dim}")
    M = sp.Matrix(dim, dim, syms)
    eqs = []
    for c in conds: eqs += list(c(M))
    A, _ = sp.linear_eq_to_matrix(eqs, syms)
    return [sp.Matrix(dim, dim, list(v)) for v in A.nullspace()]

# ---------------------------------------------------------------- G1
anti = solve_space(2, [lambda M: X * M + M * X, lambda M: Z * M + M * Z])
comm1 = solve_space(2, [lambda M: X * M - M * X, lambda M: Z * M - M * Z])
check("G1 on one site, the operators anticommuting with X and Z are the multiples of Y, those commuting with both the multiples of 1",
      len(anti) == 1 and sp.simplify(anti[0] - anti[0][0, 1] / Y[0, 1] * Y) == sp.zeros(2) and len(comm1) == 1
      and sp.simplify(comm1[0] - comm1[0][0, 0] * I2) == sp.zeros(2))
a0, a1, a2, a3 = sp.symbols("a0:4", real=True)
s = a0 * I2 + a1 * X + a2 * Y + a3 * Z
sq = sp.expand(s * s - I2)
sols = sp.solve([sq[0, 0], sq[1, 1], sq[0, 1], sq[1, 0]], [a0, a1, a2, a3], dict=True)
# s^2 = 1 reads a0^2 + |a|^2 = 1 and a0 a = 0 over the reals: the solutions are a = 0, a0 = +-1, or a0 = 0, |a| = 1
def kind(sol):
    v = [sol.get(x, x) for x in (a0, a1, a2, a3)]
    return ("scalar" if v[0] in (1, -1) and all(x == 0 for x in v[1:]) else
            "vector" if v[0] == 0 and sp.simplify(v[1] ** 2 + v[2] ** 2 + v[3] ** 2 - 1) == 0 else "other")
kinds = [kind(x) for x in sols]
check("G1 s^2 - 1 = (a0^2 + |a|^2 - 1) 1 + 2 a0 a.sigma, so s^2 = 1 forces a0 a = 0: s = +-1 or n.sigma with |n| = 1; "
      "sympy's real solutions read the same two kinds",
      sp.expand(s * s - I2) == sp.expand((a0 ** 2 + a1 ** 2 + a2 ** 2 + a3 ** 2 - 1) * I2 + 2 * a0 * (a1 * X + a2 * Y + a3 * Z))
      and len(sols) > 0 and "other" not in kinds and {"scalar", "vector"} <= set(kinds), f"{kinds}")

# ---------------------------------------------------------------- G2
PRIMES = [998244353, 469762049]
def _sqrtm1(p):
    for g in range(2, 400):
        x = pow(g, (p - 1) // 4, p)
        if x * x % p == p - 1: return x
def _red(z, p, ip):
    """a Gaussian rational modulo p, i read as a square root of -1"""
    q = lambda r: sp.Rational(r).p * pow(sp.Rational(r).q, p - 2, p)
    return (q(sp.re(z)) + q(sp.im(z)) * ip) % p
def _rank(rows, ncols, p):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]; inv = pow(rows[r][c], p - 2, p); rows[r] = [x * inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]; rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
    return r
def _nullity(H, jumps, dim, sign, p):
    """dim of {M : [H, M] = 0, A M - sign M A = 0 for every jump} modulo p (an upper bound on the rational dimension)"""
    ip = _sqrtm1(p)
    ops = [(H, -1)] + [(A, -sign) for A in jumps]
    rows = []
    for O, sg in ops:
        Om = [[_red(O[i, j], p, ip) for j in range(dim)] for i in range(dim)]
        for i in range(dim):
            for j in range(dim):
                # (O M + sg M O)[i, j] = sum_k O[i,k] M[k,j] + sg M[i,k] O[k,j]
                row = [0] * (dim * dim)
                for k in range(dim):
                    row[k * dim + j] = (row[k * dim + j] + Om[i][k]) % p
                    row[i * dim + k] = (row[i * dim + k] + sg * Om[k][j]) % p
                if any(row): rows.append(row)
    return dim * dim - _rank(rows, dim * dim, p)
def spaces(H, jumps, dim):
    """(dim of the commutant, dim of W), read modulo two primes that must agree"""
    reads = [(_nullity(H, jumps, dim, 1, p), _nullity(H, jumps, dim, -1, p)) for p in PRIMES]
    if reads[0] != reads[1]: check("the two primes agree on the nullities", False, str(reads))
    return reads[0]
def _zero(M): return M.applyfunc(sp.expand) == sp.zeros(M.rows, M.cols)
def in_W(H, jumps, M): return _zero(H * M - M * H) and all(_zero(A * M + M * A) for A in jumps)
def in_N(H, jumps, M): return _zero(H * M - M * H) and all(_zero(A * M - M * A) for A in jumps)
rnd = random.Random(20261006)
def rand_herm(dim):
    H = sp.zeros(dim)
    for i in range(dim):
        H[i, i] = rnd.randint(-3, 3)
        for j in range(i + 1, dim):
            c = rnd.randint(-3, 3) + sp.I * rnd.randint(-3, 3); H[i, j] = c; H[j, i] = sp.conjugate(c)
    return H
jumps2 = [kron(X, I2), kron(Z, I2)]
g2_ok = True; g2_ctrl = True; rows = 0
for n in [(sp.Rational(3, 5), sp.Rational(4, 5), 0), (0, sp.Rational(3, 5), sp.Rational(4, 5)), (sp.Rational(2, 3), sp.Rational(1, 3), sp.Rational(2, 3)),
          (sp.Rational(6, 7), sp.Rational(-2, 7), sp.Rational(3, 7)), (1, 0, 0)]:
    ns = n[0] * X + n[1] * Y + n[2] * Z
    Q = kron(Y, ns)
    H0 = rand_herm(4); H = (H0 + Q * H0 * Q).applyfunc(sp.expand)
    dN, dW = spaces(H, jumps2, 4)
    A = Z + ns
    carried = sp.simplify(A * Z * A - 2 * (1 + n[2]) * ns) == sp.zeros(2)
    HV = kron(I2, A) * H * kron(I2, A)
    F = kron(Y, Z)
    framed = sp.simplify(HV * F - F * HV) == sp.zeros(4)
    # exact: 1 in the commutant and Q in W give the lower bounds, the modular nullities the upper ones
    g2_ok &= dN == dW == 1 and in_N(H, jumps2, sp.eye(4)) and in_W(H, jumps2, Q) and carried and framed
    N0, W0 = spaces(H0, jumps2, 4)
    g2_ctrl &= W0 == 0  # the identity in the commutant gives dim >= 1 exactly, W read 0 modulo p is exact
    rows += 1
check(f"G2 m = 1, {rows} rational axes n: the row pairs, dim commutant = dim W = 1 with W = span(Y (x) n.sigma), and A = Z + n.sigma carries Y (x) Z onto it, "
      "A Z A = 2(1 + n_z) n.sigma and [(1 (x) A) H (1 (x) A), Y (x) Z] = 0 (exact)", g2_ok)
check("G2 control: H0 before symmetrisation does not pair (W = 0)", g2_ctrl)

# ---------------------------------------------------------------- G3
e = sp.eye(4)
B1 = e[:, 0] * e[:, 3].T + e[:, 3] * e[:, 0].T
B2 = (e[:, 1] + e[:, 2]) * e[:, 3].T + e[:, 3] * (e[:, 1] + e[:, 2]).T
Hc = kron(I2, sp.diag(1, 2, 3, 4)) + kron(X, B1) + kron(Z, B2)
jumps3 = [kron(X, I2, I2), kron(Z, I2, I2)]
dNc, dWc = spaces(Hc, jumps3, 8)
w = sp.diag(1, 1, 1, -1)
S = kron(Y, w)
check("G3 the page's row: the commutant is span(1), W = span(Y (x) diag(1, 1, 1, -1)); its Hermitian unitaries are +-S, "
      "trace of s = +-2, neither 0 nor +-4: no frame",
      dNc == dWc == 1 and in_N(Hc, jumps3, sp.eye(8)) and in_W(Hc, jumps3, S) and w.trace() == 2 and (S * S) == sp.eye(8))
# a rational unitary on U by the Cayley transform v0 = (1 - iK)(1 + iK)^-1 of a rational Hermitian K
K = sp.Matrix([[1, 2, 0, -1], [2, 0, 1, 1], [0, 1, -2, 1], [-1, 1, 1, 3]])
v0 = ((sp.eye(4) - sp.I * K) * (sp.eye(4) + sp.I * K).inv()).applyfunc(sp.simplify)
unitary = _zero(v0.H * v0 - sp.eye(4))
s_bal = (v0.H * kron(Z, I2) * v0).applyfunc(sp.simplify)
Q3 = kron(Y, s_bal)
H0 = rand_herm(8); H3 = (H0 + Q3 * H0 * Q3).applyfunc(sp.expand)
dN3, dW3 = spaces(H3, jumps3, 8)
V = kron(I2, v0)
HV = (V * H3 * V.H).applyfunc(sp.expand)
F3 = kron(Y, Z, I2)
paulis2 = [kron(P, R) for P in (I2, X, Y, Z) for R in (I2, X, Y, Z)]
s_not_pauli = all(not _zero(s_bal - P) and not _zero(s_bal + P) for P in paulis2)
check("G3 s = v0^dagger (Z (x) 1) v0 is no Pauli string up to sign, so with W' = span(s) no Clifford frame exists", s_not_pauli)
check("G3 a row on s = v0^dagger (Z (x) 1) v0, v0 a rational unitary (trace 0), pairs with dim commutant "
      "= dim W = 1, the frame V = 1 (x) v0 makes the lit string Y (x) Z (x) 1 commute with V H V^dagger (exact), and G = Y (x) 1 (x) 1 is not in W",
      unitary and dN3 == dW3 == 1 and in_W(H3, jumps3, Q3) and _zero(HV * F3 - F3 * HV) and not in_W(H3, jumps3, F3)
      and not in_W(H3, jumps3, kron(Y, I2, I2)),
      f"dim commutant = {dN3}, dim W = {dW3} (modulo two primes)")
s_unb = (v0.H * w * v0).applyfunc(sp.simplify)
Q4 = kron(Y, s_unb)
H4 = (H0 + Q4 * H0 * Q4).applyfunc(sp.expand)
dN4, dW4 = spaces(H4, jumps3, 8)
check("G3 the same construction on s = v0^dagger diag(1, 1, 1, -1) v0 (trace 2) pairs with dim commutant = dim W = 1, "
      "so its only Hermitian unitaries are +-s and the criterion says no frame",
      dN4 == dW4 == 1 and in_W(H4, jumps3, Q4) and s_unb.trace() == 2, f"dim commutant = {dN4}, dim W = {dW4}")

# ---------------------------------------------------------------- G4
def cayley(Kh):
    n = Kh.rows
    return ((sp.eye(n) - sp.I * Kh) * (sp.eye(n) + sp.I * Kh).inv()).applyfunc(sp.simplify)
a_ = cayley(sp.Matrix([[1, 2 - sp.I], [2 + sp.I, -1]])); b_ = cayley(sp.Matrix([[0, 1 + sp.I], [1 - sp.I, 2]]))
P0 = sp.Matrix([[1, 0], [0, 0]]); P1 = sp.Matrix([[0, 0], [0, 1]])
v4 = (kron(P0, a_) + kron(P1, b_)).applyfunc(sp.simplify)
s4 = (v4.H * kron(X, I2) * v4).applyfunc(sp.simplify)
S4 = kron(Y, s4)
jumps4 = [kron(X, I2, I2), kron(Z, I2, I2), kron(I2, Z, I2)]
H0 = rand_herm(8); H4 = (H0 + S4 * H0 * S4).applyfunc(sp.expand)
dN, dW = spaces(H4, jumps4, 8)
lit = [kron(Y, P, R) for P in (X, Y) for R in (I2, X, Y, Z)]
Pp = (sp.eye(8) + jumps4[2]) / 2
ok4 = _zero(v4.H * v4 - sp.eye(4)) and _zero(S4 * S4 - sp.eye(8)) and in_W(H4, jumps4, S4) and in_N(H4, jumps4, sp.eye(8))
no_colouring = all(not in_W(H4, jumps4, F) for F in lit)
for F in lit:
    V = (Pp + (sp.eye(8) - Pp) * F * S4).applyfunc(sp.expand)
    HV = (V * H4 * V.H).applyfunc(sp.expand)
    ok4 &= _zero(V.H * V - sp.eye(8)) and all(_zero(V * A - A * V) for A in jumps4) and _zero(HV * F - F * HV)
check("G4 the extension: jumps X, Z on site 0 and Z on site 1 pair, no lit string commutes with H, and V = P + (1 - P) F S with "
      "P = (1 + Z_1)/2 is unitary, commutes with every jump and frames each of the eight lit strings (exact)",
      ok4 and no_colouring and dN == dW == 1, f"dim commutant = {dN}, dim W = {dW}")
Pbad = (sp.eye(8) + jumps4[1]) / 2
Vbad = (Pbad + (sp.eye(8) - Pbad) * lit[0] * S4).applyfunc(sp.expand)
check("G4 control: with P built from Z_0, which X_0 anticommutes with, V no longer commutes with X_0",
      not _zero(Vbad * jumps4[0] - jumps4[0] * Vbad))

print("\nALL PASS" if not FAILS else f"\n{len(FAILS)} FAIL")
os.makedirs("simulations/results", exist_ok=True)
open("simulations/results/two_axis_frame_gate.txt", "w", encoding="utf-8").write("\n".join(OUT) + "\n")
sys.exit(1 if FAILS else 0)
