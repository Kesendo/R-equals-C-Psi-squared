"""The centre rung at N = 4 (docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, "The centre rung at N = 4"): where the block
(2, 2) holds the real eigenvalue -4 gamma.

Pauli book, J = 1 the unit, H = sum_bonds J_b (XX + YY + Delta ZZ) real symmetric on the six popcount-2
states of N = 4. In the block (2, 2) write B(gamma) = -i A + gamma D, A = H (x) I - I (x) H^T, D = -2 Hamming, and
T4 = B + 4 gamma (the rung k = 2). The global flip X (x -> x~ on ket and bra) commutes with B. With a_i, b_i =
(x_i +- x~_i)/sqrt 2 over the three complementary pairs, E = span{a_i} and O = span{b_i} are the X-even and X-odd states,
h_E, h_O the restrictions of H, and the X-odd half of the block is |E><O| (+) |O><E|. On it H acts by the Sylvester map
S(U) = h_E U - U h_O, and the centred dephasing touches only the six pair cells (a_i, b_i), (b_i, a_i), where it swaps
them with weight 4 gamma. When h_E and h_O share no eigenvalue,

    det T4_odd = (det S)^2 * det(I - 16 gamma^2 G^2),    G_ji = [S^-1(|a_i><b_i|)]_(a_j, b_j),

G real symmetric 3 x 3, so the X-odd half reaches -4 gamma exactly at gamma = 1/(4 |omega|), one point per modulus
|omega| > 0 among the eigenvalues omega of G, of order the number of eigenvalues of that modulus.

Rows (all exact, sympy over Q(i) and Q(i)[gamma] / Q(Delta, gamma); a row compares to 0 or to a closed form):
  C1  X commutes with B; in the basis |e><o|, |o><e| the centred dephasing is the pair swap 4 gamma (U_ii <-> V_ii) and
      zero on every other cell of the X-odd half (all other cells there have Hamming 2).
  C2  the determinant identity, polynomial in gamma, constant (det S)^2 included, on the chain at Delta = 1, 1/2, 2, 3,
      the chain (1, 3/5, 7/5), five random bonds, and the star; G symmetric in each. The chain (1, 3/5, 7/5) and the
      random bonds at Delta = 1 have no doubled modulus (the characteristic polynomial of G^2 has nonzero
      discriminant); the star at Delta = 1/2 has charpoly(G) = (4u - 1)^2 (8u + 1)/128. Control: a dephasing shift on an
      X-symmetric quadruple of off-pair cells breaks the identity. The N = 4 ring at Delta = 1, 1/2 is outside the
      hypothesis: h_E and h_O share the eigenvalue 0 and S is singular.
  C3  the chain at symbolic Delta: b = |0110> - |1001> is an eigenvector of H at -Delta; G splits under the chain
      reflection into G+ = (1 + Delta^2)/(8 Delta) and G- = (Delta/8) [[-1, 1], [1, 4 Delta^2 - 5]], so the points are
      gamma* = 2 Delta/(1 + Delta^2) (at most 1, equal only at Delta = 1) and the roots of
      Delta^4 (Delta^2 - 1)^2 u^2 - Delta^2 (4 Delta^4 - 10 Delta^2 + 7) u + 1, u = gamma^2; the X-odd determinant
      equals the closed form up to a constant, identically in Delta and gamma; det S = 32768 Delta, so the hypothesis
      holds on the chain at every Delta != 0.
  C4  the chain's X-even half: det T4_even = c * Delta^4 (Delta^2 - 1)^2 (Delta^2 + 1)^2 gamma^6 with c a nonzero
      constant, identically in Delta and gamma: nonsingular at every gamma > 0 for Delta not in {0, +-1}, singular at
      every gamma for Delta = +-1.
  C5  exact nullities at the points: chain Delta = 1, gamma = 1: nullity 2 in the X-odd half; Delta = 1/2,
      gamma = 4/5: nullity 1; Delta = 1/2, gamma = 1: nullity 0 (control).
  C6  the doubling at Delta = 1 is not forced by SU(2) and the reflection: on the chain (1, r, 1) at Delta = 1 the
      eigenvalues of G are 0, k1 = r/(4(2 - r)) and k2 = -(r^2 - 2r + 2)/(4r); their moduli agree at two isolated
      points of r > 0, r = 1 (k1 + k2 = 0, numerator (r - 1)(r^2 - 2r + 4)) and the real root r0 = 3.658967 of
      r^3 - 5r^2 + 6r - 4 (k1 = k2), and nowhere else: the doubling is a crossing of two branches.
  C7  the chain at Delta = 1: the orbit sum T2 of the Hamming-2 cells is in the kernel of T4 at symbolic gamma
      (PROOF_UNIFORM_LAW Lemma A1); the X-even half has rank 16 over Q(i)(gamma), so nullity 2 at generic gamma; exact
      at gamma = 2/3, 7/5 the block holds -4 gamma with nullity 2, one mode of pure weight 2 (T2) and one that mixes
      weights, and the X-odd half holds none, so both sit in the X-even half that C4 finds singular at every gamma.
  C8  the X-even half is the dephasing compressed onto the commutant: with C_E, C_O the commutants of h_E, h_O
      (dimensions c_E, c_O), their pair diagonals (read in the bases a_i, b_i) spanning L_E, L_O in R^3, and r the rank
      of the pairing <delta(C), delta(C')> between them, the X-even half holds -4 gamma with nullity
      (c_E - r) + (c_O - r), of which (c_E - dim L_E) + (c_O - dim L_O) are of pure weight 2; exact at gamma = 2/3, 7/5
      on the chain at Delta = 1, 1/2, 0 and on the chain (1, 3/5, 7/5), six random bonds, the star, the ring and the
      complete graph at Delta = 1, 1/2; ring (c_E, c_O, r) = (3, 5, 2); on the chain the pair diagonals of 1, h, h^2 have
      determinant 16 Delta (Delta^2 - 1) on the E side and 16 Delta (Delta^2 + 1) on the O side. The ring's whole
      block (2, 2) at Delta = 1 has nullity 4, three of pure weight 2, none X-odd. On the chain at Delta = 1 the partner
      of T2 is C = P_O - 3|b3><b3|, whose pair diagonal (1, 1, -2) is orthogonal to L_E, and C + gamma C' with C' on
      Hamming-2 cells is in the kernel at symbolic gamma.
  C9  the N = 8 chain at Delta = 1, block (4, 4): the restrictions of H to the 35 X-even and the 35 X-odd states are
      cyclic (their powers span 35 dimensions mod two primes, so the commutants are the powers), share no eigenvalue and
      have simple spectra (exact characteristic polynomials), and the E-O matrix of D + 8 between the two power bases has
      rank 34 mod both primes, so at most 2 over Q; the orbit sum of the block's Hamming-4 cells commutes with H
      (exact integers) and lies in the kernel, so exactly 2.

Run: python simulations/f50_centre_rung_green.py   (under a minute; prints "ALL GATES PASS" or a FAIL line)
"""
import itertools
import sympy as sp
from sympy import Rational as Q
from sympy.polys.matrices import DomainMatrix

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{(': ' + detail) if detail else ''}", flush=True)
    if not ok:
        FAILS.append(name)


N = 4
STATES = [x for x in range(2 ** N) if bin(x).count("1") == 2]
FLIP = 2 ** N - 1
PAIRS = [(x, x ^ FLIP) for x in STATES if x < x ^ FLIP]
CELLS = list(itertools.product(range(6), range(6)))
HAM = [bin(STATES[k] ^ STATES[l]).count("1") for k, l in CELLS]
g, D_, r, u = sp.symbols('gamma Delta r u')
CHAIN = {(0, 1): 1, (1, 2): 1, (2, 3): 1}


def bit(x, l):
    return (x >> l) & 1


def hamiltonian(bonds, delta):
    H = sp.zeros(6, 6)
    for c, x in enumerate(STATES):
        for (i, j), J in bonds.items():
            H[c, c] += J * delta * (1 if bit(x, i) == bit(x, j) else -1)
            if bit(x, i) != bit(x, j):
                H[STATES.index(x ^ (1 << i) ^ (1 << j)), c] += 2 * J
    return H


def block_T(H, gam, extra=None):
    """T = -i A + gamma D + 4 gamma on the cells |x><y| (column-stacked index (k, l))."""
    T = sp.zeros(36, 36)
    for c, (k, l) in enumerate(CELLS):
        for rr, (k2, l2) in enumerate(CELLS):
            if l2 == l:
                T[rr, c] += -sp.I * H[k2, k]
            if k2 == k:
                T[rr, c] -= -sp.I * H[l, l2]
        T[c, c] += gam * (4 - 2 * HAM[c])
    if extra:
        for c, v in extra.items():
            T[c, c] += v
    return T


def vec(coeffs):
    v = sp.zeros(6, 1)
    for x, cf in coeffs.items():
        v[STATES.index(x)] += cf
    return v


EV = [vec({x: 1, y: 1}) / sp.sqrt(2) for x, y in PAIRS]
OD = [vec({x: 1, y: -1}) / sp.sqrt(2) for x, y in PAIRS]
Eb, Ob = sp.Matrix.hstack(*EV), sp.Matrix.hstack(*OD)


def op_cell(v, w):
    """|v><w| as a 36-vector on the cells (k, l)."""
    return sp.Matrix([v[k] * w[l] for k, l in CELLS])


# X-odd basis: Y_ij = |a_i><b_j| (9), Z_ij = |b_i><a_j| (9); X-even basis: |a_i><a_j|, |b_i><b_j| (18)
ODD = [op_cell(EV[i], OD[j]) for i in range(3) for j in range(3)] + [op_cell(OD[i], EV[j]) for i in range(3) for j in range(3)]
EVEN = [op_cell(EV[i], EV[j]) for i in range(3) for j in range(3)] + [op_cell(OD[i], OD[j]) for i in range(3) for j in range(3)]
BASIS = sp.Matrix.hstack(*(ODD + EVEN))      # orthonormal (real), so BASIS^T = BASIS^-1
assert sp.simplify(BASIS.T * BASIS - sp.eye(36)).is_zero_matrix


def split(T):
    Tp = sp.simplify(BASIS.T * T * BASIS)
    return Tp[:18, :18], Tp[18:, 18:], Tp[:18, 18:], Tp[18:, :18]


def sylvester(H):
    hE, hO = sp.simplify(Eb.T * H * Eb), sp.simplify(Ob.T * H * Ob)
    assert sp.simplify(Eb.T * H * Ob).is_zero_matrix
    return sp.kronecker_product(hE, sp.eye(3)) - sp.kronecker_product(sp.eye(3), hO.T)   # Y[e, o] at e*3 + o


def green(H):
    S = sylvester(H)
    if S.det() == 0:
        return S, None
    Sinv = S.inv()
    K = sp.Matrix(3, 3, lambda j, i: sp.simplify(Sinv[j * 3 + j, i * 3 + i]))
    return S, K


def det_poly(M):
    Dm = DomainMatrix.from_Matrix(M.applyfunc(sp.expand))
    return sp.expand(Dm.domain.to_sympy(Dm.det()) if not hasattr(Dm.det(), "as_expr") else Dm.det().as_expr())


def nullity(M):
    Dm = DomainMatrix.from_Matrix(M.applyfunc(sp.nsimplify)).convert_to(sp.QQ_I)
    return M.shape[0] - Dm.rank()


print("C1  the X flip and the shape of the X-odd half")
Xop = sp.zeros(36, 36)
for c, (k, l) in enumerate(CELLS):
    Xop[CELLS.index((STATES.index(STATES[k] ^ FLIP), STATES.index(STATES[l] ^ FLIP))), c] = 1
H1 = hamiltonian(CHAIN, 1)
T1 = block_T(H1, g)
check("X commutes with T (chain, Delta = 1, symbolic gamma)", (Xop * T1 - T1 * Xop).is_zero_matrix)
Todd, Teven, Toe, Teo = split(T1)
check("T is block diagonal in X-odd (+) X-even", Toe.is_zero_matrix and Teo.is_zero_matrix)
Dodd = sp.simplify(Todd - Todd.subs(g, 0))
want = sp.zeros(18, 18)
for i in range(3):
    want[i * 3 + i, 9 + i * 3 + i] = want[9 + i * 3 + i, i * 3 + i] = 4 * g
check("centred dephasing on the X-odd half = the pair swap 4 gamma, zero elsewhere", (Dodd - want).is_zero_matrix)

print("C2  det T_odd = (det S)^2 det(I - 16 gamma^2 G^2)")
graphs = [("chain Delta=1", CHAIN, 1), ("chain Delta=1/2", CHAIN, Q(1, 2)), ("chain Delta=2", CHAIN, 2),
          ("chain Delta=3", CHAIN, 3), ("chain (1,3/5,7/5)", {(0, 1): 1, (1, 2): Q(3, 5), (2, 3): Q(7, 5)}, 1),
          ("five random bonds", {(0, 1): Q(2, 3), (1, 2): Q(5, 4), (0, 2): Q(1, 3), (2, 3): Q(3, 2), (1, 3): Q(4, 7)}, 1),
          ("star", {(0, 1): 1, (0, 2): 1, (0, 3): 1}, Q(1, 2))]
for name, bonds, delta in graphs:
    H = hamiltonian(bonds, delta)
    S, K = green(H)
    if K is None:
        check(f"{name}: h_E and h_O disjoint", False, "Sylvester map singular")
        continue
    lhs = det_poly(split(block_T(H, g))[0])
    rhs = sp.expand(S.det() ** 2 * (sp.eye(3) - 16 * g ** 2 * K * K).det())
    check(f"{name}: identity in gamma, G symmetric", sp.expand(lhs - rhs) == 0 and K == K.T,
          f"charpoly G = {sp.factor(K.charpoly(u).as_expr())}")
    if name in ("chain (1,3/5,7/5)", "five random bonds"):
        check(f"{name}: no doubled modulus of G", sp.discriminant((K * K).charpoly(u).as_expr(), u) != 0)
    if name == "star":
        check("star: charpoly G = (4u - 1)^2 (8u + 1)/128",
              sp.expand(K.charpoly(u).as_expr() - (4 * u - 1) ** 2 * (8 * u + 1) / 128) == 0)
H = hamiltonian(CHAIN, Q(1, 2))
S, K = green(H)
RING = {(0, 1): 1, (1, 2): 1, (2, 3): 1, (3, 0): 1}
for delta in (1, Q(1, 2)):
    Hr = hamiltonian(RING, delta)
    hE, hO = Eb.T * Hr * Eb, Ob.T * Hr * Ob
    shared = [e for e in hE.eigenvals() if hO.charpoly(u).as_expr().subs(u, e) == 0]
    check(f"ring Delta={delta}: outside the hypothesis, h_E and h_O share the eigenvalue 0 only",
          sylvester(Hr).det() == 0 and shared == [0], f"shared {shared}")
quad = [(0b0011, 0b0101), (0b0101, 0b0011), (0b1100, 0b1010), (0b1010, 0b1100)]   # Hamming 2, closed under X
bad_T = block_T(H, g, extra={CELLS.index((STATES.index(x), STATES.index(y))): g for x, y in quad})
assert (Xop * bad_T - bad_T * Xop).is_zero_matrix
bad = det_poly(split(bad_T)[0])
check("control: a dephasing shift on an X-symmetric quadruple of off-pair cells breaks the identity",
      sp.expand(bad - S.det() ** 2 * (sp.eye(3) - 16 * g ** 2 * K * K).det()) != 0)

print("C3  the chain at symbolic Delta")
HD = hamiltonian(CHAIN, D_)
bvec = vec({0b0110: 1, 0b1001: -1})
check("b = |0110> - |1001> is an eigenvector of H at -Delta", sp.simplify(HD * bvec + D_ * bvec).is_zero_matrix)
S, K = green(HD)
order = [PAIRS.index((0b0110, 0b1001))] + [i for i in range(3) if PAIRS[i] != (0b0110, 0b1001)]
Kp = sp.simplify(K.extract(order, order))
check("G splits 1 + 2 under the reflection", all(Kp[0, j] == 0 and Kp[j, 0] == 0 for j in (1, 2)))
check("G+ = (1 + Delta^2)/(8 Delta)", sp.simplify(Kp[0, 0] - (1 + D_ ** 2) / (8 * D_)) == 0)
Km = Kp.extract([1, 2], [1, 2])
target = (D_ / 8) * sp.Matrix([[-1, 1], [1, 4 * D_ ** 2 - 5]])
check("G- = (Delta/8)[[-1, 1], [1, 4 Delta^2 - 5]] (up to the order of the two pairs)",
      sp.simplify(Km - target).is_zero_matrix or sp.simplify(Km.extract([1, 0], [1, 0]) - target).is_zero_matrix)
closed = sp.expand(((1 + D_ ** 2) ** 2 * g ** 2 - 4 * D_ ** 2) *
                   (D_ ** 4 * (D_ ** 2 - 1) ** 2 * g ** 4 - D_ ** 2 * (4 * D_ ** 4 - 10 * D_ ** 2 + 7) * g ** 2 + 1))
ratio = sp.factor(sp.cancel(sp.expand(S.det() ** 2 * (sp.eye(3) - 16 * g ** 2 * K * K).det()) / closed))
check("det(I - 16 gamma^2 G^2) = closed form x a gamma-free factor", not ratio.has(g), f"factor {ratio}")
lhs = sp.factor(det_poly(split(block_T(HD, g))[0]))
check("det T_odd equals (det S)^2 det(I - 16 gamma^2 G^2) identically in Delta, gamma",
      sp.cancel(lhs / (ratio * closed)) == 1, f"det T_odd = {lhs}")
gs = 2 * D_ / (1 + D_ ** 2)
check("det S = 32768 Delta (the hypothesis holds at every Delta != 0)", sp.expand(S.det() - 32768 * D_) == 0)
check("gamma* = 2 Delta/(1 + Delta^2) <= 1 with equality only at Delta = 1",
      sp.factor(1 - gs) == sp.factor((D_ - 1) ** 2 / (D_ ** 2 + 1)))

print("C4  the chain's X-even half")
ev = sp.factor(det_poly(split(block_T(HD, g))[1]))
cst = sp.cancel(ev / (D_ ** 4 * (D_ ** 2 - 1) ** 2 * (D_ ** 2 + 1) ** 2 * g ** 6))
check("det T_even = c Delta^4 (Delta^2-1)^2 (Delta^2+1)^2 gamma^6, c a nonzero constant",
      cst != 0 and not cst.has(g) and not cst.has(D_), f"c = {cst}")

print("C5  exact nullities at the points")
for delta, gam, want_n in ((1, 1, 2), (Q(1, 2), Q(4, 5), 1), (Q(1, 2), 1, 0)):
    nl = nullity(split(block_T(hamiltonian(CHAIN, delta), gam))[0])
    check(f"Delta = {delta}, gamma = {gam}: nullity of T_odd = {want_n}", nl == want_n, f"read {nl}")

print("C6  the chain (1, r, 1) at Delta = 1")
S, K = green(hamiltonian({(0, 1): 1, (1, 2): r, (2, 3): 1}, 1))
cp = sp.factor(K.charpoly(u).as_expr())
k1, k2 = r / (4 * (2 - r)), -(r ** 2 - 2 * r + 2) / (4 * r)
check("eigenvalues of G are 0, r/(4(2-r)), -(r^2-2r+2)/(4r)",
      sp.simplify(cp - u * (u - k1) * (u - k2)) == 0, f"charpoly {cp}")
n_plus = sp.factor(sp.numer(sp.together(k1 + k2)))
n_minus = sp.factor(sp.numer(sp.together(k1 - k2)))
roots_plus = [x for x in sp.Poly(n_plus, r).real_roots() if x > 0]
roots_minus = [x for x in sp.Poly(n_minus, r).real_roots() if x > 0]
check("k1 + k2 = 0 for r > 0 only at r = 1",
      sp.expand(n_plus ** 2 - ((r - 1) * (r ** 2 - 2 * r + 4)) ** 2) == 0 and roots_plus == [1])
check("k1 = k2 for r > 0 only at the real root of r^3 - 5r^2 + 6r - 4",
      sp.expand(n_minus ** 2 - (r ** 3 - 5 * r ** 2 + 6 * r - 4) ** 2) == 0 and len(roots_minus) == 1,
      f"r0 = {sp.N(roots_minus[0], 7) if roots_minus else None}, gamma there = "
      f"{sp.N(1 / (4 * abs(k1.subs(r, roots_minus[0]))), 7) if roots_minus else None}")

print("C7  the permanent modes of the chain's block at Delta = 1")
T2 = sp.Matrix([1 if HAM[c] == 2 else 0 for c in range(36)])
check("the orbit sum T2 of the Hamming-2 cells is in the kernel of T4 at symbolic gamma (Lemma A1)",
      (block_T(hamiltonian(CHAIN, 1), g) * T2).is_zero_matrix)
Tev = split(block_T(hamiltonian(CHAIN, 1), g))[1]
Dm = DomainMatrix.from_Matrix(Tev.applyfunc(sp.expand))
rk = Dm.to_field().rank()
check("rank of the X-even half over Q(i)(gamma) is 16 (generic nullity 2)", rk == 16, f"rank {rk}")
for gam in (Q(2, 3), Q(7, 5)):
    T = block_T(hamiltonian(CHAIN, 1), gam)
    null = nullity(T)
    PF = [c for c in range(36) if HAM[c] != 2]
    pin = sp.zeros(len(PF), 36)
    for i, c in enumerate(PF):
        pin[i, c] = 1
    pure = nullity(sp.Matrix.vstack(T, pin)) - len(PF)
    n_odd = nullity(split(T)[0])
    check(f"gamma = {gam}: nullity 2 = one pure weight-2 + one weight-mixing, none in the X-odd half",
          null == 2 and pure == 1 and n_odd == 0, f"nullity {null}, pure weight 2: {pure}, X-odd: {n_odd}")

print("C8  the X-even half: the dephasing compressed onto the commutant")


def commutant_diagonals(h):
    """A basis of the commutant of the 3 x 3 matrix h: the 3 x c matrix of the elements' diagonals, and c."""
    X = sp.Matrix(3, 3, sp.symbols('x0:9'))
    A_, _ = sp.linear_eq_to_matrix(list(h * X - X * h), list(X))
    basis = A_.nullspace()
    return sp.Matrix.hstack(*[sp.Matrix([v[0], v[4], v[8]]) for v in basis]), len(basis)


PAIR_DIAG = [i * 3 + i for i in range(3)] + [9 + i * 3 + i for i in range(3)]   # |a_i><a_i|, |b_i><b_i| in the X-even half
PIN = sp.zeros(6, 18)
for i, c in enumerate(PAIR_DIAG):
    PIN[i, c] = 1
RING = {(0, 1): 1, (1, 2): 1, (2, 3): 1, (0, 3): 1}
C8_RUNS = [("chain", CHAIN, d) for d in (1, Q(1, 2), 0)]
C8_RUNS += [(name, bonds, d) for name, bonds in (
    ("chain (1, 3/5, 7/5)", {(0, 1): 1, (1, 2): Q(3, 5), (2, 3): Q(7, 5)}),
    ("six random bonds", {(0, 1): Q(5, 7), (1, 2): Q(9, 4), (2, 3): Q(2, 3), (0, 2): Q(1, 3), (1, 3): Q(6, 5), (0, 3): Q(3, 8)}),
    ("star", {(0, 1): 1, (0, 2): 1, (0, 3): 1}), ("ring", RING),
    ("complete", {(i, j): 1 for i in range(4) for j in range(i + 1, 4)})) for d in (1, Q(1, 2))]
for name, bonds, delta in C8_RUNS:
    H = hamiltonian(bonds, delta)
    dE, cE = commutant_diagonals(Eb.T * H * Eb)
    dO, cO = commutant_diagonals(Ob.T * H * Ob)
    rK = (dE.T * dO).rank()
    pred, pred_pure = cE + cO - 2 * rK, cE + cO - dE.rank() - dO.rank()
    got = []
    for gam in (Q(2, 3), Q(7, 5)):
        Tev = split(block_T(H, gam))[1]
        got.append((nullity(Tev), nullity(sp.Matrix.vstack(Tev, PIN)) - 6))
    check(f"{name}, Delta = {delta}: X-even nullity (c_E - r) + (c_O - r) = {pred}, of pure weight 2: {pred_pure}",
          all(n == pred and q == pred_pure for n, q in got), f"c_E {cE}, c_O {cO}, r {rK}; exact (nullity, pure) {got}")
    if name == "ring":
        check(f"ring, Delta = {delta}: (c_E, c_O, r) = (3, 5, 2)", (cE, cO, rK) == (3, 5, 2), f"{(cE, cO, rK)}")
H = hamiltonian(CHAIN, D_)
for h, want in ((Eb.T * H * Eb, 16 * D_ * (D_ ** 2 - 1)), (Ob.T * H * Ob, 16 * D_ * (D_ ** 2 + 1))):
    h = h.applyfunc(sp.simplify)
    M = sp.Matrix([[(h ** k)[i, i] for k in range(3)] for i in range(3)])
    check(f"chain: det of the pair diagonals of 1, h, h^2 = {want}", sp.expand(M.det() - want) == 0, f"{sp.factor(M.det())}")
for gam in (Q(2, 3), Q(7, 5)):
    T = block_T(hamiltonian(RING, 1), gam)
    PF = [c for c in range(36) if HAM[c] != 2]
    pin = sp.zeros(len(PF), 36)
    for i, c in enumerate(PF):
        pin[i, c] = 1
    null, pure, n_odd = nullity(T), nullity(sp.Matrix.vstack(T, pin)) - len(PF), nullity(split(T)[0])
    check(f"ring, Delta = 1, gamma = {gam}: the whole block holds -4 gamma with nullity 4, three of pure weight 2, none X-odd",
          (null, pure, n_odd) == (4, 3, 0), f"nullity {null}, pure {pure}, X-odd {n_odd}")
H = hamiltonian(CHAIN, 1)
Cm = Ob * sp.diag(1, 1, -2) * Ob.T
check("chain, Delta = 1: C = P_O - 3|b3><b3|, pair diagonal (1, 1, -2), commutes with H", (H * Cm - Cm * H).is_zero_matrix)
dE, _ = commutant_diagonals(Eb.T * H * Eb)
check("(1, 1, -2) is orthogonal to every pair diagonal the commutant of h_E reaches", (sp.Matrix([[1, 1, -2]]) * dE).is_zero_matrix)
Cv = sp.Matrix([Cm[k, l] for k, l in CELLS])
H2 = [c for c in range(36) if HAM[c] == 2]
vs = sp.symbols('v0:%d' % len(H2))
R = sp.zeros(36, 1)
for s_, c in zip(vs, H2):
    R[c] = s_
sol = sp.solve(list(block_T(H, 0) * R + (block_T(H, 1) - block_T(H, 0)) * Cv), vs, dict=True)   # -i[H, R] + (D + 4) C = 0
ok = bool(sol) and (block_T(H, g) * (Cv + g * R.subs(sol[0]).subs({s_: 0 for s_ in vs}))).applyfunc(sp.expand).is_zero_matrix
check("the partner C + gamma C', C' on Hamming-2 cells, is in the kernel of T4 at symbolic gamma", ok)

print("C9  the centre block (4, 4) of the N = 8 chain at Delta = 1: the same compression, an upper bound")


def centre_compression_n8(prime):
    """Rank mod prime of the E-O matrix of D + 8 on the power bases of h_E, h_O (block (4, 4), N = 8 chain, Delta = 1)."""
    n8 = 8
    st = [x for x in range(2 ** n8) if bin(x).count("1") == 4]
    ix = {x: i for i, x in enumerate(st)}
    full = 2 ** n8 - 1
    prs = [(x, x ^ full) for x in st if x < x ^ full]
    d = len(prs)
    Hm = [[0] * len(st) for _ in st]
    for c, x in enumerate(st):
        for i in range(n8 - 1):
            bi, bj = (x >> i) & 1, (x >> (i + 1)) & 1
            Hm[c][c] += 1 if bi == bj else -1
            if bi != bj:
                Hm[ix[x ^ (1 << i) ^ (1 << (i + 1))]][c] += 2

    def restrict(sign):   # H on span{x + sign x~}, the basis vectors orthogonal and of one norm
        return [[Hm[ix[u]][ix[a]] + sign * Hm[ix[u]][ix[b]] for (a, b) in prs] for (u, v) in prs]

    def mul(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(d)) % prime for j in range(d)] for i in range(d)]

    def rank(M):
        M = [r[:] for r in M]
        rk = 0
        for c in range(len(M[0])):
            piv = next((r for r in range(rk, len(M)) if M[r][c] % prime), None)
            if piv is None:
                continue
            M[rk], M[piv] = M[piv], M[rk]
            inv = pow(M[rk][c], prime - 2, prime)
            M[rk] = [v * inv % prime for v in M[rk]]
            for r in range(len(M)):
                if r != rk and M[r][c] % prime:
                    f = M[r][c]
                    M[r] = [(a - f * b) % prime for a, b in zip(M[r], M[rk])]
            rk += 1
        return rk

    def cells(m, s):   # sum_lk m[l][k] |u_l + s v_l><u_k + s v_k| on the cells
        X = {}
        for l, (u, v) in enumerate(prs):
            for k, (a, b) in enumerate(prs):
                if m[l][k]:
                    for r, sr in ((u, 1), (v, s)):
                        for c, sc in ((a, 1), (b, s)):
                            X[(r, c)] = (X.get((r, c), 0) + m[l][k] * sr * sc) % prime
        return X

    hs = {s: restrict(s) for s in (1, -1)}
    pw = {}
    for s in (1, -1):
        out = [[[int(i == j) for j in range(d)] for i in range(d)]]
        for _ in range(d - 1):
            out.append(mul(out[-1], hs[s]))
        pw[s] = out
    cyc = [rank([[m[i][j] for i in range(d) for j in range(d)] for m in pw[s]]) for s in (1, -1)]
    CE, CO = [cells(m, 1) for m in pw[1]], [cells(m, -1) for m in pw[-1]]
    dep = lambda r, c: n8 - 2 * bin(r ^ c).count("1")
    K = [[sum(dep(r, c) * v * co[(r, c)] for (r, c), v in ce.items() if (r, c) in co) % prime for co in CO] for ce in CE]
    return d, cyc, rank(K), [[sp.Integer(x) for x in row] for row in hs[1]], [[sp.Integer(x) for x in row] for row in hs[-1]]


lam = sp.Symbol('lam')
for prime in (1000003, 998244353):
    d8, cyc, rK8, hE8, hO8 = centre_compression_n8(prime)
    check(f"p = {prime}: h_E, h_O cyclic (commutant = powers, dim {d8} each), rank of the E-O compression 34, first-order nullity at most 2",
          cyc == [d8, d8] and rK8 == 34 and d8 == 35, f"cyclic {cyc}, rank {rK8}")
cE8, cO8 = sp.Matrix(hE8).charpoly(lam).as_expr(), sp.Matrix(hO8).charpoly(lam).as_expr()
st8 = [x for x in range(256) if bin(x).count("1") == 4]
H8 = sp.zeros(70, 70)
for c, x in enumerate(st8):
    for i in range(7):
        bi, bj = (x >> i) & 1, (x >> (i + 1)) & 1
        H8[c, c] += 1 if bi == bj else -1
        if bi != bj:
            H8[st8.index(x ^ (1 << i) ^ (1 << (i + 1))), c] += 2
T4o = sp.Matrix(70, 70, lambda k, l: 1 if bin(st8[k] ^ st8[l]).count("1") == 4 else 0)
check("the orbit sum of the Hamming-4 cells of (4, 4) commutes with H (exact), so it is in the kernel and the rank is 34 over Q",
      (H8 * T4o - T4o * H8).is_zero_matrix)
check("h_E and h_O share no eigenvalue and both spectra are simple (exact characteristic polynomials)",
      sp.gcd(cE8, cO8) == 1 and sp.gcd(cE8, sp.diff(cE8, lam)) == 1 and sp.gcd(cO8, sp.diff(cO8, lam)) == 1)

print("\nALL GATES PASS" if not FAILS else "\nFAILED: " + ", ".join(FAILS))
