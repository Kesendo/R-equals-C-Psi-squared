"""Why the crack axis and the anisotropy axis carry the same seat-blindness locus.

Companion proof: docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md
Companion page:  experiments/THE_BLIND_SEAT_ON_THE_ROAD.md  (whose section (c) this closes)

WHAT IS GATED, in the order the proof needs it.

(L1) THE SECTOR REDUCTION. On the reflection-EVEN subspace the wrap bond
     P = E_{0,N-1} + E_{N-1,0} acts exactly as the end-pair diagonal D = E_{00} + E_{N-1,N-1},
     and on the reflection-ODD subspace exactly as -D, because psi_{N-1} = +-psi_0 there. So

         (A + u P)|even = (A + u D)|even        (A + u P)|odd = (A - u D)|odd

     as matrices, at every N and symbolically in u. The folded half of this is already owned
     and gated: PROOF_CRACKED_RING_EXACT_CURVE Corollary B folds the cracked ring along R into
     two Jacobi blocks carrying +-u at the crack end (its gate P9, helper folded_blocks:
     e[0, 0] = u, o[0, 0] = -u). What L1 adds is the one comparison line: the SAME fold applied
     to the end-pair-anisotropic OPEN chain gives the same two blocks at +-Delta.

(L2) THE LOCUS DECOMPOSITION, and it is where the seat's striking is disposed of. F157's count
     is defined by striking row and column j, and striking a non-central seat destroys the
     reflection, so the sector reduction does not transfer to the struck matrix. It does not
     have to. By the fence-free Cramer argument of THE_SEAT_THAT_CUTS section 7 a seat is blind
     at lambda exactly when some eigenvector at lambda vanishes at j, a statement about
     ker(H - lambda), which IS reflection-invariant however j sits. Splitting that kernel gives

         locus(j) = E(j)  u  O(j)  u  C

     with E, O the sector-internal vanishing loci and C the DEGENERACY set, the knob values at
     which H has any eigenvalue of multiplicity at least two. C carries no seat: a 2-dimensional
     kernel contains a nonzero vector vanishing at ANY single prescribed site.

(C)  THE TWO DEGENERACY SETS, and the whole difference between the axes sits here. The
     anisotropic open chain is tridiagonal with every off-diagonal entry 1, an unreduced Jacobi
     matrix, so its spectrum is simple (Lemma J of the node-lemma proof) and C_Delta is EMPTY.
     The crack's wrap entry lifts it out of that form and C_u = {+1, -1} exactly. Both are read
     off the same solved object, the real root set of disc_x(chi(H(t))) in t. That object counts
     ALGEBRAIC multiplicity while the proof's C is the GEOMETRIC one, dim ker >= 2; for real t
     both families are real symmetric, so the two coincide and the discriminant is entitled to
     certify a kernel statement. Over the complex numbers they part (at N = 3, u = 2*sqrt(2)*i
     the discriminant vanishes with dim ker = 1, an exceptional point), which is why the reality
     of t is named here rather than assumed.

(B)  WHY C_u IS THE RING ENDS, from the boundary system alone. Writing the interior recursion's
     solution as psi_l = psi_1*U_{l-1} - psi_0*U_{l-2} at U = Chebyshev-U(lambda/2), the two
     boundary rows become a 2x2 system in (psi_0, psi_1),

         [ -(lambda + u*U_{N-3})     1 + u*U_{N-2} ]
         [    U_{N-2} + u              -U_{N-1}    ]

     and dim ker(H_u - lambda) >= 2 iff all four entries vanish. Then
     u*(1 + u*U_{N-2}) - (U_{N-2} + u) = U_{N-2}*(u^2 - 1) forces u^2 = 1 or U_{N-2} = 0, and
     U_{N-2} = 0 forces u = 0 from one entry and 1 = 0 from the other. Conversely at
     lambda = 2cos(m*pi/N) one has U_{N-2} = (-1)^(m+1), so u = +1 degenerates at the even m and
     u = -1 at the odd, in counts floor((N-1)/2) and ceil((N-1)/2), both at least 1 from N = 3.
     This replaces an earlier route through F160's simplicity clause, which is fenced to u >= 0
     while the whole u axis here includes the negative half.

(P)  PARITY, which is the only place N's parity enters. With Sigma = diag((-1)^l),
     Sigma H Sigma = -H(-t) holds for the anisotropy at every N (the companion gate's C5b) and
     for the crack at ODD N only (its C5), and Sigma R = (-1)^(N-1) R Sigma. At ODD N Sigma
     preserves each reflection sector, so E and O are EACH negation-closed; at EVEN N it SWAPS
     them, so O = -E and neither need be closed.

(T)  THE THEOREM. From L1, E_u = E_Delta and O_u(t) = O_Delta(-t). Hence
       odd N:   u-locus = E u O u {+-1} = Delta-locus u {+-1}, at every N and every seat.
       even N:  u-locus = E u {+-1} and Delta-locus = E u (-E), so the relation holds exactly
                when E u {+-1} is negation-closed -- which turns the companion page's MEASURED
                biconditional (C8: every break is a non-closure) into a consequence, and
                predicts WHICH seats break. T4 checks that prediction against the companion
                gate's committed 22-seat list, taken as a LITERAL, and T5 checks that the
                ring-ends clause in the criterion is load-bearing rather than decorative.

Every locus is SOLVED, never sampled: each is the real root set of a resultant or a
discriminant, compared as the set of its irreducible rational factors that carry a real root,
decided by Sturm counting. No root object is built and no trigonometric value is ever handed to
a simplifier: sp.roots on an irreducible cubic reports .is_real as None and sp.simplify does not
reduce P(2*cos(3*pi/7)) to 0, two defects docs/CAUGHT_ERRORS.md records for the companion gate,
and a third was met while writing this one (sp.solve returned an empty solution set at N = 11
for a system whose solutions the discriminant route exhibits). What is NOT gated here:
multiplicity. Every locus below is a SET, F157's blind COUNT at a given Delta is a multiplicity,
and the proof claims nothing about it.
"""
import sympy as sp

t, x, lam, u = sp.symbols('t x lam u')
_fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else '**FAIL**'}] {name}  {detail}")
    if not ok:
        _fails.append(name)


# ----------------------------------------------------------------- the families
def path_sym(N):
    return sp.Matrix(N, N, lambda a, b: 1 if abs(a - b) == 1 else 0)


def H_crack(N, knob=t):
    """The cracked ring: the wrap bond carries u. F160's road, u=0 open chain, u=1 ring."""
    M = path_sym(N)
    M[0, N - 1] += knob
    M[N - 1, 0] += knob
    return M


def H_aniso(N, knob=t):
    """The end-pair-anisotropic open chain, in the scale THE_BLIND_SEAT_ON_THE_ROAD uses."""
    M = path_sym(N)
    M[0, 0] += knob
    M[N - 1, N - 1] += knob
    return M


def H_oneend(N, knob=t):
    """CONTROL: the same diagonal on ONE end only, so it does not commute with the reflection."""
    M = path_sym(N)
    M[0, 0] += knob
    return M


def reflection(N):
    return sp.Matrix(N, N, lambda a, b: 1 if a + b == N - 1 else 0)


# ----------------------------------------------------------------- sector algebra
def sector_basis(N, parity):
    """Columns of the R-eigenbasis: pair p <-> {p, N-1-p}, plus the fixed centre when odd."""
    cols, norms = [], []
    for p in range(N // 2):
        v = sp.zeros(N, 1)
        v[p] = 1
        v[N - 1 - p] = 1 if parity == +1 else -1
        cols.append(v)
        norms.append(2)
    if N % 2 == 1 and parity == +1:
        v = sp.zeros(N, 1)
        v[(N - 1) // 2] = 1
        cols.append(v)
        norms.append(1)
    if not cols:
        return None, None
    return sp.Matrix.hstack(*cols), norms


def block(H, N, parity):
    """H restricted to a reflection sector, in the pair basis. Rational entries.

    Not symmetric: the pair basis is orthogonal but not orthonormal, so the block is the
    conjugate of a real symmetric matrix by the positive diagonal diag(sqrt(norm)). That
    conjugation changes neither the characteristic polynomial nor which coordinate of an
    eigenvector vanishes, which is all the criterion below reads.
    """
    S, norms = sector_basis(N, parity)
    if S is None:
        return None
    return sp.expand(sp.diag(*[sp.Rational(1, n) for n in norms]) * S.T * H * S)


def strike(M, k):
    idx = [i for i in range(M.rows) if i != k]
    return M[idx, idx]


def charpoly_expr(M):
    return sp.Poly(M.charpoly(x).as_expr(), x)


# ----------------------------------------------------------------- exact locus sets
def real_factors(poly):
    """Monic irreducible rational factors of poly that carry at least one real root.

    Two distinct irreducible rational polynomials share no root, so equality of these sets IS
    equality of the real root sets. Sturm counting only (sp.Poly.count_roots); no root built.
    """
    if poly is None or poly.is_zero:
        return 'ALL'
    out = set()
    for f, _m in sp.factor_list(poly.as_expr(), t)[1]:
        pf = sp.Poly(f, t)
        if pf.degree() >= 1 and pf.count_roots() > 0:
            out.add(sp.factor(sp.monic(pf).as_expr()))
    return frozenset(out)


def resultant_locus(M1, M2):
    """The knob values at which two matrix families share an eigenvalue, as a factor set."""
    R = sp.expand(sp.resultant(charpoly_expr(M1).as_expr(), charpoly_expr(M2).as_expr(), x))
    if R == 0:
        return 'ALL'
    return real_factors(sp.Poly(R, t))


def locus_full(Hf, N, j):
    H = Hf(N)
    return resultant_locus(H, strike(H, j))


def locus_sector(Hf, N, j, parity):
    """The seat's blindness locus INSIDE one reflection sector.

    A sector vector's amplitude at site j is its pair-basis coordinate min(j, N-1-j), so the
    same strike-and-resultant criterion applies inside the block. The centre seat of an odd
    chain has no odd-sector coordinate: every odd mode vanishes there, hence ALL.
    """
    B = block(Hf(N), N, parity)
    if B is None:
        return frozenset()
    if N % 2 == 1 and j == (N - 1) // 2:
        if parity == -1:
            return 'ALL'
        jr = N // 2
    else:
        jr = min(j, N - 1 - j)
    return resultant_locus(B, strike(B, jr))


def degeneracy(Hf, N):
    """C: the knob values at which H(t) has an eigenvalue of multiplicity >= 2. Seat-free."""
    d = sp.discriminant(charpoly_expr(Hf(N)).as_expr(), x)
    if sp.expand(d) == 0:
        return 'ALL'
    return real_factors(sp.Poly(sp.expand(d), t))


RING = frozenset({sp.factor(t - 1), sp.factor(t + 1)})


def union(*sets):
    if any(s == 'ALL' for s in sets):
        return 'ALL'
    out = set()
    for s in sets:
        out |= set(s)
    return frozenset(out)


def negate(fs):
    if fs == 'ALL':
        return fs
    return frozenset(sp.factor(sp.monic(sp.Poly(sp.expand(f.subs(t, -t)), t)).as_expr())
                     for f in fs)


# ----------------------------------------------------------------- exact blind count
def blind_count(Hnum, j):
    """F157's definition on a numeric matrix: deg gcd(chi(H), chi(H_j)), over the rationals."""
    a = sp.Poly(Hnum.charpoly(x).as_expr(), x)
    b = sp.Poly(strike(Hnum, j).charpoly(x).as_expr(), x)
    return sp.Poly(sp.gcd(a, b), x).degree()


# ----------------------------------------------------------------- the boundary system
def cheb_U(n):
    """U_n(lambda/2), with the two conventional negative-index values the recursion needs."""
    if n == -1:
        return sp.Integer(0)
    if n == -2:
        return sp.Integer(-1)
    return sp.chebyshevu(n, lam / 2)


def boundary_matrix(N):
    """The 2x2 system the crack's two boundary rows impose on (psi_0, psi_1).

    psi_l = psi_1*U_{l-1} - psi_0*U_{l-2} solves the interior recursion; substituting it into
    row 0 and row N-1 of (H_u - lambda)psi = 0 leaves exactly this.
    """
    a, b = sp.symbols('a b')

    def psi(l):
        return b * cheb_U(l - 1) - a * cheb_U(l - 2)

    rows = [sp.expand(psi(1) + u * psi(N - 1) - lam * psi(0)),
            sp.expand(psi(N - 2) + u * psi(0) - lam * psi(N - 1))]
    return sp.Matrix([[sp.expand(r.coeff(v)) for v in (a, b)] for r in rows])


# ================================================================== the gate
# The 22 even-N seats at which THE_BLIND_SEAT_ON_THE_ROAD's relation breaks, copied as a
# LITERAL from that page's gate C4 run
# (simulations/results/blind_seat_on_the_road/blind_seat_on_the_road_run.txt). T4 predicts this
# set from the proof and must reproduce it; nothing below derives it from a measurement.
COMMITTED_BREAKS = {(8, 1), (8, 6), (10, 1), (10, 8),
                    (12, 1), (12, 2), (12, 3), (12, 8), (12, 9), (12, 10),
                    (14, 1), (14, 2), (14, 11), (14, 12),
                    (16, 1), (16, 2), (16, 3), (16, 4), (16, 11), (16, 12), (16, 13), (16, 14)}

# T5's expectation: the eight seats the criterion adds when the ring-ends clause is dropped.
# A literal, so that a criterion which silently stopped depending on the clause would go red.
RING_CLAUSE_EXTRAS = {(6, 1), (6, 4), (10, 2), (10, 7), (12, 4), (12, 7), (14, 3), (14, 10)}

# L2b's expectations: how many of the 60 crack seats notice each dropped summand. Literals for
# the same reason; a threshold would pass at one seat while the proof quotes three numbers.
L2B_EXPECTED = {"C": 54, "O": 10, "E": 14}

L1_N = range(3, 15)
L2_N = range(4, 12)
C_N = range(3, 15)
ODD = [5, 7, 9, 11, 13]
EVEN = [6, 8, 10, 12, 14, 16]


def gate_L1():
    print("\n(L1) the sector reduction")

    bad = []
    for N in L1_N:
        for parity, sign in ((+1, +1), (-1, -1)):
            lhs = block(H_crack(N), N, parity)
            if lhs is None:
                continue
            rhs = block(H_aniso(N, sign * t), N, parity)
            if sp.expand(lhs - rhs) != sp.zeros(*lhs.shape):
                bad.append((N, parity))
    check("L1a (A+uP)|even == (A+uD)|even and (A+uP)|odd == (A-uD)|odd, symbolic in u",
          not bad, f"exact zero matrix at every N = {L1_N.start}..{L1_N.stop - 1}, both sectors"
                   if not bad else f"FAILS at {bad}")

    # Anti-vacuity: the same code path with the sign the reduction does NOT have must be nonzero
    # wherever the odd sector exists. A zero-assert alone would pass on an empty comparison.
    wrong = []
    for N in L1_N:
        lhs = block(H_crack(N), N, -1)
        if lhs is None:
            continue
        if sp.expand(lhs - block(H_aniso(N, +t), N, -1)) == sp.zeros(*lhs.shape):
            wrong.append(N)
    check("L1b the WRONG sign in the odd sector disagrees, at every N",
          not wrong, "(A+uP)|odd != (A+uD)|odd, so L1a is not comparing a matrix with itself"
                     if not wrong else f"agrees at {wrong}, L1a is vacuous there")

    # The MECHANISM, not merely a difference of values: both P and D commute with the reflection
    # and a one-end diagonal does not, which is why only the first two have a sector reduction
    # at all. Checked as the commutator on the full matrix, through the same door.
    mech = []
    for N in L1_N:
        R = reflection(N)
        for name, Hf, want_commute in (("crack", H_crack, True),
                                       ("aniso", H_aniso, True),
                                       ("one-end", H_oneend, False)):
            H = Hf(N)
            commutes = sp.expand(R * H * R - H) == sp.zeros(N, N)
            if commutes != want_commute:
                mech.append((N, name))
    check("L1c the reflection is what does it: R commutes with P and D and NOT with a one-end "
          "diagonal", not mech,
          f"both directions, every N = {L1_N.start}..{L1_N.stop - 1}"
          if not mech else f"FAILS at {mech}")


def gate_L2():
    print("\n(L2) the locus decomposition, and (C) the two degeneracy sets")

    for label, Hf in (("crack", H_crack), ("aniso", H_aniso)):
        bad, seats = [], 0
        for N in L2_N:
            C = degeneracy(Hf, N)
            for j in range(N):
                seats += 1
                if locus_full(Hf, N, j) != union(locus_sector(Hf, N, j, +1),
                                                 locus_sector(Hf, N, j, -1), C):
                    bad.append((N, j))
        check(f"L2a {label}: locus(j) == E(j) u O(j) u C, every seat",
              not bad, f"{seats} seats over N = {L2_N.start}..{L2_N.stop - 1}"
                       if not bad else f"FAILS at {bad}")

    # Mutation: drop a summand and count how many seats notice. The counts are LITERALS, so a
    # regression that left one seat disagreeing would go red rather than pass a threshold.
    for drop in ("C", "O", "E"):
        still = []
        for N in L2_N:
            C = degeneracy(H_crack, N)
            for j in range(N):
                E = locus_sector(H_crack, N, j, +1)
                O = locus_sector(H_crack, N, j, -1)
                parts = {"C": (E, O), "O": (E, C), "E": (O, C)}[drop]
                if locus_full(H_crack, N, j) == union(*parts):
                    still.append((N, j))
        got = sum(L2_N) - len(still)
        check(f"L2b mutation: dropping {drop} breaks the identity at exactly "
              f"{L2B_EXPECTED[drop]} of {sum(L2_N)} crack seats",
              got == L2B_EXPECTED[drop], f"disagrees at {got}")

    # C: the anisotropy is unreduced Jacobi, so no degeneracy at any Delta; the crack is not,
    # and degenerates exactly at the ring ends. Both read off the SAME solved discriminant, so
    # the empty answer is not a different code path from the nonempty one.
    empty, ends = [], []
    for N in C_N:
        if degeneracy(H_aniso, N) != frozenset():
            empty.append(N)
        if degeneracy(H_crack, N) != RING:
            ends.append(N)
    check("C1 C_Delta is EMPTY at every N", not empty,
          "the anisotropic chain is tridiagonal with nonzero off-diagonals, an unreduced Jacobi "
          f"matrix, so its spectrum is simple; N = {C_N.start}..{C_N.stop - 1}"
          if not empty else f"nonempty at {empty}")
    check("C2 C_u is exactly {+1, -1} at every N", not ends,
          f"the crack's wrap entry leaves the Jacobi form; N = {C_N.start}..{C_N.stop - 1}"
          if not ends else f"differs at {ends}")

    # The taxonomy entry the degeneracy set needs: a doubly degenerate eigenvalue blinds EVERY
    # seat, because a 2-dimensional eigenspace contains a vector vanishing at any single site.
    # Read on the ring (u=1, every level paired) against a generic u on the same code path.
    ring_bad, generic_all_blind = [], []
    for N in range(5, 12):
        pairs = (N - 1) // 2 if N % 2 else (N - 2) // 2
        Hring = H_crack(N, sp.Integer(1))
        if any(blind_count(Hring, j) < pairs for j in range(N)):
            ring_bad.append(N)
        Hgen = H_crack(N, sp.Rational(1, 3))
        if all(blind_count(Hgen, j) > 0 for j in range(N)):
            generic_all_blind.append(N)
    check("C3 at u = 1 every seat is blind by at least the number of degenerate pairs",
          not ring_bad, "N = 5..11, exact deg gcd" if not ring_bad else f"FAILS at {ring_bad}")
    check("C3b anti-vacuity: at u = 1/3 the seats are not all blind",
          not generic_all_blind,
          "so C3 reads the degeneracy and not the code path"
          if not generic_all_blind else f"all blind at {generic_all_blind}")


def gate_B():
    print("\n(B) the boundary system, which is why C_u is the ring ends")

    # The reduction from the N x N matrix to the 2 x 2 is what the whole lemma rests on, and
    # boundary_matrix() never touches H_crack. One determinant ties them: it pins both boundary
    # rows, the Chebyshev ansatz and the kernel isomorphism at once, against the real matrix.
    tie = []
    for N in C_N:
        chi = charpoly_expr(H_crack(N)).as_expr().subs({t: u, x: lam})
        if sp.expand(sp.det(boundary_matrix(N)) - chi) != 0:
            tie.append(N)
    check("B0 det(the 2x2) == chi(H_u), so the reduction is tied to the N x N matrix",
          not tie, f"exact, symbolic in lambda and u, N = {C_N.start}..{C_N.stop - 1}"
                   if not tie else f"FAILS at {tie}")

    bad = []
    for N in C_N:
        M = boundary_matrix(N)
        want = sp.expand(cheb_U(N - 2) * (u**2 - 1))
        if sp.expand(sp.expand(u * M[0, 1]) - sp.expand(M[1, 0]) - want) != 0:
            bad.append(N)
    check("B1 u*(entry 01) - (entry 10) == U_{N-2}*(u^2 - 1), symbolic in lambda and u",
          not bad, f"so a common zero of the two forces u^2 = 1 or U_{{N-2}} = 0; N = "
                   f"{C_N.start}..{C_N.stop - 1}" if not bad else f"FAILS at {bad}")

    # The U_{N-2} = 0 branch is not a second solution: entry 10 then gives u = 0 and entry 01
    # gives 1 = 0. Checked as the substitution rather than argued.
    branch = []
    for N in C_N:
        M = boundary_matrix(N)
        # On the branch U_{N-2} = 0 entry 01 reads 1 and entry 10 reads u, checked by removing
        # the U_{N-2} term from each rather than by substituting into a simplifier.
        red01 = sp.expand(M[0, 1] - u * cheb_U(N - 2))
        red10 = sp.expand(M[1, 0] - cheb_U(N - 2))
        if red01 != 1 or red10 != u:
            branch.append(N)
    check("B1b on the branch U_{N-2} = 0 the two entries read 1 and u, so the branch is empty",
          not branch, f"N = {C_N.start}..{C_N.stop - 1}" if not branch else f"FAILS at {branch}")

    # The converse, as an exact multiplicity count. At lambda = 2cos(m*pi/N) one has
    # U_{N-2} = (-1)^(m+1), so u = +1 degenerates at the even m and u = -1 at the odd. The
    # expected multiplicities are combinatorial in N and are not read off the discriminant.
    mult_bad = []
    for N in C_N:
        d = sp.Poly(sp.expand(sp.discriminant(charpoly_expr(H_crack(N)).as_expr(), x)), t)
        got = {str(f): m for f, m in sp.factor_list(d.as_expr(), t)[1]}
        want_plus, want_minus = 2 * ((N - 1) // 2), 2 * ((N // 2))
        if got.get('t - 1') != want_plus or got.get('t + 1') != want_minus:
            mult_bad.append((N, got.get('t - 1'), want_plus, got.get('t + 1'), want_minus))
    check("B2 the ring ends are attained, with multiplicity 2*floor((N-1)/2) at u = +1 and "
          "2*ceil((N-1)/2) at u = -1",
          not mult_bad, f"in disc_x(chi(H_u)); N = {C_N.start}..{C_N.stop - 1}"
                        if not mult_bad else f"FAILS at {mult_bad}")

    # B2 reads the conclusion off the discriminant; B3 runs the converse in the boundary system
    # itself, which is what the proof argues, and it is the only gate that reads entry (0,0).
    # Exact polynomial reduction modulo U_{N-1}: no trigonometric value meets a simplifier.
    conv = []
    for N in C_N:
        M = boundary_matrix(N).subs(u, -cheb_U(N - 2))
        Un1 = sp.Poly(sp.expand(cheb_U(N - 1)), lam)
        if any(sp.rem(sp.Poly(sp.expand(M[i, j]), lam), Un1).as_expr() != 0
               for i in range(2) for j in range(2)):
            conv.append(N)
    check("B3 the converse in the boundary system: at u = -U_{N-2} all FOUR entries vanish "
          "modulo U_{N-1}", not conv,
          f"exact reduction, no root object and no simplifier; N = "
          f"{C_N.start}..{C_N.stop - 1}" if not conv else f"FAILS at {conv}")


def gate_P():
    print("\n(P) parity: what Sigma does to the sectors")

    odd_bad, even_bad, swap_bad = [], [], []
    for N in ODD + EVEN:
        for j in range(N):
            E = locus_sector(H_aniso, N, j, +1)
            O = locus_sector(H_aniso, N, j, -1)
            if E == 'ALL' or O == 'ALL':
                continue
            if N % 2:
                if E != negate(E) or O != negate(O):
                    odd_bad.append((N, j))
            else:
                if O != negate(E):
                    swap_bad.append((N, j))
                if E != negate(E):
                    even_bad.append((N, j))
    check("P1 at ODD N each sector locus is negation-closed on its own",
          not odd_bad, f"N in {ODD}, every NON-CENTRE seat (40 of 45; at the centre O is all of "
                       f"R and there is nothing to close)" if not odd_bad
                       else f"FAILS at {odd_bad}")
    check("P2 at EVEN N Sigma SWAPS the sectors: O == -E",
          not swap_bad, f"N in {EVEN}, every seat" if not swap_bad else f"FAILS at {swap_bad}")
    check("P3 and at EVEN N the sectors are NOT separately closed",
          bool(even_bad),
          f"E != -E at {len(even_bad)} even-N seats, so P1 is a statement about odd N and not "
          f"about every N")


def gate_T():
    print("\n(T) the theorem")

    # T1 is a CONSISTENCY check, not independent evidence. E_u = E_Delta is L1a: the two even
    # blocks are the same matrix, so this asserts nothing L1a has not already fixed. What it
    # does test is the odd leg, that the locus of a family evaluated at -t is the negation of
    # its locus, i.e. that real_factors and negate commute with t -> -t.
    bad = []
    for N in ODD + EVEN:
        for j in range(N):
            Eu, Ed = locus_sector(H_crack, N, j, +1), locus_sector(H_aniso, N, j, +1)
            Ou, Od = locus_sector(H_crack, N, j, -1), locus_sector(H_aniso, N, j, -1)
            if Eu != Ed or Ou != negate(Od):
                bad.append((N, j))
    check("T1 consistency: E_u == E_Delta (this is L1a) and O_u(t) == O_Delta(-t)",
          not bad, f"N in {ODD + EVEN}; the second leg is the informative one"
                   if not bad else f"FAILS at {bad}")

    odd_bad = []
    for N in ODD:
        for j in range(N):
            uu, dd = locus_full(H_crack, N, j), locus_full(H_aniso, N, j)
            if uu == 'ALL' or dd == 'ALL':
                if uu != dd:
                    odd_bad.append((N, j))
                continue
            if uu != union(dd, RING):
                odd_bad.append((N, j))
    check("T2 odd N: u-locus == Delta-locus u {+-1}, every seat",
          not odd_bad, f"N in {ODD}, {sum(ODD)} seats" if not odd_bad else f"FAILS at {odd_bad}")

    predicted, no_ring, struct = set(), set(), []
    for N in EVEN:
        for j in range(N):
            E = locus_sector(H_aniso, N, j, +1)
            O = locus_sector(H_aniso, N, j, -1)
            if E == 'ALL' or O == 'ALL':
                continue
            if locus_full(H_crack, N, j) != union(E, RING):
                struct.append(("u", N, j))
            if locus_full(H_aniso, N, j) != union(E, O):
                struct.append(("D", N, j))
            if union(E, RING) != negate(union(E, RING)):
                predicted.add((N, j))
            if E != negate(E):
                no_ring.add((N, j))
    check("T3 even N: u-locus == E u {+-1} and Delta-locus == E u O",
          not struct, f"N in {EVEN}, every seat" if not struct else f"FAILS at {struct}")
    check("T4 even N: the seats the proof predicts to break ARE the committed 22",
          predicted == COMMITTED_BREAKS,
          f"{len(predicted)} predicted, {len(COMMITTED_BREAKS)} committed in "
          f"blind_seat_on_the_road_run.txt gate C4"
          + ("" if predicted == COMMITTED_BREAKS
             else f"; only predicted {sorted(predicted - COMMITTED_BREAKS)}, "
                  f"only committed {sorted(COMMITTED_BREAKS - predicted)}"))
    check("T5 the ring-ends clause in the criterion is load-bearing: dropping it adds exactly "
          "the 8 seats where the ring ends restore closure",
          no_ring - predicted == RING_CLAUSE_EXTRAS and predicted < no_ring,
          f"{len(no_ring)} predicted without the clause against {len(predicted)} with it, the "
          f"extras {sorted(no_ring - predicted)}")


def main():
    print("=" * 86)
    print("PROOF_BLIND_SEAT_TWO_AXES: the crack axis and the anisotropy axis, sector by sector")
    print("=" * 86)
    gate_L1()
    gate_L2()
    gate_B()
    gate_P()
    gate_T()
    print()
    print("=" * 86)
    print("VERDICT:", "ALL GREEN" if not _fails else f"{len(_fails)} FAILED: {_fails}")
    print("=" * 86)
    return 1 if _fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
