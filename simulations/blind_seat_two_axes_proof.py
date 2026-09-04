"""Why the crack axis and the anisotropy axis carry one seat-blindness locus, and at odd N one count.

Companion proof: docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md
Companion page:  experiments/THE_BLIND_SEAT_ON_THE_ROAD.md  (whose section (c) and whose open
                 item 1 this closes: the locus in blocks L to T, the count in block K, the
                 index on that count in block S, and the two constants of the proof's sections
                 (g) and (h) in block W)

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

(K)  THE COUNT, which is what the loci do not settle. Striking the seat splits chi(H_j) along
     the same two sectors,

         chi(H_j) = 1/2 * ( chi(E_jr)*chi(O) + chi(E)*chi(O_jr) ),   jr = min(j, N-1-j)

     (and chi(E_jr)*chi(O) at the reflection-fixed centre seat, where e_j has no odd part), so
     off the degeneracy set the count is additive, blind(j) = b_E + b_O, and the two sector
     terms are exactly what L1 and P transport. Hence
       odd N, |t| != 1:   blind_u(j, t) = blind_Delta(j, t) at every seat: the two axes carry
                the same points AND the same count.
       even N, |t| != 1:  blind_u = 2*b_E while blind_Delta = b_E(t) + b_E(-t), so the ring pays
                the even sector twice and the counts agree only where that sector pays alike at
                +-t. Where the two LOCI break the two COUNTS break with them, provably; the
                converse is measured over N = 6..14 and not proved.
       |t| = 1: additivity picks up its one correction. A shared eigenvalue has a
                two-dimensional eigenspace and the seat misses it only when BOTH sector
                eigenvectors vanish there, so
                    blind(j) = b_E + b_O + #{shared lambda that neither sector misses at j}
                and the ring pays floor/ceil((N-1)/2) at EVERY seat, the degenerate-pair count.
     Two by-products. The same split says WHERE F157's P_j gets squared: inside a sector the
     seat DISCONNECTS the tridiagonal block into two halves, the sector resultant is a t-free
     constant times its halves' resultant SQUARED, and on the CHAIN those two sector halves'
     resultants multiply to F157's own P_j. So the full resultant carries the count DOUBLED
     rather than not at all, and halving it is exact off the degeneracy set and away from the
     forced centre seat. And the centre seat's (N-1)/2, which the companion page could only
     measure at six couplings, follows at every coupling because a Jacobi eigenvector cannot
     vanish at an END of its chain, which is (J1) of the node-lemma proof and not new here.

(S)  WHICH k SITS IN WHICH SECTOR, which is the index block K leaves off b_E and b_O. F157
     indexes the locus by k through the node angle theta_k = k*pi/N_node. The blind mode there
     is one sine across the whole chain, psi_l = sin((j-l)*theta_k), the two halves F157 solves
     separately written as one; and since N_node*theta_k = +-k*pi,

         psi_{N-1-l} = (-1)^(k+1) * psi_l

     so ODD k sits in the reflection-EVEN sector and EVEN k in the odd one. Hence b_E counts
     the odd k landing on a Delta value and b_O the even ones, and block K's even-N law closes
     to blind_u(j, t) = 2*#{k odd : Delta_k = t} off the ring ends. The sign law is NOT new
     here: it is F71's, owned in six places, five of them on the chain's own Dirichlet comb at
     modulus N+1 and the sixth a bare parity with no modulus in it. What
     is new is that it survives the move to the node modulus on a chain whose two ends are
     DETUNED, where the six owners all read the uniform chain. Free with it: Delta_{N_node-k} =
     -Delta_k, which with N_node = N-1 mod 2 reproduces P's Lemma 5 from the index alone.

(H)  THE HOP, which is what BlindSeat.cs's 128 is. Two books. F157's committed one at uniform
     J = 1 with the ZZ coupling carrying Delta: hop 2 on every bond, and a ZZ diagonal paying
     -Delta at a bond's own two ends and +Delta elsewhere, bond by bond, i.e. Delta*(N-5) in the
     interior and Delta*(N-3) at the ends, which is SeatBlindnessDeltaLocusWitness.cs's
     (b^2 = 4, diagonal Delta*(N-5), end shift 2*Delta) and at Delta = 1 is the integer matrix of
     seat_cut_blindness.py se_hamiltonian_int (H0 compares the transcription entry for entry; the
     C# BlindSeat.H() is the same recipe, READ and not run from here). And this file's: hop 1, a
     bare t on the end pair. The two differ by a common shift Delta*(N-5)*I and the factor 2. The shift drops out of the halves-resultant outright (Res_x(f(x-c), g(x-c)) =
     Res_x(f, g)); the factor does not: a resultant of bidegree (a, b) is homogeneous of degree
     a*b under a common scaling of the two matrices, so

         Res_F157(Delta) = 2^(j(N-1-j)) * Res_here(t = Delta)        (J = 1)

     exactly at every interior seat, and at N = 9 seat 1 that is 2^(1*7) = 128. At a uniform
     J the factor is (2J)^(j(N-1-j)) (H1c, at J = 2 and 3); the roots do not move, the constant
     does. Res_F157 here is the HALVES-resultant in F157's matrix normalisation and not F157's
     generator P_j = Res(S_n, Delta*S_j - S_{j+1}), a different object whose constant depends on
     the Chebyshev book (1 in the monic 2cos book, 2^5 in the cos book at this seat) and which
     the witness returns primitive; in the monic book 128*P_j IS the right-first halves-resultant
     at N = 9 seat 1 (H5). The quotient is
     primitive in Delta at every seat read (N = 4..12), so the hop's power is the whole of the
     un-normalisation up to a sign, and the sign is an ORDER (bidegree 1*7 is odd): left-first
     the Sylvester determinant gives -128 and right-first +128, and sympy.resultant returns +128
     in BOTH orders, which is section (i)'s rule ((-1)^(deg f * deg g) times the Sylvester
     determinant when deg f < deg g) and not a third value. One more thing the block separates:
     the 2 in the end pair's 2*Delta is a BOND COUNT (an end site has one bond fewer, and the ZZ
     diagonal jumps by 2J*Delta per lost bond) while the hop's 2 is XX+YY on one bond; they agree as
     numbers under this normalisation at every Delta, and t = Delta needs nothing but h = 2. With
     hop h and the ZZ unchanged the reading is t = 2*Delta/h (H3, at h = 3).

Every locus is SOLVED, never sampled: each is the real root set of a resultant or a
discriminant, compared as the set of its irreducible rational factors that carry a real root,
decided by Sturm counting. No root object is built anywhere except in block W's W0, which
compares this file's Sylvester determinant against the DEFINITION of a resultant and therefore
has to evaluate at roots; and no trigonometric value is ever handed to a simplifier: sp.roots on an irreducible cubic reports .is_real as None and sp.simplify does not
reduce P(2*cos(3*pi/7)) to 0, two defects docs/CAUGHT_ERRORS.md records for the companion gate,
and a third was met while writing this one (sp.solve returned an empty solution set at N = 11
for a system whose solutions the discriminant route exhibits).

Blocks L, C, B, P and T compare SETS, block H two books of the same polynomial. Blocks K and S compare COUNTS, and a count at an irrational
locus point is a degree, so it is read over the field Q[t]/(mu) for mu the point's minimal
polynomial. That is exact and builds no root object, which is what keeps the three traps above
out of the count as well. Block W compares POLYNOMIALS: the two constants the count leaves
unidentified are a sign and a leading coefficient, and reading either needs a named argument
order for the resultant, which is what its first three checks establish.

Block S puts an INDEX on the two sector counts block K leaves as measurements, by reading F157's
node mode and its reflection parity. Its arithmetic is block K's, moved from the count to the
mode: every value lives in the quotient by an integer minimal polynomial, this time the one of
2*cos(pi/N_node), so the node angle never reaches a simplifier either.
"""
from math import comb

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


# ----------------------------------------------------------------- exact arithmetic in Q[t]/(mu)
# The count at an IRRATIONAL locus point is a degree, so it needs a field, not a sample. The
# field is Q[t]/(mu) for mu the point's minimal polynomial: exact, and it never builds a root
# object, which is what makes it immune to the three sympy traps this file's docstring names.

def _red(e, mu):
    return sp.rem(sp.expand(e), mu, t)


def _trim(a):
    while a and a[0] == 0:
        a = a[1:]
    return a


def _inv_mod_mu(c, mu):
    if c.is_number:
        return sp.Rational(1) / c
    return sp.invert(sp.Poly(c, t), sp.Poly(mu, t)).as_expr()


def _rem_mod_mu(a, b, mu):
    """Remainder of a by b, x-polynomials high-degree-first, coefficients in Q[t]/(mu)."""
    a, b = _trim(list(a)), _trim(list(b))
    inv = _inv_mod_mu(b[0], mu)
    b = [_red(c * inv, mu) for c in b]
    while True:
        a = _trim(a)
        if len(a) < len(b):
            return a
        f = a[0]
        for i in range(len(b)):
            a[i] = _red(a[i] - f * b[i], mu)
        a = a[1:]


def _gcd_mod_mu(a, b, mu):
    """gcd of two x-polynomials over the field Q[t]/(mu). Exact; mu must be irreducible."""
    a, b = _trim(list(a)), _trim(list(b))
    while b:
        a, b = b, _rem_mod_mu(a, b, mu)
    return a


def gcd_deg_mod_mu(a, b, mu):
    """Its degree, which is F157's count when the two arguments are chi(H) and chi(H_j)."""
    g = _gcd_mod_mu(a, b, mu)
    return len(g) - 1 if g else -1


_CP_CACHE = {}


def cp_coeffs(M, key):
    if key not in _CP_CACHE:
        _CP_CACHE[key] = sp.Poly(M.charpoly(x).as_expr(), x).all_coeffs()
    return _CP_CACHE[key]


def _fold(N, j):
    return N // 2 if (N % 2 == 1 and j == (N - 1) // 2) else min(j, N - 1 - j)


def _is_centre(N, j):
    return N % 2 == 1 and j == (N - 1) // 2


def blind_at(Hf, N, j, mu, sign=+1):
    """F157's count at a root of mu, on family Hf, exact over Q[t]/(mu). sign=-1 reads it at -t."""
    H = Hf(N, sign * t)
    a = [_red(c, mu) for c in cp_coeffs(H, (Hf.__name__, N, sign, 'H'))]
    b = [_red(c, mu) for c in cp_coeffs(strike(H, j), (Hf.__name__, N, sign, 'Hj', j))]
    return gcd_deg_mod_mu(a, b, mu)


def sector_counts_at(Hf, N, j, mu, sign=+1):
    """(b_E, b_O, cross) at a root of mu: the two sector counts and the cross-sector gcd degree."""
    H = Hf(N, sign * t)
    E, O = block(H, N, +1), block(H, N, -1)
    jr = _fold(N, j)
    aE = [_red(c, mu) for c in cp_coeffs(E, (Hf.__name__, N, sign, 'E'))]
    bE = gcd_deg_mod_mu(
        aE, [_red(c, mu) for c in cp_coeffs(strike(E, jr), (Hf.__name__, N, sign, 'Ej', jr))], mu)
    if O is None:
        return bE, 0, 0
    aO = [_red(c, mu) for c in cp_coeffs(O, (Hf.__name__, N, sign, 'O'))]
    # at the reflection-fixed centre seat EVERY odd mode vanishes, so the count is the sector
    bO = O.rows if _is_centre(N, j) else gcd_deg_mod_mu(
        aO, [_red(c, mu) for c in cp_coeffs(strike(O, jr), (Hf.__name__, N, sign, 'Oj', jr))], mu)
    return bE, bO, gcd_deg_mod_mu(aE, aO, mu)


def shared_neither_at(Hf, N, j, mu, parts=False):
    """#{lambda shared by the two sectors at which NEITHER sector eigenvector vanishes at j}.

    An eigenvalue carried by both sectors has a two-dimensional eigenspace, and each such level
    adds 1 + [both sector eigenvectors vanish at j] to the count while b_E + b_O already counts
    it once per vanishing sector; so the correction to additivity is
    #shared - #{even vanishes} - #{odd vanishes} + #{both vanish}, by inclusion-exclusion.

    The LAST term is 0 on both families here, and provably: off the ring ends neither family has
    a shared level at all, and AT the ring ends the two sector eigenvectors of a shared level are
    the cosine and the sine of one mode about the reflection centre, whose squares sum to 1, so
    they cannot both vanish at a site. Gate K2c reads that, because a term that is always zero
    can carry a sign error forever. Read as degrees of gcds, exactly, over Q[t]/(mu).
    """
    H = Hf(N, t)
    E, O = block(H, N, +1), block(H, N, -1)
    if O is None:
        return 0
    jr = _fold(N, j)
    cE = [_red(c, mu) for c in cp_coeffs(E, (Hf.__name__, N, 1, 'E'))]
    cO = [_red(c, mu) for c in cp_coeffs(O, (Hf.__name__, N, 1, 'O'))]
    cEj = [_red(c, mu) for c in cp_coeffs(strike(E, jr), (Hf.__name__, N, 1, 'Ej', jr))]
    cOj = ([sp.Integer(0)] if _is_centre(N, j)
           else [_red(c, mu) for c in cp_coeffs(strike(O, jr), (Hf.__name__, N, 1, 'Oj', jr))])
    g = _gcd_mod_mu(cE, cO, mu)
    a = gcd_deg_mod_mu(g, cEj, mu)
    b = gcd_deg_mod_mu(g, cOj, mu)
    c = gcd_deg_mod_mu(_gcd_mod_mu(g, cEj, mu), cOj, mu)
    shared = len(_trim(g)) - 1
    return (shared, a, b, c) if parts else shared - a - b + c


def locus_mus(N, j):
    """This seat's locus points on EITHER axis as minimal polynomials, the ring ends removed."""
    out = set()
    for Hf in (H_aniso, H_crack):
        L = locus_full(Hf, N, j)
        if L != 'ALL':
            out |= set(L)
    return sorted(out - RING, key=str)


# ------------------------------------------------------- the block's literals, by provenance
#
# PINNED MEASUREMENTS, read off a run and frozen here so that a change in the objects moves the
# measurement and not the expectation. They are regression pins, not predictions, and each says
# what would move it.
#
# K2B_CORRECTED_SEATS: at how many (seat, ring end) pairs the correction term is strictly below
# the shared-level count, that is where the naive sum b_E + b_O + shared overcounts. Moves with
# the sector split, the fold coordinate or the field arithmetic.
K2B_CORRECTED_SEATS = {5: 2, 6: 2, 7: 2, 8: 0, 9: 6, 10: 2, 11: 2, 12: 4}
# K1D_SQUARE_SURVIVORS: how many fold seats still give a constant ratio once the sector block
# carries one next-nearest entry. K3C_DISAGREEMENTS: over how many of K3's own triples the
# R-breaking one-end family fails to track the crack.
K1D_SQUARE_SURVIVORS = 0
K1E2_SURVIVORS = 0
# K2C_READINGS / K2C_SHARED: the population of the both-vanish check and the number of shared
# levels it walks past, so that "the term is always 0" is not a statement about an empty loop.
K2C_READINGS = 544
K2C_SHARED = 552
K3C_DISAGREEMENTS = 24
# K1B_SIGNS: the constant in K1b is not merely nonzero, it is +-1. The two counts are pinned so
# that a change moves the measurement and not the expectation, and so that the doc cannot say
# "always -1" again. WHICH sign is section (i)'s, and this literal is not the law: it is a reading
# of sympy's argument order AND of the seat order below, each of which moves six of the fifty
# seats at N = 4..12 (block W). The law itself needs a named order and is gate W3.
K1B_SIGNS = {-1: 22, 1: 10}
# K3_AT_ORIGIN: how many of K3's readings sit at mu = t, that is at the knob value 0, where the
# crack and the anisotropy are the SAME matrix and the comparison is a number against itself.
# K3's informative population is the difference.
K3_AT_ORIGIN = 20
# K3D_ODD / K3D_EVEN: the two populations of Lemma 5-prime read on the counts.
K6B_EXACT = 76
K6B_BROKEN = [(6, 1, 't - 1'), (6, 4, 't - 1'), (9, 1, 't - 1'), (9, 1, 't + 1'),
              (9, 7, 't - 1'), (9, 7, 't + 1'), (10, 2, 't - 1'), (10, 7, 't - 1')]
K3D_ODD = 44
K3D_EVEN = 22
# K1E_CROSS_CONSTANT_N: over how many N the cross-sector resultant Res(chi_E, chi_O) is t-free on
# the chain. It is measured, not derived: having no real roots would not make it constant.
K1E_CROSS_CONSTANT_N = 12
# K_EXPECTED: the population each check must reach, pinned. K4C_AGREEING_SEATS: the even-N
# seat-ends at which blind_u happens to equal 2*b_E anyway, which is what K4c's fence allows.
K_EXPECTED = {'K0': 98, 'K1': 92, 'K1b': 32, 'K1b2': 20, 'K1b3': 32, 'K1c': 184, 'K1e': 32,
              'K2': 96, 'K3': 44, 'K4': 22, 'K6': 72}
# K1B3_STANDING: the 12 ring seats of j = 1..N-2 at which K1b2 finds the chain identity standing,
# as a literal, so that the divisibility law (N - 1 - 2*jr) | jr is tested against a list and not
# against the measurement it explains. 32 - 20 of K1b2 is the same number reached from the other
# side. (Not "survivors": that noun is MirrorWorld's Survivor and F123's, fenced in this arc.)
K1B3_STANDING = [(4, 1), (4, 2), (6, 2), (6, 3), (7, 2), (7, 4), (8, 3), (8, 4),
                  (10, 3), (10, 4), (10, 5), (10, 6)]
K4C_AGREEING_SEATS = [(6, 1, 't - 1'), (6, 4, 't - 1'), (10, 2, 't - 1'), (10, 7, 't - 1')]
# EVEN_N_DIFFERENCES: the (N, seat, point) triples where the two axes pay different counts, with
# both values. The values FOLLOW from the even-N law once b_E is known (2*b_E against
# b_E(t) + b_E(-t)); as written they are read off a run.
EVEN_N_DIFFERENCES = {
    (8, 1, 't**2 - t - 1'): (1, 2), (8, 1, 't**2 + t - 1'): (1, 0),
    (10, 1, 't**3 - t**2 - 2*t + 1'): (1, 2), (10, 1, 't**3 + t**2 - 2*t - 1'): (1, 0),
}
#
# DERIVED here rather than read off THIS run, though both columns also appear in other committed
# runs (the chain values in seat_cut_blindness_run.txt, the ring constant in gate B3 of
# blind_seat_on_the_road_run.txt). (N, seat): (chain, ring) at u = +1. The ring value is
# Lemma 4's degenerate-pair count floor((N-1)/2). The chain value is the multiplicity of
# (Delta - 1) in F157's locus polynomial P_j. Two of the four keys are committed rows of F157's
# table: (9, 1) is Delta^5 - 4Delta^3 + 3Delta, carrying Delta = 1 simply, and (11, 1) is
# Delta^7 - 6Delta^5 + 10Delta^3 - 4Delta, with no root at 1. The other two come from F157's
# committed GENERATOR P_j = Res_x(U_{N_node-1}, Delta*U_{j-1} - U_j): at (11, 3) it gives Delta^3,
# again no root at 1, and (9, 7) is the reflection mirror of (9, 1), the same polynomial.
ODD_N_RING_END_DIFFERENCES = {(9, 1): (1, 4), (9, 7): (1, 4), (11, 1): (0, 5), (11, 3): (0, 5)}


# ----------------------------------------------------------------- (K) the count
def gate_K():
    print("\n(K) THE COUNT: the split of chi(H_j), and what each axis pays at a shared point")

    # -- K0: the split identity. Two genuinely different routes to one polynomial: the N x N
    # -- struck characteristic polynomial, and the two sector blocks. Symbolic in t AND x. The
    # -- second clause is what K1's factorization rests on and nothing else gated.
    bad, seen = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 11):
            H = Hf(N)
            E, O = block(H, N, +1), block(H, N, -1)
            chE = sp.expand(E.charpoly(x).as_expr())
            chO = sp.expand(O.charpoly(x).as_expr())
            if sp.expand(H.charpoly(x).as_expr() - chE * chO) != 0:
                bad.append((name, N, "chi(H) != chi(E)chi(O)"))
            for j in range(N):
                jr, seen = _fold(N, j), seen + 1
                lhs = sp.expand(strike(H, j).charpoly(x).as_expr())
                chEs = sp.expand(strike(E, jr).charpoly(x).as_expr())
                if _is_centre(N, j):
                    rhs = sp.expand(chEs * chO)
                else:
                    chOs = sp.expand(strike(O, jr).charpoly(x).as_expr())
                    rhs = sp.expand(sp.Rational(1, 2) * (chEs * chO + chE * chOs))
                if sp.expand(lhs - rhs) != 0:
                    bad.append((name, N, j))
    check("K0  chi(H) = chi(E)chi(O) and chi(H_j) = 1/2(chi(E_j)chi(O) + chi(E)chi(O_j))",
          not bad and seen == K_EXPECTED.get('K0'),
          f"symbolic in t and x, exact zero polynomial; {seen} seats over N = 4..10, both axes, "
          f"expected {K_EXPECTED.get('K0')}; failures {bad}")

    # -- K0b: BOTH sectors are load-bearing. At a seat the reflection does not fix, NEITHER
    # -- summand alone is the struck polynomial; at the centre seat, where e_j has no odd part,
    # -- exactly one of them IS. Two verdicts through one door, so neither is a code path.
    wrong, centres = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 11):
            H = Hf(N)
            E, O = block(H, N, +1), block(H, N, -1)
            chE, chO = sp.expand(E.charpoly(x).as_expr()), sp.expand(O.charpoly(x).as_expr())
            for j in range(N):
                jr = _fold(N, j)
                lhs = sp.expand(strike(H, j).charpoly(x).as_expr())
                even_only = sp.expand(lhs - sp.expand(strike(E, jr).charpoly(x).as_expr()) * chO)
                odd_only = sp.expand(lhs - chE * sp.expand(strike(O, jr).charpoly(x).as_expr()))
                if _is_centre(N, j):
                    centres += 1
                    if even_only != 0 or odd_only == 0:
                        wrong.append((name, N, j, "centre"))
                elif even_only == 0 or odd_only == 0:
                    wrong.append((name, N, j, "one summand sufficed"))
    check("K0b both sectors are load-bearing off the centre, and exactly the even one at it",
          not wrong and centres == 6,
          f"neither summand alone reproduces chi(H_j) at any of the {K_EXPECTED.get('K0') - centres} "
          f"non-centre seats, and the even one alone does at all {centres} centre seats, expected "
          f"6; failures {wrong}")

    # -- K0c: the mutation is on the OBJECT. H_oneend does not commute with R, so it has no
    # -- sector split; fed through the SAME door, K0 must break.
    broke = []
    for N in range(4, 11):
        H = H_oneend(N)
        E, O = block(H, N, +1), block(H, N, -1)
        chE, chO = sp.expand(E.charpoly(x).as_expr()), sp.expand(O.charpoly(x).as_expr())
        jr = _fold(N, 1)
        lhs = sp.expand(strike(H, 1).charpoly(x).as_expr())
        rhs = sp.expand(sp.Rational(1, 2) * (
            sp.expand(strike(E, jr).charpoly(x).as_expr()) * chO
            + chE * sp.expand(strike(O, jr).charpoly(x).as_expr())))
        broke.append(sp.expand(lhs - rhs) != 0)
    check("K0c mutation on the OBJECT: the one-end diagonal, which breaks R, breaks K0",
          all(broke) and len(broke) == 7,
          f"through the same door at {sum(broke)} of {len(broke)} N in 4..10, expected 7 of 7")

    # -- K1: the resultant factorization. It is a CONSEQUENCE of K0 and the multiplicativity of
    # -- the resultant, constant and sign included, and fed random polynomials with no matrix
    # -- behind them it holds; so it is a consistency check on the algebra, not independent
    # -- evidence. What it earns is that the displayed constant is the one the algebra gives.
    bad, seen = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 11):
            H = Hf(N)
            E, O = block(H, N, +1), block(H, N, -1)
            cE, cO = charpoly_expr(E).as_expr(), charpoly_expr(O).as_expr()
            dE, dO = E.rows, O.rows
            CX = sp.expand(sp.resultant(cE, cO, x))
            for j in range(N):
                if _is_centre(N, j):
                    continue
                jr, seen = _fold(N, j), seen + 1
                full = sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                              charpoly_expr(strike(H, j)).as_expr(), x))
                RE = sp.expand(sp.resultant(cE, charpoly_expr(strike(E, jr)).as_expr(), x))
                RO = sp.expand(sp.resultant(cO, charpoly_expr(strike(O, jr)).as_expr(), x))
                pred = sp.expand((-1) ** (dE * dO) * RE * RO * CX ** 2)
                if sp.expand(2 ** N * full - pred) != 0:
                    bad.append((name, N, j))
    check("K1  2^N Res(chi_H, chi_H_j) = (-1)^(dE dO) R_E R_O Res(chi_E, chi_O)^2",
          not bad and seen == K_EXPECTED.get('K1'),
          f"exact polynomial identity in t, {seen} seats over N = 4..10, both axes, expected "
          f"{K_EXPECTED.get('K1')}; failures {bad}")

    # -- K1b: the bridge to F157. Its definition-route polynomial is the resultant of the two
    # -- halves the seat cuts the FULL chain into; that polynomial is, up to sign, the PRODUCT of
    # -- the same resultants taken inside the two sectors.
    bad, seen, signs = [], 0, {}
    for N in range(4, 11):
        H = H_aniso(N)
        for j in range(N):
            if _is_centre(N, j) or j in (0, N - 1):
                continue
            seen += 1
            Pfull = sp.expand(sp.resultant(charpoly_expr(H[:j, :j]).as_expr(),
                                           charpoly_expr(H[j + 1:, j + 1:]).as_expr(), x))
            prod, jr = sp.Integer(1), _fold(N, j)
            for par in (+1, -1):
                S = block(H, N, par)
                SL, SR = S[:jr, :jr], S[jr + 1:, jr + 1:]
                prod *= (sp.expand(sp.resultant(charpoly_expr(SL).as_expr(),
                                                charpoly_expr(SR).as_expr(), x))
                         if SL.rows and SR.rows else sp.Integer(1))
            q = sp.cancel(Pfull / prod) if prod != 0 else sp.nan
            if q not in (sp.Integer(1), sp.Integer(-1)):
                bad.append((N, j, str(q)))
            else:
                signs[int(q)] = signs.get(int(q), 0) + 1
    check("K1b on the CHAIN, F157's P_j is +-1 times the two SECTOR halves' resultants",
          not bad and seen == K_EXPECTED.get('K1b') and signs == K1B_SIGNS,
          f"the open chain only, where striking DISCONNECTS: {seen} seats over N = 4..10, "
          f"expected {K_EXPECTED.get('K1b')}; the sign splits {signs}, expected {K1B_SIGNS}; the "
          f"split is this READING's, sympy's order and the seat order, and the law is W3; "
          f"failures {bad}")

    # -- K1b2: and it is a chain statement. On the ring the strike leaves ONE path, so the two
    # -- principal submatrices are not the pieces the seat cuts; the discriminator the check
    # -- actually reads is that on the ring the knob sits only in the wrap entry, so BOTH
    # -- submatrices are knob-free while the sector halves' product is not, and the ratio stops
    # -- being a constant. It fails at 20 of the 32 seats and survives at the other 12, which is
    # -- why the count is pinned rather than the property asserted universally.
    ring_fail, ring_seen = 0, 0
    for N in range(4, 11):
        H = H_crack(N)
        for j in range(N):
            if _is_centre(N, j) or j in (0, N - 1):
                continue
            ring_seen += 1
            Pfull = sp.expand(sp.resultant(charpoly_expr(H[:j, :j]).as_expr(),
                                           charpoly_expr(H[j + 1:, j + 1:]).as_expr(), x))
            prod, jr = sp.Integer(1), _fold(N, j)
            for par in (+1, -1):
                S = block(H, N, par)
                SL, SR = S[:jr, :jr], S[jr + 1:, jr + 1:]
                prod *= (sp.expand(sp.resultant(charpoly_expr(SL).as_expr(),
                                                charpoly_expr(SR).as_expr(), x))
                         if SL.rows and SR.rows else sp.Integer(1))
            q = sp.cancel(Pfull / prod) if prod != 0 else sp.nan
            if not (q.is_number and q != 0):
                ring_fail += 1
    check("K1b2 control: on the RING the identity fails at most seats, the knob sitting in the wrap alone",
          ring_fail == K_EXPECTED.get('K1b2') and ring_seen == K_EXPECTED.get('K1b'),
          f"{ring_fail} of {ring_seen} ring seats fail through the same door, expected "
          f"{K_EXPECTED.get('K1b2')} of {K_EXPECTED.get('K1b')}; it SURVIVES at the other "
          f"{ring_seen - ring_fail}, so this pins a count and not a universal")

    # -- K1b3: what stands where K1b stood, on the ring. Striking cuts nothing there, so there is
    # -- no P_j; the object that carries the u-locus is Res(chi_H, chi_H_j) itself, and on the
    # -- ring it is, exactly and at every interior non-centre seat, a seat-independent sign times
    # -- the sector halves' product SQUARED times the two ring ends, each raised to twice the
    # -- ring-end blind count of the count theorem. The cross factor Res(chi_E, chi_O) is those
    # -- ring ends and nothing else: no other root, real or complex. Second half: the twelve seats
    # -- where K1b2's identity stands are the seats where Q_E*Q_O is knob-free, a divisibility: the
    # -- distance from the seat to its mirror seat divides the seat's distance to the near end.
    # -- Two mutations through the same door: the ring-end exponents swapped (invisible at odd N,
    # -- where they agree, so it must redden at every even-N seat and only there), and the
    # -- modulus N - 2*jr in place of N - 1 - 2*jr. The end seats j = 0, N-1 are gated too, with
    # -- Q_E*Q_O = 1 (the seat strikes coordinate 0 of each sector, no inner half).
    bad, seen, standing, mut_law_red, mut_law_seen, mut_div_red, ends = [], 0, [], 0, 0, 0, 0
    for N in range(4, 11):
        H = H_crack(N)
        E, O = block(H, N, +1), block(H, N, -1)
        a, b = (N - 1) // 2, N // 2                     # the ring-end counts, u = +1 and u = -1
        cross2 = sp.expand(sp.resultant(charpoly_expr(E).as_expr(), charpoly_expr(O).as_expr(), x) ** 2)
        law_cross2 = sp.expand(4 ** (N // 2) * (t - 1) ** (2 * a) * (t + 1) ** (2 * b))
        if sp.expand(cross2 - law_cross2) != 0:
            bad.append(('cross', N))
        for j in (0, N - 1):
            ends += 1
            lhs_end = sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                             charpoly_expr(strike(H, j)).as_expr(), x))
            rhs_end = sp.expand((-1) ** (N // 2) * (t - 1) ** (2 * a) * (t + 1) ** (2 * b))
            if sp.expand(lhs_end - rhs_end) != 0:
                bad.append(('end seat', N, j))
        for j in range(1, N - 1):
            if _is_centre(N, j):
                continue
            seen += 1
            jr = _fold(N, j)
            prod = sp.Integer(1)
            for S in (E, O):
                SL, SR = S[:jr, :jr], S[jr + 1:, jr + 1:]
                prod *= sp.expand(sp.resultant(charpoly_expr(SL).as_expr(),
                                               charpoly_expr(SR).as_expr(), x))
            lhs = sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                         charpoly_expr(strike(H, j)).as_expr(), x))
            rhs = sp.expand((-1) ** (N // 2) * prod ** 2 * (t - 1) ** (2 * a) * (t + 1) ** (2 * b))
            if sp.expand(lhs - rhs) != 0:
                bad.append(('ring-law', N, j))
            mut = sp.expand((-1) ** (N // 2) * prod ** 2 * (t - 1) ** (2 * b) * (t + 1) ** (2 * a))
            if N % 2 == 0:
                mut_law_seen += 1
                mut_law_red += int(sp.expand(lhs - mut) != 0)
            elif sp.expand(lhs - mut) != 0:
                bad.append(('mutation visible at odd N', N, j))
            if prod == 0:
                bad.append(('Q_E Q_O vanishes, u-free reading undefined', N, j))
            knobfree = sp.Poly(prod, t).degree() == 0
            if knobfree:
                standing.append((N, j))
            if knobfree != (jr % (N - 1 - 2 * jr) == 0):
                bad.append(('divisibility law', N, j))
            mut_div_red += int((jr % (N - 2 * jr) == 0) != knobfree)
    check("K1b3 on the RING, Res(chi_H, chi_H_j) = (-1)^floor(N/2) (Q_E Q_O)^2 (u-1)^2a (u+1)^2b, and the cross factor is the ring ends alone",
          not bad and seen == K_EXPECTED.get('K1b3') and ends == 14 and standing == K1B3_STANDING
          and mut_law_red == mut_law_seen and mut_law_seen > 0 and mut_div_red > 0,
          f"exact in u at {seen} seats over N = 4..10, expected {K_EXPECTED.get('K1b3')}, and at "
          f"{ends} end seats with Q_E Q_O = 1, expected 14; the "
          f"seats where K1b2's identity stands are exactly those with (N-1-2jr) | jr: {standing}; exponents "
          f"swapped reddens {mut_law_red} of {mut_law_seen} even-N seats, modulus N-2jr reddens "
          f"{mut_div_red}; failures {bad}")

    # -- K1c: and inside a sector the seat DISCONNECTS the block, so the sector resultant is a
    # -- t-free constant times the sector halves' resultant SQUARED. With K1b that is where the
    # -- squaring of P_j comes from, and K1e is the consequence.
    bad, seen = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 11):
            H = Hf(N)
            for par in (+1, -1):
                S = block(H, N, par)
                for j in range(N):
                    if _is_centre(N, j):
                        continue
                    jr = _fold(N, j)
                    R = sp.expand(sp.resultant(charpoly_expr(S).as_expr(),
                                               charpoly_expr(strike(S, jr)).as_expr(), x))
                    if R == 0:
                        continue
                    seen += 1
                    L, Rt = S[:jr, :jr], S[jr + 1:, jr + 1:]
                    P = sp.Integer(1) if L.rows == 0 or Rt.rows == 0 else sp.expand(
                        sp.resultant(charpoly_expr(L).as_expr(), charpoly_expr(Rt).as_expr(), x))
                    q = sp.cancel(R / P ** 2)
                    if not (q.is_number and q != 0):
                        bad.append((name, N, par, j, str(q)))
    check("K1c the sector resultant is a t-free constant times the sector halves' resultant SQUARED",
          not bad and seen == K_EXPECTED.get('K1c'),
          f"{seen} (axis, N, sector, seat) readings over N = 4..10, expected {K_EXPECTED.get('K1c')}; "
          f"failures {bad}")

    # -- K1d: the square is the DISCONNECTION, so it must go when the block stops being
    # -- tridiagonal. One next-nearest entry through the same door, and the ratio stops being a
    # -- constant.
    survived, tried = 0, 0
    for N in range(5, 11):
        H = sp.Matrix(H_aniso(N))
        H[0, 2] = H[0, 2] + 1
        H[2, 0] = H[2, 0] + 1
        H[N - 1, N - 3] = H[N - 1, N - 3] + 1
        H[N - 3, N - 1] = H[N - 3, N - 1] + 1
        S = block(H, N, +1)
        d = S.rows
        if d < 4:
            continue
        for jr in range(1, d - 1):
            R = sp.expand(sp.resultant(charpoly_expr(S).as_expr(),
                                       charpoly_expr(strike(S, jr)).as_expr(), x))
            L, Rt = S[:jr, :jr], S[jr + 1:, jr + 1:]
            if R == 0 or L.rows == 0 or Rt.rows == 0:
                continue
            P = sp.expand(sp.resultant(charpoly_expr(L).as_expr(), charpoly_expr(Rt).as_expr(), x))
            q = sp.cancel(R / P ** 2) if P != 0 else sp.nan
            tried += 1
            if q.is_number and q != 0:
                survived += 1
    check("K1d control: one next-nearest entry, and the square is gone",
          survived == K1D_SQUARE_SURVIVORS and tried == 10,
          f"{survived} of {tried} fold seats over N = 5..10 still give a constant ratio, "
          f"expected {K1D_SQUARE_SURVIVORS} of 10; so K1c reads the disconnection, not the "
          f"construction")

    # -- K1e: the headline four repo passages state. Off the centre seat the FULL resultant is a
    # -- t-free rational constant times F157's P_j SQUARED. Like K1 it is a CONSEQUENCE of what
    # -- sits above it, K1 with K1b and K1c and the t-freeness of the cross factor, on exactly
    # -- this population; it earns its place by computing the composed statement directly, from
    # -- the N x N matrix rather than from the pieces, so a slip in the composition would show.
    bad, seen, cross_const = [], 0, 0
    for N in range(3, 15):
        H = H_aniso(N)
        E, O = block(H, N, +1), block(H, N, -1)
        cx = sp.expand(sp.resultant(charpoly_expr(E).as_expr(), charpoly_expr(O).as_expr(), x))
        if cx.is_number and abs(cx) == 2 ** (N // 2):
            cross_const += 1
        else:
            bad.append((N, "cross factor not +-2^floor(N/2)", str(cx)))
    for N in range(4, 11):
        H = H_aniso(N)
        for j in range(N):
            if _is_centre(N, j) or j in (0, N - 1):
                continue
            full = sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                          charpoly_expr(strike(H, j)).as_expr(), x))
            Pj = sp.expand(sp.resultant(charpoly_expr(H[:j, :j]).as_expr(),
                                        charpoly_expr(H[j + 1:, j + 1:]).as_expr(), x))
            if Pj == 0:
                continue
            seen += 1
            q = sp.cancel(full / Pj ** 2)
            if not (q.is_number and q != 0):
                bad.append((N, j, str(q)))
    check("K1e on the CHAIN, Res(chi_H, chi_H_j) is a t-free constant times F157's P_j SQUARED",
          not bad and seen == K_EXPECTED.get('K1e') and cross_const == K1E_CROSS_CONSTANT_N,
          f"the claim four repo passages stated and none derived before the companion proof's "
          f"section (g); {seen} readings over N = 4..10, expected {K_EXPECTED.get('K1e')}; the "
          f"cross-sector factor is t-free at {cross_const} of {K1E_CROSS_CONSTANT_N} N here, and "
          f"that proof derives it for every N from the geometric multiplicity of a tridiagonal "
          f"block; failures {bad}")

    # -- K1e2: K1e reads the DISCONNECTION too, not the chain's shape: moving the diagonal leaves
    # -- it green, and only losing tridiagonality breaks it. One next-nearest entry, same door.
    survived, tried = 0, 0
    for N in range(5, 11):
        H = sp.Matrix(H_aniso(N))
        H[0, 2] = H[0, 2] + 1
        H[2, 0] = H[2, 0] + 1
        for j in range(1, N - 1):
            full = sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                          charpoly_expr(strike(H, j)).as_expr(), x))
            Pj = sp.expand(sp.resultant(charpoly_expr(H[:j, :j]).as_expr(),
                                        charpoly_expr(H[j + 1:, j + 1:]).as_expr(), x))
            if full == 0 or Pj == 0:
                continue
            tried += 1
            q = sp.cancel(full / Pj ** 2)
            if q.is_number and q != 0:
                survived += 1
    check("K1e2 control: one next-nearest entry on the chain, and K1e's constant is gone",
          survived == K1E2_SURVIVORS and tried == 32,
          f"{survived} of {tried} seats over N = 5..10 still give a t-free ratio, expected "
          f"{K1E2_SURVIVORS} of 32; so K1e reads the disconnection. That it does NOT read the "
          f"diagonal's POSITION is a separate statement and no check here makes it")

    # -- K2a: the field arithmetic against two independent oracles. On a LINEAR mu the committed
    # -- rational blind_count; on a QUADRATIC one sympy's own gcd over the number field, which is
    # -- the only check here that exercises the polynomial branch of the modular inverse.
    bad = []
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 12):
            for j in range(N):
                for val in (sp.Integer(0), sp.Rational(1, 3), sp.Integer(1)):
                    if blind_at(Hf, N, j, t - val) != blind_count(Hf(N, val), j):
                        bad.append((name, N, j, str(val)))
    # The DEGREE-4 mu is not decoration: the modular inverse's polynomial branch is reached only
    # from degree 3 up. Counted over THIS check's own loop: 0 calls at t^2-3, 0 at t^2-1/2, and
    # 258 at t^4-4t^2+2. So a quadratic oracle alone would leave that branch unread here.
    nf, poly_branch = [], [0]
    _orig_inv = _inv_mod_mu

    def _counting_inv(c, mu, _f=_orig_inv, _n=poly_branch):
        if not c.is_number:
            _n[0] += 1
        return _f(c, mu)

    globals()['_inv_mod_mu'] = _counting_inv
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in (7, 9, 11):
            for alpha, mu in ((sp.sqrt(3), t ** 2 - 3),
                              (sp.sqrt(2) / 2, t ** 2 - sp.Rational(1, 2)),
                              (sp.sqrt(2 + sp.sqrt(2)), t ** 4 - 4 * t ** 2 + 2)):
                K = sp.QQ.algebraic_field(alpha)
                for j in range(N):
                    H = Hf(N, alpha)
                    a = sp.Poly(H.charpoly(x).as_expr(), x, domain=K)
                    b = sp.Poly(strike(H, j).charpoly(x).as_expr(), x, domain=K)
                    if blind_at(Hf, N, j, mu) != sp.gcd(a, b).degree():
                        nf.append((name, N, j, str(mu)))
    globals()['_inv_mod_mu'] = _orig_inv
    check("K2a the Q[t]/(mu) route against two independent oracles, rational and number-field",
          not bad and not nf and poly_branch[0] > 0,
          f"the committed rational blind_count at 3 knob values over N = 4..11, and sympy's gcd "
          f"over Q(alpha) at two quadratic and one QUARTIC mu over N = 7, 9, 11; the modular "
          f"inverse's polynomial branch was entered {poly_branch[0]} times, COUNTED rather than "
          f"inferred; the two quadratic mu contribute 0 of those calls, which is why the quartic "
          f"is in the loop; failures {bad + nf}")

    # -- K2: the count decomposition. Off the degeneracy set, blind(j) = b_E + b_O exactly.
    bad, seen = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 13):
            for j in range(N):
                for mu in locus_mus(N, j):
                    bE, bO, cross = sector_counts_at(Hf, N, j, mu)
                    if cross != 0:
                        continue
                    seen += 1
                    if blind_at(Hf, N, j, mu) != bE + bO:
                        bad.append((name, N, j, str(mu)))
    check("K2  blind(j) = b_E + b_O at every locus point off the degeneracy set",
          not bad and seen == K_EXPECTED.get('K2'),
          f"{seen} (axis, N, seat, point) readings over N = 4..12, expected {K_EXPECTED.get('K2')}; "
          f"failures {bad}")

    # -- K2b: the correction term, where the additivity is replaced. blind = b_E + b_O + the
    # -- shared levels NEITHER sector eigenvector misses at the seat, at both ring ends.
    wrong, bitten = [], {}
    for N in range(5, 13):
        n_bitten = 0
        for j in range(N):
            for end in (t - 1, t + 1):
                bE, bO, cross = sector_counts_at(H_crack, N, j, end)
                neither = shared_neither_at(H_crack, N, j, end)
                if blind_at(H_crack, N, j, end) != bE + bO + neither:
                    wrong.append((N, j, str(end)))
                if neither != cross:
                    n_bitten += 1
        bitten[N] = n_bitten
    check("K2b the correction term at u = +-1: blind = b_E + b_O + #{shared lambda neither misses}",
          not wrong and bitten == K2B_CORRECTED_SEATS,
          f"every seat over N = 5..12, both ring ends; the correction is strictly below the "
          f"shared count at {bitten} seat-ends per N, expected {K2B_CORRECTED_SEATS}; "
          f"arithmetic failures {wrong}")

    # -- K2c: the third inclusion-exclusion term is ZERO at every reading, on both families and
    # -- at both ring ends, and it is zero for a reason rather than by luck: no site is a node of
    # -- both the cosine and the sine of one mode. Without this the term could carry a sign error
    # -- forever, since flipping it changes nothing the other checks read.
    nonzero, shared_seen, seen = [], 0, 0
    for Hf in (H_aniso, H_crack):
        for N in range(5, 13):
            for j in range(N):
                for end in (t - 1, t + 1, t, t ** 2 - 3):
                    shared, a, b, c = shared_neither_at(Hf, N, j, end, parts=True)
                    seen += 1
                    shared_seen += shared
                    if c != 0:
                        nonzero.append((Hf.__name__, N, j, str(end), c))
    check("K2c the both-vanish term of the correction is 0 at every reading, and not vacuously",
          not nonzero and seen == K2C_READINGS and shared_seen == K2C_SHARED,
          f"{seen} readings over both families and N = 5..12, expected {K2C_READINGS}, carrying "
          f"{shared_seen} shared levels in total, expected {K2C_SHARED}, so the term had somewhere "
          f"to be nonzero; readings with a nonzero term {nonzero}")

    # -- K3: the count theorem at ODD N. Same points AND same count, off the ring ends.
    bad, seen, at_origin = [], 0, 0
    for N in (5, 7, 9, 11, 13):
        for j in range(N):
            for mu in locus_mus(N, j):
                seen += 1
                if mu == t:
                    at_origin += 1
                if blind_at(H_crack, N, j, mu) != blind_at(H_aniso, N, j, mu):
                    bad.append((N, j, str(mu)))
    check("K3  ODD N: blind_u(j, t) = blind_Delta(j, t) at every locus point off the ring ends",
          not bad and seen == K_EXPECTED.get('K3') and at_origin == K3_AT_ORIGIN,
          f"{seen} (N, seat, point) readings over N = 5..13, expected {K_EXPECTED.get('K3')}; "
          f"{at_origin} of them sit at mu = t, where the two families are the SAME matrix and the "
          f"comparison is a number against itself, so the informative population is "
          f"{seen - at_origin}; failures {bad}")

    # -- K3b: and at the ring ends they DIFFER, against literals derived from F157's committed
    # -- P_j and from Lemma 4. Its ring column repeats what K5 asserts over a wider range; its
    # -- content is the four CHAIN numbers.
    got = {}
    for (N, j) in ODD_N_RING_END_DIFFERENCES:
        got[(N, j)] = (blind_at(H_aniso, N, j, t - 1), blind_at(H_crack, N, j, t - 1))
    check("K3b at u = +1 the two axes pay DIFFERENT counts, at the four seats named as literals",
          got == ODD_N_RING_END_DIFFERENCES,
          f"(chain, ring) = {got}, expected {ODD_N_RING_END_DIFFERENCES}")

    # -- K3c: K3 is an equality between two readings of one routine, so a defect that moves BOTH
    # -- axes alike survives it, a uniform miscount included. This is the control against the
    # -- other failure: that the equality is a property of any two families rather than of these.
    disagree, same, same_off_origin = 0, 0, 0
    for N in (5, 7, 9, 11, 13):
        for j in range(N):
            for mu in locus_mus(N, j):
                if blind_at(H_oneend, N, j, mu) != blind_at(H_crack, N, j, mu):
                    disagree += 1
                else:
                    same += 1
                    if mu != t:
                        same_off_origin += 1
    check("K3c control: a THIRD family tracks the crack only at the origin, where all three are one",
          disagree == K3C_DISAGREEMENTS and same + disagree == K_EXPECTED.get('K3')
          and same_off_origin == 0,
          f"{disagree} of {disagree + same} triples disagree, expected {K3C_DISAGREEMENTS}, and "
          f"the {same} that do not are ALL at mu = t, where all three families are one matrix "
          f"({same_off_origin} agreements off the origin); so the separation off the origin is "
          f"total and K3 reads these two families and not any two")

    # -- K3d: Lemma 5-prime itself, which K3 uses and never tests. At ODD N the staggering bijection
    # -- preserves each sector, so each sector count is negation-invariant; at EVEN N it swaps
    # -- them, so the invariance must FAIL and the two counts cross over instead. Both halves
    # -- through the same door, so neither is a code path.
    bad, odd_seen, even_seen = [], 0, 0
    for N in (5, 7, 9, 11, 13, 6, 8, 10, 12):
        for j in range(N):
            for mu in locus_mus(N, j):
                bE, bO, _ = sector_counts_at(H_aniso, N, j, mu)
                bEm, bOm, _ = sector_counts_at(H_aniso, N, j, mu, sign=-1)
                if N % 2:
                    odd_seen += 1
                    if (bE, bO) != (bEm, bOm):
                        bad.append(("odd sector not negation-invariant", N, j, str(mu)))
                else:
                    even_seen += 1
                    if (bE, bO) != (bOm, bEm):
                        bad.append(("even sectors do not cross over", N, j, str(mu)))
    check("K3d Lemma 5-prime on the COUNTS: odd N preserves each sector, even N swaps them",
          not bad and odd_seen == K3D_ODD and even_seen == K3D_EVEN,
          f"{odd_seen} odd-N and {even_seen} even-N readings, expected {K3D_ODD} and {K3D_EVEN}; "
          f"failures {bad}")

    # -- K4: EVEN N, off the ring ends. The ring pays twice the even sector; the chain pays the
    # -- even sector at t and at -t.
    bad, seen = [], 0
    for N in (6, 8, 10, 12):
        for j in range(N):
            for mu in locus_mus(N, j):
                bE, bO, cross = sector_counts_at(H_aniso, N, j, mu)
                if cross != 0:
                    continue
                seen += 1
                bEm = sector_counts_at(H_aniso, N, j, mu, sign=-1)[0]
                if blind_at(H_crack, N, j, mu) != 2 * bE:
                    bad.append(("ring != 2 b_E", N, j, str(mu)))
                if blind_at(H_aniso, N, j, mu) != bE + bEm:
                    bad.append(("chain != b_E(t) + b_E(-t)", N, j, str(mu)))
                # the crossover b_O(t) = b_E(-t) is K3d's even branch, on exactly this
                # population; it is not repeated here
    check("K4  EVEN N off the ring ends: blind_u = 2 b_E, blind_Delta = b_E(t) + b_E(-t)",
          not bad and seen == K_EXPECTED.get('K4'),
          f"{seen} readings over N = 6..12, expected {K_EXPECTED.get('K4')}; failures {bad}")

    # -- K4b: the even-N counts genuinely PART, at four (N, seat, point) triples named with both
    # -- values. A locus that moves must redden this rather than crash it.
    got, missing = {}, []
    for (N, j, sname) in EVEN_N_DIFFERENCES:
        mu = next((m for m in locus_mus(N, j) if str(m) == sname), None)
        if mu is None:
            missing.append((N, j, sname))
            continue
        got[(N, j, sname)] = (blind_at(H_aniso, N, j, mu), blind_at(H_crack, N, j, mu))
    check("K4b EVEN N: the counts differ at the four literal triples, two of them the ring paying 0",
          got == EVEN_N_DIFFERENCES and not missing,
          f"(chain, ring) = {got}, expected {EVEN_N_DIFFERENCES}; points no longer in the locus: "
          f"{missing}")

    # -- K4c: THE FENCE on K4, and the reason the even-N law needs one. At the ring ends the ring
    # -- pays Lemma 4's count at every seat while 2*b_E is something else entirely, so K4's
    # -- exclusion of {+-1} is load-bearing rather than an artefact of locus_mus.
    agree, seen = [], 0
    for N in (6, 8, 10, 12):
        for j in range(N):
            for end in (t - 1, t + 1):
                seen += 1
                b = blind_at(H_crack, N, j, end)
                two = 2 * sector_counts_at(H_aniso, N, j, end)[0]
                if b == two:
                    agree.append((N, j, str(end)))
    check("K4c the |t| != 1 fence on K4 is load-bearing: at the ring ends blind_u is mostly NOT 2 b_E",
          sorted(agree) == sorted(K4C_AGREEING_SEATS) and seen == 72,
          f"the two readings coincide at {len(agree)} of {seen} even-N seat-ends, expected "
          f"{len(K4C_AGREEING_SEATS)} of 72: {agree}")

    # -- K4d: WHERE the even-N counts part. If a seat's two loci differ, some point of the
    # -- symmetric difference has b_E(t) >= 1 and b_E(-t) = 0, so the counts differ there too:
    # -- that direction is a theorem. The converse is NOT, and this reads it: no seat carries
    # -- the same locus and a different count anywhere in the range.
    parting, same_locus_diff = set(), []
    for N in (6, 8, 10, 12, 14):
        for j in range(N):
            lu, ld = locus_full(H_crack, N, j), locus_full(H_aniso, N, j)
            if lu == 'ALL' or ld == 'ALL':
                continue
            same_locus = (lu == (ld | RING))
            for mu in locus_mus(N, j):
                if blind_at(H_aniso, N, j, mu) != blind_at(H_crack, N, j, mu):
                    parting.add((N, j))
                    if same_locus:
                        same_locus_diff.append((N, j, str(mu)))
    breaks = {(N, j) for (N, j) in COMMITTED_BREAKS if N <= 14}
    check("K4d EVEN N: the counts part at exactly the seats where the loci break, and nowhere else",
          parting == breaks and not same_locus_diff,
          f"{len(parting)} parting seats over N = 6..14 against the {len(breaks)} committed break "
          f"seats in that range; seats with the SAME locus and a different count: "
          f"{same_locus_diff}")

    # -- K5b: at the ring ends, WHERE the two axes still agree. Exactly the reflection-fixed
    # -- centre seat of an odd chain, and no seat at all at even N.
    wrong, seen = [], 0
    for N in range(5, 12):
        for end in (t - 1, t + 1):
            seen += 1
            agree = [j for j in range(N)
                     if blind_at(H_aniso, N, j, end) == blind_at(H_crack, N, j, end)]
            want = [(N - 1) // 2] if N % 2 else []
            if agree != want:
                wrong.append((N, str(end), agree, want))
    check("K5b at the ring ends the axes agree at the odd-N centre seat ALONE, at even N nowhere",
          not wrong and seen == 14,
          f"{seen} (N, end) rows over N = 5..11, expected 14; the centre half follows from K5 and "
          f"K7, the 'and nowhere else' half is this check's own content; disagreements {wrong}")

    # -- K5c: WHICH spectrum the two ring ends carry, which the third bullet of the count theorem
    # -- reads. Decided by an exact integer identity rather than by divisibility or a simplifier:
    # -- the wrap bond +1 is a periodic boundary condition and -1 an ANTIperiodic one, so
    # -- chi(H_u(+-1)) = 2*T_N(x/2) -+ 2, whose roots are 2cos(2*pi*m/N) and 2cos((2m+1)*pi/N).
    bad, seen = [], 0
    for N in range(3, 15):
        T = sp.expand(sp.chebyshevt(N, x / 2))
        for uval, offset in ((sp.Integer(1), -2), (sp.Integer(-1), 2)):
            seen += 1
            ch = sp.expand(H_crack(N, uval).charpoly(x).as_expr())
            if sp.expand(ch - sp.expand(2 * T + offset)) != 0:
                bad.append((N, int(uval)))
    check("K5c the ring ends are the periodic and the ANTIperiodic comb: chi = 2*T_N(x/2) -+ 2",
          not bad and seen == 24,
          f"exact integer polynomial identity, {seen} readings over N = 3..14 at both ends, "
          f"expected 24; failures {bad}")

    # -- K5: the ring ends pay the number of degenerate pairs, at EVERY seat, from Lemma 4's
    # -- combinatorial counts rather than from the measurement.
    bad, seen = [], 0
    for N in range(5, 13):
        for j in range(N):
            seen += 2
            if blind_at(H_crack, N, j, t - 1) != (N - 1) // 2:
                bad.append(("u = +1", N, j))
            if blind_at(H_crack, N, j, t + 1) != -(-(N - 1) // 2):
                bad.append(("u = -1", N, j))
    check("K5  at u = +-1 EVERY seat pays floor/ceil((N-1)/2), Lemma 4's degenerate-pair count",
          not bad and seen == 136,
          f"{seen} seat-end readings over N = 5..12, expected 136; failures {bad}")

    # -- K6: the multiplicity IS in the resultant, doubled. Off C the order of vanishing of the
    # -- full resultant at mu is exactly twice the blind count.
    bad, seen = [], 0
    for name, Hf in (("aniso", H_aniso), ("crack", H_crack)):
        for N in range(4, 12):
            for j in range(N):
                if _is_centre(N, j):
                    continue
                H = Hf(N)
                full = sp.Poly(sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                                      charpoly_expr(strike(H, j)).as_expr(), x)), t)
                if full.is_zero:
                    continue
                mult = {sp.factor(sp.monic(sp.Poly(f, t)).as_expr()): m
                        for f, m in sp.factor_list(full.as_expr(), t)[1]}
                for mu in locus_mus(N, j):
                    if sector_counts_at(Hf, N, j, mu)[2] != 0:
                        continue
                    seen += 1
                    if mult.get(mu, 0) != 2 * blind_at(Hf, N, j, mu):
                        bad.append((name, N, j, str(mu)))
    check("K6  off C the order of vanishing of Res(chi_H, chi_H_j) at mu is exactly 2 * blind(j)",
          not bad and seen == K_EXPECTED.get('K6'),
          f"{seen} readings over N = 4..11, both axes, expected {K_EXPECTED.get('K6')}; failures {bad}")

    # -- K6b: and AT the ring ends halving is mostly still exact, which is the opposite of what a
    # -- reader would guess from K6's fence. It fails at exactly the seat-ends where some shared
    # -- level has ONE sector eigenvector vanishing at the seat; the cross-sector factor vanishes
    # -- at every ring-end seat and is therefore NOT what breaks it.
    exact, broken = 0, []
    for N in range(5, 11):
        H = H_crack(N)
        for j in range(N):
            if _is_centre(N, j):
                continue
            full = sp.Poly(sp.expand(sp.resultant(charpoly_expr(H).as_expr(),
                                                  charpoly_expr(strike(H, j)).as_expr(), x)), t)
            if full.is_zero:
                continue
            mult = {sp.factor(sp.monic(sp.Poly(f, t)).as_expr()): m
                    for f, m in sp.factor_list(full.as_expr(), t)[1]}
            for end in (t - 1, t + 1):
                if mult.get(sp.factor(end), 0) == 2 * blind_at(H_crack, N, j, end):
                    exact += 1
                else:
                    bE, bO, _ = sector_counts_at(H_crack, N, j, end)
                    broken.append((N, j, str(end), bE + bO))
    check("K6b at the ring ends halving stays exact except where a sector eigenvector vanishes",
          exact == K6B_EXACT and [b[:3] for b in broken] == K6B_BROKEN
          and all(b[3] > 0 for b in broken),
          f"exact at {exact} of {exact + len(broken)} ring-end seat-ends, expected {K6B_EXACT}; "
          f"the exceptions are {[b[:3] for b in broken]} and every one of them has b_E + b_O > 0")

    # -- K7: the centre seat. b_O is the whole odd sector by reflection and b_E is 0 because the
    # -- centre is an END of the folded even block, so the count is (N-1)/2 at EVERY knob value,
    # -- the ring ends included.
    bad, seen = [], 0
    for N in (5, 7, 9, 11, 13):
        c = (N - 1) // 2
        for Hf in (H_aniso, H_crack):
            for mu in (t, t - 1, t + 1, t - sp.Rational(7, 5), t ** 2 - 3, t ** 3 - 3 * t - 1):
                seen += 1
                if blind_at(Hf, N, c, mu) != (N - 1) // 2:
                    bad.append((Hf.__name__, N, str(mu)))
    check("K7  the centre seat pays (N-1)/2 at EVERY knob value on BOTH axes, rational or not",
          not bad and seen == 60,
          f"{seen} readings, 6 knob values incl. two irrational and both ring ends, N = 5..13, "
          f"both axes, expected 60; failures {bad}")


# ------------------------------------------------- (S) which k sits in which reflection sector
#
# F157 indexes the Delta-locus by k in 1..N_node-1 through the node angle theta_k = k*pi/N_node.
# Lemma 8 reads the mode at that angle, psi_l = sin((j-l)*theta_k), and its reflection parity,
# (-1)^(k+1). Everything below is a polynomial identity in v = 2*cos(theta_k) reduced modulo an
# integer minimal polynomial: no trigonometric value is ever handed to a simplifier and no root
# object is built, which is block K's discipline applied to the mode instead of to the count.
#
# The field is L = Q[z]/(Phi_n) with z = 2*cos(pi/n), n = N_node. Every 2*cos(k*pi/n) lives in it
# as 2*T_k(z/2), so two node indices can be compared for an EQUAL Delta exactly, which grouping
# by a simplified trigonometric value could not do.
#
# The sign law itself is NOT new here and the file does not claim it: R*psi_k = (-1)^(k+1)*psi_k
# is owned for the uniform chain's own comb in PROOF_COLLISION_GAP_ODD_ORDERS (b), in
# experiments/SLOW_MODE_R_PARITY.md, and typed as MirrorWorld's Formulas.F71_ReflectionParity.
# What is read here is that it TRANSPORTS to the node comb, whose modulus N_node = |N-1-2j| is at
# most N-1 and so never the chain's N+1 (arithmetic, not a gate: an earlier build asserted it as
# one and no input could have reddened it), and what that transport buys is the k-index on b_E and
# b_O. The matrix is A + Delta_k*D and NOT the uniform chain the six owners read, which is why the
# single sine is a claim and not a restatement.

z = sp.symbols('z')


def n_node(N, j):
    return abs(N - 1 - 2 * j)


def S_cheb(n, v):
    """sin(n*theta)/sin(theta) as an integer polynomial in v = 2*cos(theta); S_{-n} = -S_n."""
    if n < 0:
        return -S_cheb(-n, v)
    return sp.expand(sp.chebyshevu(n - 1, v / 2))


def T_cheb(n, v):
    """cos(n*theta) as an integer polynomial in v = 2*cos(theta); T_{-n} = T_n."""
    return sp.expand(sp.chebyshevt(abs(n), v / 2))


_PHI_CACHE = {}


def phi_2cos(n):
    """The minimal polynomial of 2*cos(pi/n) over Q, as a Poly in z. Verified by S0."""
    if n not in _PHI_CACHE:
        _PHI_CACHE[n] = sp.Poly(sp.minimal_polynomial(2 * sp.cos(sp.pi / n), z), z, domain='QQ')
    return _PHI_CACHE[n]


def _pz(e):
    return sp.Poly(e, z, domain='QQ')


def _redz(e, phi):
    """e mod phi, as a Poly in z over Q. Accepts an expression or a Poly."""
    return (e if isinstance(e, sp.Poly) else _pz(e)).rem(phi)


def _invz(p, phi):
    return sp.invert(p, phi)


_XK_CACHE, _SC_CACHE = {}, {}


def _x_k(n, k):
    """2*cos(k*pi/n) as an element of L, that is 2*T_k(z/2) reduced."""
    if (n, k) not in _XK_CACHE:
        _XK_CACHE[(n, k)] = _redz(2 * T_cheb(k, z), phi_2cos(n))
    return _XK_CACHE[(n, k)]


def _S_at(m, n, k):
    """S_m evaluated at x_k, in L, by the recursion S_{m+1} = x*S_m - S_{m-1}."""
    key = (m, n, k)
    if key in _SC_CACHE:
        return _SC_CACHE[key]
    if m < 0:
        out = -_S_at(-m, n, k)
    elif m == 0:
        out = _pz(0)
    elif m == 1:
        out = _pz(1)
    else:
        out = _redz(_x_k(n, k) * _S_at(m - 1, n, k) - _S_at(m - 2, n, k), phi_2cos(n))
    _SC_CACHE[key] = out
    return out


def _T_at(m, n, k):
    """T_m evaluated at x_k, in L. T_0 = 1, T_1 = x/2, T_{m+1} = x*T_m - T_{m-1}."""
    m = abs(m)
    a, b = _pz(1), _redz(_x_k(n, k) * sp.Rational(1, 2), phi_2cos(n))
    if m == 0:
        return a
    for _ in range(m - 1):
        a, b = b, _redz(_x_k(n, k) * b - a, phi_2cos(n))
    return b


def mode_data(N, j, k):
    """(phi, x_k, Delta_k, psi) in L for seat j at node index k, or None at one of F157's poles."""
    n = n_node(N, j)
    phi = phi_2cos(n)
    Sj = _S_at(j, n, k)
    if Sj.is_zero:
        return None
    dk = _redz(_S_at(j + 1, n, k) * _invz(Sj, phi), phi)
    psi = [_S_at(j - l, n, k) for l in range(N)]
    return phi, _x_k(n, k), dk, psi


def _coeff_vec(p, d):
    c = p.all_coeffs()[::-1]
    return [c[i] if i < len(c) else sp.Integer(0) for i in range(d)]


def min_poly_of(el, phi, var):
    """Minimal polynomial over Q of an element of L, by exact linear algebra on its powers.

    Not an algebraic-number oracle: S3 checks that the answer annihilates the element and is
    irreducible, so a wrong dependency cannot pass for a right one.
    """
    d = phi.degree()
    powers, cur = [], _pz(1)
    for r in range(d + 1):
        powers.append(_coeff_vec(cur, d))
        ns = sp.Matrix(powers).T.nullspace()
        if ns:
            w = ns[0] / ns[0][r]
            return sp.expand(sum(w[i] * var ** i for i in range(r + 1)))
        cur = _redz(cur * el, phi)
    raise RuntimeError("no dependency within the degree of the field")


def _sector_split_by_k(N, j):
    """{Delta as an element of L: [#odd k, #even k]} at seat j, poles dropped."""
    out = {}
    for k in range(1, n_node(N, j)):
        md = mode_data(N, j, k)
        if md is None:
            continue
        key = md[2].as_expr()
        out.setdefault(key, [0, 0])[0 if k % 2 else 1] += 1
    return out


# The sweep ranges. S3's and S4's are shorter because a sector count is a gcd over Q[t]/(mu).
S_N = range(5, 21)
S3_N = range(5, 15)
S4_N = [6, 8, 10, 12, 14]

# PINNED MEASUREMENTS, frozen from a run so that a shrinking sweep goes red instead of quiet.
# Every check in this block carries one: the ranges above are bare literals, so without a pinned
# population a narrowed sweep would pass in silence, which is the defect the 2026-09-02 entry of
# docs/CAUGHT_ERRORS.md logs for block K and which came back here.
S1_MODES, S1_POLES = 812, 412          # non-pole modes, and poles, over S_N
S2_ANGLES = 1224                       # (N, seat, k) angles over S_N, poles included
S2C_PAIRS = 812                        # non-pole k for which {k, N_node-k} is read; each
                                       # unordered pair is therefore read twice
S2D_TRIED, S2D_SURVIVORS = 851, 0      # wrong-modulus constructions tried, and how many keep
                                       # a parity anyway; the population is pinned too, or a
                                       # narrowed sweep would satisfy 'none survive' by
                                       # trying none
S3_READINGS = 196                      # (seat, Delta) readings over S3_N
S3_ASYMMETRIC = 188                    # of those, the ones with #odd k != #even k
S4_READINGS = 36                       # (N, seat, locus point) readings over S4_N


def _entry_times(e, p, dk, phi):
    """The matrix entry e, an element of Z[t], applied to a field element p at knob dk.

    This is what lets S1 read the FAMILY OBJECT H_aniso rather than a second transcription of
    it: an earlier build wrote the tridiagonal rows out by hand, and moving the anisotropy in
    H_aniso left S1 green. Entries of both families are affine in t, which the assert pins.
    """
    q = sp.Poly(e, t)
    assert q.degree() <= 1, "a family entry is not affine in the knob"
    c0, c1 = q.nth(0), q.nth(1)
    out = _pz(0)
    if c0:
        out = out + int(c0) * p
    if c1:
        out = out + int(c1) * _redz(dk * p, phi)
    return out


def gate_S():
    print("\n(S) THE SECTOR OF k: the node mode, its parity, and the k-index on the count")

    # -- S0: the field is verified, not trusted. Phi_n is irreducible, carries the degree the
    # -- theory of Q(2cos(pi/n)) forces, and divides the integer polynomial whose roots are
    # -- exactly the 2cos(k*pi/n); an irreducible divisor of U_{n-1}(z/2) of that degree can only
    # -- be the minimal polynomial of a primitive one, and those are all conjugate, so the three
    # -- together pin Phi_n. Every reduction below happens modulo it, and the range is DERIVED
    # -- from the sweep rather than written as a literal, so widening S_N cannot leave a modulus
    # -- unverified.
    n_max = max(n_node(N, j) for N in S_N for j in range(N))
    bad, seen = [], 0
    for n in range(2, n_max + 1):
        phi = phi_2cos(n)
        seen += 1
        deg_ok = phi.degree() == sp.totient(2 * n) // 2
        irr_ok = phi.is_irreducible
        div_ok = _pz(sp.expand(sp.chebyshevu(n - 1, z / 2))).rem(phi).is_zero
        if not (deg_ok and irr_ok and div_ok):
            bad.append((n, deg_ok, irr_ok, div_ok))
    check("S0  the node field Q(2cos(pi/n)): irreducible, degree phi(2n)/2, divides U_{n-1}",
          not bad and seen == n_max - 1,
          f"{seen} moduli n = 2..{n_max}, the range derived from the sweep and not pinned as a "
          f"literal; failures {bad}")

    # -- S1: the mode itself. psi_l = sin((j-l)theta_k)/sin(theta_k) solves every row of
    # -- (A + Delta_k*D)psi = 2cos(theta_k)*psi exactly in L, the two END rows included, which is
    # -- where Delta_k comes from. The matrix is H_aniso, the file's own family object, read
    # -- entry by entry. What this check does NOT assert is psi_j = 0 or psi != 0: psi_j is
    # -- _S_at(0, ...) and psi_0 is the pole test's own quantity, so both would be true of any
    # -- input and an earlier build named them as red-makers when nothing could redden them.
    bad, seen, poles = [], 0, 0
    for N in S_N:
        Hsym = H_aniso(N)
        for j in range(N):
            n = n_node(N, j)
            if n < 2:
                continue                              # no k in range at all
            for k in range(1, n):
                md = mode_data(N, j, k)
                if md is None:
                    poles += 1
                    continue
                phi, xk, dk, psi = md
                seen += 1
                for a in range(N):
                    r = -_redz(xk * psi[a], phi)
                    for b in range(N):
                        if Hsym[a, b] != 0:
                            r += _entry_times(Hsym[a, b], psi[b], dk, phi)
                    if not _redz(r, phi).is_zero:
                        bad.append((N, j, k, a))
    check("S1  psi_l = sin((j-l)theta_k) is an EXACT eigenvector of H_aniso at Delta_k",
          not bad and seen == S1_MODES and poles == S1_POLES,
          f"{seen} (N, seat, k) modes over N = 5..20, expected {S1_MODES}, every row the zero "
          f"element of Q[z]/(Phi), read off the family object H_aniso; {poles} poles skipped, "
          f"expected {S1_POLES}; failures {bad}")

    # -- S2: the two ingredients of the parity, separately, because the parity is their product
    # -- and a compensating pair of sign errors would be invisible in it. S_{N_node}(x_k) = 0
    # -- says the node angle closes at a multiple of pi; T_{N_node}(x_k) = (-1)^k says which.
    bad, seen = [], 0
    for N in S_N:
        for j in range(N):
            n = n_node(N, j)
            if n < 2:
                continue
            for k in range(1, n):
                seen += 1
                if (not _S_at(n, n, k).is_zero
                        or not (_T_at(n, n, k) - _pz((-1) ** k)).is_zero):
                    bad.append((N, j, k))
    check("S2  the node angle closes: S_{N_node}(x_k) = 0 and T_{N_node}(x_k) = (-1)^k",
          not bad and seen == S2_ANGLES,
          f"{seen} (N, seat, k) angles over N = 5..20, expected {S2_ANGLES}, exact in "
          f"Q[z]/(Phi); failures {bad}")

    # -- S2b: Lemma 8's parity, psi_{N-1-l} = (-1)^(k+1) psi_l at EVERY site. An earlier build
    # -- guarded it with "and the opposite sign at none", which cannot fire: both signs at once
    # -- force 2*eps*psi = 0 in characteristic zero, and psi = 0 is what mode_data's pole test
    # -- already excludes. The live control is S2d, which reads a DIFFERENT modulus.
    bad, seen, poles = [], 0, 0
    for N in S_N:
        for j in range(N):
            n = n_node(N, j)
            if n < 2:
                continue
            for k in range(1, n):
                md = mode_data(N, j, k)
                if md is None:
                    poles += 1
                    continue
                phi, xk, dk, psi = md
                seen += 1
                eps = (-1) ** (k + 1)
                if any(not _redz(psi[N - 1 - l] - eps * psi[l], phi).is_zero for l in range(N)):
                    bad.append((N, j, k))
    check("S2b psi_{N-1-l} = (-1)^(k+1) psi_l at every site",
          not bad and seen == S1_MODES and poles == S1_POLES,
          f"{seen} modes over N = 5..20, expected {S1_MODES}, {poles} poles skipped; "
          f"failures {bad}")

    # -- S2d: the control S2b's dead fence was pretending to be. The SAME construction on the
    # -- WRONG modulus, N_node + 1, keeps neither sign at any site pattern; so S2b reads the node
    # -- angle and not the shape of the formula. It asserts that something BREAKS, so no mutation
    # -- of the objects tests it: what does is feeding it the RIGHT modulus, under which it goes
    # -- red. It PASSES in a healthy block, like every other check here.
    survivors, seen = [], 0
    for N in S_N:
        for j in range(N):
            n = n_node(N, j)
            if n < 2:
                continue
            phi = phi_2cos(n + 1)
            for k in range(1, n):
                Sj = _S_at(j, n + 1, k)
                if Sj.is_zero:
                    continue
                psi = [_S_at(j - l, n + 1, k) for l in range(N)]
                seen += 1
                for eps in (+1, -1):
                    if all(_redz(psi[N - 1 - l] - eps * psi[l], phi).is_zero for l in range(N)):
                        survivors.append((N, j, k, eps))
    check("S2d control: the same mode built on the WRONG modulus carries NO reflection parity",
          len(survivors) == S2D_SURVIVORS and seen == S2D_TRIED,
          f"{seen} constructions at modulus N_node+1 over N = 5..20, expected {S2D_TRIED}; "
          f"{len(survivors)} of them carry either sign, expected {S2D_SURVIVORS}: {survivors[:6]}")

    # -- S2c: the Remark. Delta_{n-k} = -Delta_k and the pole condition shared by k and n-k are
    # -- its two readings. Its third clause in an earlier build, that k and n-k share a parity
    # -- exactly when n is even and that n is even exactly when N is odd, was deleted: both are
    # -- identities of the integers with no reference to any object here.
    bad, seen = [], 0
    for N in S_N:
        for j in range(N):
            n = n_node(N, j)
            if n < 2:
                continue
            phi = phi_2cos(n)
            for k in range(1, n):
                a, b = mode_data(N, j, k), mode_data(N, j, n - k)
                if (a is None) != (b is None):
                    bad.append((N, j, k, "pole condition not shared"))
                    continue
                if a is None:
                    continue
                seen += 1
                if not _redz(a[2] + b[2], phi).is_zero:
                    bad.append((N, j, k, "Delta_{n-k} != -Delta_k"))
    check("S2c Delta_{N_node-k} = -Delta_k, and k is a pole exactly when N_node-k is",
          not bad and seen == S2C_PAIRS,
          f"{seen} index pairs over N = 5..20, expected {S2C_PAIRS}; failures {bad}")

    # -- S3: the corroboration, and it is a second instrument and not a restatement. b_E and b_O
    # -- here come from the SECTOR BLOCKS' characteristic polynomials by a gcd over Q[t]/(mu),
    # -- block K's route, which never sees a mode or a k. The prediction comes from Lemma 8's
    # -- closed form. Two roads to one pair of integers.
    # --
    # -- The SWAPPED assignment needs no clause of its own and gets none. An earlier build put it
    # -- in a check of its own, S3b, and the repair then folded it into this predicate; both were
    # -- inert, because (b_E, b_O) = (n_odd, n_even) at a reading with n_odd != n_even already
    # -- excludes (n_even, n_odd). What refutes the swap is this check's own equality together
    # -- with the pinned count of readings that HAVE n_odd != n_even, and that count is live.
    bad, asym, mp_bad, balanced_off_zero, seen = [], 0, [], [], 0
    for N in S3_N:
        for j in range(N):
            if n_node(N, j) < 2:
                continue
            phi = phi_2cos(n_node(N, j))
            for dexpr, (n_odd, n_even) in _sector_split_by_k(N, j).items():
                mu = min_poly_of(_pz(dexpr), phi, t)
                muP = sp.Poly(mu, t)
                # the minimal polynomial is checked and not trusted: it annihilates the element
                # and it is irreducible, which is what the field Q[t]/(mu) below requires.
                cs = muP.all_coeffs()[::-1]
                ann = _redz(sum(cs[i] * _pz(dexpr) ** i for i in range(muP.degree() + 1)), phi)
                if not ann.is_zero or not muP.is_irreducible:
                    mp_bad.append((N, j, str(mu)))
                    continue
                bE, bO, _ = sector_counts_at(H_aniso, N, j, mu)
                seen += 1
                if n_odd != n_even:
                    asym += 1
                elif not sp.Poly(mu - t, t).is_zero:
                    balanced_off_zero.append((N, j, str(mu)))
                if (bE, bO) != (n_odd, n_even):
                    bad.append((N, j, str(mu), (bE, bO), (n_odd, n_even)))
    check("S3  b_E = #{k odd on that Delta} and b_O = #{k even}, against the sector blocks",
          not bad and not mp_bad and not balanced_off_zero and seen == S3_READINGS
          and asym == S3_ASYMMETRIC,
          f"{seen} (seat, Delta) readings over N = 5..14, expected {S3_READINGS}, of which "
          f"{asym} have #odd != #even and so exclude the swapped assignment, expected "
          f"{S3_ASYMMETRIC}; the other {S3_READINGS - S3_ASYMMETRIC} carry #odd = #even and "
          f"are asserted, not narrated, to sit at Delta = 0, exceptions {balanced_off_zero}; "
          f"minimal-polynomial failures {mp_bad}; count failures {bad}")

    # -- S4: the payoff, and the reason the item was open. The Theorem (count) gives the crack
    # -- 2*b_E at even N off the ring ends; Lemma 8 turns that into an integer read off F157's
    # -- own index set. Read against blind_at, which is block K's own reader on the full N x N
    # -- matrices and knows nothing about sectors or modes. A locus point carrying NO node index
    # -- predicts 0 and is READ rather than skipped, since that is the shape of the counterexample
    # -- this check exists to catch.
    bad, seen = [], 0
    for N in S4_N:
        for j in range(N):
            if n_node(N, j) < 2:
                continue
            phi = phi_2cos(n_node(N, j))
            split = _sector_split_by_k(N, j)
            for mu in locus_mus(N, j):
                # Several Delta_k can be CONJUGATE and share one minimal polynomial, which is
                # the rule rather than the exception at even N: at N = 8 seat 1 the two roots of
                # t^2 - t - 1 are Delta_1 and Delta_3. b_E is a gcd degree over Q[t]/(mu) and so
                # is Galois-invariant, which forces those to agree; a disagreement would be a
                # defect and is read as one rather than resolved by taking the last match. S3
                # would already have caught such a disagreement over a wider range, so this is a
                # defensive read and not an independent red-maker.
                hits = {n_odd for dexpr, (n_odd, _e) in split.items()
                        if sp.Poly(min_poly_of(_pz(dexpr), phi, t) - mu, t).is_zero}
                if len(hits) > 1:
                    bad.append((N, j, str(mu), "conjugate Delta_k disagree on #odd k", hits))
                    continue
                seen += 1
                got, pred = blind_at(H_crack, N, j, mu), 2 * (hits.pop() if hits else 0)
                if got != pred:
                    bad.append((N, j, str(mu), got, pred))
    check("S4  EVEN N off the ring ends: the crack pays 2*#{k odd on that Delta}, closed form",
          not bad and seen == S4_READINGS,
          f"{seen} (N, seat, locus point) readings over N = 6..14, expected {S4_READINGS}, "
          f"against blind_at on the full crack matrix; failures {bad}")




# ------------------------------------------- (W) the two constants section (g) and (h) leave open
#
# Section (g) reads P_j = +-Q_E*Q_O off a run and says WHICH sign it does not identify; section (h)
# names the roots of each factor and says the factorisation does not follow from root sets alone.
# Block W closes both, and it does so on a DIFFERENT ROAD from block K: every resultant that
# carries a LAW here is the Sylvester determinant built in this file, and sympy's own routine
# appears in W0b, W0c and the sympy leg of W4/W4b, where it is the object under test rather than
# an instrument. That distinction is not fastidiousness. It is
# what made the question answerable: sympy.resultant does not keep the order it is given, agreeing
# with the Sylvester determinant when deg f >= deg g and being (-1)^(deg f * deg g) times it below,
# so a sign law cannot be read off it without saying which order was meant (W0 to W0c). Block K is not wrong; it reads the pair, and W4 pins its
# committed literal from here and shows what it carries.
#
# The mechanism, in one congruence. Write S_m for sin(m*theta)/sin(theta) as a monic integer
# polynomial in x = 2*cos(theta), which is S_cheb above and block S's own S, and put
#
#     alpha_p = S_{p+1} - t*S_p     the charpoly of the path of p sites carrying t at coordinate 0
#     p = jr = min(j, N-1-j)        the fold coordinate,  n = N_node = |N-1-2j|
#
# Then S_p * alpha_{N-1-p} = -S_n (mod alpha_p): the addition formula twice, plus Cassini
# S_{p+1}^2 - S_p*S_{p+2} = 1, which is PROOF_CRACKED_RING_EXACT_CURVE's own gate P2. The two
# sector left halves are literally alpha_p, and the two right halves are the reflection blocks of
# the MIDDLE ROUTE, the uniform path the seat leaves between its two mirror images, whose comb is
# F71's on its own Dirichlet modulus. So Q_E*Q_O = Res(alpha_p, S_n) and the whole question is the
# common factor Res(alpha_p, S_p), which is a sign.

def res_std(f, g):
    """The Sylvester determinant of f and g in x, built here.

    This is the textbook resultant, lc(f)^deg(g) * prod over roots of f of g(root), and W0 checks
    it against that DEFINITION rather than against another resultant. sympy.resultant is not this
    when deg f < deg g; see W0b. Matrices are refused outright, because feeding one silently
    returns a number.
    """
    for e in (f, g):
        assert not isinstance(e, sp.MatrixBase), "res_std takes polynomials; matrices go to res_m"
    pf, pg = sp.Poly(f, x), sp.Poly(g, x)
    m, n = pf.degree(), pg.degree()
    if n <= 0:
        return sp.expand(pg.LC() ** max(m, 0))
    if m <= 0:
        return sp.expand(pf.LC() ** n)
    a, b = pf.all_coeffs(), pg.all_coeffs()
    M = sp.zeros(m + n, m + n)
    for i in range(n):
        for k, c in enumerate(a):
            M[i, i + k] = c
    for i in range(m):
        for k, c in enumerate(b):
            M[n + i, i + k] = c
    return sp.expand(M.det())


def res_m(M1, M2):
    """res_std of two matrices' characteristic polynomials; an empty block contributes Res(f,1)=1."""
    if M1.rows == 0 or M2.rows == 0:
        return sp.Integer(1)
    return res_std(charpoly_expr(M1).as_expr(), charpoly_expr(M2).as_expr())


def alpha_p(p):
    return sp.expand(S_cheb(p + 1, x) - t * S_cheb(p, x))


def beta_closed(n, parity):
    """The middle route's sector charpolys in closed form: beta_E carries odd k, beta_O even k."""
    if n <= 1:
        return sp.Integer(1)
    if n % 2 == 0:
        return sp.expand(sp.cancel(S_cheb(n, x) / S_cheb(n // 2, x))) if parity == +1 \
            else S_cheb(n // 2, x)
    r = (n - 1) // 2
    return sp.expand(S_cheb(r + 1, x) - S_cheb(r, x)) if parity == +1 \
        else sp.expand(S_cheb(r + 1, x) + S_cheb(r, x))


def eval_at_node(beta, n, k):
    """beta(x) evaluated at x_k = 2*cos(k*pi/n), as an element of Q[z]/Phi, by Horner.

    Block S's field arithmetic, reused: no root object is built and no trigonometric value ever
    reaches a simplifier, so "beta_E vanishes exactly at the odd k" is decided exactly.
    """
    phi, xk = phi_2cos(n), _x_k(n, k)
    acc = _pz(0)
    for c in sp.Poly(beta, x).all_coeffs():
        acc = _redz(acc * xk + _pz(int(c)), phi)
    return acc


def w_seats(N):
    """Interior, non-R-fixed seats. W3c says why the centre seat is out, rather than assuming it."""
    return [j for j in range(1, N - 1) if not _is_centre(N, j)]


def pole_split(beta, p):
    """(g, h, n_S, r_S): beta split against S_p into its pole part and its pole-free part.

    F157's pole indices are exactly the k with S_p(x_k) = 0, so the gcd IS the pole part; nothing
    is enumerated and no node index is built.
    """
    g = sp.gcd(sp.Poly(beta, x), sp.Poly(S_cheb(p, x), x))
    h = sp.Poly(sp.cancel(beta / g.as_expr()), x)
    return g, h, h.degree(), sp.Poly(beta, x).degree()


def w_constant(g, h, n_S, r_S, p):
    return ((-1) ** ((p * r_S + n_S) % 2)
            * res_std(h.as_expr(), S_cheb(p, x)) * res_std(g.as_expr(), S_cheb(p + 1, x)))


def w_monic(h, n_S, p):
    """prod over non-pole k of that parity of (t - Delta_k), exactly and with no root object.

    Delta_k = S_{p+1}(x_k)/S_p(x_k) and alpha_p = S_{p+1} - t*S_p, so the product of (t - Delta_k)
    is Res(h, -alpha_p)/Res(h, S_p). This is what lets W6 read the ROOTS; degree and leading
    coefficient cannot see them, which is how "every factor simple" survived a first build.
    """
    return sp.cancel((-1) ** (n_S % 2) * res_std(h.as_expr(), alpha_p(p))
                     / res_std(h.as_expr(), S_cheb(p, x)))


def w_Q(N, j, par, order="fold-first", conv="sylvester"):
    """Q_S. Both the CONVENTION (W0b) and the ORDER (W8b) are load-bearing and both are read."""
    p = _fold(N, j)
    S = block(H_aniso(N), N, par)
    L, R = S[:p, :p], S[p + 1:, p + 1:]
    a, b = (L, R) if order == "fold-first" else (R, L)
    if conv == "sympy":
        return (sp.expand(sp.resultant(charpoly_expr(a).as_expr(),
                                       charpoly_expr(b).as_expr(), x))
                if a.rows and b.rows else sp.Integer(1))
    return res_m(a, b)


def w_P(H, N, j, order="fold-first", conv="sylvester"):
    lo, hi = H[:j, :j], H[j + 1:, j + 1:]
    if order == "fold-first" and j > N - 1 - j:
        lo, hi = hi, lo
    if conv == "sympy":
        return sp.expand(sp.resultant(charpoly_expr(lo).as_expr(),
                                      charpoly_expr(hi).as_expr(), x))
    return res_m(lo, hi)


def w_q(N, j, order="fold-first", conv="sylvester", H=None):
    """The ratio P_j / (Q_E*Q_O).

    order 'fold-first' is the law's, 'seat' is block K1b's outer order, 'sector-reversed' is W8b's
    control. K1b already takes the fold half first INSIDE a sector, so 'seat' differs from the law
    in the outer argument order alone, which is what lets W4b separate the two causes.
    """
    H = H_aniso(N) if H is None else H
    outer = "seat" if order == "seat" else "fold-first"
    inner = "right-first" if order == "sector-reversed" else "fold-first"
    p, pr = _fold(N, j), sp.Integer(1)
    for par in (+1, -1):
        S = block(H, N, par)
        L, R = S[:p, :p], S[p + 1:, p + 1:]
        a, b = (L, R) if inner == "fold-first" else (R, L)
        if conv == "sympy":
            pr *= (sp.expand(sp.resultant(charpoly_expr(a).as_expr(),
                                          charpoly_expr(b).as_expr(), x))
                   if a.rows and b.rows else sp.Integer(1))
        else:
            pr *= res_m(a, b)
    return sp.cancel(w_P(H, N, j, outer, conv) / pr)


# PINNED MEASUREMENTS for block W, frozen from a run. The three sign readings are pinned APART
# because the total hides the per-seat difference: over N = 4..10 the sympy reading and the
# Sylvester reading at the same seat order agree in TOTAL while disagreeing at four seats, so the
# committed K1B_SIGNS is accidentally convention-stable in exactly the committed range.
W_SIGN = {(10, 'sympy'): {-1: 22, 1: 10}, (10, 'sylv'): {-1: 22, 1: 10},
          (10, 'canon'): {-1: 24, 1: 8},
          (12, 'sympy'): {-1: 34, 1: 16}, (12, 'sylv'): {-1: 32, 1: 18},
          (12, 'canon'): {-1: 34, 1: 16}}
W_FLIP_CONV = [(5, 1), (7, 1), (9, 3), (9, 7)]      # convention alone, N = 4..10
W_FLIP_BOTH = [(5, 1), (5, 3), (7, 1), (7, 5), (9, 3), (9, 5)]   # K1b's road against the law's
W_CAUSE_COUNTS = (6, 6)             # of the 50 seats at N = 4..12: convention alone, order alone
W_SECTOR_ORDER_BREAKS = 12          # of 50 seats, N = 4..12, taking Q_S right-half-first
W_ONEEND_SURVIVORS = 3              # of 24, N = 4..9: the R-breaking control
W_SEATS_PER_P = {1: 22, 2: 18, 3: 14, 4: 10, 5: 6, 6: 2}         # over N = 4..14
# W_REPEATED: the readings whose Q_S carries a REPEATED factor. Without it the words "with
# repetition" would be untested: over N = 4..9 this population is EMPTY, which is exactly how an
# earlier build of this block came to write "every factor simple" and stay green.
W_REPEATED = [(10, 2, +1), (10, 2, -1), (10, 7, +1), (10, 7, -1),
              (11, 2, -1), (11, 3, +1), (11, 7, +1), (11, 8, -1),
              (14, 3, +1), (14, 3, -1), (14, 4, +1), (14, 4, -1),
              (14, 9, +1), (14, 9, -1), (14, 10, +1), (14, 10, -1)]
W_N = range(4, 15)
W_SEATS = 72                        # non-centre interior seats over W_N
W_F157_SCALES = {(9, 1): -1, (9, 2): -1, (11, 1): -1, (11, 2): 1}
W1B_READINGS = (1260, 1014)         # (n, p, k) over n = 2..15, p = 0..11, and how many of
                                    # those are not the trivial 0 = 0 at S_p(x_k) = 0 or p = 0
W5_READINGS = (144, 110)            # (N, seat, sector) readings, and how many have a nonempty
                                    # right half rather than the trivial 1 = 1
W6_SHORTFALL = 20                   # readings where a pole root leaves deg_t Q_S below deg beta_S
W6C_DRAWS = (54, 0, 26, 4)            # measured below: squarefree draws, their failures, then
                                    # the non-squarefree draws and theirs
W0D_ASYM = (12, 55)                 # measured below: asymmetric pairs, and equal-degree pairs
W0C_READ = (33, 4)                  # W0c's own population: pairs with a nonzero resultant, and
                                    # how many of those sit at EQUAL degrees
W5B_READINGS = 500                  # (N, seat, sector, k) over N = 4..14
W5C_READINGS = 113                  # non-pole k at the 36 seats where j > N-1-j
W0_NONMONIC, W0_DEG0, W0_DIFF = 31, 15, 5   # of the 60 random pairs: lc(f) != 1; a degree-0
                                    # argument; and where sympy parts from the Sylvester value


def gate_W():
    print("\n(W) the two constants: the sign of section (g) and the factorisation of section (h)")

    # -- W0: the road. res_std is the Sylvester determinant, checked against the DEFINITION.
    import random
    random.seed(20260903)
    pairs = []
    for _ in range(60):
        m, n = random.randint(0, 5), random.randint(0, 5)
        lf, lg = random.choice([1, 1, 2, -3]), random.choice([1, 1, -2, 5])
        pairs.append((sp.expand(lf * sp.prod([x - random.randint(-3, 3) for _ in range(m)])),
                      sp.expand(lg * sp.prod([x - random.randint(-3, 3) for _ in range(n)])), m, n))
    bad, diff, nonmonic, deg0 = [], 0, 0, 0
    for f, g, m, n in pairs:
        pf, pg = sp.Poly(f, x), sp.Poly(g, x)
        prod = pf.LC() ** n                      # the factor a monic-only sweep never exercises
        for r in sp.roots(pf, multiple=True):
            prod *= pg.eval(r)
        if sp.expand(res_std(f, g) - prod) != 0:
            bad.append((m, n))
        if sp.expand(sp.resultant(f, g, x) - prod) != 0:
            diff += 1
        nonmonic += (pf.LC() != 1)
        deg0 += (m == 0 or n == 0)
    check("W0  res_std IS lc(f)^deg(g)*prod g(alpha), against the definition and not another routine",
          not bad and (nonmonic, deg0) == (W0_NONMONIC, W0_DEG0),
          f"60 random pairs of degrees 0..5, {nonmonic} with lc(f) != 1 and {deg0} with a "
          f"degree-0 argument, expected {W0_NONMONIC} and {W0_DEG0}: without those the "
          f"lc(f)^deg(g) factor and both early returns go untested; failures {bad}")
    check("W0b sympy.resultant DIFFERS from it, at a pinned population", diff == W0_DIFF,
          f"{diff} of 60 pairs at this seed, expected {W0_DIFF}; every resultant that carries a "
          f"LAW in block W is res_std, sympy's appearing only where it is the object under test")
    bad, read, equal = [], 0, 0
    for f, g, m, n in pairs:
        if sp.resultant(f, g, x) == 0:
            continue
        read += 1
        equal += (m == n)
        want = res_std(f, g) * (-1) ** ((m * n) % 2) if m < n else res_std(f, g)
        if sp.expand(sp.resultant(f, g, x) - want) != 0:
            bad.append((m, n))
    check("W0c the difference IS the swap: sympy = res_std times (-1)^(mn) exactly when deg f < deg g",
          not bad and (read, equal) == W0C_READ,
          f"{read} of 60 pairs are read, {equal} of them at EQUAL degrees where the rule's else "
          f"branch asserts plain agreement, expected {W0C_READ}; sympy is order-dependent there "
          f"too, correctly so, which is why the rule is stated by degree and not by symmetry; "
          f"failures {bad}")

    # -- W1: the congruence. This is the mechanism; everything below is its consequence.
    bad, seen = [], 0
    for N in W_N:
        for j in w_seats(N):
            p, n = _fold(N, j), abs(N - 1 - 2 * j)
            seen += 1
            lhs = sp.expand(S_cheb(p, x) * alpha_p(N - 1 - p) + S_cheb(n, x))
            if sp.rem(sp.Poly(lhs, x), sp.Poly(alpha_p(p), x)) != 0:
                bad.append((N, j))
    check("W1  the congruence S_p * alpha_{N-1-p} = -S_{N_node} (mod alpha_p)",
          not bad and seen == W_SEATS,
          f"{seen} seats over N = 4..14, expected {W_SEATS}; exact remainder in Z[t][x]; "
          f"failures {bad}")

    # -- W2: the common factor, in closed form and t-free.
    bad = []
    for p in range(1, 14):
        got = res_std(alpha_p(p), S_cheb(p, x))
        if sp.expand(got - (-1) ** (comb(p, 2) % 2)) != 0:
            bad.append((p, str(got)))
    check("W2  Res(alpha_p, S_p) = (-1)^binom(p,2), and it carries no t",
          not bad, f"p = 1..13; failures {bad}")

    # -- W3: the sign law. Fold half first, in the sectors as well as outside.
    bad, split, per_p = [], {}, {}
    for N in W_N:
        for j in w_seats(N):
            p = _fold(N, j)
            q = w_q(N, j)
            if sp.expand(q - (-1) ** (comb(p + 1, 2) % 2)) != 0:
                bad.append((N, j, str(q)))
            split[int(q)] = split.get(int(q), 0) + 1
            per_p[p] = per_p.get(p, 0) + 1
    check("W3  the sign law: Res(alpha_p, alpha_{N-1-p}) = (-1)^binom(p+1,2) * Q_E * Q_O",
          not bad and sum(split.values()) == W_SEATS,
          f"{sum(split.values())} seats over N = 4..14, split {split}; the exponent reads the fold "
          f"coordinate and N does not appear in it; failures {bad}")
    check("W3b the N-independence is thinnest at the largest fold coordinate, and that is pinned",
          per_p == W_SEATS_PER_P,
          f"seats per p: {per_p}, expected {W_SEATS_PER_P}; p = 6 rests on 2 seats over this "
          f"range, so the claim that N does not enter is thin exactly there and nowhere else")

    # -- W3c: the fence, MEASURED. At the R-fixed centre seat there is no ratio to have a sign.
    bad, seen = [], 0
    for N in W_N:
        if N % 2 == 0:
            continue
        j = (N - 1) // 2
        seen += 1
        if w_P(H_aniso(N), N, j) != 0:
            bad.append(N)
    check("W3c the centre seat is fenced BECAUSE P_j vanishes identically there, not by fiat",
          not bad and seen == 5,
          f"{seen} odd chains over N = 4..14; the seat cuts the chain into two EQUAL matrices, so "
          f"the outer resultant is 0 in t and section (a)'s third convention is why; failures {bad}")

    # -- W4 / W4b: what K1b's committed literal is a reading of.
    counts = {k: {} for k in W_SIGN}
    flip_conv, flip_both, conv_far, order_only = [], [], [], []
    for N in range(4, 13):
        for j in w_seats(N):
            a = int(w_q(N, j, order="seat", conv="sympy"))
            b = int(w_q(N, j, order="seat"))
            c = int(w_q(N, j))
            for hi in (10, 12):
                if N <= hi:
                    for key, v in ((( hi, 'sympy'), a), ((hi, 'sylv'), b), ((hi, 'canon'), c)):
                        counts[key][v] = counts[key].get(v, 0) + 1
            if a != b:
                (flip_conv if N <= 10 else conv_far).append((N, j))
            if b != c:
                order_only.append((N, j))
            if N <= 10 and a != c:
                flip_both.append((N, j))
    check("W4  K1b's committed K1B_SIGNS is reproduced, and the three readings are pinned apart",
          counts == W_SIGN,
          f"{ {k: counts[k] for k in sorted(counts)} }; the (10, 'sympy') row IS K1B_SIGNS")
    check("W4b the committed split carries BOTH causes, named seat by seat and not by count",
          flip_conv == W_FLIP_CONV and flip_both == W_FLIP_BOTH
          and (len(flip_conv) + len(conv_far), len(order_only)) == W_CAUSE_COUNTS,
          f"the convention alone flips {flip_conv} (expected {W_FLIP_CONV}); K1b's road against "
          f"the law's flips {flip_both} (expected {W_FLIP_BOTH}). At N = 4..10 the two TOTALS "
          f"agree while four seats do not, so K1B_SIGNS is accidentally convention-stable in "
          f"exactly the range it is pinned over, and is not so at N = 4..12. Over the whole "
          f"sweep N = 4..12 the CONVENTION alone moves {len(flip_conv) + len(conv_far)} of 50 "
          f"seats and the SEAT ORDER alone moves {len(order_only)}, expected {W_CAUSE_COUNTS}: "
          f"attributing the committed split to either cause alone would be wrong")

    # -- W5: the sector halves are the middle route's two combs.
    bad, seen, live = [], 0, 0
    for N in W_N:
        for j in w_seats(N):
            p, n = _fold(N, j), abs(N - 1 - 2 * j)
            for par in (+1, -1):
                R = block(H_aniso(N), N, par)[p + 1:, p + 1:]
                got = charpoly_expr(R).as_expr() if R.rows else sp.Integer(1)
                if sp.expand(got - beta_closed(n, par)) != 0:
                    bad.append((N, j, par))
                seen += 1
                live += bool(R.rows)
    check("W5  the sector right halves ARE the middle route's combs: beta_E odd k, beta_O even k",
          not bad and (seen, live) == W5_READINGS,
          f"{seen} (N, seat, sector) readings, closed form against the block, of which {live} "
          f"have a nonempty right half and are not 1 = 1, expected {W5_READINGS}; failures {bad}")

    # -- W6 / W6b: the factorisation of section (h), READ ON ITS ROOTS.
    bad, seen, repeated, shortfall = [], 0, [], 0
    for N in W_N:
        for j in w_seats(N):
            p = _fold(N, j)
            for par in (+1, -1):
                R = block(H_aniso(N), N, par)[p + 1:, p + 1:]
                beta = charpoly_expr(R).as_expr() if R.rows else sp.Integer(1)
                g, h, nS, rS = pole_split(beta, p)
                Q, c = w_Q(N, j, par), w_constant(g, h, nS, rS, p)
                if sp.expand(Q - c * w_monic(h, nS, p)) != 0:
                    bad.append((N, j, par, "product"))
                Qp = sp.Poly(Q, t)
                if Qp.degree() != nS or sp.expand(Qp.LC() - c) != 0:
                    bad.append((N, j, par, "degree or leading coefficient"))
                if nS >= 2 and sp.degree(sp.gcd(Qp, Qp.diff(t)), t) > 0:
                    repeated.append((N, j, par))
                shortfall += (nS != rS)
                seen += 1
    check("W6  Q_S = c_S * prod over the NON-POLE roots of beta_S of (t - Delta), one per root",
          not bad and seen == 2 * W_SEATS and shortfall == W6_SHORTFALL,
          f"{seen} readings, on the PRODUCT first and then on degree and leading coefficient. A "
          f"pole root contributes a constant and no factor, so deg_t Q_S is below deg beta_S at "
          f"{shortfall} of these, pinned. The identity holds for any monic SQUAREFREE beta (W6c), "
          f"so it reads the roots and NOT which node indices they are; the parity is W5b's and "
          f"W5's; failures {bad}")
    check("W6b the factors REPEAT, at a pinned set, so 'per index' is not idle wording",
          repeated == W_REPEATED,
          f"{len(repeated)} readings carry a repeated factor, expected {len(W_REPEATED)}: "
          f"{repeated}. The smallest is N = 10 seat 2, where the odd k = 1 and k = 3 give one "
          f"Delta and Q_E = -(t-1)^2; at N = 14 seat 3 the multiplicity is three. Over N = 4..9 "
          f"this population is EMPTY, which is why the sweep may not be narrowed")

    # -- W7 / W7b: the two indeterminacies compose, and the product is F157's own polynomial.
    D = sp.symbols("D")
    rows = {(9, 1): D**5 - 4*D**3 + 3*D, (9, 2): 2*D**2 - 1,
            (11, 1): D**7 - 6*D**5 + 10*D**3 - 4*D, (11, 2): 3*D**4 - 4*D**2}
    bad, seen, scales = [], 0, {}
    for N in W_N:
        for j in w_seats(N):
            p = _fold(N, j)
            P = w_P(H_aniso(N), N, j)
            cs, ns = {}, {}
            for par in (+1, -1):
                R = block(H_aniso(N), N, par)[p + 1:, p + 1:]
                beta = charpoly_expr(R).as_expr() if R.rows else sp.Integer(1)
                g, h, nS, rS = pole_split(beta, p)
                cs[par], ns[par] = w_constant(g, h, nS, rS, p), nS
            Pp = sp.Poly(P, t)
            want = (-1) ** (comb(p + 1, 2) % 2) * cs[+1] * cs[-1]
            if Pp.degree() != ns[+1] + ns[-1] or sp.expand(Pp.LC() - want) != 0:
                bad.append((N, j, Pp.degree(), ns[+1] + ns[-1], str(Pp.LC()), str(want)))
            seen += 1
            if (N, j) in rows:
                scales[(N, j)] = sp.cancel(P.subs(t, D) / rows[(N, j)])
    check("W7  the two constants compose: lc_t(P_j) = (-1)^binom(p+1,2) * c_E * c_O, deg = n_E+n_O",
          not bad and seen == W_SEATS, f"{seen} seats over N = 4..14; failures {bad}")
    check("W7b P_j built here IS F157's four committed rows, times a pinned sign",
          scales == W_F157_SCALES,
          f"scales {{{', '.join(f'{k}: {v}' for k, v in sorted(scales.items()))}}}, expected "
          f"{W_F157_SCALES}; F157 normalises its rows primitive with a positive leading "
          f"coefficient, so the sign is exactly what that normalisation discards")

    # -- W1b: Lemma 11, the node identity that carries the seat index onto the fold coordinate.
    # -- It rests on the two readings gate S2 already pins, S_n(x_k) = 0 and S_{n+1}(x_k) = (-1)^k,
    # -- and it is what makes F157's SEAT-indexed Delta_k this section's FOLD-indexed one.
    bad, seen, live = [], 0, 0
    for n in range(2, 16):
        phi = phi_2cos(n)
        for p in range(0, 12):
            for k in range(1, n):
                seen += 1
                sp_val = _S_at(p, n, k)
                live += not sp_val.is_zero
                if (_S_at(p + n, n, k) - _redz((-1) ** (k % 2) * sp_val, phi)).rem(phi) != 0:
                    bad.append((n, p, k))
    check("W1b Lemma 11: S_{p+n}(x_k) = (-1)^k * S_p(x_k) at every node of the modulus-n comb",
          not bad and (seen, live) == W1B_READINGS,
          f"{seen} (n, p, k) readings over n = 2..15 and p = 0..11, of which {live} have "
          f"S_p(x_k) != 0 and so are not 0 = 0, expected {W1B_READINGS}; exact in "
          f"Q(2cos(pi/n)); failures {bad}")

    # -- W5b: the content Corollary 10b actually needs and W6 does NOT read. W6's product identity
    # -- is true of any monic SQUAREFREE beta (W6c), so it says nothing about which indices sit in
    # -- sector; that is this check, and it is read in the node field rather than asserted.
    bad, seen = [], 0
    for N in W_N:
        for j in w_seats(N):
            n = abs(N - 1 - 2 * j)
            for par in (+1, -1):
                beta = beta_closed(n, par)
                for k in range(1, n):
                    seen += 1
                    vanishes = eval_at_node(beta, n, k).is_zero
                    if vanishes != ((k % 2 == 1) == (par == +1)):
                        bad.append((N, j, par, k))
    check("W5b beta_E vanishes at x_k for ODD k and beta_O for EVEN k, read in the node field",
          not bad and seen == W5B_READINGS,
          f"{seen} (N, seat, sector, k) readings over N = 4..14, expected {W5B_READINGS}; this is "
          f"the parity W6 is blind to, W6's identity holding for any monic SQUAREFREE beta; it "
          f"carries that jointly with W5 and neither alone; failures {bad}")

    # -- W5c: and the same lemma makes F157's own SEAT form of Delta_k equal the FOLD form this
    # -- section computes in. It is read only at the seats where the two differ as written.
    bad, seen = [], 0
    for N in W_N:
        for j in w_seats(N):
            if j <= N - 1 - j:
                continue
            n, p = abs(N - 1 - 2 * j), _fold(N, j)
            phi = phi_2cos(n)
            for k in range(1, n):
                Sj, Sp = _S_at(j, n, k), _S_at(p, n, k)
                if Sj.is_zero or Sp.is_zero:
                    if Sj.is_zero != Sp.is_zero:
                        bad.append(("pole", N, j, k))
                    continue
                seen += 1
                a = _redz(_S_at(j + 1, n, k) * _invz(Sj, phi), phi)
                b = _redz(_S_at(p + 1, n, k) * _invz(Sp, phi), phi)
                if (a - b).rem(phi) != 0:
                    bad.append((N, j, k))
    check("W5c F157's seat-indexed Delta_k IS this section's fold-indexed one, and the poles agree",
          not bad and seen == W5C_READINGS,
          f"{seen} non-pole k at the {sum(1 for N in W_N for j in w_seats(N) if j > N - 1 - j)} "
          f"seats where the two are different expressions, expected {W5C_READINGS}; failures {bad}")

    # -- W7c: and the object section (i) factorises IS F157's own definition route, with the sign
    # -- in closed form (Corollary 11b) rather than a reading at the four committed rows.
    bad, seen = [], 0
    for N in W_N:
        for j in w_seats(N):
            n, p = abs(N - 1 - 2 * j), _fold(N, j)
            G = res_std(S_cheb(n, x), sp.expand(t * S_cheb(j, x) - S_cheb(j + 1, x)))
            e = (n - 1) * (p + 1) + p + comb(p, 2) + (comb(n, 2) if j > N - 1 - j else 0)
            seen += 1
            if sp.expand(sp.cancel(w_P(H_aniso(N), N, j) / G) - (-1) ** (e % 2)) != 0:
                bad.append((N, j))
    check("W7c what is factorised here IS F157's generator Res(S_n, Delta*S_j - S_{j+1}), signed",
          not bad and seen == W_SEATS,
          f"{seen} seats over N = 4..14; the sign is Corollary 11b's closed form and the "
          f"[j > N-1-j] term is Lemma 11's price; failures {bad}")

    # -- W6c: squarefreeness is a HYPOTHESIS of W6's identity, not a convenience. beta_S has it
    # -- because beta_S divides S_n; a repeated pole root leaves one in h_S, Res(h_S, S_p) is then
    # -- 0 and the monic product divides by zero. Without this the section's "any monic beta"
    # -- would be a wider claim than anything measured, which is the shape the ledger records.
    random.seed(11)
    sf_bad = nsf_bad = sf = nsf = 0
    for _ in range(80):
        p = random.randint(1, 6)
        roots = [random.randint(-4, 4) for _ in range(random.randint(1, 5))]
        beta = sp.expand(sp.prod([x - r for r in roots]))
        squarefree = len(set(roots)) == len(roots)
        g, h, nS, rS = pole_split(beta, p)
        try:
            ok = sp.expand(res_std(alpha_p(p), beta) - w_constant(g, h, nS, rS, p)
                           * w_monic(h, nS, p)) == 0
        except ZeroDivisionError:
            ok = False
        if squarefree:
            sf += 1
            sf_bad += not ok
        else:
            nsf += 1
            nsf_bad += not ok
    check("W6c CONTROL: W6's identity NEEDS beta squarefree, and breaks without it",
          (sf, sf_bad, nsf, nsf_bad) == W6C_DRAWS,
          f"{sf} squarefree draws with {sf_bad} failures and {nsf} non-squarefree with {nsf_bad}, "
          f"expected {W6C_DRAWS}; so 'any monic beta' is false and 'any monic SQUAREFREE beta' is "
          f"what the section may say")

    # -- W0d: and sympy's routine is order-dependent at EQUAL degrees, correctly so. W0c compares
    # -- one fixed order against res_std and cannot see this, so the sentence that says it needs
    # -- its own reading.
    random.seed(11)
    asym, equal = 0, 0
    for _ in range(200):
        d1, d2 = random.randint(1, 4), random.randint(1, 4)
        a = sp.expand(random.choice([1, 2, -3, 5])
                      * sp.prod([x - random.randint(-3, 3) for _ in range(d1)]))
        b = sp.expand(random.choice([1, 2, -3, 5])
                      * sp.prod([x - random.randint(-3, 3) for _ in range(d2)]))
        equal += (d1 == d2)
        asym += sp.resultant(a, b, x) != sp.resultant(b, a, x)
    check("W0d sympy is ORDER-DEPENDENT at equal degrees, so it is not a function of the pair",
          (asym, equal) == W0D_ASYM,
          f"{asym} of 200 random pairs have sympy(f,g) != sympy(g,f), {equal} of the 200 sitting "
          f"at equal degrees, expected {W0D_ASYM}; the antisymmetry there is CORRECT, which is "
          f"why the rule W0c reads is stated by degree and never by symmetry")

    # -- W8: CONTROL. Break the reflection and the law must go, through the same door.
    surv, tried = [], 0
    for N in range(4, 10):
        M = path_sym(N)
        M[0, 0] += t                      # ONE end only: this does not commute with R
        for j in w_seats(N):
            tried += 1
            if w_q(N, j, H=M) == (-1) ** (comb(_fold(N, j) + 1, 2) % 2):
                surv.append((N, j))
    check("W8  CONTROL: on the one-end family the law goes, at a pinned survivor count",
          len(surv) == W_ONEEND_SURVIVORS and tried == 24,
          f"{len(surv)} of {tried} seats over N = 4..9 keep the value anyway, expected "
          f"{W_ONEEND_SURVIVORS}: {surv}; it is fed through the same door the real family uses")

    # -- W8b: CONTROL. The order inside the sectors is load-bearing too, not only the outer one.
    brk, tried = 0, 0
    for N in range(4, 13):
        for j in w_seats(N):
            tried += 1
            if w_q(N, j, order="sector-reversed") != (-1) ** (comb(_fold(N, j) + 1, 2) % 2):
                brk += 1
    check("W8b CONTROL: taking Q_S right-half-first breaks the law, at a pinned count",
          brk == W_SECTOR_ORDER_BREAKS and tried == 50,
          f"{brk} of {tried} seats over N = 4..12, expected {W_SECTOR_ORDER_BREAKS}; so 'fold half "
          f"first' is a statement about the SECTORS as well and not only about P_j")


# ------------------------------------------------- (H) the hop, and what the 128 in BlindSeat.cs is
H_SEATS_4_TO_11 = 40                 # interior non-centre seats over N = 4..11 (44 interior seats less 4 odd-N centres)
H_CONTENT_ROWS = 50                  # N = 4..12: interior non-centre seats, every quotient primitive


def H_f157(N, knob=t, hop=2):
    """F157's committed book (seat_cut_blindness.py se_hamiltonian_int, BlindSeat.H()), the ZZ
    coupling carrying Delta: hop on every bond, and a diagonal that pays -Delta at a bond's own two
    ends and +Delta elsewhere, bond by bond. That diagonal is Delta*(N-5) in the interior and
    Delta*(N-3) at the two ends, i.e. a common shift Delta*(N-5)*I plus 2*Delta on the end pair."""
    M = sp.zeros(N, N)
    for b in range(N - 1):
        M[b, b + 1] += hop
        M[b + 1, b] += hop
    for s_ in range(N):
        for b in range(N - 1):
            M[s_, s_] += (-knob if (s_ == b or s_ == b + 1) else knob)
    return M


def _halves_res(M, j):
    """Res_x(chi_L, chi_R) of the two principal submatrices seat j leaves behind, left first."""
    L, R = M[:j, :j], M[j + 1:, j + 1:]
    cL = (x * sp.eye(j) - L).det()
    cR = (x * sp.eye(M.rows - 1 - j) - R).det()
    return sp.expand(res_std(sp.expand(cL), sp.expand(cR)))


def gate_H():
    print("\n(H) the hop: F157's book is this file's times 2^(j(N-1-j)), the shift falling out")

    # -- H0: the transcription. H_f157 is built here, so before anything is read off it, it is
    #        compared entry for entry with the committed integer builder at Delta = 1, uniform J = 1,
    #        both books. Without this every check below could pass on a mis-copied matrix.
    import importlib.util as _ilu
    import os as _os
    _spec = _ilu.spec_from_file_location(
        "seat_cut_blindness", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "seat_cut_blindness.py"))
    _scb = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_scb)
    off = []
    for N in range(3, 12):
        bonds = [(b, b + 1, 1) for b in range(N - 1)]
        theirs = sp.Matrix(_scb.se_hamiltonian_int(N, bonds, zz=True))
        if H_f157(N, knob=1) != theirs:
            off.append(N)
    check("H0  H_f157 at Delta = 1 is se_hamiltonian_int's matrix entry for entry, N = 3..11",
          not off, f"off at {off}")

    # -- H1: the identity, exact in Delta. Both sides are the SAME left-first resultant of the two
    #        half-blocks, one in F157's book and one in this file's, so the only things between them
    #        are the hop and the shift; the claim is that the shift leaves no trace and the hop
    #        leaves exactly its power at the resultant's bidegree.
    bad, seats, zero_both = [], 0, 0
    for N in range(4, 12):
        for j in range(1, N - 1):
            a, b = j, N - 1 - j
            lhs = _halves_res(H_f157(N), j)
            rhs = sp.expand(2 ** (a * b) * _halves_res(H_aniso(N), j))
            if lhs == 0 and rhs == 0:
                zero_both += 1          # the centre seat: both books' resultants vanish identically
                continue
            seats += 1
            if sp.expand(lhs - rhs) != 0:
                bad.append((N, j))
    check("H1  Res_F157(Delta) = 2^(j(N-1-j)) * Res_here(t = Delta), exact at every non-centre seat",
          not bad and seats == H_SEATS_4_TO_11 and zero_both == 4,
          f"{seats} seats over N = 4..11 (expected {H_SEATS_4_TO_11}), {zero_both} centre seats "
          f"vanishing in both books (expected 4, the odd N); off at {bad}")

    # -- H1b: MUTATION. One more power of 2 reddens every seat, so H1 reads the exponent.
    brk = 0
    for N in range(4, 12):
        for j in range(1, N - 1):
            a, b = j, N - 1 - j
            rhs = _halves_res(H_aniso(N), j)
            if rhs == 0:
                continue
            if sp.expand(_halves_res(H_f157(N), j) - 2 ** (a * b + 1) * rhs) != 0:
                brk += 1
    check("H1b MUTATION: exponent j(N-1-j)+1 breaks H1 at every non-centre seat",
          brk == H_SEATS_4_TO_11, f"{brk} of {H_SEATS_4_TO_11} broken")

    # -- H1c: the J half of the hop. At a uniform J the hop is 2J and the ZZ diagonal J*Delta, so
    #         the common factor is 2J and the constant is (2J)^(ab) while the roots stay at Delta.
    #         Read at J = 2 and 3, N = 4..8, the two surfaces that say "(2J)^(j(N-1-j)) in general"
    #         (the proof's Remark, BlindSeat.cs) being otherwise ungated.
    bad, tried = [], 0
    for J in (2, 3):
        for N in range(4, 9):
            for j in range(1, N - 1):
                a, b = j, N - 1 - j
                rhs = _halves_res(H_aniso(N), j)
                if rhs == 0:
                    continue
                tried += 1
                MJ = sp.zeros(N, N)
                for bb in range(N - 1):
                    MJ[bb, bb + 1] += 2 * J
                    MJ[bb + 1, bb] += 2 * J
                for s_ in range(N):
                    for bb in range(N - 1):
                        MJ[s_, s_] += (-J * t if (s_ == bb or s_ == bb + 1) else J * t)
                if sp.expand(_halves_res(MJ, j) - (2 * J) ** (a * b) * rhs) != 0:
                    bad.append((J, N, j))
    check("H1c uniform J: the constant is (2J)^(j(N-1-j)) and the roots stay at Delta, J = 2 and 3 over N = 4..8",
          not bad and tried == 36, f"{tried} seats (expected 36); off at {bad}")

    # -- H2: the shift half. A COMMON shift c*I leaves the halves-resultant untouched, symbolically
    #        in c; the same c on ONE half only does not. Without the second half this check would
    #        pass on any matrix whose resultant ignores its diagonal altogether.
    c = sp.symbols('c')
    bad_common, bad_one, tried = [], 0, 0
    for N in range(4, 10):
        for j in range(1, N - 1):
            M = H_aniso(N)
            r0 = _halves_res(M, j)
            if r0 == 0:
                continue
            tried += 1
            if sp.expand(_halves_res(M + c * sp.eye(N), j) - r0) != 0:
                bad_common.append((N, j))
            Mo = M.copy()
            for i_ in range(j + 1, N):
                Mo[i_, i_] += c
            if sp.expand(_halves_res(Mo, j) - r0) == 0:
                bad_one += 1
    check("H2  a COMMON shift c*I drops out of the halves-resultant, symbolic in c; a one-sided shift never does",
          not bad_common and bad_one == 0 and tried == 24,
          f"{tried} seats over N = 4..9; common-shift residue at {bad_common}; one-sided shift invisible at {bad_one}")

    # -- H3: the hop half is a LAW, not a 2, and it shows which 2 is which. F157's end pair carries
    #        2*Delta because an end site has ONE bond fewer than an interior site and the ZZ diagonal
    #        pays -J at a bond's own ends and +J elsewhere, a jump of 2J per lost bond; the hop is the
    #        OTHER 2J, XX+YY on one bond. So with hop h and the ZZ diagonal unchanged the factor is
    #        h^(ab) and this file's t reads 2*Delta/h, which is Delta at every Delta exactly when
    #        h = 2. An earlier build of this check asked hop 3 for t = Delta
    #        and went red at 16 of 24 seats; that was the check reading the two 2s as one.
    bad = []
    for N in range(4, 10):
        for j in range(1, N - 1):
            a, b = j, N - 1 - j
            rhs = _halves_res(H_aniso(N), j)
            if rhs == 0:
                continue
            pred = sp.expand(3 ** (a * b) * rhs.subs(t, 2 * t / 3))
            if sp.expand(_halves_res(H_f157(N, hop=3), j) - pred) != 0:
                bad.append((N, j))
    check("H3  hop 3, ZZ unchanged: Res = 3^(ab) * Res_here(t = 2*Delta/3) at every non-centre seat over N = 4..9",
          not bad, f"off at {bad}")

    # -- H4: a READING, not a claim. Over N = 4..12 the quotient Res_F157 / 2^(ab) is primitive in
    #        Delta (content 1), so the hop's power is the WHOLE of the un-normalisation up to sign;
    #        nothing else hides in BlindSeat.cs's 128. Pinned as a count because no lemma here says
    #        it must hold at every N.
    rows, prim = 0, 0
    for N in range(4, 13):
        for j in range(1, N - 1):
            a, b = j, N - 1 - j
            r = _halves_res(H_f157(N), j)
            if r == 0:
                continue
            rows += 1
            q = sp.Poly(sp.expand(r / 2 ** (a * b)), t)
            if abs(sp.gcd_list(q.all_coeffs())) == 1:
                prim += 1
    check("H4  READING: the quotient is primitive in Delta at every non-centre seat over N = 4..12",
          rows == H_CONTENT_ROWS and prim == rows, f"{prim} of {rows} primitive, expected {H_CONTENT_ROWS} rows")

    # -- H5: the literal BlindSeat.cs and Formulas.cs quote, in its own book, and its sign. The
    #        bidegree 1*7 is odd, so Res(chi_L, chi_R) = -Res(chi_R, chi_L) and the literal's sign is
    #        a statement about an ORDER, as section (i) says of every sign here. Left-first, the order
    #        this whole block names, the Sylvester determinant gives -128 and right-first +128. The
    #        +128 the two MirrorWorld comments carry is what sympy.resultant returns, and it returns
    #        it in BOTH orders here, which is section (i)'s rule at deg 1 < deg 7 (W0c) and is gated
    #        as such. The 128 is identified either way; the sign was never a number of the physics.
    #        (2^7 == 128 is not asserted: a literal identity cannot fail.)
    # The literal is READ from BlindSeat.cs, not retyped, so a drift in the comment reddens this.
    import re as _re
    _bs = open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "compute", "MirrorWorld",
                             "BlindSeat.cs"), encoding="utf-8").read()
    _m = _re.search(r"//\s*(128\*D\*\(D-1\)\*\(D\+1\)\*\(D\^2-3\))", _bs)
    lit = sp.expand(sp.sympify(_m.group(1).replace("^", "**").replace("D", "t"))) if _m else None
    left_first = _halves_res(H_f157(9), 1)
    M = H_f157(9)
    cL = sp.expand((x * sp.eye(1) - M[:1, :1]).det())
    cR = sp.expand((x * sp.eye(7) - M[2:, 2:]).det())
    right_first = sp.expand(res_std(cR, cL))
    sym_lr = sp.expand(sp.resultant(cL, cR, x))
    sym_rl = sp.expand(sp.resultant(cR, cL, x))
    # F157's generator in the monic 2cos book: S_m = U_m(x/2) monic, P_1 = Res_x(S_5, t*S_0 - S_1).
    S_ = lambda m: sp.expand(sp.chebyshevu(m, x / 2))
    P_1 = sp.expand(res_std(S_(5), sp.expand(t * S_(0) - S_(1))))
    check("H5  N = 9 seat 1: BlindSeat.cs's literal read from the file is right-first +128*(...), left-first -128, sympy +128 both orders, and 128*P_1 (monic book) = right-first",
          lit is not None
          and sp.expand(left_first + lit) == 0 and sp.expand(right_first - lit) == 0
          and sp.expand(sym_lr - lit) == 0 and sp.expand(sym_rl - lit) == 0
          and sp.expand(128 * P_1 - right_first) == 0,
          f"literal {sp.factor(lit) if lit is not None else 'NOT FOUND in BlindSeat.cs'}; left-first {sp.factor(left_first)}; "
          f"right-first {sp.factor(right_first)}; sympy {sp.factor(sym_lr)} / {sp.factor(sym_rl)}; P_1 {sp.factor(P_1)}")


def main():
    print("=" * 86)
    print("PROOF_BLIND_SEAT_TWO_AXES: the crack axis and the anisotropy axis, sector by sector")
    print("=" * 86)
    gate_L1()
    gate_L2()
    gate_B()
    gate_P()
    gate_T()
    gate_K()
    gate_S()
    gate_W()
    gate_H()
    print()
    print("=" * 86)
    print("VERDICT:", "ALL GREEN" if not _fails else f"{len(_fails)} FAILED: {_fails}")
    print("=" * 86)
    return 1 if _fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
