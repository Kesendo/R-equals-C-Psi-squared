"""F157's blind count split into a forced half and an accidental half, and its Delta locus
carried onto the crack's u axis.

WHAT IS GATED, in the order the work found it.

(a) ON THE UNIFORM CHAIN THE COUNT IS AN ALIASING DEGREE. F157's uniform XY law
    blind(j) = gcd(j+1, N+1) - 1 is the aliasing degree of multiplication-by-(j+1) on Z/n,
    n = N+1, minus one: n divided by the size of that map's image. This is the WHOLE count,
    not one part of it: at the centre seat of an odd chain it returns the same (N-1)/2 that the
    reflection forces, because the arithmetic finds the reflection-odd modes too and does not
    know they are forced. Forced versus accidental is not a decomposition of this number; it is
    a statement about which SEATS keep their count when the chain is spoiled. The left-hand side is
    computed from F157's DEFINITION (N minus the seat's Krylov rank, exact elimination) and
    no gcd is called on that side. For ODD multipliers the seat is sighted exactly when the
    multiplier is a Galois automorphism of Q(zeta_2n), which is the condition a remark in
    F161's section (e) requires on the same modulus (Theorem D's own Statement is the X ladder). At EVEN multipliers the two differ and the gate
    exhibits that rather than skipping it.

(b) WHICH SEATS KEEP THE COUNT. Where a reflection of H fixes a seat, the reflection-odd
    modes vanish there at every coupling, so that seat's blindness is a representation count
    and not an arithmetic one. Gated as the identically-vanishing resultant at the centre
    seat of an odd chain, and as the constant ring-end value at BOTH parities of N.

(c) WHERE THE ACCIDENTAL BLINDNESS GOES, AND IT IS ALREADY OWNED. Off the comb the
    accidental blindness survives only on the real roots of a polynomial. On the Delta axis
    that polynomial is committed: F157's P_j(Delta) = Res_x(U_{Nnode-1}(x),
    Delta*U_{j-1}(x) - U_j(x)) with Nnode = |N-1-2j| (docs/ANALYTICAL_FORMULAS.md, live
    witness `inspect --root blindlocus`). This file gates that the u axis carries the SAME
    locus together with the ring ends {+1,-1}, at odd N, which hands the u axis that closed
    form; and that the parity restriction is exactly the staggering identity
    Sigma*H(x)*Sigma = -H(-x), which holds for the crack only at odd N and for the
    anisotropy at every N.

(d) WHICH PERTURBATIONS CANNOT MOVE A BLIND SEAT AT ALL. Block (c)'s CONTAINS column counts two
    things under one number, and at many seats the perturbation leaves the resultant identically
    zero in its own knob: the seat is blind at EVERY knob value and contains every locus for free.
    That branch is read by ONE criterion for both perturbation kinds, the overlap of a blind
    eigenvector of the unperturbed chain with the perturbation direction, a sum of squares for a
    diagonal and a product for a bond. Corollary D of
    docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md, gates D0 to D8b.

Every locus is SOLVED (the real roots of a resultant), never sampled: gate B1 evaluates the
blindness ON the solved locus, because a sampled interior point misses it by construction.
Gates C1a and C1b take their expectations as LITERALS from F157's committed table.

Companion page: experiments/THE_BLIND_SEAT_ON_THE_ROAD.md
"""
from fractions import Fraction
import itertools
import math

import sympy as sp

# ------------------------------------------------- part (a): the accidental half

def _path_int(N):
    H = [[0] * N for _ in range(N)]
    for i in range(N - 1):
        H[i][i + 1] = H[i + 1][i] = 1
    return H

def _rank_exact(rows):
    M = [[Fraction(v) for v in r] for r in rows]
    nr, nc = len(M), len(M[0])
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [v / pv for v in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == nr:
            break
    return r

def blind_by_krylov(N, j):
    """F157's DEFINITION: N minus the seat's Krylov dimension. Exact. No gcd."""
    H = _path_int(N)
    v = [0] * N
    v[j] = 1
    rows = []
    for _ in range(N):
        rows.append(list(v))
        v = [sum(H[i][k] * v[k] for k in range(N)) for i in range(N)]
    return N - _rank_exact(rows)

def aliasing_degree(m, mod):
    """|Z/mod| divided by |image of multiplication by m|. Brute force. No gcd."""
    return mod // len({(m * k) % mod for k in range(mod)})

def is_automorphism(m, mod):
    """Is k -> k*m a permutation of Z/mod? Brute force. No gcd."""
    return len({(m * k) % mod for k in range(mod)}) == mod

# ----------------------------------------------- parts (b) and (c): the two axes

t, x = sp.symbols('t x')

def _path_sym(N):
    H = sp.zeros(N, N)
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = 1
    return H

def H_u(N):
    """The crack: the wrap bond carries t. The end pair's OFF-diagonal."""
    H = _path_sym(N)
    H[N - 1, 0] = H[0, N - 1] = t
    return H

def H_delta(N):
    """XXZ anisotropy on the single-excitation block, in F157's scale. The end pair's
    DIAGONAL. F157's committed convention (docs/ANALYTICAL_FORMULAS.md) writes hop 2 and
    +2*Delta per end; the Delta VALUES are unchanged by that common scaling, and gates C1a
    and C1b check this convention against F157's committed roots."""
    H = _path_sym(N)
    H[0, 0] = t
    H[N - 1, N - 1] = t
    return H

def H_midbond(N):
    """CONTROL: a detuned bond that is not the wrap bond."""
    H = _path_sym(N)
    c = N // 2
    H[c, c + 1] = H[c + 1, c] = t
    return H

def H_symdiag(k):
    """CONTROL: a diagonal shift on the REFLECTION-SYMMETRIC interior pair (k, N-1-k).
    H_u and H_delta both commute with the chain reflection at every parameter value, so a
    control must too; otherwise it differs from the object in more ways than the one under
    test and isolates nothing."""
    def build(N):
        if k == N - 1 - k:
            return None            # the pair degenerates to the reflection-FIXED centre site
        H = _path_sym(N)
        H[k, k] = t
        H[N - 1 - k, N - 1 - k] = t
        return H
    return build

def H_onesite(k):
    """CONTROL for part (d): the SAME diagonal shift on the single interior site k alone.
    It breaks the reflection the other controls keep, which is exactly what gate D7 reads."""
    def build(N):
        if k <= 0 or k >= N - 1:
            return None
        H = _path_sym(N)
        H[k, k] = t
        return H
    return build

def H_symbond(k):
    """CONTROL: the reflection-symmetric interior BOND pair (k,k+1) and (N-2-k,N-1-k)."""
    def build(N):
        H = _path_sym(N)
        H[k, k + 1] = H[k + 1, k] = t
        H[N - 2 - k, N - 1 - k] = H[N - 1 - k, N - 2 - k] = t
        return H
    return build

def _strike(M, j):
    idx = [i for i in range(M.shape[0]) if i != j]
    return M[idx, idx]

def blind_fence_free(H, j):
    """F157's fence-free form: deg gcd(chi(H), chi(H with row and column j struck))."""
    a = sp.Poly(H.charpoly(x).as_expr(), x)
    b = sp.Poly(_strike(H, j).charpoly(x).as_expr(), x)
    return sp.gcd(a, b).degree()

def _resultant(Hf, N, j):
    """Res_x(chi(H), chi(H with row and column j struck)) as a polynomial in the knob.

    Returns the string 'SKIP' when the builder declines this N (a symmetric pair that would
    degenerate onto the reflection-fixed centre site).
    """
    H = Hf(N)
    if H is None:
        return 'SKIP'
    R = sp.resultant(sp.Poly(H.charpoly(x).as_expr(), x),
                     sp.Poly(_strike(H, j).charpoly(x).as_expr(), x), x)
    return None if sp.expand(R) == 0 else sp.Poly(sp.expand(R), t)

def _squarefree(P):
    return sp.Poly(sp.quo(P, sp.gcd(P, P.diff(t))), t)

def _n_real(P):
    """Number of DISTINCT real roots, by Sturm's theorem. No root object is ever built."""
    return _squarefree(P).count_roots()

def loci_equal(P, Q):
    """Do two polynomials have the same real root set? Exact.

    Deliberately NOT done by extracting roots and comparing them. sp.roots returns an
    irreducible cubic's roots in casus irreducibilis form, where `.is_real` is None rather
    than True; a truthiness filter then DROPS them silently. That defect stood in this file
    and made every even-N count wrong. Squarefree parts plus Sturm counting cannot express it.
    """
    A, B = _squarefree(P), _squarefree(Q)
    g = sp.Poly(sp.gcd(A, B), t)
    if g.degree() < 1:
        return _n_real(A) == 0 and _n_real(B) == 0
    return _n_real(g) == _n_real(A) == _n_real(B)

def locus_is_identical(Hf, N, j):
    return _resultant(Hf, N, j) is None


# ------------------------- part (d): blindness at every knob value, one overlap law

# Block (c)'s CONTAINS column counts two different things under one number. At some seats the
# perturbation leaves the resultant identically zero in the knob: the seat is blind at EVERY KNOB
# VALUE, its locus is the whole line, and it therefore contains every locus for free. Part (d)
# separates that branch from the rest and reads it as ONE criterion for both perturbation kinds.
#
# The criterion is the overlap. A perturbation direction V cannot move a level whose eigenvector
# does not overlap it, and for the two kinds the overlap is a different shape:
#     a diagonal on the site set S:  v^T V v = sum over m in S of v_m^2   (a sum of SQUARES)
#     a bond (b, b+1):               v^T V v = 2 * v_b * v_(b+1)          (a PRODUCT)
# so the diagonal needs a node on EVERY site of S and the bond needs one at EITHER end.
#
# h = gcd(j+1, N+1) is F157's own letter for this integer, and the entry that owns it
# reserves g for gcd(2j+1, N), the OTHER book; h's aliasing degree minus one is
# the uniform XY count (block (a) above reads the aliasing on the modulus n = N+1, which is a
# different number and gate A4 mutates it), and F157 also gives h as the node-count
# reduction gcd(j+1, N_node). The BLIND EIGENVECTORS at seat j of the unperturbed path,
# the eigenvectors that vanish there and span what the seat cannot touch, are
# v_l = sin((j-l)*c*pi/h) for c = 1..h-1, derived in Corollary D of
# PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA, so v vanishes at site m exactly when h divides c*(j-m).

def _cp_tridiag(diag, sym, bonds=None):
    """Characteristic polynomial of a tridiagonal matrix, by the Jacobi three-term recursion.

    `bonds[i]` is the entry between sites i and i+1, default 1 throughout.
    """
    p, pm1 = sp.Integer(1), sp.Integer(0)
    for i, a in enumerate(diag):
        b = 1 if bonds is None or i == 0 else bonds[i - 1]
        p, pm1 = sp.expand((sym - a) * p - b**2 * pm1), p
    return p

def _bond_halves(N, j, bond, knob):
    """The two halves when the single bond (bond, bond+1) carries the knob instead of 1."""
    left = [1] * max(j - 1, 0)
    right = [1] * max(N - j - 2, 0)
    if bond < j - 1:
        left[bond] = knob
    elif bond > j:
        right[bond - j - 1] = knob
    # a bond incident on the seat is struck away with it and carries nothing
    return left, right

def _bond_shares_root_at(N, j, bond, tv):
    bl, br = _bond_halves(N, j, bond, tv)
    dl, dr = [0] * j, [0] * (N - j - 1)
    if not dl or not dr:
        return False
    return sp.gcd(sp.Poly(_cp_tridiag(dl, x, bl), x),
                  sp.Poly(_cp_tridiag(dr, x, br), x)).degree() >= 1

def bond_blind_at_every_knob(N, j, bond):
    """The bond twin of `blind_at_every_knob`, same two routes and the same honest status:
    the integer probes decide every cell in the swept range and the symbolic resultant is
    never the one that answers."""
    for tv in _PROBES:
        if not _bond_shares_root_at(N, j, bond, tv):
            return False
    bl, br = _bond_halves(N, j, bond, t)
    dl, dr = [0] * j, [0] * (N - j - 1)
    return sp.expand(sp.resultant(sp.Poly(_cp_tridiag(dl, x, bl), x),
                                  sp.Poly(_cp_tridiag(dr, x, br), x), x)) == 0

def _halves_diag(N, j, S, knob):
    """The two principal submatrices' diagonals when the knob sits on the site set S."""
    return ([knob if i in S else 0 for i in range(0, j)],
            [knob if i in S else 0 for i in range(j + 1, N)])

def _shares_root_at(N, j, S, tv):
    dl, dr = _halves_diag(N, j, S, tv)
    if not dl or not dr:
        return False                       # an end seat: Lemma J1 forbids blindness outright
    return sp.gcd(sp.Poly(_cp_tridiag(dl, x), x),
                  sp.Poly(_cp_tridiag(dr, x), x)).degree() >= 1

_PROBES = (1, 2, 3, 5, 7, 11, 13, 17, 19, 23)

def blind_at_every_knob(N, j, S):
    """Is seat j blind at EVERY value of the knob carried by the diagonal site set S? Exact.

    A single integer knob value with no shared root PROVES the resultant is not identically
    zero, which is why the probes come first; only the survivors pay for the symbolic
    resultant. Both routes are exact, neither is a sample of a continuum.

    HONEST STATUS OF THE TWO ROUTES: over the populations swept below they never disagree,
    and no cell exists where the ten probes agree and the symbolic resultant is nonzero. So
    the symbolic route is unexercised here and the redundancy is a design property rather
    than a measured one; removing either route leaves every number in this file unchanged.
    It is kept because the probes alone would make the verdict a sample, and the absence of
    a separating cell in range is itself the finding rather than a reason to drop it.
    """
    for tv in _PROBES:
        if not _shares_root_at(N, j, S, tv):
            return False
    dl, dr = _halves_diag(N, j, S, t)
    return sp.expand(sp.resultant(sp.Poly(_cp_tridiag(dl, x), x),
                                  sp.Poly(_cp_tridiag(dr, x), x), x)) == 0

def overlap_is_zero(N, j, V, modulus=None):
    """Does a blind eigenvector at seat j have zero overlap with the perturbation direction V?

    V is given as ('diag', site set) or ('bond', b). The two shapes of v^T V v are a sum of
    squares and a product, so the diagonal asks for a node on every site and the bond for a
    node at either end. One criterion, two faces.

    Basis and not span, deliberately: A's spectrum is simple, so every eigenvector vanishing
    at the seat is a multiple of one v^(c) and the two readings agree.
    """
    h = math.gcd(j + 1, N + 1) if modulus is None else modulus(N, j)
    if h < 2:
        return False
    kind, arg = V
    for c in range(1, h):
        if kind == 'diag':
            if all((c * (j - m)) % h == 0 for m in arg):
                return True
        else:
            if (c * (j - arg)) % h == 0 or (c * (j - arg - 1)) % h == 0:
                return True
    return False

def _interior_sets(N, j, size):
    """Site sets for the diagonal sweeps: interior sites only, the watched seat excluded.

    Two exclusions, both deliberate and both narrowing what the sweeps certify. The seat's own
    site is struck away with the seat and changes neither block, so it is trivial. The two chain
    ends are dropped because the perturbations this file compares are interior ones. Nothing
    below is measured at an end site or at the seat.
    """
    sites = [i for i in range(1, N - 1) if i != j]
    return itertools.combinations(sites, size)


# ------------------------------------------------------------------------ gates

_fails = []

def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else '**FAIL**'}] {name}  {detail}")
    if not ok:
        _fails.append(name)

def locus_contains(Pbig, Psmall):
    """Is every real root of Psmall also a real root of Pbig? Exact, same route as loci_equal."""
    A, B = _squarefree(Pbig), _squarefree(Psmall)
    if _n_real(B) == 0:
        return True
    g = sp.Poly(sp.gcd(A, B), t)
    return g.degree() >= 1 and _n_real(g) == _n_real(B)

def _plus_ring_ends(P):
    """The polynomial whose real roots are P's together with +1 and -1."""
    return sp.Poly(P.as_expr() * (t**2 - 1), t)

def _matches(Pa, Pb):
    """Does axis A's locus equal axis B's locus together with the ring ends?"""
    if Pa == 'SKIP' or Pb == 'SKIP':
        return False
    if Pa is None or Pb is None:
        return Pa is None and Pb is None
    return loci_equal(Pa, _plus_ring_ends(Pb))

def main():
    A_RANGE = range(3, 31)          # one range for the whole of block (a)
    ODD = [5, 7, 9, 11, 13, 15, 17]     # block (c)'s odd range
    EVEN = [6, 8, 10, 12, 14, 16]       # and its even one
    ODD_SEATS = sum(ODD)
    lo, hi = A_RANGE.start, A_RANGE.stop - 1

    print("=" * 86)
    print(f"(a) ON THE UNIFORM CHAIN THE COUNT IS AN ALIASING DEGREE   [N={lo}..{hi}]")
    print("=" * 86)

    bad, nontrivial = [], 0
    for N in A_RANGE:
        for j in range(N):
            d = aliasing_degree(j + 1, N + 1)
            if blind_by_krylov(N, j) + 1 != d:
                bad.append((N, j))
            if d > 1:
                nontrivial += 1
    check("A1 blind(j)+1 == aliasing degree of mu_{j+1} on Z/n, every seat",
          not bad, f"({nontrivial} seats alias, so the identity is not 1 == 1)")

    bad, notsighted = [], 0
    for N in A_RANGE:
        for j in range(N):
            if (j + 1) % 2 == 0:
                continue
            sighted = blind_by_krylov(N, j) == 0
            if sighted != is_automorphism(j + 1, 2 * (N + 1)):
                bad.append((N, j))
            if not sighted:
                notsighted += 1
    check("A2 ODD multiplier: seat sighted <=> mu_{j+1} is a Galois automorphism of Q(zeta_2n)",
          not bad, f"({notsighted} of the tested seats are NOT sighted)")

    even_fail = [(N, j) for N in A_RANGE for j in range(N)
                 if (j + 1) % 2 == 0
                 and (blind_by_krylov(N, j) == 0) != is_automorphism(j + 1, 2 * (N + 1))]
    # Not a discovery: for even m, gcd(m, 2n) >= 2 always, so mu_m is NEVER an automorphism of
    # Q(zeta_2n) and the Galois criterion is simply silent there. What the count says is how many
    # of those silent seats are sighted anyway, i.e. how much A2 would over-claim without its fence.
    check("A3 THE FENCE: the Galois criterion is silent at EVEN multipliers, and misses seats there",
          len(even_fail) > 0,
          f"{len(even_fail)} even-multiplier seats are sighted with no automorphism, smallest "
          f"{even_fail[0]}; A2's odd restriction is a fence and not a convenience")

    print("\n  NEGATIVE CONTROLS on the same range (each must report disagreement)")
    for mod_of, label in [(lambda N: 2 * (N + 1), "2n"), (lambda N: N, "N"),
                          (lambda N: 2 * N, "2N")]:
        d = sum(1 for N in A_RANGE for j in range(N)
                if (blind_by_krylov(N, j) + 1) != aliasing_degree(j + 1, mod_of(N)))
        check(f"A4 reading the aliasing on modulus {label} breaks A1", d > 0, f"{d} disagreements")
    d = sum(1 for N in A_RANGE for j in range(1, N)
            if (blind_by_krylov(N, j) == 0) != is_automorphism(j, 2 * (N + 1)))
    check("A5 an off-by-one seat<->multiplier map breaks A2", d > 0, f"{d} disagreements")

    print()
    print("=" * 86)
    print("(b) WHICH SEATS KEEP THE COUNT WHEN THE CHAIN IS SPOILED")
    print("=" * 86)

    bad = []
    for N in range(5, 14, 2):
        c = N // 2
        if not (locus_is_identical(H_u, N, c) and locus_is_identical(H_delta, N, c)):
            bad.append(N)
    check("B1 odd N=5..13: the centre seat's resultant vanishes IDENTICALLY on BOTH axes", not bad)

    # B1 says only "blind >= 1 everywhere". The prose says the VALUE does not move.
    bad = []
    for N in range(5, 12, 2):
        c = N // 2
        for val in [sp.Integer(0), sp.Rational(1, 3), sp.Rational(7, 5), sp.Integer(1),
                    sp.Integer(2), sp.Rational(-1, 3)]:
            if blind_fence_free(H_u(N).subs(t, val), c) != (N - 1) // 2:
                bad.append((N, "u", str(val)))
            if blind_fence_free(H_delta(N).subs(t, val), c) != (N - 1) // 2:
                bad.append((N, "D", str(val)))
    check("B2 the centre seat's VALUE is (N-1)/2 at every tested coupling on both axes",
          not bad, "(not merely nonzero)")

    bad = []
    for N in range(5, 14):
        row = [blind_fence_free(H_u(N).subs(t, 1), j) for j in range(N)]
        expect = (N - 1) // 2 if N % 2 else (N - 2) // 2
        if any(v != expect for v in row):
            bad.append((N, row))
    check("B3 the ring end u=1 is that same constant at EVERY seat, both parities",
          not bad, "((N-1)/2 at odd N, (N-2)/2 at even N; a ring fixes every seat)")

    acc = [(N, j) for N in range(5, 12) for j in range(N)
           if not (N % 2 and j == N // 2) and blind_by_krylov(N, j) > 0]
    stuck = [(N, j) for (N, j) in acc if locus_is_identical(H_u, N, j)]
    check("B4 every ACCIDENTALLY blind seat has a FINITE locus, i.e. it moves",
          not stuck and bool(acc), f"{len(acc)} accidental seats, {len(stuck)} stuck")

    both = [(N, j) for N in range(5, 12, 2) for j in range(N)
            if j == N // 2 and blind_by_krylov(N, j) != (N - 1) // 2]
    check("B5 AT u=0 the forced seat's whole count is the forced one: blind(centre) == (N-1)/2",
          not both, "(given A1 this is integer arithmetic, not a second measurement; and it is a "
                    "u=0 statement only, since at u=1 a ring fixes EVERY seat)")

    print()
    print("=" * 86)
    print("(c) WHERE THE ACCIDENTAL BLINDNESS GOES")
    print("=" * 86)

    for N, j, poly, label in [(9, 1, t**5 - 4*t**3 + 3*t, "D^5 - 4D^3 + 3D"),
                              (9, 2, 2*t**2 - 1, "2D^2 - 1"),
                              (11, 1, t**7 - 6*t**5 + 10*t**3 - 4*t, "D^7 - 6D^5 + 10D^3 - 4D"),
                              (11, 2, 3*t**4 - 4*t**2, "3D^4 - 4D^2")]:
        check(f"C1 N={N} seat {j}: Delta-locus == real roots of F157's committed {label}",
              loci_equal(_resultant(H_delta, N, j), sp.Poly(poly, t)))

    bad, tested, empty = [], 0, 0
    for N in range(5, 12):
        for j in range(N):
            nn = abs(N - 1 - 2 * j)
            P = _resultant(H_delta, N, j)
            if nn < 1 or P is None:
                continue
            # Membership is decided by MINIMAL-POLYNOMIAL DIVISIBILITY, never by evaluating P
            # at a nested trigonometric expression: sp.simplify does not reduce
            # P(2*cos(3*pi/7)) to 0 and would report a true root as a counterexample. A
            # heuristic simplifier must not judge an exact question; that is the same defect
            # class as the is_real filter this file used to carry.
            sq = _squarefree(P)
            vals = {sp.sin((j + 1) * k * sp.pi / nn) / sp.sin(j * k * sp.pi / nn)
                    for k in range(1, nn) if (j * k) % nn != 0}
            minpolys = {sp.Poly(sp.minimal_polynomial(v, t), t) for v in vals}
            if (any(sp.rem(sq, mp) != 0 for mp in minpolys)
                    or sum(mp.count_roots() for mp in minpolys) != _n_real(sq)):
                bad.append((N, j))
            tested += 1
            if not vals:
                empty += 1
    check("C2 the Delta-locus is F157's committed Delta_k = sin((j+1)k*pi/Nn)/sin(jk*pi/Nn)",
          not bad, f"({tested} non-centre seats over N=5..11, {tested - empty} with a nonempty value set and {empty} the empty-locus case N_node | j; the continuation stays in comb country and "
                   f"only changes modulus)")

    bad, matched, idz = [], 0, 0
    for N in ODD:
        for j in range(N):
            Pu, Pd = _resultant(H_u, N, j), _resultant(H_delta, N, j)
            if Pu is None or Pd is None:
                if (Pu is None) != (Pd is None):
                    bad.append((N, j, "IDZERO mismatch"))
                else:
                    idz += 1
            elif _matches(Pu, Pd):
                matched += 1
            else:
                bad.append((N, j))
    content = sum(1 for N in ODD for j in range(N)
                  if _resultant(H_delta, N, j) is not None
                  and _n_real(_resultant(H_delta, N, j)) > 0)
    check(f"C3 odd N={ODD[0]}..{ODD[-1]}: u-locus == Delta-locus + the ring ends, all {ODD_SEATS} seats",
          not bad, f"({matched} non-IDZERO, {idz} identically-blind centres; the relation has "
                   f"CONTENT at the {content} seats whose Delta-locus is nonempty; {ODD_SEATS - content - idz} "
                   f"more read u-locus == the ring ends alone, and the {idz} forced centres are "
                   f"identically blind on BOTH axes and have no locus at all)")

    breaks, even_seats = [], 0
    for N in EVEN:
        for j in range(N):
            even_seats += 1
            Pu, Pd = _resultant(H_u, N, j), _resultant(H_delta, N, j)
            if Pu is not None and Pd is not None and not _matches(Pu, Pd):
                breaks.append((N, j))
    check("C4 THE FENCE: even N breaks, so C3's odd-N scope is measured, not decorative",
          len(breaks) > 0, f"{len(breaks)} of {even_seats} even-N seats: {breaks}")

    S = lambda N: sp.diag(*[(-1)**l for l in range(N)])
    ident = lambda Hf, N: sp.expand(S(N) * Hf(N) * S(N) + Hf(N).subs(t, -t)) == sp.zeros(N, N)
    check("C5 Sigma*H*Sigma == -H(-t) holds for the crack at ODD N and FAILS at every even N",
          all(ident(H_u, N) for N in ODD)
          and not any(ident(H_u, N) for N in EVEN),
          "(the wrap entry picks up (-1)^(N-1))")
    check("C5b the same identity holds for the anisotropy at EVERY N",
          all(ident(H_delta, N) for N in range(ODD[0], EVEN[-1] + 1)),
          "(which is why the Delta-locus is negation-closed at both parities)")

    print()
    print("  WHAT DISTINGUISHES THE CRACK?  Two earlier versions of this block got the question")
    print("  wrong in opposite directions. The first asked whether a control satisfies the CRACK's")
    print("  relation, locus == Delta-locus + {+-1}; the ring ends belong to the u axis and a")
    print("  diagonal has none, so every control failed before its locus was looked at. The second")
    print("  scored the controls on EQUALITY with the Delta-locus while scoring the crack on the")
    print("  relation with the ring ends: two predicates, two denominators. Asked the same way of")
    print("  everything, the predicate that means 'carries the whole Delta-locus' is CONTAINMENT.")
    all_seats = [(N, j) for N in ODD for j in range(N)]
    nonempty = [(N, j) for (N, j) in all_seats
                if _resultant(H_delta, N, j) not in (None, 'SKIP')
                and _n_real(_resultant(H_delta, N, j)) > 0]

    PERTURBATIONS = [(H_u, "the crack (the wrap bond)"),
                     (H_midbond, "an interior bond, reflection-ASYMMETRIC"),
                     (H_symbond(1), "a reflection-SYMMETRIC interior bond pair"),
                     (H_symdiag(1), "a reflection-SYMMETRIC interior diagonal pair (1, N-2)"),
                     (H_symdiag(2), "a reflection-SYMMETRIC interior diagonal pair (2, N-3)")]

    print(f"    READ, all five perturbations on the same {len(nonempty)} seats and the same two")
    print("    questions.  EQUAL: locus == Delta-locus.  CONTAINS: Delta-locus is a subset of it.")
    print("    A row shows a smaller denominator where its builder is not defined at some N (the")
    print("    (2, N-3) pair degenerates onto the reflection-fixed centre at N = 5).")
    contains_rate, equal_rate, denom = {}, {}, {}
    for Hf, label in PERTURBATIONS:
        eq = ct = seats = 0
        for (N, j) in nonempty:
            P, Pd = _resultant(Hf, N, j), _resultant(H_delta, N, j)
            if P == 'SKIP':
                continue               # this builder is not defined at this N; not its denominator
            seats += 1
            if P is None:                      # identically blind: contains everything
                ct += 1
                continue
            if loci_equal(P, Pd):
                eq += 1
            if locus_contains(P, Pd):
                ct += 1
        equal_rate[label], contains_rate[label], denom[label] = eq, ct, seats
        print(f"      {label}:  EQUAL {eq}/{seats}   CONTAINS {ct}/{seats}")

    crack = "the crack (the wrap bond)"
    # The crack CANNOT win on equality: it always carries the ring ends as well, so on that
    # predicate it scores worse than a diagonal. Containment is the question the sentence
    # "carries the whole Delta-locus" actually asks, and it is asked of all five identically.
    # C6 is C3 restated, not a second measurement: C3 already asserts equality with the Delta-locus
    # plus the ring ends over these seats, and containment follows from that. It is kept because the
    # controls are scored on containment and the crack's row of the table has to be the same column.
    check("C6 (C3 restated) the crack CONTAINS the whole Delta-locus at every seat with one",
          contains_rate[crack] == denom[crack], f"{contains_rate[crack]}/{denom[crack]}")
    others = {k: v for k, v in contains_rate.items() if k != crack}
    check("C7 none of the four interior perturbations tried does the same",
          all(v < denom[k] for k, v in others.items()),
          "best of them " + max(f"{v}/{denom[k]}" for k, v in others.items())
          + "; 'the end pair is special' is REFUTED for this locus. How wide the remaining gap "
          + "is, this line does NOT say: part (d) splits the column and 32 of that 44 are seats "
          + "blind at every knob value, so the genuine comparison is 48 against 12")

    # The break list and the negation-closure failure list, as SETS. Note which half is content:
    # C5b makes the Delta-locus negation-closed at every N and {+-1} is closed, so the right-hand
    # side of C3 is closed always and a non-closed u-locus CANNOT equal it. That direction is a
    # theorem and this check cannot report it false. The converse is what is measured.
    breaks, nonclosed, holds = set(), set(), 0
    for N in EVEN:
        for j in range(N):
            Pu, Pd = _resultant(H_u, N, j), _resultant(H_delta, N, j)
            if Pu is None or Pd is None:
                continue
            if _matches(Pu, Pd):
                holds += 1
            else:
                breaks.add((N, j))
            if not loci_equal(Pu, sp.Poly(Pu.as_expr().subs(t, -t), t)):
                nonclosed.add((N, j))
    check(f"C8 even N={EVEN[0]}..{EVEN[-1]}: every break IS a non-closure (the converse is forced)",
          breaks == nonclosed and bool(breaks),
          f"both sets are the same {len(breaks)} seats")
    print(f"    READ: the relation still HOLDS at {holds} of the {holds + len(breaks)} even-N seats.")
    print("    Parity decides whether the staggering identity FORCES it, not whether it holds.")

    # Does the u-locus, WITH THE RING ENDS REMOVED, take exactly one member of each +- pair of
    # the Delta-locus? An earlier version of the page asserted this at every break and gated
    # nothing. It is checked here, and it is not universal.
    def strip_ring_ends(P):
        for e in (t - 1, t + 1):
            while sp.rem(P, sp.Poly(e, t)) == 0:
                P = sp.Poly(sp.quo(P, sp.Poly(e, t)), t)
        return P

    good, bad = [], []
    for (N, j) in sorted(breaks):
        Pd = _squarefree(_resultant(H_delta, N, j))
        Pu = strip_ring_ends(_squarefree(_resultant(H_u, N, j)))
        ok = True
        for r in {abs(v) for v in Pd.real_roots() if v != 0}:
            mp_p = sp.Poly(sp.minimal_polynomial(r, t), t)
            mp_m = sp.Poly(sp.minimal_polynomial(-r, t), t)
            got = (Pu.degree() >= 1 and sp.rem(Pu, mp_p) == 0) +                   (Pu.degree() >= 1 and sp.rem(Pu, mp_m) == 0 and mp_m != mp_p)
            if got != 1:
                ok = False
        (good if ok else bad).append((N, j))
    # What is GATED here is non-universality. The 20-of-22 rate beside it is a READ.
    check("C9 THE SECTOR-SPLIT SHAPE IS NOT UNIVERSAL at a break",
          bool(bad),
          f"u minus the ring ends takes exactly one of each +- pair at {len(good)} of "
          f"{len(breaks)} breaks and NOT at {bad}, where the Delta-locus contains +-1 itself "
          f"and u carries those two only as ring ends")

    # ------------------------------------- part (d): blindness at every knob value, one law
    print()
    print("  " + "-" * 84)
    print("  (d) WHAT THE CONTAINS COLUMN WAS COUNTING, AND THE ONE LAW UNDER IT.")
    print("  At some seats the perturbation leaves the resultant identically zero in its own knob:")
    print("  the seat is blind at EVERY KNOB VALUE, its locus is the whole line, and it contains")
    print("  every locus for free. That branch is not about the Delta-locus at all, and it closes")
    print("  to one criterion for both kinds, the OVERLAP v^T V v of a BLIND EIGENVECTOR of")
    print("  the UNPERTURBED chain with the perturbation direction: a sum of squares for a")
    print("  diagonal, so a node on every perturbed site; a product for a bond, so a node at")
    print("  either end. Corollary D of docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md.")

    # The split of block (c)'s own table. A READ; the gate on it is D0.
    print()
    split, fail_seats = {}, {}
    for Hf, label in PERTURBATIONS:
        every_knob = genuine = 0
        fails = []
        for (N, j) in nonempty:
            P, Pd = _resultant(Hf, N, j), _resultant(H_delta, N, j)
            if P == 'SKIP':
                continue
            if P is None:
                every_knob += 1
            elif locus_contains(P, Pd):
                genuine += 1
            else:
                fails.append((N, j))
        split[label] = (every_knob, genuine, len(fails))
        fail_seats[label] = fails
        print(f"      {label}:  CONTAINS {every_knob + genuine}/{every_knob + genuine + len(fails)}"
              f"  =  {every_knob} blind at every knob + {genuine} genuine,  {len(fails)} fail")
    best = "a reflection-SYMMETRIC interior diagonal pair (1, N-2)"
    check("D0 the best interior row's CONTAINS column is mostly seats blind at EVERY knob value",
          split[best][0] > split[best][1],
          f"{split[best][0]} of its {split[best][0] + split[best][1]}; the {split[best][2]} it "
          f"fails outright are {fail_seats[best]}, which is a READ printed here so the page's "
          f"list of them has a producer, and it is reflection-closed: "
          f"{sorted((N, N - 1 - j) for (N, j) in fail_seats[best]) == sorted(fail_seats[best])}")

    # ------------------------------------------------ the law against the table it explains
    # Block (c) runs on ODD N = 5..17; the sweeps below stop earlier. Without this join nothing
    # would connect D0's split to the criterion, which is the whole motivation of part (d).
    # TWO populations, deliberately. The table's own seats are the 48 with a nonempty Delta-locus,
    # and that list EXCLUDES every centre seat, whose Delta-resultant vanishes identically. The
    # centre seat is exactly where D3 shows the criterion's converse failing on this same (1, N-2)
    # row, so certifying only the table's seats would take the green from the counterexamples
    # being outside the loop. The full seat list is therefore read as well, and its disagreements
    # are gated to sit AT centre seats, which is the one direction and not a set equality.
    join_bad, join_cells, join_pos = [], 0, 0
    full_bad, full_cells = [], 0
    for k in (1, 2):
        for N in ODD:
            if k >= N - 1 - k:
                continue
            for j in range(N):
                meas = _resultant(H_symdiag(k), N, j) is None
                pred = overlap_is_zero(N, j, ('diag', {k, N - 1 - k}))
                full_cells += 1
                if meas != pred:
                    full_bad.append((N, j, k))
                if (N, j) not in nonempty:
                    continue
                join_cells += 1
                join_pos += pred
                if meas != pred:
                    join_bad.append((N, j, k))
    check("D8 THE JOIN: the criterion decides block (c)'s own two diagonal rows, seat by seat, "
          "on the 48 seats that table runs on, over ODD N = 5..17",
          not join_bad and 0 < join_pos < join_cells,
          f"{join_cells} cells, {join_pos} blind at every knob and {join_cells - join_pos} not, "
          f"0 disagreements")
    check("D8b AND THE SEATS THAT LIST DROPS ARE NOT A FREE PASS: read over ALL seats of the same "
          "N, the criterion DOES disagree, and exactly at the centre seats D3 names",
          bool(full_bad) and all(N - 1 - 2 * j == 0 for (N, j, _) in full_bad),
          f"{len(full_bad)} disagreements of {full_cells} cells, every one at a centre seat: "
          f"{full_bad}; the 48-seat list excludes them because their Delta-resultant vanishes "
          f"identically, so D8's green is scoped to that list and not to the N range")

    # ------------------------------------------------------------- the law, both halves, diagonal
    D_RANGE = range(5, 15)
    onesided, straddling = [], []
    for N in D_RANGE:
        for j in range(1, N - 1):
            for size in (1, 2, 3):
                for S in _interior_sets(N, j, size):
                    (onesided if (max(S) < j or min(S) > j) else straddling).append(
                        (N, j, frozenset(S)))
    one_truth = {k: blind_at_every_knob(k[0], k[1], set(k[2])) for k in onesided}
    one_law = {k: overlap_is_zero(k[0], k[1], ('diag', set(k[2]))) for k in onesided}
    fwd = [k for k in onesided if one_law[k] and not one_truth[k]]
    cnv = [k for k in onesided if one_truth[k] and not one_law[k]]
    n_law, n_true = sum(one_law.values()), sum(one_truth.values())
    print()
    print(f"    {len(onesided)} one-sided triples (N, seat, diagonal site set) over "
          f"N = {D_RANGE.start}..{D_RANGE.stop - 1}, EVERY seat, the")
    print(f"    centre seat included, site sets of size 1, 2 and 3. {n_true} are blind at every")
    print(f"    knob value and {len(onesided) - n_true} are not.")
    check("D1 THE FORWARD HALF (a theorem): a blind eigenvector with zero overlap leaves the "
          "seat blind at every knob value",
          not fwd and n_law > 0,
          f"0 exceptions in the {n_law} triples where the LAW fires; the two counts below are "
          f"equal BECAUSE both directions hold here, and each check guards the population that "
          f"could refute IT, so killing either route reddens exactly one of D1 and D2")
    check("D2 THE CONVERSE when the perturbed sites lie on ONE side of the seat (a theorem)",
          not cnv and n_true > 0,
          f"0 exceptions in the {n_true} triples that ARE blind at every knob value")

    # Straddling is where the converse genuinely fails, and the centre seat is where it fails
    # first: there the two halves have EQUAL length, so a mirrored knob pattern makes them the
    # same matrix. This is what the earlier version of part (d) hid by dropping the centre seat.
    str_fail = [k for k in straddling
                if blind_at_every_knob(k[0], k[1], set(k[2]))
                and not overlap_is_zero(k[0], k[1], ('diag', set(k[2])))]
    small2 = sorted((N, j, tuple(sorted(S))) for (N, j, S) in str_fail if len(S) == 2)
    L5, R5 = _halves_diag(5, 2, {1, 3}, t)
    check("D3 THE CONVERSE IS FALSE WHEN S STRADDLES THE SEAT, and it fails first at the CENTRE "
          "seat, at |S| = 2, on the road page's own (1, N-2) pair",
          bool(small2) and all(N - 1 - 2 * j == 0 for (N, j, _) in small2)
          and sp.Poly(_cp_tridiag(L5, x), x) == sp.Poly(_cp_tridiag(R5, x), x),
          f"{len(small2)} failures at |S| = 2 over N = {D_RANGE.start}..{D_RANGE.stop - 1}, "
          f"EVERY one at the centre seat; the smallest is N = 5 seat 2 with S = (1, 3), where "
          f"the two halves are literally the same polynomial")
    away = sorted((N, j, tuple(sorted(S))) for (N, j, S) in str_fail if N - 1 - 2 * j != 0)
    print(f"    READ, not a census: away from the centre seat the smallest failures need "
          f"|S| = 3; there are {len(away)} of them in this range,")
    print(f"    {sorted({(N, j) for (N, j, _) in away})}, and the count GROWS with N, so no "
          f"total here is a closed number.")

    # Are the straddling failures the criterion in disguise? They are not: their shared factor
    # still carries the knob. NOTE this cannot go red on a data input, since a knob-free shared
    # level would give a blind eigenvector with zero overlap and the triple would not be in the list.
    moving = []
    for (N, j, S) in away:
        dl, dr = _halves_diag(N, j, set(S), t)
        gshared = sp.gcd(sp.Poly(_cp_tridiag(dl, x), x), sp.Poly(_cp_tridiag(dr, x), x))
        facs = sp.factor_list(gshared.as_expr())[1]
        moving.append(all(t in sp.Poly(f, x).as_expr().free_symbols for f, _ in facs))
    check("D3b FORCED, and computed rather than asserted: every straddling failure's shared "
          "factor carries the knob in EVERY one of its factors",
          bool(moving) and all(moving),
          f"{sum(moving)} of {len(moving)}; a knob-free factor would BE a blind eigenvector with zero "
          f"overlap, so this check cannot go red on a data input and is here as the arithmetic "
          f"face of that argument")

    # Mutations of the modulus, scored against the one-sided truth table.
    MUT = {"the node modulus |N-1-2j|": lambda N, j: abs(N - 1 - 2 * j),
           "gcd(2j+1, N)": lambda N, j: math.gcd(2 * j + 1, N),
           "gcd(j+1, N)": lambda N, j: math.gcd(j + 1, N),
           "gcd(j, N+1)": lambda N, j: math.gcd(j, N + 1),
           "N+1": lambda N, j: N + 1}
    mut_scores = {name: sum(1 for k in onesided
                            if overlap_is_zero(k[0], k[1], ('diag', set(k[2])), modulus=mod)
                            != one_truth[k])
                  for name, mod in MUT.items()}
    check("D4 five wrong moduli in the same law all redden",
          all(v > 0 for v in mut_scores.values()),
          ", ".join(f"{k}: {v}" for k, v in mut_scores.items())
          + f" disagreements of {len(onesided)}")

    shiftable = [(N, j, S) for (N, j, S) in onesided
                 if max(S) + 1 <= N - 2 and j not in {m + 1 for m in S}]
    shifted = sum(1 for (N, j, S) in shiftable
                  if overlap_is_zero(N, j, ('diag', {m + 1 for m in S})) != one_truth[(N, j, S)])
    check("D5 the law reads the SITES too: shifting every perturbed site by one reddens it",
          shifted > 0,
          f"{shifted} disagreements of the {len(shiftable)} triples whose sites can all be "
          f"shifted without hitting the seat or the chain end")

    # ------------------------------------------------------------------------ the bond face
    # An earlier version of this block claimed a bond could not take this route at all, because
    # a blind eigenvector would need two consecutive zeros and Lemma J1 forbids that. That transplanted
    # the DIAGONAL's condition. The bond's overlap is a PRODUCT, so one node suffices, and J1
    # never applied. Read per cell, not as a total. D6b shifts the knob to the neighbouring bond
    # and compares against the UNSHIFTED measurement: the criterion then disagrees at many cells,
    # which is what pins the bond INDEX rather than only the bond count.
    bond_cells, bond_bad, bond_pos = [], [], 0
    for N in D_RANGE:
        for j in range(1, N - 1):
            for b in range(0, N - 1):
                if b in (j - 1, j):
                    continue                  # incident on the seat; struck away with it
                meas = bond_blind_at_every_knob(N, j, b)
                pred = overlap_is_zero(N, j, ('bond', b))
                bond_cells.append((N, j, b))
                bond_pos += meas
                if meas != pred:
                    bond_bad.append((N, j, b))
    shiftable_bonds = [(N, j, b) for (N, j, b) in bond_cells if b > 0]
    shift_bad = sum(1 for (N, j, b) in shiftable_bonds
                    if overlap_is_zero(N, j, ('bond', b - 1)) != bond_blind_at_every_knob(N, j, b))
    shift_fires = sum(1 for (N, j, b) in shiftable_bonds if overlap_is_zero(N, j, ('bond', b - 1)))
    true_fires = sum(1 for (N, j, b) in shiftable_bonds if overlap_is_zero(N, j, ('bond', b)))
    check("D6 THE BOND IS THE SAME LAW, cell by cell: a node at EITHER end of the knob-bearing "
          "bond, the product where the diagonal has a sum of squares",
          not bond_bad and 0 < bond_pos < len(bond_cells),
          f"{len(bond_cells)} cells over N = {D_RANGE.start}..{D_RANGE.stop - 1}, every seat, "
          f"{bond_pos} blind at every knob value and {len(bond_cells) - bond_pos} not, "
          f"0 disagreements")
    check("D6b the criterion reads the bond INDEX and not only how many bonds qualify: moving "
          "the knob to the neighbouring bond disagrees with the unshifted measurement",
          shift_bad > 0,
          f"{shift_bad} disagreements of the {len(shiftable_bonds)} cells with a bond to its left; "
          f"the shifted criterion also fires a different NUMBER of times ({shift_fires} against "
          f"{true_fires}), so the off-by-one is not count-preserving either")

    # The reflection buys nothing on this branch and everything on the other.
    pair_vs_single = pair_seats = 0
    for N in ODD:
        for j in range(1, N - 1):
            for m in range(1, (N - 1) // 2):
                if m == j or N - 1 - m == j:
                    continue
                pair_seats += 1
                if (overlap_is_zero(N, j, ('diag', {m, N - 1 - m}))
                        != overlap_is_zero(N, j, ('diag', {m}))):
                    pair_vs_single += 1
    eq_pair = eq_single = live_pair = live_single = 0
    for (N, j) in nonempty:
        Pd = _resultant(H_delta, N, j)
        P2 = _resultant(H_symdiag(1), N, j)
        if P2 not in (None, 'SKIP'):
            live_pair += 1
            eq_pair += loci_equal(P2, Pd)
        P1 = _resultant(H_onesite(1), N, j)
        if P1 not in (None, 'SKIP'):
            live_single += 1
            eq_single += loci_equal(P1, Pd)
    check("D7 the reflection buys nothing on this branch and IS what buys equality",
          eq_pair > eq_single and pair_vs_single == 0,
          f"on EQUALITY with the Delta-locus the pair scores {eq_pair} of {live_pair} live seats "
          f"and the single site {eq_single} of {live_single}; the two agree on blindness at every "
          f"knob value at all {pair_seats} pairs tried BECAUSE gcd(j+1, N+1) divides N-1-2j, not "
          f"because it was measured")

    print()
    print("=" * 86)
    print("VERDICT:", "ALL GREEN" if not _fails else f"{len(_fails)} FAILED: {_fails}")
    print("=" * 86)
    return 1 if _fails else 0

if __name__ == "__main__":
    raise SystemExit(main())
