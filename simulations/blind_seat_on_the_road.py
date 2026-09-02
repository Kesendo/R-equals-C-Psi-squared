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

Every locus is SOLVED (the real roots of a resultant), never sampled: gate B1 evaluates the
blindness ON the solved locus, because a sampled interior point misses it by construction.
Gates C1a and C1b take their expectations as LITERALS from F157's committed table.

Companion page: experiments/THE_BLIND_SEAT_ON_THE_ROAD.md
"""
from fractions import Fraction
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
          + "; the margin is narrow, and 'the end pair is special' is REFUTED for this locus")

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

    print()
    print("=" * 86)
    print("VERDICT:", "ALL GREEN" if not _fails else f"{len(_fails)} FAILED: {_fails}")
    print("=" * 86)
    return 1 if _fails else 0

if __name__ == "__main__":
    raise SystemExit(main())
