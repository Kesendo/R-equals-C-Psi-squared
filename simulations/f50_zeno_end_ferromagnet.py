"""F50's Zeno end: the population generator is the ferromagnet, and on every tree its gap is lambda_1(T).

Gate for the Zeno-end paragraphs of docs/proofs/PROOF_WEIGHT1_DEGENERACY.md, section "The count as a plane
crossing": "The Zeno end and the dispersion", Theorem E with its proofs ("By the lift the token gap is at most" and "The slow levels are the eigenvalues"), "Which filling lasts longest", and the Zeno half of "The Zeno end and the Hamiltonian end side by side" (its overlaps are read by
f50_isotropy_and_energy_mode.py, rows E3 and E4).

Objects (Pauli book, H = sum_b J_b (X_a X_b + Y_a Y_b + Delta Z_a Z_b), uniform Z-dephasing, x = J/gamma):
in the diagonal joint-popcount block (p, p) the populations P see, as gamma -> infinity, the generator
-(A^2)_PP / (4 gamma) with A = H (x) I - I (x) H^T; the exclusion graph of block p has the popcount-p
configurations as vertices and an edge of weight c_b for every bond b whose two sites disagree.
Site l is bit N-1-l throughout, as in f50_plane_crossing_count.py (whose helpers are copied, not imported:
importing a gate runs it).

Rows (an exact row compares integers, rationals or Gaussian rationals with ==; a float row states its error model and prints
the ratio to it; a reading prints and gates nothing):
  Z1  The exclusion Laplacian is the ferromagnet: 2 Lap(block p) == block p of sum_b c_b (I - X X - Y Y - Z Z)
      built from Kronecker products, integer weights, on the chain, ring, star, complete graph and random graphs,
      every p.  Control: the same operator built without its Y Y term goes red.
  Z2  (A^2)_PP is 8 times that Laplacian with the bond weights |t_b|^2/4, t_b the bond's hop amplitude, so the Zeno
      generator -(A^2)_PP/(4 gamma) is -(2/gamma) times it:
      (A^2)_PP == 8 Lap(J_b^2) exactly at Delta = 0, 1/2, 1, 2, -1 on chains with integer bond profiles and
      Delta = 3/2 on a random graph; unchanged by terms diagonal in the configurations (the longitudinal
      fields (2, -1, 3, 0, -2) and next-nearest ZZ couplings); with a Dzyaloshinskii-Moriya term D_b the weight is J_b^2 + D_b^2.
      The rows in Delta, fields and next-nearest ZZ are structural (the diagonal of H cancels in A_QP, so they
      check the construction); the Liouvillian rows of Z6 test the blindness.  The assembled generator commutes
      with S^- exactly on three of these rows.  Controls: the weights J_b in
      place of J_b^2 miss the first identity; a correlated hop, its amplitude doubled when a third site is
      occupied, breaks the commutation (SU(2) gone).
  Z2b Under a rational profile of rates gamma_l (real H: Delta = 1/2 with fields and a next-nearest ZZ on a chain;
      Delta = 3/2 with fields on a graph with a triangle, since on a graph without triangles the next term vanishes for every
      H: its diagonal parts cancel pair by pair and the rest needs three configurations pairwise one hop apart) minus the leading term of the Schur complement is exact: A_PQ Gamma_Q^-1 A_QP == Lap with the
      weights 4 J_b^2/(gamma_i + gamma_j), Gamma the cells' rates 2 sum gamma over the disagreeing sites, and the
      next term A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP == 0 exactly (the transpose makes every odd order vanish for a
      real H).  Controls: the uniform-rate weights 2 J_b^2/gamma_mean miss the first identity; a
      Dzyaloshinskii-Moriya term on a bond of the triangle (a flux, H no longer real) makes the next term nonzero,
      computed exactly in Gaussian rationals; that term is purely imaginary and antisymmetric there, so i M3 moves
      no real part (row Z6c reads what it does to a degenerate level); and on a chain and a 4-ring that carry a
      Dzyaloshinskii-Moriya term and no triangle the next term vanishes, exactly, H not real.
  Z3  For the uniform Heisenberg graph Lap(block p) == ((#bonds) I - H_p)/2 exactly: the Hamiltonian read from
      the top of its spectrum.  Control: Delta = 2 goes red.
  Z4  The lift is SU(2)'s lowering: (S^-)^(p-1) sum_i f(i)|i> == (p-1)! sum_S F(S)|S>, F(S) = sum_{i in S} f(i),
      and [sum_b c_b (1 - SWAP_b), S^-] == 0, and Lap F == lift of (L_site f), all exact on integer f.
      Controls: S^- with one site missing breaks the lift; the anisotropic sum_b c_b (I - XX - YY - 2 ZZ) breaks
      the commutation with S^-; a site Laplacian with one diagonal entry raised breaks Lap F == lift(L_site f).
  Z5  The theorem, without an eigensolver: for every tree (paths N = 2..9, stars N = 3..9, random integer-weighted
      trees N = 3..9, two chains with random integer bond profiles) and every block 1 <= p <= N-1, the number of
      eigenvalues below a rational q is read exactly (Jacobi's rule on Bareiss minors of b*Lap - a*I over the
      integers; an eigensolver only places the rationals, which the counts then verify).  Below q_lo < lambda_1(T)
      only the zero mode; below q_hi in (lambda_1, lambda_1 + 1e-6 max(1, lambda_1)) exactly as
      many as the tree's Laplacian has; below q* placed 1 to 2 ppm under max_leaf lambda_1(T - v), again exactly the tree's count
      (the new multiplets lie at or above lambda_1(T - v)).  Each bracket is verified exactly on the site Laplacians.
      On the stars with N >= 4 max_leaf lambda_1(T - v) = lambda_1(T), so q* sits below lambda_1 there and adds nothing to the
      q_lo bracket; the row's label says so.
      Controls, through the same row, each counting only with exact brackets and a definite count that differs:
      the blocks of the N = 7 chain with its bond (3, 4) doubled, and with it halved, counted against the
      unchanged chain's brackets, miss (through q_hi and through q_lo); and, on what the theorem adds, the term
      -c (1 - SWAP_01)(1 - SWAP_23), which vanishes on the ferromagnetic multiplet and the one-magnon band and
      moves only multiplets of lower spin, added to every block of the N = 6 chain: at c = 1 a new multiplet
      falls below lambda_1 (every bracket misses), at bond weights 100 and c = 95 one lands between lambda_1 and
      lambda_1(T - v) (q* alone misses).  A cross-check compares the exact counter with a float count far from
      the eigenvalues.
  Z6  The Zeno limit of the full Liouvillian is blind to Delta: on the chains N = 4, 5 the C(N,p) slowest
      eigenvalues of B(x) are -2 x^2 ell_j + O(x^4) with ell_j the Laplacian's spectrum at Delta = 0, 1/2, 1, 2,
      and at Delta = 1 with the fields (2, -1, 3, 0) and a next-nearest ZZ on N = 4.  The error model is the
      expansion in even powers (for a real H, as on every row here, every odd order vanishes by the transpose of
      Lemma 1): deviation/x^4 = a + b x^2 + ...
      A correction in x^k makes the successive differences over x = 0.005, 0.01, 0.02 stand in the ratio 2^k;
      the law is k = 2, and the gate decides it against k = 1 and k = 3 at the half-integer powers,
      2^(3/2) < ratio < 2^(5/2) (measured within 3 % of 4, the float route's rounding included).  A leading generator off by a relative
      delta adds 2 delta ell/x^2, k = -2.  Control on every row: ell (1 + 1e-7), fed through the same
      comparison, is not read as the law.  Under a random profile of rates gamma_l and couplings J_b at
      Delta = 7/10 (chain N = 5, every block) the slow eigenvalues are minus the spectrum of the Laplacian with
      weights 4 J_b^2/(gamma_i + gamma_j); the deviation times the cube of the overall rate scale s, over
      s = 50, 100, 200, follows s^k with the law k = -2, decided the same way, 2^(-5/2) < ratio < 2^(-3/2).
      Control: the weights 2 J_b^2/gamma_mean of a uniform profile, fed through the same comparison, give
      k = 2, a wrong leading term, in every block.
  Z6c What the next term does where it does not vanish: on the triangle with the same Dzyaloshinskii-Moriya D
      on every bond along it (a flux 3 arctan(D/J) through it), block p = 1, the leading pair is degenerate and
      the gamma^-2 term splits it into a complex pair, |Im lambda| falling as gamma^-2 (ratio 2^2 per doubling of
      gamma, decided at the half-integer powers), the real parts of the pair printed; on the ring N = 5 with the
      same D, no triangle, |Im lambda| falls as gamma^-4 (ratio 2^4); each row is read against the other's law,
      which it must not meet, as its control.  The exact counterpart without a flux, the next
      term zero, is row Z2b's.
  Z7  Reading, not gated: rings N = 4..9, the complete graph and random graphs (where no leaf removal is available)
      under the same exact counter, each bracket's lower end checked to have only the zero mode below it: whether
      any block has an eigenvalue below lambda_1, and the multiplicity at lambda_1 per block (on the N = 4 ring a
      singlet joins the two lifted modes at lambda_1 = 2 in block p = 2).  The equality holds there by the route
      the proof recognizes afterwards (the exclusion process a quotient of the interchange process, whose gap
      Caputo, Liggett and Richthammer proved equal to lambda_1); this row reads it exactly.
  Z8  Reading, not gated: the chain prefactors that experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC.md and
      hypotheses/F1_DISSIPATION_GAP_PATTERN.md quote, in their book H = (J/4) sum sigma.sigma, gamma (Z rho Z - rho):
      the slowest rate over the diagonal blocks (every gap printed lies below 2 gamma, so by Theorem A it is a real
      mode of a diagonal block), as <n_XY> N^2/Q^2 and gap N^2/(gamma Q^2), for the XXX and the XX chain at Q = 2,
      N = 4, 5, 6, the N = 5 sweep Q = 0.5 to 2.5 and the Q = 0.5 mean over N = 4, 5, 6, beside the Zeno values
      N^2 (1 - cos(pi/N))/8 and /4.

  Z9  Beyond the leading order (Theorem E (d)): the slow levels are those of -M1 + i M3 + M4 + lambda K up to O(s^-4),
      M1 = A_PQ G^-1 A_QP, K = A_PQ G^-2 A_QP, M4 = A_PQ G^-1 A_QQ G^-1 A_QQ G^-1 A_QP.  The population ladder E'
      (add a particle anywhere, (E'v)(B) = sum over l in B of v(B - l)) commutes with M1, K and M4 of the hops alone
      in every block, exact over the rationals, on a chain with a bond and a rate profile, the chain N = 6, the
      rings N = 5, 6, the stars N = 4..7 (one with a rate profile) and two random trees; its M1 is the leading term,
      Lap(4 J_b^2/(gamma_i + gamma_j)), on the profile chain.  On the chain this is the
      Jordan-Wigner ladder of PROOF_CODIM1_BY_ADDITIVITY at every order; on the rings N >= 5 the wrapped string
      enters only after N hops; on the stars and trees the identity is read here and not proved.  Controls: on the
      ring N = 4 (a path of four hops winds around it) and with the ZZ term present (Heisenberg chain N = 5) E'
      commutes with M1 and K and not with M4; under a hop doubled when a third site is occupied (SU(2) gone) with
      none of the three.  The graphs with squares K2,3 and K2,4, hops only: E' commutes with M1 and K and not with M4
      on both (Z13 reads K2,3's fillings split at x^4 and K2,4's tied there).
  Z10 The diagonal terms enter M4 only as the detuning form D_V, the exclusion graph's Laplacian with the weight
      2 |t_e|^2 dV_e^2 / G_e^3 on each hop e, dV_e the change of the diagonal of H across it: M4(H) - M4(hops of H)
      == D_V exactly in every block, on a chain with Delta = 1/2, fields, a next-nearest ZZ and a rate profile, the
      ring N = 5, the star N = 5 with fields and a random tree N = 6.  Control: on the triangle graph the cross terms
      of hop and diagonal survive.
  Z11 The detuning weights W_b(p), read off F^T D_V F with F the lift (a bond sum by construction, every hop crossing
      one bond): on the uniform XXZ chain, ring and star (J = gamma = 1) W_b(p) is 4 Delta^2 C(N-4, p-2) on a bond with two outer
      neighbours (a hop there is detuned, by 4 Delta, exactly when they disagree), Delta^2 C(N-2, p-1)/2 on an end bond
      of the chain (detuned by 2 Delta always) and Delta^2 (N - 2p)^2 C(N-2, p-1)/2 on every bond of the star
      (detuned by 2 Delta |N - 2p|, the other arms' imbalance), exact for the chains and stars N = 4..7 and the rings
      N = 5..7, Delta = 1, 1/2, every block.  Controls through the same comparison: the form with every bulk hop
      detuned misses the chain N = 6 in every block, and the chain N = 6 with ZZ couplings (+1, +1, -1, -1, +1), whose
      bulk hops are detuned when their outer neighbours agree, misses the uniform form in blocks 1, 3 and 5 and meets it
      in blocks 2 and 4, where the two counts agree (2 C(2, p-1) + 2 C(2, p-3) against 4 C(2, p-2)).
  Z11c The shape of the split on a chain with nearest-neighbour ZZ couplings K_b of any sign and no fields (J = gamma
      = 1): W_b(p)/C(N-2, p-1) == ((K1 - K2)^2 (1 - phi(p)) + (K1 + K2)^2 phi(p))/2 on a bulk bond, K1 and K2 the
      couplings beside it, phi(p) = 2 (p-1)(N-p-1)/((N-2)(N-3)), and K^2/2 on an end bond, K its neighbouring
      coupling, exact in every block on four chains N = 6, 7 (one uniform, one with fractional couplings), so that
      only 2 K1 K2 phi(p) depends on the filling.  Control: a next-nearest ZZ term leaves that form.  At any range,
      ZZ couplings K_kl and no fields: W_b(p)/C(N-2, p-1) == ((sum c)^2 (1 - 2 phi(p)) + 2 phi(p) sum c^2)/2 on every
      bond (i, i+1), c_k = K_ik - K_(i+1)k over the other sites, which hold p - 1 particles in every arrangement
      alike, exact on the chain N = 6 with the next-nearest term (0, 2) and on two chains N = 7 with next-nearest
      and longer-range terms.  Control: a field on one site of the Heisenberg chain N = 6 misses it.
  Z11b The displayed delta_p of PROOF_WEIGHT1_DEGENERACY ("Which filling lasts longest"), the chain's in sines, the
      ring's with its own lambda_1 = 4 sin^2(pi/N) and the star's, against the quotient on the closed-form weights
      (which Z11 checks against D_V up to N = 7) with the site Laplacian's lowest nonconstant mode in closed form,
      N = 4..12 (rings from 5), every block, at 30 and at 60 digits: the residual falls with the working precision.
      Control: the ring's form with the chain's lambda_1 misses.
  Z12 The filling split at order x^4: each block's slowest level less block 1's, over x^4, against
      the difference of the detuning form's quotients on the lifted wave, on the chain N = 5 at Delta = 1, the chain
      N = 6 at Delta = 1 and 1/2, the ring N = 6, the stars N = 5, 6, and the chain N = 6 with ZZ couplings
      (1, 2, -1, -3, 1) and the uniform chain N = 6 with a next-nearest coupling 3/2, whose quotients use the exact
      weights; and delta_p itself, each block's slowest level less the same block's without the ZZ term, over x^4,
      against the quotient, on the chain N = 5 at Delta = 1, the chain N = 6 at Delta = 1/2, the ring N = 6 and the
      star N = 5 (the hops' own correction, the same with and without the ZZ term, cancels); each deviation's law is
      x^2 (the next even order), decided at the half-integer powers over x = 0.01, 0.02, 0.04 (ratio 4 per doubling),
      the prediction scaled by 1 + 1e-2 fed through the same comparison as the control, read outside the window at
      both doublings.
  Z13 Reading, not gated: the block of the slowest non-stationary rate at x = 0.02 (blocks within 64 eps times the
      largest block norm of it tie) and the growth of the spread across the blocks per doubling of x (x 16 for a split
      at x^4, x 64 at x^6), for the Heisenberg chains N = 4, 5, 6, the
      XXZ chain at Delta = 1/2, the XY chain, the Heisenberg ring N = 6, the Heisenberg star N = 6, the XY stars
      N = 4..7, the XY
      rings N = 4, 5, 6 and 7, the XY graphs K2,3 and K2,4, the XY and Heisenberg ladders of 2 x 3 sites (squares: the hops' own correction differs
      between the fillings at x^4) and the XY star N = 6 with the Jordan-Wigner sign on its hops, free fermions, where
      every filling ties, so that the spin star's split at x^6 is consistent with its exchange statistics, and the
      chain N = 6 with ZZ couplings (+1, +1, -1, -1, +1), whose hops have flanking couplings of opposite sign and whose
      split is the uniform chain's reversed, the boundary slowest; each row prints the split
      (block p's slowest level less block 1's)/x^4; and the Heisenberg chain N = 5 under experiments/GAMMA_AS_BINDING.md's profile of rates (IBM
      Torino T2, four of five sites far below J), at its own rates and with the same shape at the Zeno end, with the
      slowest rate over every block (p, q) at those rates, the (0, 1) band edge's, and the next block's, and the same
      chain at uniform rates 0.05, 0.5, 2; and
      the chain N = 6 with fields, XY (every filling ties, the fields quadratic in the fermions) and Heisenberg
      with a field 4 on both end sites (a third filling, block 4, picked).
  Z14 The odd ring: the sublattice sign K = prod_{l odd} Z_l flips every bond of an odd ring but the wrapped one, so
      K P_T A(ring) P_T K == A(the ring with its wrapped bond negated, a pi flux) entry by entry in every block at
      N = 3, 5, 7, P_T the transpose |a><b| -> |b><a| and A the commutator with H (integer entries, signed
      permutations: exact).  Control: on the even rings N = 4, 6 the same map returns A(ring) itself and misses the
      flux.  Reading: the XY ring's slowest rate per block at x = 0.5, 1, 3 (N = 7 at 1), the odd rings tied to
      rounding, the even rings split into the even fillings and the odd ones, each class tied within itself, and the
      ring N = 8 at x = 0.3 (blocks 1..4, the rest by X^N), the classes' agreement read against the rounding, 64 eps
      times the largest block norm.
  Z15 Reading, not gated: the stars' survivor in block (1, 1) at Q = J/gamma = 1 to 5 and 20, N = 4..8, Heisenberg
      and XY: the lifted wave (hub empty, arms carrying populations), the hub-against-arms mode or a mode without
      populations, with its <n_XY> at Q = 20 against 4/N and 4/(N - 1) and its largest |Im| over those Q; and the Heisenberg ring N = 4's block (2, 2) at Q = 0.3 to 100: how many of
      its modes lie below the plane <n_XY> = 1, the height of its slowest mode and that mode's populations against
      the domain-wall count.

Run: python simulations/f50_zeno_end_ferromagnet.py   (about two minutes; prints "ALL GATES PASS" or a FAIL line)
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
FAILS = []
RNG = np.random.default_rng(20261008)


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{(': ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


# ----------------------------------------------------------------------------------------------------------
# builders (site l is bit N-1-l)
# ----------------------------------------------------------------------------------------------------------
def bit(N, x, l):
    return (x >> (N - 1 - l)) & 1


def confs_of(N, p):
    return [x for x in range(2 ** N) if bin(x).count("1") == p]


def chain(N, w=None):
    return [(i, i + 1, 1 if w is None else w[i]) for i in range(N - 1)]


def ring(N):
    return [(i, (i + 1) % N, 1) for i in range(N)]


def star(N):
    return [(0, i, 1) for i in range(1, N)]


def complete(N):
    return [(i, j, 1) for i in range(N) for j in range(i + 1, N)]


def random_tree(N):
    return [(int(RNG.integers(0, v)), v, int(RNG.integers(1, 6))) for v in range(1, N)]


def random_graph(N, extra):
    """a random spanning tree plus `extra` random chords, integer weights"""
    E = {(min(a, b), max(a, b)): c for a, b, c in random_tree(N)}
    while len(E) < N - 1 + extra:
        a, b = sorted(RNG.choice(N, 2, replace=False).tolist())
        E.setdefault((a, b), int(RNG.integers(1, 6)))
    return [(a, b, c) for (a, b), c in E.items()]


def exclusion_laplacian(N, p, edges):
    """Integer Laplacian of the exclusion graph of block p: an edge of weight c per bond whose sites disagree."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    L = [[0] * len(cf) for _ in cf]
    for a, b, c in edges:
        for x in cf:
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                L[pos[x]][pos[x]] += c
                L[pos[x]][pos[y]] -= c
    return L


def site_laplacian(N, edges):
    L = [[0] * N for _ in range(N)]
    for a, b, c in edges:
        L[a][a] += c; L[b][b] += c; L[a][b] -= c; L[b][a] -= c
    return L


I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_MINUS = np.array([[0, 0], [1, 0]], dtype=complex)   # |0> -> |1>: adds one excitation


def kron_op(N, site_ops):
    out = np.ones((1, 1), dtype=complex)
    for l in range(N):
        out = np.kron(out, site_ops.get(l, I2))
    return out


def ferromagnet_times_two(N, edges, with_yy=True):
    """sum_b c_b (I - XX - YY - ZZ) = 2 sum_b c_b (1 - SWAP_b), from Kronecker products (integer entries)."""
    d = 2 ** N
    F = np.zeros((d, d), dtype=complex)
    for a, b, c in edges:
        F += c * (np.eye(d) - kron_op(N, {a: PX, b: PX}) - kron_op(N, {a: PZ, b: PZ})
                  - (kron_op(N, {a: PY, b: PY}) if with_yy else 0))
    return F


def ham_matrix(N, edges, delta=Fraction(1), dm=0):
    """H = sum_b c_b (XX + YY + delta ZZ) + dm (XY - YX) on the first bond; entries exact (Fraction) in a dict-of-dicts."""
    d = 2 ** N
    H = [dict() for _ in range(d)]

    def add(r, s, v):
        H[r][s] = H[r].get(s, 0) + v
    for k, (a, b, c) in enumerate(edges):
        for x in range(d):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                add(y, x, 2 * c)
                if dm and k == 0:
                    add(y, x, 2j * dm * (1 if bit(N, x, a) else -1))
                add(x, x, -delta * c)
            else:
                add(x, x, delta * c)
    return H


def add_diagonal(N, H, fields=(), zz=()):
    """Add sum_l h_l Z_l + sum K Z_a Z_b (terms diagonal in the configurations) to H in place."""
    for x in range(2 ** N):
        z = [1 - 2 * bit(N, x, l) for l in range(N)]
        v = sum(Fraction(h) * z[l] for l, h in enumerate(fields)) + sum(Fraction(K) * z[a] * z[b] for a, b, K in zz)
        if v:
            H[x][x] = H[x].get(x, 0) + v
    return H


def ham_correlated_hop(N, edges, k):
    """XX + YY + ZZ, the hop on the first bond doubled when site k is occupied (number conserving, not SU(2))."""
    d = 2 ** N
    H = [dict() for _ in range(d)]
    for idx_b, (a, b, c) in enumerate(edges):
        for x in range(d):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                amp = 2 * c * (2 if (idx_b == 0 and bit(N, x, k)) else 1)
                H[y][x] = H[y].get(x, 0) + amp
                H[x][x] = H[x].get(x, 0) - c
            else:
                H[x][x] = H[x].get(x, 0) + c
    return H


class GQ:
    """An exact Gaussian rational a + b i with Fraction parts."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=0):
        if isinstance(re, complex):
            re, im = Fraction(re.real), Fraction(re.imag)   # integer-valued floats only (asserted below)
        self.re, self.im = Fraction(re), Fraction(im)

    @staticmethod
    def of(v):
        if isinstance(v, GQ):
            return v
        if isinstance(v, complex):
            assert v.real == int(v.real) and v.imag == int(v.imag), v
        return GQ(v)

    def __add__(self, o):
        o = GQ.of(o); return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        o = GQ.of(o); return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ.of(o) - self

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __mul__(self, o):
        o = GQ.of(o); return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = GQ.of(o)
        n = o.re * o.re + o.im * o.im
        return GQ((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def __eq__(self, o):
        o = GQ.of(o); return self.re == o.re and self.im == o.im

    def __ne__(self, o):
        return not self.__eq__(o)

    __hash__ = None


def to_gq(H):
    return [{k: GQ.of(v) for k, v in row.items()} for row in H]


def zeno_terms(N, H, p, gam):
    """Exact A_PQ Gamma^-1 A_QP and A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP on block (p, p); Gamma_xy = 2 sum gamma over
    the sites where x and y disagree (gam rational)."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    H_cols = [dict() for _ in range(2 ** N)]
    for r in range(2 ** N):
        for c_, v in H[r].items():
            H_cols[c_][r] = v

    def A_vec(vec):
        out = {}
        for (x, y), a in vec.items():
            for xp, v in H_cols[x].items():
                if xp in pos:
                    out[(xp, y)] = out.get((xp, y), 0) + v * a
            for yp, v in H[y].items():
                if yp in pos:
                    out[(x, yp)] = out.get((x, yp), 0) - v * a
        return {k_: v for k_, v in out.items() if v != 0}

    def rate(x, y):
        return 2 * sum(gam[l] for l in range(N) if bit(N, x, l) != bit(N, y, l))

    n = len(cf)
    M1 = [[0] * n for _ in range(n)]
    M3 = [[0] * n for _ in range(n)]
    for j, x in enumerate(cf):
        first = {k_: v / rate(*k_) for k_, v in A_vec({(x, x): 1}).items()}           # Gamma^-1 A_QP, all in Q
        second = A_vec(first)
        for (u, w), v in second.items():
            if u == w:
                M1[pos[u]][j] += v
        mid = {k_: v / rate(*k_) for k_, v in second.items() if k_[0] != k_[1]}     # Gamma^-1 A_QQ Gamma^-1 A_QP
        third = A_vec(mid)
        for (u, w), v in third.items():
            if u == w:
                M3[pos[u]][j] += v
    return M1, M3


def ad_block_squared_pp(N, H, p):
    """(A^2)_PP for A = H (x) I - I (x) H^T on block (p, p), exact (Fraction / complex Fraction pairs)."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}

    def A_col(x, y):
        out = {}
        for xp, v in H_cols[x].items():
            if xp in pos:
                out[(xp, y)] = out.get((xp, y), 0) + v
        for yp, v in H[y].items():            # (rho H)_{x, yp} = sum_y rho_{x y} H_{y yp}
            if yp in pos:
                out[(x, yp)] = out.get((x, yp), 0) - v
        return out
    H_cols = [dict() for _ in range(2 ** N)]  # H_cols[x][xp] = H[xp][x]
    for r in range(2 ** N):
        for s, v in H[r].items():
            H_cols[s][r] = v
    n = len(cf)
    M = [[0] * n for _ in range(n)]
    for j, x in enumerate(cf):
        first = A_col(x, x)
        second = {}
        for (u, v), a1 in first.items():
            for cell, a2 in A_col(u, v).items():
                second[cell] = second.get(cell, 0) + a2 * a1
        for i, z in enumerate(cf):
            M[i][j] = second.get((z, z), 0)
    return M


def block_of(Mfull, N, p):
    cf = confs_of(N, p)
    return Mfull[np.ix_(cf, cf)]


# ----------------------------------------------------------------------------------------------------------
# exact eigenvalue counting: Jacobi's rule on the leading principal minors, Bareiss over the integers
# ----------------------------------------------------------------------------------------------------------
def leading_minors(M):
    A = [row[:] for row in M]; n = len(A); minors = []; prev = 1
    for k in range(n):
        if A[k][k] == 0:
            return None
        minors.append(A[k][k])
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev   # exact division (Bareiss)
        prev = A[k][k]
    return minors


def count_below(L, q):
    """# eigenvalues of the integer symmetric matrix L strictly below the rational q (exact; None on a zero minor)."""
    a, b = q.numerator, q.denominator
    n = len(L)
    m = leading_minors([[b * L[i][j] - (a if i == j else 0) for j in range(n)] for i in range(n)])
    if m is None:
        return None
    seq = [1] + m
    return sum(1 for k in range(n) if (seq[k] > 0) != (seq[k + 1] > 0))


def rational_in(lo, hi):
    """a rational in the middle half of (lo, hi); floats only place it, the counts verify it exactly"""
    for den in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10, 10 ** 12):
        q = Fraction(round((lo + hi) / 2 * den), den)
        if lo + (hi - lo) / 4 < float(q) < hi - (hi - lo) / 4:
            return q
    raise RuntimeError("window too narrow")


def remove_vertex(N, edges, v):
    keep = [k for k in range(N) if k != v]; ren = {k: n for n, k in enumerate(keep)}
    return N - 1, [(ren[a], ren[b], c) for a, b, c in edges if v not in (a, b)]


def leaves(N, edges):
    deg = [0] * N
    for a, b, c in edges:
        deg[a] += 1; deg[b] += 1
    return [v for v in range(N) if deg[v] == 1]


def site_spectrum(N, edges):
    return np.linalg.eigvalsh(np.array(site_laplacian(N, edges), float))


# ----------------------------------------------------------------------------------------------------------
# Z1: the exclusion Laplacian is the ferromagnet
# ----------------------------------------------------------------------------------------------------------
print("Z1  2 Lap(exclusion graph, block p) == block p of sum_b c_b (I - XX - YY - ZZ), exact")
graphs_z1 = [(f"chain N={N}", N, chain(N)) for N in (3, 4, 5, 6)] + [(f"ring N={N}", N, ring(N)) for N in (4, 5)] + \
            [(f"star N={N}", N, star(N)) for N in (4, 5)] + [("K5", 5, complete(5))] + \
            [(f"random graph N=6 #{t}", 6, random_graph(6, 3)) for t in range(3)]
for label, N, E in graphs_z1:
    F2 = ferromagnet_times_two(N, E)
    ok = np.array_equal(F2.imag, np.zeros_like(F2.imag))
    for p in range(N + 1):
        ok &= np.array_equal(block_of(F2.real, N, p), 2 * np.array(exclusion_laplacian(N, p, E), float))
    check(f"{label}: every block p = 0..{N}", ok)
N, E = 5, random_graph(5, 2)
F2bad = ferromagnet_times_two(N, E, with_yy=False)
bad = any(not np.array_equal(block_of(F2bad.real, N, p), 2 * np.array(exclusion_laplacian(N, p, E), float)) for p in range(1, N))
check("control: the operator built without its YY term differs from 2 Lap in some block 1 <= p <= N-1", bad)

# ----------------------------------------------------------------------------------------------------------
# Z2: the Zeno generator of every XXZ is the Laplacian with weights J_b^2
# ----------------------------------------------------------------------------------------------------------
print("Z2  (A^2)_PP == 8 Lap(weights J_b^2) exactly, every Delta")
for N, p, prof in ((4, 1, [1, 1, 1]), (4, 2, [1, 1, 1]), (5, 2, [1, 2, 3, 1]), (6, 3, [2, 1, 1, 3, 1])):
    E = chain(N, prof)
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    for delta in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(-1)):
        M = ad_block_squared_pp(N, ham_matrix(N, E, delta=delta), p)
        same = all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M)))
        check(f"chain N={N} profile {prof} ({p},{p}) Delta = {delta}", same)
for N, p in ((5, 2),):
    E = random_graph(N, 2)
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    M = ad_block_squared_pp(N, ham_matrix(N, E, delta=Fraction(3, 2)), p)
    check(f"random graph N={N} ({p},{p}) Delta = 3/2", all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M))))
N = 5
E = chain(N, [1, 2, 3, 1])
H_diag = add_diagonal(N, ham_matrix(N, E, delta=Fraction(1, 2)), fields=(2, -1, 3, 0, -2), zz=((0, 2, 2), (1, 3, -1), (2, 4, 3)))
ok = True
for p in (1, 2):
    Lsq = exclusion_laplacian(N, p, [(a, b, c * c) for a, b, c in E])
    M = ad_block_squared_pp(N, H_diag, p)
    ok &= all(M[i][j] == 8 * Lsq[i][j] for i in range(len(M)) for j in range(len(M)))
check("chain N=5, Delta = 1/2 with the fields (2, -1, 3, 0, -2) and next-nearest ZZ: (A^2)_PP == 8 Lap(J_b^2), blocks 1, 2", ok)
N = 4
E = chain(N)
ok = True
for p in (1, 2):
    Lw = exclusion_laplacian(N, p, [(a, b, c * c + (1 if a == 0 else 0)) for a, b, c in E])   # D_0 = 1 on bond 0
    M = ad_block_squared_pp(N, ham_matrix(N, E, dm=1), p)
    ok &= all(M[i][j] == 8 * Lw[i][j] for i in range(len(M)) for j in range(len(M)))
check("chain N=4 with a DM term D = 1 on bond 0: (A^2)_PP == 8 Lap(J_b^2 + D_b^2), blocks 1, 2", ok)
E = chain(5, [1, 2, 3, 1])
Lw = exclusion_laplacian(5, 2, E)
M = ad_block_squared_pp(5, ham_matrix(5, E), 2)
check("control: the weights J_b in place of J_b^2 miss (A^2)_PP on the chain N=5 with profile [1, 2, 3, 1]",
      any(M[i][j] != 8 * Lw[i][j] for i in range(len(M)) for j in range(len(M))))


def zeno_full(N, H):
    """The Zeno generator (A^2)_PP of every block assembled on the 2^N configurations (exact lists)."""
    Z = [[0] * 2 ** N for _ in range(2 ** N)]
    for p in range(N + 1):
        cf = confs_of(N, p)
        M = ad_block_squared_pp(N, H, p)
        for i, x in enumerate(cf):
            for j, y in enumerate(cf):
                Z[x][y] = M[i][j]
    return Z


def commutes_with_lowering(N, Z):
    Sm = np.real(sum(kron_op(N, {l: SIGMA_MINUS}) for l in range(N))).astype(int).tolist()
    d = 2 ** N
    ZS = [[sum(Z[i][k] * Sm[k][j] for k in range(d) if Sm[k][j]) for j in range(d)] for i in range(d)]
    SZ = [[sum(Sm[i][k] * Z[k][j] for k in range(d) if Sm[i][k]) for j in range(d)] for i in range(d)]
    return all(ZS[i][j] == SZ[i][j] for i in range(d) for j in range(d))


rows = [("chain N=5, Delta = 1/2, fields and next-nearest ZZ", 5, H_diag),
        ("chain N=4, DM D = 1 on bond 0", 4, ham_matrix(4, chain(4), dm=1)),
        ("random graph N=5, Delta = 2", 5, ham_matrix(5, random_graph(5, 2), delta=Fraction(2)))]
for label, N, H in rows:
    check(f"{label}: the assembled Zeno generator commutes with S^- exactly", commutes_with_lowering(N, zeno_full(N, H)))
check("control: a correlated hop (doubled when site 2 is occupied) breaks the commutation with S^-",
      not commutes_with_lowering(5, zeno_full(5, ham_correlated_hop(5, chain(5), 2))))

# ----------------------------------------------------------------------------------------------------------
# Z2b: a rational profile of rates, the leading Schur term exact and the next one zero
# ----------------------------------------------------------------------------------------------------------
print("Z2b A_PQ Gamma^-1 A_QP == Lap(4 J_b^2/(gamma_i + gamma_j)) and A_PQ Gamma^-1 A_QQ Gamma^-1 A_QP == 0, exact")
gam_q = [Fraction(int(v), 7) for v in RNG.integers(3, 15, size=5)]
TRI = [(0, 1, 1), (1, 2, 2), (0, 2, 1), (2, 3, 1), (3, 4, 2)]
for label, N, E, H in (
        ("chain N=5, Delta = 1/2, fields, next-nearest ZZ", 5, chain(5, [1, 2, 1, 3]),
         add_diagonal(5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(1, 2)), fields=(1, 0, -2, 1, 3), zz=((0, 2, 1),))),
        ("triangle graph N=5 (bonds 01, 12, 02, 23, 34), Delta = 3/2, fields", 5, TRI,
         add_diagonal(5, ham_matrix(5, TRI, delta=Fraction(3, 2)), fields=(1, -2, 0, 3, 1)))):
    ok1, ok3 = True, True
    for p in (1, 2):
        M1, M3 = zeno_terms(N, H, p, gam_q)
        Lw = exclusion_laplacian(N, p, [(a, b, Fraction(4) * Fraction(c) ** 2 / (gam_q[a] + gam_q[b])) for a, b, c in E])
        ok1 &= all(M1[i][j] == Lw[i][j] for i in range(len(M1)) for j in range(len(M1)))
        ok3 &= all(v == 0 for row in M3 for v in row)
    check(f"{label}, gamma_l = {[str(g) for g in gam_q]}: leading term == Lap(4 J_b^2/(gamma_i + gamma_j)), blocks 1, 2", ok1)
    check(f"{label}: the next term vanishes, blocks 1, 2", ok3)
E = chain(5, [1, 2, 1, 3])
M1, _ = zeno_terms(5, ham_matrix(5, E), 2, gam_q)
gbar = sum(gam_q) / 5
Lwrong = exclusion_laplacian(5, 2, [(a, b, Fraction(2) * Fraction(c) ** 2 / gbar) for a, b, c in E])
check("control: the uniform-rate weights 2 J_b^2/gamma_mean miss the leading term",
      any(M1[i][j] != Lwrong[i][j] for i in range(len(M1)) for j in range(len(M1))))
_, M3f = zeno_terms(5, to_gq(ham_matrix(5, TRI, delta=Fraction(3, 2), dm=1)), 2, gam_q)
check("control: a DM flux through the triangle (H not real) makes the next term nonzero (exact, Gaussian rationals)",
      any(v != 0 for row in M3f for v in row))
# that term is i M3 with M3 purely imaginary and antisymmetric: real and antisymmetric on the populations, so it moves
# no real part (L commutes with X -> X^dagger, which conjugates the populations)
check("the flux term M3 is purely imaginary and antisymmetric, so i M3 is real and antisymmetric (exact)",
      all(GQ.of(v).re == 0 for row in M3f for v in row)
      and all(GQ.of(M3f[r][c]) == -GQ.of(M3f[c][r]) for r in range(len(M3f)) for c in range(len(M3f))))
# a complex H on a graph without a triangle: the next term vanishes all the same
for label, Nn, Hn in (("chain N=5 with a DM term", 5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(1, 2), dm=1)),
                      ("ring N=4 with a DM term (a flux)", 4, ham_matrix(4, ring(4), delta=Fraction(1), dm=1))):
    zero = all(v == 0 for p in (1, 2) for row in zeno_terms(Nn, to_gq(Hn), p, gam_q[:Nn])[1] for v in row)
    check(f"{label}, H not real, no triangle: the next term vanishes, blocks 1, 2 (exact, Gaussian rationals)", zero)

# ----------------------------------------------------------------------------------------------------------
# Z3: the uniform Heisenberg graph: Lap = (#bonds - H)/2, the Hamiltonian read from the top
# ----------------------------------------------------------------------------------------------------------
print("Z3  uniform Heisenberg: 2 Lap(block p) == #bonds I - H_p exactly")


def dense_block(H, N, p):
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    n = len(cf); M = [[0] * n for _ in range(n)]
    for x in cf:
        for y, v in H[x].items():
            if y in pos:
                M[pos[x]][pos[y]] = v
    return M


for label, N, E in ((f"chain N=6", 6, chain(6)), ("ring N=5", 5, ring(5)), ("star N=5", 5, star(5)), ("K4", 4, complete(4))):
    ok = True
    for p in range(N + 1):
        Hp = dense_block(ham_matrix(N, E, delta=Fraction(1)), N, p)
        L = exclusion_laplacian(N, p, E)
        n = len(L)
        ok &= all(2 * L[i][j] == (len(E) if i == j else 0) - Hp[i][j] for i in range(n) for j in range(n))
    check(f"{label}: every block", ok)
N = 6; E = chain(N); p = 3
Hp = dense_block(ham_matrix(N, E, delta=Fraction(2)), N, p); L = exclusion_laplacian(N, p, E); n = len(L)
check("control: Delta = 2 breaks the identity", any(2 * L[i][j] != (len(E) if i == j else 0) - Hp[i][j] for i in range(n) for j in range(n)))

# ----------------------------------------------------------------------------------------------------------
# Z4: the lift is SU(2)'s lowering
# ----------------------------------------------------------------------------------------------------------
print("Z4  (S^-)^(p-1) sum_i f(i)|i> == (p-1)! sum_S F(S)|S>;  [ferromagnet, S^-] == 0;  Lap F == lift(L_site f)")


def lowering(N, missing=None):
    return sum(kron_op(N, {l: SIGMA_MINUS}) for l in range(N) if l != missing)


for label, N, E in (("chain N=6", 6, chain(6, [1, 2, 1, 3, 1])), ("random tree N=7", 7, random_tree(7)), ("ring N=5", 5, ring(5))):
    f = RNG.integers(-5, 6, size=N)
    v = np.zeros(2 ** N, dtype=complex)
    for i in range(N):
        v[1 << (N - 1 - i)] = f[i]
    Sm = lowering(N)
    F2 = ferromagnet_times_two(N, E)
    ok_comm = np.array_equal(F2 @ Sm - Sm @ F2, np.zeros_like(F2))
    ok_lift, ok_eig = True, True
    Lsite = np.array(site_laplacian(N, E))
    w = v.copy()
    for p in range(1, N):
        cf = confs_of(N, p)
        F = np.array([sum(f[l] for l in range(N) if bit(N, x, l)) for x in cf])
        ok_lift &= np.array_equal(w[cf], factorial(p - 1) * F.astype(complex)) and \
            np.array_equal(np.delete(w, cf), np.zeros(2 ** N - len(cf)))
        Lf = Lsite @ f
        liftLf = np.array([sum(Lf[l] for l in range(N) if bit(N, x, l)) for x in cf])
        ok_eig &= np.array_equal(np.array(exclusion_laplacian(N, p, E)) @ F, liftLf)
        w = Sm @ w
    check(f"{label}: [2 sum c_b (1 - SWAP_b), S^-] == 0", ok_comm)
    check(f"{label}: (S^-)^(p-1) f == (p-1)! F in every block 1..N-1", ok_lift)
    check(f"{label}: Lap F == lift(L_site f) in every block 1..N-1", ok_eig)
N = 5; f = np.array([3, -1, 4, 1, -5])
v = np.zeros(2 ** N, dtype=complex)
for i in range(N):
    v[1 << (N - 1 - i)] = f[i]
w = lowering(N, missing=2) @ v
cf = confs_of(N, 2)
F = np.array([sum(f[l] for l in range(N) if bit(N, x, l)) for x in cf])
check("control: S^- with site 2 missing does not give the lift at p = 2", not np.array_equal(w[cf], F.astype(complex)))
# controls on the commutation and on the lift identity, through the same doors
Ec = chain(6, [1, 2, 1, 3, 1])
Sm6 = lowering(6)
F2a = ferromagnet_times_two(6, Ec) - sum(c * kron_op(6, {a: PZ, b: PZ}) for a, b, c in Ec)
check("control: the anisotropic sum_b c_b (I - XX - YY - 2 ZZ) does not commute with S^-",
      not np.array_equal(F2a @ Sm6 - Sm6 @ F2a, np.zeros_like(F2a)))
f6 = np.array([2, -3, 1, 4, -1, -3])
Lbad = np.array(site_laplacian(6, Ec)); Lbad[2, 2] += 1
Lf6 = Lbad @ f6
F6 = np.array([sum(f6[l] for l in range(6) if bit(6, x, l)) for x in confs_of(6, 2)])
liftbad = np.array([sum(Lf6[l] for l in range(6) if bit(6, x, l)) for x in confs_of(6, 2)])
check("control: a site Laplacian with one diagonal entry raised breaks Lap F == lift(L_site f) at p = 2",
      not np.array_equal(np.array(exclusion_laplacian(6, 2, Ec)) @ F6, liftbad))

# ----------------------------------------------------------------------------------------------------------
# Z5: the theorem on trees, by exact counts
# ----------------------------------------------------------------------------------------------------------
print("Z5  every tree, every block: only the zero mode below q_lo < lambda_1; the tree's count below q_hi and below q*")


def tree_rows(label, N, E, token_edges=None, extra=None, expect=True, miss_only_at=None):
    ev = site_spectrum(N, E)
    l1 = ev[1]
    above = ev[ev > l1 * (1 + 1e-9) + 1e-12]
    l2 = above[0] if len(above) else l1 + 1.0
    w = 1e-6 * max(1.0, l1)
    q_lo = rational_in(l1 - w, l1)
    q_hi = rational_in(l1, min(l1 + w, (l1 + l2) / 2))
    LT = site_laplacian(N, E)
    m_lo, m_hi = count_below(LT, q_lo), count_below(LT, q_hi)
    bracket_ok = (m_lo == 1 and m_hi is not None and m_hi >= 2)
    # q*: just under max over leaves of lambda_1(T - v), verified exactly on that leaf's Laplacian (N >= 3)
    qs, m_s = None, None
    if N >= 3:
        best = max(leaves(N, E), key=lambda v: site_spectrum(*remove_vertex(N, E, v))[1])
        Nm, Em = remove_vertex(N, E, best)
        lv = site_spectrum(Nm, Em)[1]
        qs = rational_in(lv * (1 - 2e-6), lv * (1 - 1e-6))
        bracket_ok &= count_below(site_laplacian(Nm, Em), qs) == 1
        m_s = count_below(LT, qs)
        bracket_ok &= m_s is not None
    ok = bracket_ok
    misses = []
    for p in range(1, N):
        Lt = exclusion_laplacian(N, p, E if token_edges is None else token_edges)
        if extra is not None:
            Xp = extra(N, p)
            Lt = [[Lt[i][j] + Xp[i][j] for j in range(len(Lt))] for i in range(len(Lt))]
        got = (count_below(Lt, q_lo), count_below(Lt, q_hi), count_below(Lt, qs) if qs is not None else None)
        if None in got[:2] or (qs is not None and got[2] is None) or got != (m_lo, m_hi, m_s):
            misses.append((p, got))
    ok &= not misses
    qs_note = ""
    if qs is not None:
        qs_note = f"; below q* = {float(qs):.6f}, under lambda_1(T - leaf): {m_s}"
        if float(qs) < l1:
            qs_note += " (q* below lambda_1: the leaf bound adds nothing here)"
    if expect:
        check(f"{label}: brackets exact, every block counts as the tree (zero mode; {m_hi - 1} at lambda_1 = {l1:.6f}"
              + qs_note + ")", ok, f"misses {misses}" if misses else "")
    else:
        # a control counts only with exact brackets and a definite count that differs (a zero minor is no miss);
        # miss_only_at names the brackets that must differ, the others agreeing
        want = (m_lo, m_hi, m_s)
        definite = [(p, g) for p, g in misses if None not in g[:2] and (qs is None or g[2] is not None)]
        if miss_only_at is not None:
            definite = [(p, g) for p, g in definite if {k for k in range(3) if g[k] != want[k]} == set(miss_only_at)]
        check(f"{label}: brackets exact, some block counts differently from the tree", bracket_ok and bool(definite),
              f"misses {misses}" if misses else "")
    return ok


for N in range(2, 10):
    tree_rows(f"path N={N}", N, chain(N))
for N in range(3, 10):
    tree_rows(f"star N={N}", N, star(N))
for N in range(3, 10):
    for t in range(2):
        tree_rows(f"random tree N={N} #{t}", N, random_tree(N))
for N in (6, 8):
    prof = [int(v) for v in RNG.integers(1, 4, size=N - 1)]
    tree_rows(f"chain N={N}, Zeno weights J_b^2 for J_b = {prof}", N, chain(N, [j * j for j in prof]))
# control: the same door, the blocks of a reweighted chain counted against the original chain's brackets
E_ctrl = chain(7)
tree_rows("control: chain N=7 brackets, blocks built with the bond (3, 4) doubled", 7, E_ctrl,
          token_edges=[(a, b, 2 * c if a == 3 else c) for a, b, c in E_ctrl], expect=False)
E_ctrl2 = chain(7, [2] * 6)
tree_rows("control: chain N=7 (bonds 2) brackets, blocks built with the bond (3, 4) halved", 7, E_ctrl2,
          token_edges=[(a, b, 1 if a == 3 else c) for a, b, c in E_ctrl2], expect=False)


# controls on what the theorem adds: -c (1 - SWAP_01)(1 - SWAP_23) is SU(2)-invariant and nonzero only on states
# with a singlet on both pairs, so it leaves the ferromagnetic multiplet and the one-magnon band in place and moves
# only the multiplets of lower spin, the ones the leaf bound puts at or above lambda_1(T - v)
def four_site(c):
    def extra(N, p):
        A, B = exclusion_laplacian(N, p, [(0, 1, 1)]), exclusion_laplacian(N, p, [(2, 3, 1)])
        n = len(A)
        return [[-c * sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return extra


tree_rows("control: chain N=6 with -(1 - SWAP_01)(1 - SWAP_23) in every block, a new multiplet below lambda_1", 6,
          chain(6), extra=four_site(1), expect=False, miss_only_at=(0, 1, 2))
tree_rows("control: chain N=6, bonds 100, with -95 (1 - SWAP_01)(1 - SWAP_23), a new multiplet between lambda_1 and "
          "lambda_1(T - v): the q* count alone differs", 6, chain(6, [100] * 5), extra=four_site(95), expect=False,
          miss_only_at=(2,))
# the exact counter against a float count on a matrix with a known spectrum, far from the eigenvalues
Lp = exclusion_laplacian(6, 3, chain(6))
fl = np.linalg.eigvalsh(np.array(Lp, float))
qq = Fraction(7, 3)
check("cross-check: the exact counter agrees with a float count at q = 7/3 on the N=6 half-filling block",
      count_below(Lp, qq) == int(np.sum(fl < 7 / 3)) and np.abs(fl - 7 / 3).min() > 1e-3)

# ----------------------------------------------------------------------------------------------------------
# Z6: the full Liouvillian's Zeno limit is blind to Delta
# ----------------------------------------------------------------------------------------------------------
print("Z6  slow eigenvalues of B(x) = -2 x^2 ell_j + O(x^4) at Delta = 0, 1/2, 1, 2, the next order even (differences of "
      "deviation/x^4 in the ratio 4); a rate profile (deviation*s^3, differences in the ratio 1/4)")


def float_block(N, H, p):
    cf = confs_of(N, p)
    idx = [(x, y) for x in cf for y in cf]
    pos = {e: i for i, e in enumerate(idx)}
    n = len(idx)
    Hd = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for r in range(2 ** N):
        for s, v in H[r].items():
            Hd[r, s] = complex(v)
    A = np.zeros((n, n), dtype=complex)
    for c, (x, y) in enumerate(idx):
        for xp in cf:
            if Hd[xp, x]:
                A[pos[(xp, y)], c] += Hd[xp, x]
        for yp in cf:
            if Hd[y, yp]:
                A[pos[(x, yp)], c] -= Hd[y, yp]
    ham = np.array([bin(x ^ y).count("1") for x, y in idx])
    return A, ham


XS = (0.005, 0.01, 0.02)


def x4_law(slow, ell_used):
    """deviation/x^4 over XS and the ratio of its successive differences, which names the power of the leading
    correction: 4 for x^2, 16 for x^4, 1/2 for x^-1 (an odd order), 1/4 for the x^-2 term 2 delta ell/x^2 that a
    leading generator off by a relative delta adds"""
    rho = [np.abs(slow[x] - np.sort(-2 * x * x * ell_used)).max() / x ** 4 for x in XS]
    return rho, (rho[2] - rho[1]) / (rho[1] - rho[0])


def names_x2(ratio):
    """the x^2 law against its neighbours x^1 and x^3, decided at the half-integer powers"""
    return 2 ** 1.5 < ratio < 2 ** 2.5


for N, p in ((4, 1), (4, 2), (5, 2)):
    E = chain(N)
    ell = np.sort(np.linalg.eigvalsh(np.array(exclusion_laplacian(N, p, E), float)))
    cases = [(str(d), ham_matrix(N, E, delta=d)) for d in (Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2))]
    if N == 4:
        cases.append(("1 with fields (2, -1, 3, 0) and ZZ on (0, 2)",
                      add_diagonal(N, ham_matrix(N, E, delta=Fraction(1)), fields=(2, -1, 3, 0), zz=((0, 2, 1),))))
    for delta, Hc in cases:
        A, ham = float_block(N, Hc, p)
        slow = {x: np.sort(np.sort(np.linalg.eigvals(-1j * x * A + np.diag(-2.0 * ham)).real)[::-1][:comb(N, p)])
                for x in XS}
        rho, ratio = x4_law(slow, ell)
        _, ratio_c = x4_law(slow, ell * (1 + 1e-7))
        check(f"chain N={N} ({p},{p}) Delta = {delta}: deviation/x^4 = {rho[0]:.4f}, {rho[1]:.4f}, {rho[2]:.4f} at "
              f"x = 0.005, 0.01, 0.02, differences in the ratio {ratio:.4f} (x^2: 4, decided at the half-integer powers); "
              f"control ell (1 + 1e-7): {ratio_c:.4f}", names_x2(ratio) and not names_x2(ratio_c))

# Z6b: a profile of rates and couplings, the weights 4 J_b^2/(gamma_i + gamma_j)
N = 5
J_prof = [Fraction(int(v)) for v in RNG.integers(1, 4, size=N - 1)]
gam_prof = RNG.uniform(0.5, 2.0, size=N)
E = [(i, i + 1, J_prof[i]) for i in range(N - 1)]
H7 = ham_matrix(N, E, delta=Fraction(7, 10))


def profile_rows(weights_of, label, expect_law):
    ratios = []
    for p in range(1, N):
        A, _ = float_block(N, H7, p)
        cf = confs_of(N, p)
        rate = np.array([sum(gam_prof[l] for l in range(N) if bit(N, x, l) != bit(N, y, l)) for x in cf for y in cf])
        devs = {}
        for scale in (50.0, 100.0, 200.0):
            w = weights_of(scale)
            Lw = np.zeros((len(cf), len(cf)))
            pos = {x: k for k, x in enumerate(cf)}
            for (a, b, c), wb in zip(E, w):
                for x in cf:
                    if bit(N, x, a) != bit(N, x, b):
                        y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                        Lw[pos[x], pos[x]] += wb
                        Lw[pos[x], pos[y]] -= wb
            ev = np.linalg.eigvals(-1j * A + np.diag(-2.0 * scale * rate))
            slow = np.sort(np.sort(ev.real)[::-1][:len(cf)])
            devs[scale] = np.abs(slow - np.sort(-np.linalg.eigvalsh(Lw))).max()
        rho = [devs[sc] * sc ** 3 for sc in (50.0, 100.0, 200.0)]
        ratios.append((rho[2] - rho[1]) / (rho[1] - rho[0]))
    holds = [2 ** -2.5 < r < 2 ** -1.5 for r in ratios]    # s^-2 against s^-3 and s^-1, at the half-integer powers
    ok = all(holds) if expect_law else not any(holds)
    check(f"{label}: deviation*s^3 over s = 50, 100, 200, differences in the ratio (s^-2: 1/4, decided at the "
          f"half-integer powers; a wrong leading term: 4) per block p = 1..{N - 1}: {', '.join(f'{r:.4f}' for r in ratios)}", ok)


profile_rows(lambda scale: [4 * float(J_prof[i]) ** 2 / (scale * (gam_prof[i] + gam_prof[i + 1])) for i in range(N - 1)],
             f"chain N={N}, Delta = 7/10, J_b = {[int(j) for j in J_prof]}, gamma_l = {np.round(gam_prof, 3).tolist()}: "
             "weights 4 J_b^2/(gamma_i + gamma_j), the s^-3 law", True)
profile_rows(lambda scale: [2 * float(J_prof[i]) ** 2 / (scale * gam_prof.mean()) for i in range(N - 1)],
             "control: the uniform-profile weights 2 J_b^2/gamma_mean miss the s^-3 law", False)

# Z6c: what the next term does where it does not vanish, on a triangle with a flux
print("Z6c the triangle with a flux: the gamma^-2 term splits the degenerate pair (|Im| ~ gamma^-2); no triangle, no such term")


def ham_flux_ring(N, J, D):
    """J (XX + YY + ZZ) around the ring 0 -> 1 -> ... -> N-1 -> 0 with the same Dzyaloshinskii-Moriya D on every
    bond, oriented along the ring: a flux N arctan(D/J) through it"""
    edges = [(l, (l + 1) % N, J) for l in range(N)]
    H = ham_matrix(N, edges, delta=Fraction(1))
    for a, b, c in edges:
        for x in range(2 ** N):
            if bit(N, x, a) != bit(N, x, b):
                y = x ^ (1 << (N - 1 - a)) ^ (1 << (N - 1 - b))
                H[y][x] = H[y].get(x, 0) + 2j * D * (1 if bit(N, x, a) else -1)
    return H


for N, k_law, k_other in ((3, 2, 4), (5, 4, 2)):
    A_f, ham_f = float_block(N, ham_flux_ring(N, 1, 1), 1)
    im, re = [], []
    for g in (40.0, 80.0, 160.0):
        ev = np.linalg.eigvals(-1j * A_f + np.diag(-2.0 * g * ham_f))
        pair = sorted(ev, key=lambda z: -z.real)[1:3]
        im.append(max(abs(z.imag) for z in pair)); re.append(tuple(round(float(z.real * g), 4) for z in pair))
    ratios = [im[0] / im[1], im[1] / im[2]]
    check(f"{'triangle' if N == 3 else 'ring N=5, no triangle'} with D = J on every bond, block p = 1, the leading "
          f"pair: Re*gamma = {re} at gamma = 40, 80, 160; |Im| per doubling of gamma falls by {ratios[0]:.3f}, "
          f"{ratios[1]:.3f} (gamma^-{k_law}: {2 ** k_law}, decided at the half-integer powers; control: not read as "
          f"gamma^-{k_other})",
          all(2 ** (k_law - 0.5) < r < 2 ** (k_law + 0.5) for r in ratios)
          and not all(2 ** (k_other - 0.5) < r < 2 ** (k_other + 0.5) for r in ratios))

# ----------------------------------------------------------------------------------------------------------
# Z7: readings beyond trees
# ----------------------------------------------------------------------------------------------------------
print("Z7  reading: graphs without the leaf argument, the same exact counter")
for label, N, E in [(f"ring N={N}", N, ring(N)) for N in range(4, 10)] + [("K6", 6, complete(6))] + \
                   [(f"random graph N=7 #{t}", 7, random_graph(7, 4)) for t in range(2)]:
    ev = site_spectrum(N, E); l1 = ev[1]
    q_lo = rational_in(l1 * (1 - 1e-6), l1); q_hi = rational_in(l1, l1 * (1 + 1e-6))
    LT = site_laplacian(N, E); m_lo, m_hi = count_below(LT, q_lo), count_below(LT, q_hi)
    assert m_lo == 1, (label, m_lo)   # the bracket's lower end has only the zero mode below it
    rows = [(count_below(exclusion_laplacian(N, p, E), q_lo), count_below(exclusion_laplacian(N, p, E), q_hi)) for p in range(1, N)]
    gap_same = all(r[0] == m_lo for r in rows)
    mult = [r[1] - r[0] for r in rows]
    print(f"    {label}: lambda_1 = {l1:.6f} (multiplicity {m_hi - m_lo}); nothing below it in any block: {gap_same}; "
          f"multiplicity at lambda_1 per block p = 1..{N - 1}: {mult}")

# ----------------------------------------------------------------------------------------------------------
# Z8: the chain prefactors of CHAIN_GAP_SECTOR_DIAGNOSTIC and F1_DISSIPATION_GAP_PATTERN, in their spin book
# ----------------------------------------------------------------------------------------------------------
print("Z8  reading: the chain prefactors in the book H = (J/4) sum sigma.sigma, gamma (Z rho Z - rho), gamma = 1/2")


def spin_book_gap(N, Q, delta, g=0.5):
    """Slowest nonzero rate over the diagonal blocks of the chain, H = (J/4)(XX + YY + delta ZZ), J = Q g."""
    best = np.inf
    H = ham_matrix(N, chain(N), delta=Fraction(delta))
    for p in range(1, N):
        A, ham = float_block(N, H, p)
        ev = np.linalg.eigvals(-1j * (Q * g / 4) * A + np.diag(-2.0 * g * ham)).real
        r = -ev[-ev > 1e-9]
        best = min(best, r.min())
    return best


g = 0.5
for N in (4, 5, 6):
    z8 = N * N * (1 - np.cos(np.pi / N)) / 8
    gx, g0 = spin_book_gap(N, 2.0, 1), spin_book_gap(N, 2.0, 0)
    print(f"    Q = 2, N = {N}: <n_XY> N^2/Q^2 = {gx * N * N / (2 * g * 4):.4f} (XXX), {g0 * N * N / (2 * g * 4):.4f} (XX), "
          f"Zeno {z8:.4f};  gap N^2/(gamma Q^2) = {gx * N * N / (g * 4):.4f}, {g0 * N * N / (g * 4):.4f}, Zeno {2 * z8:.4f};  "
          f"gaps below 2 gamma: {max(gx, g0) < 2 * g}")
print("    N = 5 XXX, <n_XY> N^2/Q^2 at Q = 0.5, 1, 1.5, 2, 2.5: "
      + ", ".join(f"{spin_book_gap(5, Q, 1) * 25 / (2 * g * Q * Q):.4f}" for Q in (0.5, 1.0, 1.5, 2.0, 2.5))
      + f"  (Zeno {25 * (1 - np.cos(np.pi / 5)) / 8:.4f})")
m = np.mean([spin_book_gap(N, 0.5, 1) * N * N / (g * 0.25) for N in (4, 5, 6)])
print(f"    Q = 0.5, XXX, gap N^2/(gamma Q^2) mean over N = 4, 5, 6: {m:.4f}  "
      f"(Zeno mean {np.mean([N * N * (1 - np.cos(np.pi / N)) / 4 for N in (4, 5, 6)]):.4f})")


# ----------------------------------------------------------------------------------------------------------
# Z9-Z13: beyond the leading order, the filling the diagonal terms pick
# ----------------------------------------------------------------------------------------------------------
def zeno_terms4(N, H, p, gam):
    """Exact M1 = A_PQ G^-1 A_QP, K = A_PQ G^-2 A_QP and M4 = A_PQ G^-1 A_QQ G^-1 A_QQ G^-1 A_QP on block (p, p),
    G_xy = 2 sum gamma over the sites where x and y disagree; the slow levels are those of -M1 + i M3 + M4 + lambda K."""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    H_cols = [dict() for _ in range(2 ** N)]
    for r in range(2 ** N):
        for c_, v in H[r].items():
            H_cols[c_][r] = v

    def A_vec(vec):
        out = {}
        for (x, y), a in vec.items():
            for xp, v in H_cols[x].items():
                if xp in pos:
                    out[(xp, y)] = out.get((xp, y), 0) + v * a
            for yp, v in H[y].items():
                if yp in pos:
                    out[(x, yp)] = out.get((x, yp), 0) - v * a
        return {k_: v for k_, v in out.items() if v != 0}

    def rate(x, y):
        return 2 * sum(gam[l] for l in range(N) if bit(N, x, l) != bit(N, y, l))

    n = len(cf)
    M1 = [[0] * n for _ in range(n)]; K = [[0] * n for _ in range(n)]; M4 = [[0] * n for _ in range(n)]

    def to_pop(vec, M, j):
        for (u, w), v in A_vec(vec).items():
            if u == w:
                M[pos[u]][j] += v
    for j, x in enumerate(cf):
        first = {k_: v / rate(*k_) for k_, v in A_vec({(x, x): 1}).items()}               # G^-1 A_QP, all in Q
        to_pop(first, M1, j)
        to_pop({k_: v / rate(*k_) for k_, v in first.items()}, K, j)
        second = {k_: v / rate(*k_) for k_, v in A_vec(first).items() if k_[0] != k_[1]}
        third = {k_: v / rate(*k_) for k_, v in A_vec(second).items() if k_[0] != k_[1]}
        to_pop(third, M4, j)
    return cf, M1, K, M4


def add_particle(N, p):
    """E': block p -> block p + 1 on the populations, (E' v)(B) = sum over l in B of v(B - l)"""
    lo, hi = confs_of(N, p), confs_of(N, p + 1)
    pos = {x: k for k, x in enumerate(lo)}
    Em = [[0] * len(lo) for _ in hi]
    for r, y in enumerate(hi):
        for l in range(N):
            if bit(N, y, l):
                Em[r][pos[y ^ (1 << (N - 1 - l))]] += 1
    return Em


def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)) if A[i][k]) for j in range(len(B[0]))] for i in range(len(A))]


def ladder_commutes(N, H, gam):
    """per term (M1, K, M4): does E' M^(p) == M^(p+1) E' hold for every p = 1..N-2, exactly?"""
    T = {p: zeno_terms4(N, H, p, gam) for p in range(1, N)}
    out = []
    for idx in (1, 2, 3):
        out.append(all(mat_mul(add_particle(N, p), T[p][idx]) == mat_mul(T[p + 1][idx], add_particle(N, p))
                       for p in range(1, N - 1)))
    return tuple(out)


def off_diagonal(H):
    return [{k_: v for k_, v in row.items() if k_ != r} for r, row in enumerate(H)]


def detuning_form(N, H, p, gam):
    """D_V = sum over the hops e = {x, y} of 2 |t_e|^2 dV_e^2 / G_e^3 times the edge Laplacian, dV_e = H_yy - H_xx"""
    cf = confs_of(N, p); pos = {x: k for k, x in enumerate(cf)}
    n = len(cf); D = [[0] * n for _ in range(n)]
    for x in cf:
        for y, t in H[x].items():                 # H[x][y] = H_xy, the hop y -> x
            if y == x or y not in pos:
                continue
            dV = H[y].get(y, 0) - H[x].get(x, 0)
            G = 2 * sum(gam[l] for l in range(N) if bit(N, x, l) != bit(N, y, l))
            w = 2 * t * t.conjugate() * dV * dV / G ** 3
            D[pos[x]][pos[x]] += w; D[pos[x]][pos[y]] -= w
    return D


def lift_matrix(N, p):
    """F: site functions -> block p, F[S][i] = 1 if i in S"""
    return [[bit(N, x, i) for i in range(N)] for x in confs_of(N, p)]


def transpose(A):
    return [list(r) for r in zip(*A)]


GAM1 = [Fraction(1)] * 9
print("Z9  beyond the leading order: the population ladder E' commutes with the hopping's own terms M1, K, M4 (exact)")
gam_r = [Fraction(int(v), 5) for v in RNG.integers(3, 12, size=7)]
for label, N, E, gam in (
        ("chain N=5, bonds (1, 2, 1, 3), a rate profile", 5, chain(5, [1, 2, 1, 3]), gam_r[:5]),
        ("chain N=6", 6, chain(6), GAM1[:6]),
        ("ring N=5", 5, ring(5), GAM1[:5]), ("ring N=6", 6, ring(6), GAM1[:6]),
        ("star N=4", 4, star(4), GAM1[:4]), ("star N=5", 5, star(5), GAM1[:5]),
        ("star N=6, a rate profile", 6, star(6), gam_r[:6]), ("star N=7", 7, star(7), GAM1[:7]),
        ("random tree N=6", 6, random_tree(6), GAM1[:6]), ("random tree N=7", 7, random_tree(7), gam_r[:7])):
    res = ladder_commutes(N, ham_matrix(N, E, delta=Fraction(0)), gam)
    check(f"{label}, H = hops only: E' commutes with M1, K, M4 in every block", res == (True, True, True))
cf5, M1p, _, _ = zeno_terms4(5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(0)), 2, gam_r[:5])
Lw5 = exclusion_laplacian(5, 2, [(a, b, Fraction(4) * Fraction(c) ** 2 / (gam_r[a] + gam_r[b])) for a, b, c in chain(5, [1, 2, 1, 3])])
check("the M1 of these rows is the leading term of (a): M1 == Lap(4 J_b^2/(gamma_i + gamma_j)) on the profile chain, block 2",
      all(M1p[i][j] == Lw5[i][j] for i in range(len(M1p)) for j in range(len(M1p))))
gm5 = sum(gam_r[:5]) / 5
Lu5 = exclusion_laplacian(5, 2, [(a, b, Fraction(4) * Fraction(c) ** 2 / (2 * gm5)) for a, b, c in chain(5, [1, 2, 1, 3])])
check("control: the uniform-rate weights 4 J_b^2/(2 gamma_mean) miss that M1",
      any(M1p[i][j] != Lu5[i][j] for i in range(len(M1p)) for j in range(len(M1p))))
res = ladder_commutes(4, ham_matrix(4, ring(4), delta=Fraction(0)), GAM1[:4])
check("control: on the ring N=4 (a path of four hops winds around it) E' commutes with M1 and K and not with M4",
      res == (True, True, False))
K23 = [(a, b, 1) for a in (0, 1) for b in (2, 3, 4)]                  # the complete bipartite graphs, squares inside
K24 = [(a, b, 1) for a in (0, 1) for b in (2, 3, 4, 5)]
for label, N, E in (("K2,3", 5, K23), ("K2,4", 6, K24)):
    res = ladder_commutes(N, ham_matrix(N, E, delta=Fraction(0)), GAM1[:N])
    check(f"squares: on {label}, H = hops only, E' commutes with M1 and K and not with M4", res == (True, True, False))
res = ladder_commutes(5, ham_correlated_hop(5, chain(5), 2), GAM1[:5])
check("control: a hop doubled when a third site is occupied (number conserving, SU(2) gone) breaks E' with M1, K and M4",
      res == (False, False, False))
res = ladder_commutes(5, ham_matrix(5, chain(5), delta=Fraction(1)), GAM1[:5])
check("control: with the ZZ term (Heisenberg chain N=5) E' commutes with M1 and K and not with M4", res == (True, True, False))

print("Z10 the diagonal terms enter M4 only as the detuning form, on graphs without triangles (exact)")
for label, N, E, H, gam in (
        ("chain N=5, Delta = 1/2, fields, next-nearest ZZ, a rate profile", 5, chain(5, [1, 2, 1, 3]),
         add_diagonal(5, ham_matrix(5, chain(5, [1, 2, 1, 3]), delta=Fraction(1, 2)), fields=(1, 0, -2, 1, 3), zz=((0, 2, 1),)),
         gam_r[:5]),
        ("ring N=5, Delta = 1", 5, ring(5), ham_matrix(5, ring(5), delta=Fraction(1)), GAM1[:5]),
        ("star N=5, Delta = 3/2, fields", 5, star(5), add_diagonal(5, ham_matrix(5, star(5), delta=Fraction(3, 2)), fields=(2, -1, 0, 1, 3)),
         GAM1[:5]),
        ("random tree N=6, Delta = 2", 6, random_tree(6), None, GAM1[:6])):
    if H is None:
        H = ham_matrix(N, E, delta=Fraction(2))
    ok = True
    for p in range(1, N):
        M4 = zeno_terms4(N, H, p, gam)[3]; M4h = zeno_terms4(N, off_diagonal(H), p, gam)[3]
        D = detuning_form(N, H, p, gam)
        ok &= all(M4[i][j] - M4h[i][j] == D[i][j] for i in range(len(D)) for j in range(len(D)))
    check(f"{label}: M4(H) - M4(hops of H) == D_V in every block", ok)
H = add_diagonal(5, ham_matrix(5, TRI, delta=Fraction(3, 2)), fields=(1, -2, 0, 3, 1))
M4 = zeno_terms4(5, H, 2, GAM1[:5])[3]; M4h = zeno_terms4(5, off_diagonal(H), 2, GAM1[:5])[3]
D = detuning_form(5, H, 2, GAM1[:5])
check("control: on the triangle graph (bonds 01, 12, 02, 23, 34) the cross terms of hop and diagonal survive, "
      "M4(H) - M4(hops) != D_V in block 2", any(M4[i][j] - M4h[i][j] != D[i][j] for i in range(len(D)) for j in range(len(D))))

print("Z11 the detuning weights W_b(p), read off F^T D_V F, in closed form (exact)")


def bond_weights(N, E, H, p, gam):
    """W_b(p) = sum over the hops across bond b of 2 |t|^2 dV^2 / G^3, read off F^T D_V F (a bond sum by construction:
    every hop crosses exactly one bond, and F maps the hop's Laplacian to that bond's)"""
    D = detuning_form(N, H, p, gam); F = lift_matrix(N, p)
    Q = mat_mul(mat_mul(transpose(F), D), F)
    return {(a, b): -Q[a][b] for a, b, _ in E}


def closed_W(kind, N, p, delta, a, b):
    """uniform book J = 1, gamma = 1: a detuned hop weighs 2 * 4 * dV^2 / 64 = dV^2/8"""
    if kind == "star":
        return Fraction(delta) ** 2 * (N - 2 * p) ** 2 * comb(N - 2, p - 1) / 2
    end = kind == "chain" and (a == 0 or b == N - 1)
    if end:
        return Fraction(delta) ** 2 * comb(N - 2, p - 1) / 2
    return 4 * Fraction(delta) ** 2 * (comb(N - 4, p - 2) if p >= 2 else 0)


for kind, Ns, builder in (("chain", (4, 5, 6, 7), chain), ("ring", (5, 6, 7), ring), ("star", (4, 5, 6, 7), star)):
    for N in Ns:
        for delta in (Fraction(1), Fraction(1, 2)):
            E = builder(N); H = ham_matrix(N, E, delta=delta); ok = True
            for p in range(1, N):
                W = bond_weights(N, E, H, p, GAM1[:N])
                ok &= all(W[(a, b)] == closed_W(kind, N, p, delta, a, b) for a, b, _ in E)
            check(f"{kind} N={N}, Delta = {delta}: W_b(p) equals the closed form on every bond, every block", ok)
def naive_W(kind, N, p, delta, a, b):
    """every hop across a bulk bond taken as detuned by 4 Delta: 2 Delta^2 C(N-2, p-1); the end bonds as closed_W"""
    end = kind == "chain" and (a == 0 or b == N - 1)
    return closed_W(kind, N, p, delta, a, b) if end else 2 * Fraction(delta) ** 2 * comb(N - 2, p - 1)


def weights_match(N, E, H, p, form, kind, delta):
    """the comparison of the rows above: W_b(p) read off F^T D_V F against a form, on every bond of block p"""
    W = bond_weights(N, E, H, p, GAM1[:N])
    return all(W[(a, b)] == form(kind, N, p, delta, a, b) for a, b, _ in E)


H6 = ham_matrix(6, chain(6), delta=Fraction(1))
check("control: the form with every hop across a bulk bond detuned, 2 Delta^2 C(N-2, p-1), misses the chain N=6 in "
      "every block", not any(weights_match(6, chain(6), H6, p, naive_W, "chain", 1) for p in range(1, 6)))
H6m = add_diagonal(6, ham_matrix(6, chain(6), delta=Fraction(0)), zz=((0, 1, 1), (1, 2, 1), (2, 3, -1), (3, 4, -1), (4, 5, 1)))
meets = [p for p in range(1, 6) if weights_match(6, chain(6), H6m, p, closed_W, "chain", 1)]
check(f"control: the chain N=6 with ZZ couplings (+1, +1, -1, -1, +1), its bulk hops detuned when their outer "
      f"neighbours agree, misses the uniform closed form in blocks 1, 3, 5 and meets it in blocks {meets}, where "
      f"2 C(2, p-1) + 2 C(2, p-3) = 4 C(2, p-2)", meets == [2, 4])

print("Z11c the split's shape on a chain with ZZ couplings of any sign and range: the filling enters only through phi(p) (exact)")


def mixed_chain(N, Ks, extra=()):
    """hops of the uniform XY chain, ZZ coupling K_b on bond (b, b+1), plus any extra diagonal ZZ terms"""
    return add_diagonal(N, ham_matrix(N, chain(N), delta=Fraction(0)),
                        zz=tuple((b, b + 1, K) for b, K in enumerate(Ks)) + tuple(extra))


def weights_have_the_shape(N, H, Ks):
    """W_b(p)/C(N-2, p-1) == ((K1 - K2)^2 (1 - phi) + (K1 + K2)^2 phi)/2 on a bulk bond, K^2/2 on an end bond"""
    for p in range(1, N):
        W = bond_weights(N, chain(N), H, p, GAM1[:N])
        phi = Fraction(2 * (p - 1) * (N - p - 1), (N - 2) * (N - 3))
        for b in range(N - 1):
            w = W[(b, b + 1)] / comb(N - 2, p - 1)
            if b == 0 or b == N - 2:
                expected = Fraction(Ks[1] if b == 0 else Ks[N - 3]) ** 2 / 2
            else:
                K1, K2 = Fraction(Ks[b - 1]), Fraction(Ks[b + 1])
                expected = ((K1 - K2) ** 2 * (1 - phi) + (K1 + K2) ** 2 * phi) / 2
            if w != expected:
                return False
    return True


for N, Ks in ((6, (1, 2, -1, -3, 1)), (7, (2, -1, 1, 3, -2, 1)), (6, (1, 1, 1, 1, 1)),
              (7, (Fraction(1, 2), -2, 1, 1, Fraction(-3, 2), 1))):
    check(f"chain N={N}, ZZ couplings ({', '.join(str(Fraction(k)) for k in Ks)}): the weights have the shape in every block",
          weights_have_the_shape(N, mixed_chain(N, Ks), Ks))
check("control: with a next-nearest ZZ term (0, 2) added to the chain N=6 the weights leave the shape",
      not weights_have_the_shape(6, mixed_chain(6, (1, 2, -1, -3, 1), extra=((0, 2, 1),)), (1, 2, -1, -3, 1)))


def weights_have_the_general_shape(N, H, K):
    """ZZ couplings K[(k, l)], k < l, of any range and no fields: W_b(p)/C(N-2, p-1) == ((sum c)^2 (1 - 2 phi)
    + 2 phi sum c^2)/2 on every bond (i, i+1), c_k = K_ik - K_(i+1)k over the other sites"""
    def k_of(a, b):
        return Fraction(K.get((min(a, b), max(a, b)), 0))
    for p in range(1, N):
        W = bond_weights(N, chain(N), H, p, GAM1[:N])
        phi = Fraction(2 * (p - 1) * (N - p - 1), (N - 2) * (N - 3))
        for i in range(N - 1):
            c = [k_of(i, k) - k_of(i + 1, k) for k in range(N) if k not in (i, i + 1)]
            s1, s2 = sum(c), sum(v * v for v in c)
            if W[(i, i + 1)] / comb(N - 2, p - 1) != (s1 * s1 * (1 - 2 * phi) + 2 * phi * s2) / 2:
                return False
    return True


LONGER = (("nearest (1, 2, -1, -3, 1) and the next-nearest (0, 2) 1", 6,
           {(0, 1): 1, (1, 2): 2, (2, 3): -1, (3, 4): -3, (4, 5): 1, (0, 2): 1}),
          ("nearest 1 and next-nearest 3/2 throughout", 7,
           {**{(b, b + 1): 1 for b in range(6)}, **{(b, b + 2): Fraction(3, 2) for b in range(5)}}),
          ("eleven couplings of ranges 1 to 6", 7,
           {(0, 1): 2, (1, 2): -1, (2, 3): 1, (3, 4): 3, (4, 5): -2, (5, 6): 1, (0, 2): Fraction(1, 2), (1, 3): -1,
            (2, 5): Fraction(3, 2), (0, 6): -1, (3, 6): 2}))
for label, N, K in LONGER:
    H = add_diagonal(N, ham_matrix(N, chain(N), delta=Fraction(0)), zz=tuple((a, b, v) for (a, b), v in K.items()))
    check(f"chain N={N}, ZZ couplings of longer range, {label}: the weights have the general shape in every block",
          weights_have_the_general_shape(N, H, K))
check("control: the uniform Heisenberg chain N=6 with a field 1 on site 0 misses the general shape (a field adds a "
      "part linear in the filling)",
      not weights_have_the_general_shape(6, add_diagonal(6, mixed_chain(6, (1, 1, 1, 1, 1)), fields=(1, 0, 0, 0, 0, 0)),
                                         {(b, b + 1): 1 for b in range(5)}))

print("Z11b the displayed delta_p of the chain, the ring and the star against the quotient on the exact weights")
import mpmath


def displayed_residuals(dps):
    """max |displayed - quotient| over N = 4..12 (rings from 5), every block, at J = gamma = Delta = 1, f the site
    Laplacian's lowest nonconstant mode in closed form; the second value is the ring with the chain's lambda_1"""
    mpmath.mp.dps = dps
    mpf, cos, sin, pi = mpmath.mpf, mpmath.cos, mpmath.sin, mpmath.pi
    worst, wrong = mpf(0), mpf(0)
    for N in range(4, 13):
        fc = [cos(pi * (l + mpf(1) / 2) / N) for l in range(N)]
        fr = [cos(2 * pi * l / N) for l in range(N)]
        fs = [mpf(0), mpf(1), mpf(-1)] + [mpf(0)] * (N - 3)
        for p in range(1, N):
            C = comb(N - 2, p - 1)
            Wc = [closed_W("chain", N, p, 1, a, a + 1) for a in range(N - 1)]
            q_chain = sum(mpf(Wc[a].numerator) / Wc[a].denominator * (fc[a + 1] - fc[a]) ** 2 for a in range(N - 1)) / (sum(v * v for v in fc) * C)
            lam_c = 4 * sin(pi / (2 * N)) ** 2
            shown_chain = (2 * lam_c / N) * (2 * (N - 4 * sin(pi / N) ** 2) * (p - 1) * (N - p - 1) / mpf((N - 2) * (N - 3)) + sin(pi / N) ** 2)
            Ws = closed_W("star", N, p, 1, 0, 1)
            q_star = sum(mpf(Ws.numerator) / Ws.denominator * (fs[0] - fs[l]) ** 2 for l in range(1, N)) / (sum(v * v for v in fs) * C)
            shown_star = mpf(1) * (N - 2 * p) ** 2 / 2
            worst = max(worst, abs(q_chain - shown_chain), abs(q_star - shown_star))
            if N >= 5:
                Wr = closed_W("ring", N, p, 1, 0, 1)
                q_ring = sum(mpf(Wr.numerator) / Wr.denominator * (fr[(a + 1) % N] - fr[a]) ** 2 for a in range(N)) / (sum(v * v for v in fr) * C)
                shown_ring = 4 * 4 * sin(pi / N) ** 2 * (p - 1) * (N - p - 1) / mpf((N - 2) * (N - 3))
                wrong_ring = 4 * lam_c * (p - 1) * (N - p - 1) / mpf((N - 2) * (N - 3))
                worst = max(worst, abs(q_ring - shown_ring))
                wrong = max(wrong, abs(q_ring - wrong_ring))
    return worst, wrong


r30, w30 = displayed_residuals(30)
r60, w60 = displayed_residuals(60)
mpmath.mp.dps = 15
check(f"the displayed forms equal the quotient: residual {mpmath.nstr(r30, 3)} at 30 digits, {mpmath.nstr(r60, 3)} at 60 "
      f"(an identity: the residual is the working precision's, falling by about 30 decades)",
      r30 < mpmath.mpf(10) ** -25 and r60 < mpmath.mpf(10) ** -55)
check(f"control: the ring's form with the chain's lambda_1 = 4 sin^2(pi/(2N)) in place of 4 sin^2(pi/N) misses, by "
      f"{mpmath.nstr(w60, 3)} at 60 digits", w30 > mpmath.mpf(10) ** -2 and w60 > mpmath.mpf(10) ** -2)

print("Z12 the filling split at order x^4 is the detuning form's difference (float; the law of the deviation is x^2)")


def slowest_per_block(N, H, x, with_scale=False):
    """the slowest non-stationary level of every block; with_scale also returns the largest block 1-norm, the
    eigensolver's rounding being of order eps times it"""
    out, scale = [], 0.0
    for p in range(1, N):
        A_f, ham_f = float_block(N, H, p)
        B = -1j * x * A_f + np.diag(-2.0 * ham_f)
        scale = max(scale, np.abs(B).sum(axis=0).max())
        ev = np.linalg.eigvals(B)
        ev = ev[np.abs(ev) > 1e-11]
        out.append(ev[np.argmin(np.abs(ev.real))].real)
    return (np.array(out), scale) if with_scale else np.array(out)


def lifted_quotient(N, E, W, p):
    """sum_b W_b(p) (f_a - f_b)^2 / (|f|^2 C(N-2, p-1)) on the lowest nonconstant site mode f (any one, if degenerate)"""
    w, V = np.linalg.eigh(np.array(site_laplacian(N, E), float))
    f = V[:, 1]
    return sum(float(W[(a, b)]) * (f[a] - f[b]) ** 2 for a, b, _ in E) / (f @ f * comb(N - 2, p - 1))


XS12 = (0.01, 0.02, 0.04)
for label, kind, N, builder, delta in (("chain", "chain", 5, chain, 1), ("chain", "chain", 6, chain, 1),
                                       ("chain", "chain", 6, chain, Fraction(1, 2)), ("ring", "ring", 6, ring, 1),
                                       ("star", "star", 5, star, 1), ("star", "star", 6, star, 1)):
    E = builder(N); H = ham_matrix(N, E, delta=Fraction(delta))
    pred = np.array([lifted_quotient(N, E, {(a, b): closed_W(kind, N, p, delta, a, b) for a, b, _ in E}, p)
                     for p in range(1, N)])
    pred_split = pred - pred[0]
    devs, wrong = [], []
    for x in XS12:
        r = slowest_per_block(N, H, x)
        d = (r - r[0]) / x ** 4
        devs.append(np.max(np.abs(d - pred_split)))
        wrong.append(np.max(np.abs(d - pred_split * (1 + 1e-2))))
    ratios = [devs[1] / devs[0], devs[2] / devs[1]]
    wr = [wrong[1] / wrong[0], wrong[2] / wrong[1]]
    check(f"{kind} N={N}, Delta = {delta}: (block p's slowest level less block 1's)/x^4 -> {np.round(pred_split, 5).tolist()} "
          f"for p = 1..{N - 1}, the deviation growing by {ratios[0]:.2f}, {ratios[1]:.2f} per doubling of x (x^2: 4, decided "
          f"at the half-integer powers; control: the prediction scaled by 1 + 1e-2, {wr[0]:.2f}, {wr[1]:.2f}, read outside "
          f"the window both times)",
          all(2 ** 1.5 < q < 2 ** 2.5 for q in ratios) and not any(2 ** 1.5 < q < 2 ** 2.5 for q in wr))

for Ks, extra in (((1, 2, -1, -3, 1), ()), ((1, 1, 1, 1, 1), tuple((b, b + 2, Fraction(3, 2)) for b in range(4)))):
    N = 6
    E = chain(N); H = mixed_chain(N, Ks, extra=extra)
    pred = np.array([lifted_quotient(N, E, bond_weights(N, E, H, p, GAM1[:N]), p) for p in range(1, N)])
    pred_split = pred - pred[0]
    devs, wrong = [], []
    for x in XS12:
        r = slowest_per_block(N, H, x)
        d = (r - r[0]) / x ** 4
        devs.append(np.max(np.abs(d - pred_split)))
        wrong.append(np.max(np.abs(d - pred_split * (1 + 1e-2))))
    ratios = [devs[1] / devs[0], devs[2] / devs[1]]
    wr = [wrong[1] / wrong[0], wrong[2] / wrong[1]]
    check(f"chain N=6, ZZ couplings {Ks}{' and next-nearest 3/2' if extra else ''}, exact weights: the split -> "
          f"{np.round(pred_split, 5).tolist()}, the deviation "
          f"growing by {ratios[0]:.2f}, {ratios[1]:.2f} per doubling of x (control {wr[0]:.2f}, {wr[1]:.2f} outside the window)",
          all(2 ** 1.5 < q < 2 ** 2.5 for q in ratios) and not any(2 ** 1.5 < q < 2 ** 2.5 for q in wr))

for kind, N, builder, delta in (("chain", 5, chain, 1), ("chain", 6, chain, Fraction(1, 2)), ("ring", 6, ring, 1),
                                ("star", 5, star, 1)):
    E = builder(N); H = ham_matrix(N, E, delta=Fraction(delta)); H0 = off_diagonal(H)
    pred = np.array([lifted_quotient(N, E, {(a, b): closed_W(kind, N, p, delta, a, b) for a, b, _ in E}, p)
                     for p in range(1, N)])
    devs, wrong = [], []
    for x in XS12:
        d = (slowest_per_block(N, H, x) - slowest_per_block(N, H0, x)) / x ** 4
        devs.append(np.max(np.abs(d - pred)))
        wrong.append(np.max(np.abs(d - pred * (1 + 1e-2))))
    ratios = [devs[1] / devs[0], devs[2] / devs[1]]
    wr = [wrong[1] / wrong[0], wrong[2] / wrong[1]]
    check(f"{kind} N={N}, Delta = {delta}: delta_p itself, (block p's slowest level less the same block's without the "
          f"ZZ term)/x^4 -> {np.round(pred, 5).tolist()} for p = 1..{N - 1}, the deviation growing by {ratios[0]:.2f}, "
          f"{ratios[1]:.2f} per doubling of x (control {wr[0]:.2f}, {wr[1]:.2f} outside the window)",
          all(2 ** 1.5 < q < 2 ** 2.5 for q in ratios) and not any(2 ** 1.5 < q < 2 ** 2.5 for q in wr))

LADDER23 = [(0, 1, 1), (1, 2, 1), (3, 4, 1), (4, 5, 1), (0, 3, 1), (1, 4, 1), (2, 5, 1)]   # two rails of 3, three rungs


def ham_star_free_fermion(N):
    """the XY star with the Jordan-Wigner sign on every hop between the hub 0 and an arm l, (-1) to the number of
    occupied arms between them in the site order: free fermions on the star's graph"""
    d = 2 ** N
    H = [dict() for _ in range(d)]
    for l in range(1, N):
        for x in range(d):
            if bit(N, x, 0) != bit(N, x, l):
                y = x ^ (1 << (N - 1)) ^ (1 << (N - 1 - l))
                sign = (-1) ** sum(bit(N, x, k) for k in range(1, l))
                H[y][x] = H[y].get(x, 0) + 2 * sign
    return H


print("Z13 reading, not gated: which block holds the slowest non-stationary rate at x = J/gamma = 0.02, and the split's order")
for label, N, E, delta in (("Heisenberg chain", 4, chain(4), 1), ("Heisenberg chain", 5, chain(5), 1),
                           ("Heisenberg chain", 6, chain(6), 1), ("XXZ chain Delta = 1/2", 6, chain(6), Fraction(1, 2)),
                           ("XY chain", 6, chain(6), 0), ("Heisenberg ring", 6, ring(6), 1),
                           ("Heisenberg star", 6, star(6), 1), ("XY star", 6, star(6), 0),
                           ("XY star", 4, star(4), 0), ("XY star", 5, star(5), 0), ("XY star", 7, star(7), 0),
                           ("XY ring", 4, ring(4), 0), ("XY ring", 6, ring(6), 0),
                           ("XY ring", 5, ring(5), 0), ("XY ring", 7, ring(7), 0),
                           ("XY K2,3", 5, K23, 0), ("XY K2,4", 6, K24, 0),
                           ("XY ladder 2x3", 6, LADDER23, 0), ("Heisenberg ladder 2x3", 6, LADDER23, 1),
                           ("XY star, Jordan-Wigner-signed hops (free fermions)", 6, None, 0),
                           ("chain, ZZ couplings (+1, +1, -1, -1, +1)", 6, chain(6), None)):
    if E is None:
        H = ham_star_free_fermion(N)
    elif delta is None:
        H = add_diagonal(N, ham_matrix(N, E, delta=Fraction(0)), zz=((0, 1, 1), (1, 2, 1), (2, 3, -1), (3, 4, -1), (4, 5, 1)))
    else:
        H = ham_matrix(N, E, delta=Fraction(delta))
    (r1, s1), r2 = slowest_per_block(N, H, 0.02, with_scale=True), slowest_per_block(N, H, 0.04)
    tie = 64 * np.finfo(float).eps * s1                      # the eigensolver's rounding, eps times the block norm
    spread1, spread2 = r1.max() - r1.min(), r2.max() - r2.min()
    top = [p + 1 for p in range(N - 1) if r1[p] >= r1.max() - tie] if spread1 > tie else "every block (spread at rounding)"
    order = f"spread x {spread2 / spread1:.1f} per doubling of x" if spread1 > tie else "no split"
    split = np.round((r1 - r1[0]) / 0.02 ** 4, 3).tolist() if spread1 > tie else "at rounding"
    print(f"    {label} N={N}: slowest block(s) {top}; {order}; (block p's slowest level less block 1's)/x^4 at x = 0.02: {split}")


for label, N, delta, fields in (("XY chain with fields (1, -2, 1/2, 3, -1, 2)", 6, 0, (1, -2, Fraction(1, 2), 3, -1, 2)),
                                ("Heisenberg chain with a field 4 on both end sites", 6, 1, (4, 0, 0, 0, 0, 4))):
    H = add_diagonal(N, ham_matrix(N, chain(N), delta=Fraction(delta)), fields=fields)
    (r1, s1), r2 = slowest_per_block(N, H, 0.02, with_scale=True), slowest_per_block(N, H, 0.04)
    tie = 64 * np.finfo(float).eps * s1
    spread1, spread2 = r1.max() - r1.min(), r2.max() - r2.min()
    top = [p + 1 for p in range(N - 1) if r1[p] >= r1.max() - tie] if spread1 > tie else "every block (spread at rounding)"
    order = f"spread x {spread2 / spread1:.1f} per doubling of x" if spread1 > tie else "no split"
    split = np.round((r1 - r1[0]) / 0.02 ** 4, 3).tolist() if spread1 > tie else "at rounding"
    print(f"    {label} N={N}: slowest block(s) {top}; {order}; (block p's slowest level less block 1's)/x^4 at x = 0.02: {split}")


def slowest_rate_per_block_profile(N, H, gam, x):
    """-Re of the slowest non-stationary eigenvalue of each block (p, p), B = -i x A - diag(2 sum gamma over the
    disagreeing sites), the rates in units of the hop's J; also the largest block 1-norm"""
    out, scale = [], 0.0
    for p in range(1, N):
        A_f, _ = float_block(N, H, p)
        cf = confs_of(N, p)
        rate = np.array([2 * sum(gam[l] for l in range(N) if bit(N, a, l) != bit(N, b, l)) for a in cf for b in cf])
        B = -1j * x * A_f - np.diag(rate)
        scale = max(scale, np.abs(B).sum(axis=0).max())
        ev = np.linalg.eigvals(B)
        ev = ev[np.abs(ev) > 1e-11]
        out.append(-ev[np.argmin(np.abs(ev.real))].real)
    return np.array(out), scale


# experiments/GAMMA_AS_BINDING.md's profile: rates proportional to 1/T2 of IBM Torino, the smallest 0.05, at J = 1
T2_TORINO = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
g_ibm = (1 / (2 * T2_TORINO)) / (1 / (2 * T2_TORINO)).min() * 0.05
H5 = ham_matrix(5, chain(5), delta=Fraction(1))
own, own_scale = slowest_rate_per_block_profile(5, H5, g_ibm, 1.0)
zeno, zeno_scale = slowest_rate_per_block_profile(5, H5, g_ibm / g_ibm.min(), 0.02)


def slowest_blocks(r, scale):
    """the blocks within the eigensolver's rounding (64 eps times the block norm) of the slowest rate"""
    return [p + 1 for p in range(len(r)) if r[p] <= r.min() + 64 * np.finfo(float).eps * scale]


def slowest_over_every_block(N, H, gam, x):
    """the slowest non-stationary eigenvalue of every block (p, q), cells |a><b| with popcount a = p and b = q,
    B = -i x (H_p (x) 1 - 1 (x) H_q^T) - diag(2 sum of gamma_l over the sites where a and b differ): the list of
    (rate, |Im|, (p, q)), slowest first, and the largest block 1-norm"""
    out, scale = [], 0.0
    for p in range(N + 1):
        for q in range(N + 1):
            ca, cb = confs_of(N, p), confs_of(N, q)
            Hp = np.array([[complex(H[r].get(s, 0)) for s in ca] for r in ca])
            Hq = np.array([[complex(H[r].get(s, 0)) for s in cb] for r in cb])
            rate = [2.0 * sum(float(gam[l]) for l in range(N) if bit(N, a, l) != bit(N, b, l)) for a in ca for b in cb]
            B = -1j * x * (np.kron(Hp, np.eye(len(cb))) - np.kron(np.eye(len(ca)), Hq.T)) - np.diag(rate)
            scale = max(scale, np.abs(B).sum(axis=0).max())
            ev = np.linalg.eigvals(B)
            ev = ev[np.abs(ev) > 1e-11]
            if len(ev):
                k = np.argmax(ev.real)
                out.append((-ev[k].real, abs(ev[k].imag), (p, q)))
    return sorted(out), scale


every, every_scale = slowest_over_every_block(5, H5, g_ibm, 1.0)
tie_every = 64 * np.finfo(float).eps * every_scale
first_blocks = [b for r, _, b in every if r <= every[0][0] + tie_every]
second_rate = min(r for r, _, _ in every if r > every[0][0] + tie_every)
second_blocks = [b for r, _, b in every if abs(r - second_rate) <= tie_every]
uniform = {g: slowest_blocks(*slowest_rate_per_block_profile(5, H5, [g] * 5, 1.0)) for g in (0.05, 0.5, 2.0)}
print(f"    Heisenberg chain N=5 under GAMMA_AS_BINDING's profile, rates {np.round(g_ibm, 4).tolist()} against J = 1: "
      f"slowest rate per block p = 1..4 {np.round(own, 3).tolist()}, slowest block(s) {slowest_blocks(own, own_scale)}; "
      f"the same shape at the Zeno end, x = J/gamma_min = 0.02: rate/x^2 {np.round(zeno / 0.02 ** 2, 5).tolist()}, "
      f"slowest block(s) {slowest_blocks(zeno, zeno_scale)}; the slowest rate over every block (p, q) at these "
      f"rates {every[0][0]:.4f}, |Im| {every[0][1]:.4f}, in the blocks {first_blocks}, the next {second_rate:.4f} in "
      f"{second_blocks}; at uniform rates 0.05, 0.5, 2 the slowest diagonal block(s) "
      f"{uniform[0.05]}, {uniform[0.5]}, {uniform[2.0]}")

print("Z14 the odd ring: a pi flux through it changes no block's spectrum (exact), and its fillings at finite Q (reading)")


def ring_flux(N):
    """the ring with its wrapped bond negated: a pi flux through it"""
    return ring(N)[:-1] + [(N - 1, 0, -1)]


def flux_map(N, p):
    """K P_T A(ring) P_T K in block p, and A(ring): K = prod_{l odd} Z_l acts on |a><b| by k(a) k(b), P_T is the
    transpose |a><b| -> |b><a|; A has integer entries and both maps are signed permutations, so nothing rounds"""
    A0, _ = float_block(N, ham_matrix(N, ring(N), delta=Fraction(0)), p)
    cf = confs_of(N, p)
    m = len(cf)
    k = np.array([(-1) ** sum(bit(N, a, l) for l in range(1, N, 2)) for a in cf])
    sgn = np.outer(k, k).reshape(-1)
    perm = np.array([j * m + i for i in range(m) for j in range(m)])
    return A0[np.ix_(perm, perm)] * np.outer(sgn, sgn), A0


for N in (3, 5, 7):
    ok = True
    for p in range(1, N):
        M, _ = flux_map(N, p)
        ok &= np.array_equal(M, float_block(N, ham_matrix(N, ring_flux(N), delta=Fraction(0)), p)[0])
    check(f"ring N={N}: K and the transpose carry A(ring) onto A(ring with a pi flux) in every block, exactly", ok)
for N in (4, 6):
    ok = True
    for p in range(1, N):
        M, A0 = flux_map(N, p)
        ok &= np.array_equal(M, A0) and not np.array_equal(M, float_block(N, ham_matrix(N, ring_flux(N), delta=Fraction(0)), p)[0])
    check(f"control: on the even ring N={N} the same map returns A(ring) itself in every block and misses the flux", ok)
for N, xs in ((3, (0.5, 1.0, 3.0)), (5, (0.5, 1.0, 3.0)), (7, (1.0,)), (4, (0.5, 1.0, 3.0)), (6, (0.5, 1.0, 3.0))):
    H = ham_matrix(N, ring(N), delta=Fraction(0))
    out = []
    for x in xs:
        r, s = slowest_per_block(N, H, x, with_scale=True)
        tie = 64 * np.finfo(float).eps * s
        spread = r.max() - r.min()
        if spread <= tie:
            out.append(f"x = {x}: every block tied (spread {spread:.0e}, rounding {tie:.0e})")
        else:
            even, odd = r[1::2], r[0::2]                 # blocks 2, 4, ... and 1, 3, ...
            out.append(f"x = {x}: the even blocks within {even.max() - even.min():.0e} of each other, the odd ones within "
                       f"{odd.max() - odd.min():.0e}, the classes {abs(even.max() - odd.max()):.3g} apart, slowest block(s) "
                       f"{[p + 1 for p in range(N - 1) if r[p] >= r.max() - tie]}")
    print(f"    XY ring N={N}: " + "; ".join(out))
H8 = ham_matrix(8, ring(8), delta=Fraction(0))
H8d = np.zeros((256, 256))
for r_ in range(256):
    for s_, v_ in H8[r_].items():
        H8d[r_, s_] = float(v_)
r8, scale8 = [], 0.0
for p in range(1, 5):
    cf = confs_of(8, p)
    Hs, I_ = H8d[np.ix_(cf, cf)], np.eye(len(cf))                       # vectorised: the block (p, p) as in float_block
    A_f = np.kron(Hs, I_) - np.kron(I_, Hs.T)
    ham_f = np.array([bin(a ^ b).count("1") for a in cf for b in cf], dtype=float)
    B8 = -1j * 0.3 * A_f + np.diag(-2.0 * ham_f)
    scale8 = max(scale8, np.abs(B8).sum(axis=0).max())
    ev = np.linalg.eigvals(B8)
    ev = ev[np.abs(ev) > 1e-11]
    r8.append(ev[np.argmin(np.abs(ev.real))].real)
tie8 = 64 * np.finfo(float).eps * scale8
def same8(a, b):
    return "share one level" if abs(r8[a] - r8[b]) <= tie8 else f"differ by {abs(r8[a] - r8[b]):.1e}"
print("    XY ring N=8, x = 0.3, slowest level of blocks 1..4 (5..7 their X^N images): " + ", ".join(f"{v:.12f}" for v in r8) +
      f"; blocks 2 and 4, the half filling, {same8(1, 3)}; blocks 1 and 3 {same8(0, 2)}; the classes {abs(r8[1] - r8[0]):.2e} "
      f"apart (rounding {tie8:.0e}), the slower {'the even one' if r8[1] > r8[0] else 'the odd one'}")

print("Z15 reading, not gated: the stars' survivor in block (1, 1) across Q = J/gamma, the lifted wave or the hub-against-arms mode")
QS15 = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 20.0)
for label, delta in (("Heisenberg", 1), ("XY", 0)):
    for N in (4, 5, 6, 7, 8):
        A, ham = float_block(N, ham_matrix(N, star(N), delta=Fraction(delta)), 1)
        cf = confs_of(N, 1)
        m = len(cf)
        hub = cf.index(1 << (N - 1))
        kinds, last, im_max = [], None, 0.0
        for Q in QS15:
            ev, V = np.linalg.eig(-1j * Q * A - 2.0 * np.diag(ham.astype(float)))
            cand = [k for k in range(len(ev)) if abs(ev[k]) > 1e-9]
            k = max(cand, key=lambda k: ev[k].real)
            pops = np.array([V[i * m + i, k] for i in range(m)])
            if np.linalg.norm(pops) <= 1e-9 * np.linalg.norm(V[:, k]):
                kinds.append("no populations")
            else:
                kinds.append("wave" if abs(pops[hub]) <= 1e-9 * np.abs(pops).max() else "hub")
            last = ev[k]
            im_max = max(im_max, abs(ev[k].imag))
        print(f"    {label} star N={N}: " + ", ".join(f"Q = {Q:g} {kd}" for Q, kd in zip(QS15, kinds)) +
              f"; <n_XY> at Q = 20: {-last.real / 2:.4f} (4/N = {4 / N:.4f}, 4/(N - 1) = {4 / (N - 1):.4f}); the "
              f"slowest mode's largest |Im| over these Q {im_max:.1e}")
A4, ham4 = float_block(4, ham_matrix(4, ring(4), delta=Fraction(1)), 2)
cf4 = confs_of(4, 2)
walls = np.array([sum(bit(4, c, a) != bit(4, c, b) for a, b, _ in ring(4)) for c in cf4], float)
walls -= walls.mean()
rows4 = []
for Q in (0.3, 0.6, 1.0, 2.0, 4.0, 10.0, 30.0, 100.0):
    ev, V = np.linalg.eig(-1j * Q * A4 - 2.0 * np.diag(ham4.astype(float)))
    cand = [k for k in range(len(ev)) if abs(ev[k]) > 1e-9]
    k = max(cand, key=lambda k: ev[k].real)
    pops = np.array([V[i * len(cf4) + i, k] for i in range(len(cf4))])
    cos = abs(np.vdot(walls, pops)) / (np.linalg.norm(walls) * np.linalg.norm(pops))
    below = sum(1 for kk in cand if -ev[kk].real / 2 < 1)
    rows4.append(f"Q = {Q:g}: {below} below the plane, the slowest at height {-ev[k].real / 2:.4f}, its populations' "
                 f"|cos| with the wall count {cos:.3f}")
print("    Heisenberg ring N=4, block (2, 2): " + "; ".join(rows4))

print()
if FAILS:
    print("FAILED:", *FAILS, sep="\n  ")
    sys.exit(1)
print("ALL GATES PASS")
