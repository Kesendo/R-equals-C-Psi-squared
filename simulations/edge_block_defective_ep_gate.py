"""The edge block carries a defective EP under a per-site rate profile, and on the XY chain
the EP set has codimension one.

WHAT THIS SETTLES. The FIRST item of the (b2) inventory in the arc `site_resolved_vacuum_block`
(compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs:4118) says, verbatim, "'no defective EP
can live on an edge block' is no longer forbidden a priori ... Settle it before fencing, do not
assume either way." Settled in the positive: the EP exists, at real positive rates, at the
canonical base rate gamma_0 = 0.05, and at the sizes where the lemma that forbids it is used.
The COUPLING is NOT a canonical anchor, and the comment at CANONICAL_J says why: J = 0.075 in
this file's book is hop 0.15, carrier Q = 3.0, past every entry of docs/Q_REGIME_ANCHORS.md.
G12 carries the conversion and the canonical Q = 1.5 point beside it.

TWO OPERATORS, AND THE FIRST VERSION OF THIS GATE MEASURED THE WRONG ONE. On the (0,1)
coherence block, spanned by |0><j|:

    Delta = 1 (XXX)  M = -2i*(D_J - A_J) - 2*diag(gamma)     the Laplacian form, F152's
    Delta = 0 (XY)   M = +2i*A_J        - 2*diag(gamma)      the adjacency form

`PROOF_CODIM1_BY_ADDITIVITY.md` and F125, which carry the Edge lemma, are scoped to the **XY
chain, Delta = 0**. F152 is fenced to |Delta| = 1 and says so. A first version of this gate ran
every witness on the Laplacian form and was 81/81 green while measuring an operator the refuted
lemma is not about; on the adjacency form those same profiles are not EPs at all, splits of
order 0.1 to 1.7. Those three splits (1.54, 1.67, 0.107) were measured on the superseded
fixture, which this file no longer carries: they are the size of the miss, not numbers
anything here reproduces. Both forms are gated here, and every WITNESS is on Delta = 0.

THE THREE PARTS OF THE ANSWER.

  A. DEFECTIVENESS IS STRUCTURAL, EXACT, AND Delta-BLIND (G2, G3). On a PATH with every
     J_b nonzero, M is tridiagonal with off-diagonal entries that are nonzero for either
     Delta. Delete row 0 and column N-1 from (M - lam*I): what is left is triangular with
     those off-diagonals on its diagonal, so its determinant is their product, free of lam
     and nonzero. Hence rank(M - lam*I) >= N-1 for EVERY lam: M is non-derogatory, every
     eigenvalue has geometric multiplicity one, and ANY repeated eigenvalue is therefore
     defective. No defective-versus-diabolic instrument is invoked anywhere in this file,
     because there is no diabolic alternative left to exclude.

  B. ON THE XY CHAIN THE DISCRIMINANT IS REAL, SO THE EP SET HAS CODIMENSION ONE (G5, G6).
     A path is bipartite, so with S = diag((-1)^k) one has S*conj(M)*S^-1 = M exactly at
     Delta = 0. M is similar to its own conjugate, the characteristic polynomial has real
     coefficients, and the discriminant is real. disc = 0 is then ONE real equation in the
     rate space rather than two, so its solution set is a HYPERSURFACE of dimension N-1 in
     the rate space and turning a single rate crosses it. Existence follows from a sign
     change, by the intermediate value theorem on a real continuous function; no topological
     degree is needed and none is used. SCOPE, because the reality argument does not reach
     further: a real characteristic polynomial makes disc = 0 one equation for a coalescence
     ON THE REAL AXIS, and every EP this route produces has Im lambda = 0 (G6 asserts the
     computed residual under the bound eps*||M||_2*n, and G6c re-asserts that bound at
     every witness across four decades of scale; no raw magnitude is quoted, because a raw
     magnitude moves with the units, and no FLATNESS is claimed, because the residual is a
     cancellation and its spread beneath the bound is the object, not a defect). EPs
     between eigenvalues that are not each other's conjugates still come in conjugate pairs
     and stay codimension two. The claim is about the real-lambda stratum.
     At Delta = 1 the degree diagonal D breaks that antisymmetry (S*D*S^-1 = D, not -D) and
     the characteristic polynomial is genuinely complex, so nothing reduces the count of
     equations and the EPs sit at codimension two. G5 gates the similarity and its Delta=1
     failure; the complex-polynomial half is a consequence, not a separate measurement.

     THE IDENTITY ITSELF IS NOT NEW AND IS NOT CLAIMED AS NEW. It is ChiralKClaim (Tier 1,
     compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs), K*H*K = -H for the same
     K = diag((-1)^l) on the single-excitation site basis at any bond profile, together with
     the real rate diagonal; hypotheses/INHERITED_RULES_AND_THE_OWN.md:85 records the same
     pairing and shelves it BECAUSE the Heisenberg generator keeps the ZZ degree diagonal,
     which is the Delta=1 book. What is new is the step from the identity to codimension.

  C. IT IS THE RATE PROFILE, AND NOT ANY PROFILE (G9). Lemma A applies to the Hermitian
     Jacobi matrix A_J as well, so a path's hopping spectrum is simple for ANY bond profile.
     At uniform gamma the block therefore has a simple spectrum and nothing to be defective;
     a J profile does not change that. The per-site RATE profile is what creates the
     degeneracy, and part A then leaves it no form but a Jordan block.

WHAT SURVIVES INTACT. The window-edge lemma (PROOF_CODIM1_BY_ADDITIVITY.md:121) forbids a
defective eigenvalue only AT an edge of the rate window, and needs only that the Hermitian
part be Hermitian. Every EP located here sits strictly inside [-2*gamma_max, -2*gamma_min]
(G7), which is where that lemma permits one. Containment is a theorem, not a measurement; G7
corroborates it and does not establish it.

SCOPE. Part A is about PATHS: it uses tridiagonality. G10 is the control, and it is the
control that isolates the right hypothesis rather than a convenient one: on a STAR with a
rate PROFILE, M is not normal and yet the repeated eigenvalue is semisimple. So non-normality
does not force defectiveness. What tridiagonality forces is the CONDITIONAL, repeated implies
defective, and nothing more: at uniform gamma a path block has a simple spectrum and there is
nothing defective on it at all (part C). Note also that the star is bipartite too, so part B's
reality is not path-specific; only part A is.

Run: python simulations/edge_block_defective_ep_gate.py
"""
from __future__ import annotations

import sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.path.insert(0, "simulations")
import d10_block_closure_verify as d10

EPS = 2.220446049250313e-16

PASSED = 0
FAILED = 0


def gate(name, ok, detail=""):
    global PASSED, FAILED
    if ok:
        PASSED += 1
        print(f"  [PASS] {name}   {detail}")
    else:
        FAILED += 1
        print(f"  [FAIL] {name}   {detail}")


# --------------------------------------------------------------------------------------
# the object, with Delta as a parameter


def block_M(gamma, bonds_J, n, delta):
    """delta = 0 gives the XY (adjacency) form, delta = 1 the XXX (Laplacian) form.

    The adjacency accumulates rather than assigns, so a repeated bond adds instead of
    overwriting; a first version assigned, which would have built a non-Laplacian on a
    multigraph while no fixture in this file could have noticed."""
    adj = np.zeros((n, n))
    deg = np.zeros((n, n))
    for (u, v), j in bonds_J:
        adj[u, v] += j
        adj[v, u] += j
        deg[u, u] += j
        deg[v, v] += j
    return -2j * (delta * deg - adj) - 2 * np.diag(gamma)


def chain_bonds(n, j):
    return [((b, b + 1), j) for b in range(n - 1)]


def discriminant(M):
    w = np.linalg.eigvals(M)
    n = len(w)
    p = 1.0 + 0j
    for i in range(n):
        for j in range(i + 1, n):
            p *= (w[i] - w[j]) ** 2
    return p


def closest_pair(M):
    w = np.linalg.eigvals(M)
    n = len(w)
    i0, j0 = min(((i, j) for i in range(n) for j in range(i + 1, n)),
                 key=lambda p: abs(w[p[0]] - w[p[1]]))
    return (w[i0] + w[j0]) / 2, abs(w[i0] - w[j0])


# The witnesses. Each is a BRACKET on one rate, not a point: the certificate is a sign change
# of the real discriminant across it, and the root is found by the gate rather than quoted.
# The coupling, and WHICH BOOK it is in, because the anchor file does not pin one. This file
# writes H = sum_b J_b*(XX+YY), hop element 2J: the "doubled" book, q in docs/GLOSSARY.md's
# factor-2 table. docs/Q_REGIME_ANCHORS.md reads Q = J/gamma_0 in the OTHER book (hop J, the
# ChainSystem/carrier normalization), so J = 0.075 HERE is hop 0.15 and the carrier ratio is
# Q = 0.15/0.05 = 3.0, not the anchor's 1.5. G12 carries the conversion and the canonical
# Q = 1.5 point beside it.
CANONICAL_J = 0.075        # hop 2J = 0.15, i.e. carrier Q = 3.0 at gamma_0 = 0.05
EXACT_J = sp.Rational(75, 1000)
CANONICAL_G0 = 0.05        # the canonical base rate


def build_exact(build, x):
    """The same gamma profile over the rationals: every entry is an exact decimal."""
    return [v if isinstance(v, sp.Basic) else sp.Rational(str(v))
            for v in build(sp.Rational(str(x)))]


WITNESSES = [
    ("N=4", 4, lambda x: [0.05, x, 0.07, 0.09], (0.220915, 0.221873)),
    ("N=5 (a)", 5, lambda x: [0.05, x, 0.343164, 0.07, 0.09], (0.169192, 0.170150)),
    ("N=5 (b)", 5, lambda x: [0.05, x, 0.343164, 0.07, 0.09], (0.195054, 0.196011)),
    ("N=5 (c)", 5, lambda x: [0.05, x, 0.343164, 0.07, 0.09], (0.492938, 0.493896)),
    ("N=6 (a)", 6, lambda x: [0.05, x, 0.191302, 0.07, 0.09, 0.11], (0.229535, 0.230493)),
    ("N=6 (b)", 6, lambda x: [0.05, x, 0.191302, 0.07, 0.09, 0.11], (0.342559, 0.343517)),
    ("N=6 (c)", 6, lambda x: [0.05, x, 0.191302, 0.07, 0.09, 0.11], (0.846376, 0.847334)),
]


def disc_real(build, x, n):
    return discriminant(block_M(build(x), chain_bonds(n, CANONICAL_J), n, 0)).real


def ep_root(build, n, bracket):
    return brentq(lambda x: disc_real(build, x, n), bracket[0], bracket[1],
                  xtol=1e-15, rtol=8.9e-16)


# --------------------------------------------------------------------------------------


def g1_anchor():
    """No exact route: the reference contracts the full 4^N Liouvillian while block_M writes
    the closed form directly, so the two differ at rounding. The tolerance is a LAW, not a
    number. The model is

        residual  <=  C * eps * ||block||_max * N,

    and what is asserted is a BOUND on C, per configuration. Nothing else: an earlier version
    of this line claimed the gate also showed C does not grow, and that claim and the gate
    that carried it are both gone, for the reason below.

    WHAT IS ASSERTED IS THE BOUND ON C, AND NOTHING ABOUT THE SHAPE OF THE SIZE FACTOR.
    Three things stood in this docstring across earlier rounds and all three are gone, for
    one reason: the N=3 and N=7 maxima are each the largest of four rounding residuals, so
    every verdict built on their RATIO moves with the draw. Over 200 draws of these same 20
    configurations that ratio runs to 4.57 with a median of 1.07, and a non-growth test at
    2.0 fails on 14 of them, roughly one draw in fourteen.
    That killed, in order: a SPREAD test (max/min, the wrong statistic for an upper bound,
    since a configuration whose entries happen to be exactly representable makes the residual
    SMALLER than the model and that violates nothing); then two verdicts on competing models
    (size-free REJECTED, sqrt(N) SURVIVES), which are the same ratio against two other
    thresholds; then the non-growth test itself, which is that ratio against a third. All
    three were one number wearing three hats.
    What survives is the per-configuration bound, and its constant is MEASURED, but only the
    STABLE statistics are quoted. Over thousands of draws of these same 20 configurations the
    per-draw worst C has median 0.295 and 99th percentile near 0.47, both reproducing across
    independent sweeps; the gated 1.0 sits about twice the 99th percentile. The sample MAXIMUM
    is deliberately not quoted as headroom: it does not stabilise (0.56 at 2000 draws, 0.70 at
    2200, 0.60 at 6000), and a first version of this paragraph did quote it, which is the same
    draw-dependent species the paragraph above spends itself killing.
    An earlier version gated 0.5, which reads well against the shipped seed's 0.284 and fails
    on about one draw in 250: the passes-on-this-fixture defect one level down, and it took a
    sweep twenty times larger than the one that "confirmed" it to see.
    The linear N is the block dimension and is written for that structural reason; nothing
    here picks it out over any other shape."""
    print("G1  both forms are the block, read off the full 4^N Liouvillian")
    rng = np.random.default_rng(20260820)
    ratios = []
    for delta, letters, name in ((0, (d10.X, d10.Y), "XY  "), (1, (d10.X, d10.Y, d10.Z), "XXX ")):
        for n in range(3, 8):
            for j in (1.0, CANONICAL_J):
                gamma = rng.uniform(0.01, 2.5, n)
                ref, leak = d10.measure_block(n, d10.chain(n), gamma, j=j, letters=letters)
                mine = block_M(gamma, chain_bonds(n, j), n, delta)
                scale = np.abs(ref).max()
                resid = np.abs(ref - mine).max()
                c = resid / (EPS * scale * n)
                ratios.append((n, c, resid / (EPS * scale),
                               resid / (EPS * scale * np.sqrt(n))))
                gate(f"{name} n={n}, J={j}: closure leak is exactly zero", leak == 0.0,
                     f"leak = {leak:.1e}")
                gate(f"{name} n={n}, J={j}: residual within the noise model", c < 1.0,
                     f"C = residual/(eps*||block||*N) = {c:.3f} against a bound of 1.0, "
                     f"about twice the measured 99th percentile of 0.47")
    # NOTHING IS GATED ON THE N=3 TO N=7 RATIO, and the docstring says why: it is the ratio
    # of two maxima of four rounding residuals each, it runs to 4.57 over 200 draws, and any
    # threshold on it fails about one draw in fourteen. Three gates were built on it across
    # three rounds (a spread, two model verdicts, a non-growth test) and all three are gone.
    # The readings stay, because they are worth seeing; they are not evidence.
    def band(idx):
        return (max(r[idx] for r in ratios if r[0] == 3),
                max(r[idx] for r in ratios if r[0] == 7))
    lo1, hi1 = band(1)
    lo0, hi0 = band(2)
    los, his = band(3)
    print(f"       (context, not gated: max C at N=3 {lo1:.3f}, at N=7 {hi1:.3f}; without the "
          f"size factor {lo0:.3f} to {hi0:.3f}; with sqrt(N) {los:.3f} to {his:.3f}.")
    print("        One measurement divided three ways, and the draw moves all three; the "
          "linear N is structural.)")


def g2_non_derogatory_symbolic():
    """Exact route, so compare exactly: the minor determinant must be literally the product
    of the off-diagonal entries, with no lam in it, for BOTH Delta."""
    print("\nG2  Lemma A, symbolic and Delta-blind: the (N-1)x(N-1) minor is lam-free")
    lam = sp.Symbol("lam")
    for delta in (0, 1):
        for n in range(2, 7):
            g = sp.symbols(f"g0:{n}", positive=True)
            J = sp.symbols(f"J0:{max(n - 1, 1)}", positive=True)
            M = sp.zeros(n, n)
            for l in range(n):
                M[l, l] = -2 * g[l]
            for b in range(n - 1):
                M[b, b] += -2 * sp.I * delta * J[b]
                M[b + 1, b + 1] += -2 * sp.I * delta * J[b]
                M[b, b + 1] += 2 * sp.I * J[b]
                M[b + 1, b] += 2 * sp.I * J[b]
            det = sp.simplify(sp.expand(((M - lam * sp.eye(n))[1:, :-1]).det()))
            expected = sp.simplify(sp.prod([2 * sp.I * J[b] for b in range(n - 1)]))
            gate(f"Delta={delta}, N={n}: minor det == prod(2i*J_b), exactly",
                 sp.simplify(det - expected) == 0, f"det = {det}")
            gate(f"Delta={delta}, N={n}: it is free of lam, exactly", sp.diff(det, lam) == 0,
                 "d/dlam = 0")


def g3_nullity_where_it_bites():
    """The rank statement can only be tested where a degeneracy actually is. At a SIMPLE
    eigenvalue nullity 1 is automatic for any matrix whatsoever, so random draws prove
    nothing: the first two versions of this gate tested exactly that and could not fail. The
    minimum gap over those draws is printed below to show how far from the question they
    were. The assertion is made at the EPs, where algebraic multiplicity is 2 and the nullity
    must still be 1."""
    print("\nG3  Lemma A where it bites: nullity 1 at an eigenvalue of multiplicity 2")
    rng = np.random.default_rng(5)
    mn = np.inf
    for n in range(2, 9):
        for _ in range(40):
            M = block_M(rng.uniform(0.01, 3.0, n),
                        [((b, b + 1), v) for b, v in enumerate(rng.uniform(0.05, 2.0, n - 1))],
                        n, 0)
            mn = min(mn, closest_pair(M)[1])
    print(f"         (context: over 280 random profiles the closest eigenvalue pair is never "
          f"nearer than {mn:.4f}, so a random draw cannot test this)")
    for label, n, build, bracket in WITNESSES:
        x = ep_root(build, n, bracket)
        M = block_M(build(x), chain_bonds(n, CANONICAL_J), n, 0)
        lam, sep = closest_pair(M)
        sv = np.linalg.svd(M - lam * np.eye(n), compute_uv=False)
        gate(f"{label}: exactly one singular value collapses, so nullity 1 and not 2",
             sv[-1] / sv[-2] < 1e-6 and sep < 1e-6,
             f"split {sep:.1e}, two smallest {sv[-2]:.3e}, {sv[-1]:.3e}")


def g4_dimer_closed_form():
    """Exact route throughout."""
    print("\nG4  N=2 in closed form: the EP sits exactly at |gamma_1 - gamma_0| = 2J")
    lam = sp.Symbol("lam")
    g0, g1, J = sp.symbols("g0 g1 J", positive=True)
    M_xy = sp.Matrix([[-2 * g0, 2 * sp.I * J], [2 * sp.I * J, -2 * g1]])
    # Both books are BUILT, from the same generator the rest of the file uses, and then
    # compared. An earlier version wrote M_xxx := M_xy - 2iJ*I and then asserted that
    # difference, which is the definition restated: a sign error in the degree diagonal would
    # have left the gate green.
    D2 = sp.Matrix([[J, 0], [0, J]])                       # the N=2 degree diagonal
    A2 = sp.Matrix([[0, J], [J, 0]])                       # and its adjacency
    M_xxx = -2 * sp.I * (D2 - A2) - 2 * sp.Matrix([[g0, 0], [0, g1]])
    gate("at N=2 the two books differ by the scalar -2iJ, exactly",
         sp.simplify(M_xxx - M_xy + 2 * sp.I * J * sp.eye(2)) == sp.zeros(2, 2),
         "the Laplacian form built from D and A, then compared: so the EP condition is shared")
    for n_chk in (2, 3, 4):
        gam = np.array([0.05 + 0.01 * k for k in range(n_chk)])
        bd = chain_bonds(n_chk, CANONICAL_J)
        built = block_M(gam, bd, n_chk, 1) - block_M(gam, bd, n_chk, 0)
        deg = np.diag([sum(j for (u, v), j in bd if u == k or v == k)
                       for k in range(n_chk)])
        gate(f"N={n_chk}: the two books of block_M differ by exactly -2i*D_J",
             np.array_equal(built, -2j * deg),
             "the degree diagonal, built independently of block_M")
    for name, M in (("XY", M_xy), ("XXX", M_xxx)):
        disc = sp.simplify(sp.discriminant(sp.Poly(M.charpoly(lam).as_expr(), lam)))
        gate(f"{name}: discriminant == 4*(g1-g0)^2 - 16*J^2, exactly",
             sp.simplify(disc - (4 * (g1 - g0) ** 2 - 16 * J ** 2)) == 0, f"disc = {disc}")
    Mep = M_xy.subs(g1, g0 + 2 * J)          # the branch gamma_1 = gamma_0 + 2J
    K = sp.simplify(Mep - sp.simplify(Mep.trace() / 2) * sp.eye(2))
    gate("at the EP the nilpotent is nonzero, exactly", sp.simplify(K) != sp.zeros(2, 2),
         f"K = {list(K)}")
    gate("at the EP its square is the zero matrix, exactly",
         sp.simplify(K * K) == sp.zeros(2, 2), "K^2 == 0, a genuine 2x2 Jordan block")


def g5_bipartite_reality():
    """Exact route: the similarity is an entry-wise rearrangement with signs, so the residual
    is compared to literal zero and a nonzero one would be a finding, not a tolerance."""
    print("\nG5  a path is bipartite, so at Delta=0 the discriminant is REAL (codimension 1)")
    rng = np.random.default_rng(99)
    for n in (3, 4, 5, 6, 7):
        S = np.diag([(-1) ** k for k in range(n)])
        gamma = rng.uniform(0.01, 2.0, n)
        js = rng.uniform(0.05, 2.0, n - 1)
        bonds = [((b, b + 1), js[b]) for b in range(n - 1)]
        M0 = block_M(gamma, bonds, n, 0)
        gate(f"N={n}, Delta=0: S conj(M) S^-1 == M, exactly",
             np.array_equal(S @ M0.conj() @ S, M0),
             f"residual {np.abs(S @ M0.conj() @ S - M0).max():.1e}")
        c0 = np.poly(M0)
        gate(f"N={n}, Delta=0: the characteristic polynomial is real",
             np.abs(c0.imag).max() < 32 * EPS * np.abs(c0).max(),
             f"max |Im coeff| = {np.abs(c0.imag).max():.1e} against "
             f"|coeff| <= {np.abs(c0).max():.1e}")
        M1 = block_M(gamma, bonds, n, 1)
        gate(f"N={n}, Delta=1: the same similarity FAILS, so that book is different",
             np.abs(S @ M1.conj() @ S - M1).max() > 1e-3,
             f"residual {np.abs(S @ M1.conj() @ S - M1).max():.3f}, the degree diagonal breaks it")


def exact_block_M(gamma_rat, n):
    """The Delta=0 block over the rationals: entries are exact, so everything downstream is."""
    M = sp.zeros(n, n)
    for k in range(n):
        M[k, k] = -2 * gamma_rat[k]
    for b in range(n - 1):
        M[b, b + 1] = M[b + 1, b] = 2 * sp.I * EXACT_J
    return M


def exact_disc(M):
    lam = sp.Symbol("lam")
    return sp.expand(sp.discriminant(sp.Poly(M.charpoly(lam).as_expr(), lam), lam))


def g6_existence_by_sign_change():
    """The certificate, and it is EXACT. Every input is rational: J = 75/1000, every gamma
    entry and both bracket endpoints are exact decimals. So the discriminant has an exact
    route and the house rule puts this in case 1, compare exactly, no tolerance anywhere.
    A real continuous function with opposite signs at the ends of an interval has a zero
    inside; that is all that is used.

    TWO REASONS, and BOTH magnitude arguments tried here before were wrong, so they are
    recorded rather than quietly dropped. First, the old reality test asked whether
    |Im(disc)| was small RELATIVE TO |disc| itself, which cannot fail for an input whose
    imaginary part IS a rounding residual, and that is every input of interest here.
    Second, and this is what the exactness buys: a family whose discriminant vanishes
    identically has a float value that is pure rounding with a sign nothing determines, and
    G6b is that family.
    What is NOT a reason, twice over. The first version of this comment said the witness and
    null magnitudes were comparable, near 1e-24 apiece; they are not, the null family lands
    at 1e-46 and 1e-40. The second said the witnesses span decades, 3.65e-21 to 4.95e-10, as
    if a spread ACROSS fixtures were a noise floor WITHIN one; it is not, and every one of
    the fourteen float endpoints reproduces its exact value to ten significant figures or
    better (worst relative disagreement 1.1e-10). That figure is itself ungated, and is
    quoted here only to say what the withdrawn arguments got wrong.
    The seven witnesses are unchanged: their exact signs are the ones the float route
    reported, so the exactness here is discipline, not rescue."""
    print("\nG6  existence: the EXACT discriminant changes sign across each bracket")
    for label, n, build, (a, b) in WITNESSES:
        da = exact_disc(exact_block_M(build_exact(build, a), n))
        db = exact_disc(exact_block_M(build_exact(build, b), n))
        gate(f"{label}: the discriminant is exactly real at both ends",
             sp.im(da) == 0 and sp.im(db) == 0,
             "Im(disc) == 0 as a rational, not at a tolerance")
        gate(f"{label}: and it changes sign across [{a:.6f}, {b:.6f}]",
             sp.sign(sp.re(da)) * sp.sign(sp.re(db)) == -1,
             f"exact signs {sp.sign(sp.re(da))} -> {sp.sign(sp.re(db))}")
        x = ep_root(build, n, (a, b))
        M = block_M(build(x), chain_bonds(n, CANONICAL_J), n, 0)
        lam, sep = closest_pair(M)
        gate(f"{label}: the bracketed root is a coalescence", sep < 1e-6,
             f"gamma_1 = {x:.12f}, |split| = {sep:.2e}")
        # WHICH QUANTITY, because the two available ones obey different laws and an earlier
        # version of this comment named the wrong one. closest_pair returns the pair MEAN,
        # which is a trace over a well-separated invariant 2-space (the next eigenvalue is
        # 0.17 to 0.64 away here, gated below), so it is a stable quantity and its computed
        # imaginary part is the cancellation residual of two conjugates: bounded by
        # eps*||M||_2*n. The SPLIT is the unstable one and obeys sqrt(eps*||M||_2) instead,
        # a factor 8.5e6 to 2.3e7 above it across the seven as a ratio of the two
        # MODELS (the measured splits give a wider 3.1e6 to 3.5e7, neither band gated);
        # G12's docstring records that trap already, and calling this bound a bound on the
        # split would repeat it.
        # The bound is a BOUND, not a constant: a cancellation residual may land anywhere
        # beneath it, and across these seven it does (ratios 0.0007 to 0.25, a spread of 370
        # that is the object, not a defect). What must not degrade is the bound itself as the
        # SCALE moves, and G6c gates that over four decades and all seven witnesses.
        floor = EPS * np.linalg.norm(M, 2) * n
        # Isolation against the PERTURBATION eps*||M||_2, not against the split. An earlier
        # version compared with the split, which scales the wrong way: as the bracket tightens
        # and the split shrinks toward a better EP, that requirement relaxes toward vacuity.
        # What the cluster mean needs is that no third eigenvalue lies within reach of the
        # rounding, and the margin here is qualitative (many orders), not a tight law.
        rest = sorted(abs(np.linalg.eigvals(M) - lam))[2:]
        pert = EPS * np.linalg.norm(M, 2)
        gate(f"{label}: the coalescing pair is isolated, so its mean is a stable quantity",
             min(rest) > 1e8 * pert,
             f"next eigenvalue {min(rest):.4f} away, which is {min(rest) / pert:.1e} times "
             f"the perturbation eps*||M||_2")
        gate(f"{label}: and the meeting point is REAL, the pair MEAN under eps*||M||_2*n",
             abs(lam.imag) < floor,
             f"|Im lam| = {abs(lam.imag):.2e}, ratio = {abs(lam.imag) / floor:.4f}")


def g6c_the_reality_floor_holds_across_scales():
    """The law behind G6's Im assertion, and the reason the proof quotes a RATIO and never a
    magnitude. Scaling gamma and J together rescales the block exactly (G11), so the true
    meeting point stays real while the computed residual scales with it: the invariant is the
    ratio to eps*||M||_2*n, and the assertion is that BOUND, re-read at every witness and
    every scale. Non-degradation across scale is NOT asserted: with three scale points under
    a spread that runs to 690 it could not be tested, and an earlier version of this line
    claimed it anyway.

    A BOUND, and deliberately not a spread. An earlier version gated max/min of the ratio
    over three scales at one witness. That is the wrong statistic for exactly the reason G1's
    docstring gives, and it was green only because it ran on the single witness that clears
    it: the other six spread by 3.5 to 168, because a cancellation residual resamples freely
    beneath its bound every time the rounding changes. The spread IS the object, and no
    flatness is asserted anywhere here.

    WHAT IS ASSERTED IS THE BOUND, and nothing about the size factor's shape. Competing
    models are printed below as a reading, never gated: over two rounds three assertions
    about them stood here and all three were withdrawn as draw-dependent or tautological,
    and the comment at the print says which. G1 does not supply the shape either. The linear
    size factor is structural in both places and gated nowhere, which is the honest state."""
    print("\nG6c the reality floor holds across four decades of scale, at every witness")
    worst, worst_at = 0.0, ""
    alt = {"constant 5": 0.0, "no size factor": 0.0, "eps*max|entry|": 0.0}
    for label, n, build, bracket in WITNESSES:
        x0 = ep_root(build, n, bracket)
        for scale in (1.0, 1e2, 1e4):
            gam = [v * scale for v in build(x0)]
            M = block_M(gam, chain_bonds(n, CANONICAL_J * scale), n, 0)
            lam, _ = closest_pair(M)
            nrm = np.linalg.norm(M, 2)
            r = abs(lam.imag) / (EPS * nrm * n)
            if r > worst:
                worst, worst_at = r, f"{label} at scale {scale:.0e}"
            r_free = abs(lam.imag) / (EPS * nrm)
            alt["constant 5"] = max(alt["constant 5"], abs(lam.imag) / (EPS * nrm * 5))
            alt["no size factor"] = max(alt["no size factor"], r_free)
            alt["eps*max|entry|"] = max(alt["eps*max|entry|"],
                                        abs(lam.imag) / (EPS * np.abs(M).max()))
            gate(f"{label} at scale {scale:.0e}: the ratio stays under the bound",
                 r < 1.0, f"|Im lam| = {abs(lam.imag):.2e}, ratio {r:.4f}")
    # The counter-models, COMPUTED on the same 21 readings rather than quoted in a message.
    # An earlier version put 2.28 and 2.44 in an f-string beside an assertion that was a
    # restatement of the 21 gates above it: a tautology carrying two ungated numbers, which
    # is the defect G13's literals were added to close.
    # WHAT IS NOT GATED HERE, and why nothing was added to fix it. Three assertions stood
    # in this place across two earlier rounds and all three were withdrawn:
    #   - that the size-free model is n times this one. True, but it is IEEE division under
    #     two orders, an arithmetic identity that passes on any input at all.
    #   - that the constant 5 is indistinguishable from n. True today, and only because the
    #     largest of 21 cancellation residuals happens to land on an N=5 witness; reversing
    #     the site order (a relabelling that changes no spectrum) moves it to N=4.
    #   - that eps*max|entry| is not an independent alternative. Also draw-dependent.
    # The competing magnitudes are printed as a reading. The shape of the size factor is
    # gated nowhere in this file and the proof says so; it is structural in both places.
    print(f"       (context, not gated: worst {worst:.3f} at {worst_at}; with no size factor "
          f"{alt['no size factor']:.3f}, with the constant 5 {alt['constant 5']:.3f}, "
          f"with eps*max|entry| {alt['eps*max|entry|']:.3f}.")
    print("        N spans only 4 to 6 here, so these separate an overall factor at best, "
          "never the shape.)")


def g6b_negative_control():
    """The control the float certificate could not have failed. A bipartite STAR with THREE
    equal leaf rates carries a permanently repeated eigenvalue: the equal leaves span a
    two-dimensional space of difference modes and the repetition is semisimple, so the
    certificate has nothing here to certify, and the discriminant is identically zero.
    SCOPE, because an earlier wording said "no EP at any parameter" while gating two points:
    what is shown is that THIS repetition is semisimple at the two gated rates. A third
    eigenvalue colliding with the doubled one elsewhere in the family would be a different
    object, and disc == 0 would hide it from the certificate either way. The nearest such
    approach is at gamma_1 -> 0.09, where all four leaves coincide, the null space simply
    grows to three, and it is still semisimple.
    THREE and not two, for the reason G10's comment gives below: two equal leaves span only
    ONE difference mode, the eigenvalue is simple, and the exact discriminant of the two-leaf
    variant is NOT zero. An earlier version of this docstring said two while the code said
    three, which is the same slip G10 records fifty lines down. Semisimplicity is gated here
    rather than deferred to G10, whose fixture is a different one."""
    print("\nG6b negative control: a family with NO EP is refused by the exact route")
    n = 5
    for x in (sp.Rational(10, 100), sp.Rational(30, 100)):
        gam = [sp.Rational(5, 100), x, sp.Rational(9, 100), sp.Rational(9, 100),
               sp.Rational(9, 100)]
        M = sp.zeros(n, n)
        for k in range(n):
            M[k, k] = -2 * gam[k]
        for k in range(1, n):
            M[0, k] = M[k, 0] = 2 * sp.I * EXACT_J
        d = exact_disc(M)
        gate(f"star, gamma_1 = {x}: the exact discriminant is exactly zero", d == 0,
             "no sign to change, so the certificate cannot fire here")
        Mf = block_M(np.array([float(g) for g in gam]), [((0, k), CANONICAL_J)
                     for k in range(1, n)], n, 0)
        df = discriminant(Mf)
        gate(f"star, gamma_1 = {x}: and the FLOAT discriminant is pure noise",
             abs(df) < 1e-20,
             f"float disc = {df.real:+.3e}, which a bare sign test would read as a value")
        # SEMISIMPLE means geometric == algebraic, and a float nullity reads only the
        # geometric half: an alg 3 / geom 2 family would carry a Jordan block and pass it.
        # Both halves are taken exactly, off the same sympy matrix the exact discriminant
        # above used. An earlier version took the geometric half from a float SVD while the
        # exact route was three lines up.
        lam_sym = sp.Symbol("lam")
        roots = sp.roots(sp.Poly(M.charpoly(lam_sym).as_expr(), lam_sym))
        rep = [(r, m) for r, m in roots.items() if m > 1]
        # Sentinels chosen so a malformed root set can never satisfy alg == geo, rather
        # than relying on the len(rep) conjunct alone to save it.
        alg = rep[0][1] if len(rep) == 1 else -1
        geo = (n - (M - rep[0][0] * sp.eye(n)).rank()) if len(rep) == 1 else -2
        # Semisimple is alg == geo, and that is what is asserted; the 2 is printed, not
        # required. An earlier version hard-coded alg == 2, which would have REJECTED the
        # four-equal-leaf point this docstring names as semisimple.
        gate(f"star, gamma_1 = {x}: the repetition is SEMISIMPLE, exactly, so it is no EP",
             len(rep) == 1 and alg == geo,
             f"one repeated root, algebraic multiplicity {alg} == geometric multiplicity "
             f"{geo}, both exact over Q(i) (the off-diagonals carry 2iJ)")


def g7_window_and_positivity():
    print("\nG7  every witness has positive rates and sits STRICTLY inside the rate window")
    for label, n, build, bracket in WITNESSES:
        x = ep_root(build, n, bracket)
        gamma = build(x)
        lam, _ = closest_pair(block_M(gamma, chain_bonds(n, CANONICAL_J), n, 0))
        lo, hi = -2 * max(gamma), -2 * min(gamma)
        gate(f"{label}: every rate is positive", all(g > 0 for g in gamma),
             f"gamma = {np.round(gamma, 6).tolist()}, contrast "
             f"max/min = {max(gamma) / min(gamma):.2f}")
        gate(f"{label}: Re lam strictly inside [{lo:.4f}, {hi:.4f}]", lo < lam.real < hi,
             f"Re lam = {lam.real:.8f}")


def g8_sqrt_branch():
    """Corroboration, and it is reported as such: the claim rests on G2 (exact) and G6 (a sign
    change), neither of which needs this. The discriminator is the repo's own, the log-log
    SLOPE of the gap against the distance to the point, 1/2 for a square-root branch and 1 for
    a node. An earlier version compared the per-point exponent to 1/2 instead, which decreases
    monotonically for ANY prefactor and therefore tested arithmetic rather than the object."""
    print("\nG8  corroboration: the gap scales as a square root, slope 1/2 and not 1")
    for label, n, build, bracket in WITNESSES:
        x = ep_root(build, n, bracket)
        es = [1e-3, 1e-4, 1e-5, 1e-6]
        ss = [closest_pair(block_M(build(x + e), chain_bonds(n, CANONICAL_J), n, 0))[1]
              for e in es]
        slope = float(np.polyfit(np.log(es), np.log(ss), 1)[0])
        gate(f"{label}: log-log slope is 1/2 to three decimals", abs(slope - 0.5) < 5e-3,
             f"slope = {slope:.6f} (a node would give 1.0)")


def g9_the_rate_profile_is_what_does_it():
    """Part C. Lemma A applies to the Hermitian Jacobi matrix too, so no bond profile can make
    the hopping spectrum degenerate; a first version of this argument cited the closed form
    2J(1 - cos(m*pi/N)), which is true only at UNIFORM J and therefore did not cover the case
    it was invoked for."""
    print("\nG9  it is the RATE profile: no J profile degenerates the hopping spectrum")
    rng = np.random.default_rng(3)
    worst = np.inf
    for n in range(2, 11):
        for _ in range(300):
            js = rng.uniform(0.05, 2.0, n - 1)
            M = block_M([CANONICAL_G0] * n, [((b, b + 1), js[b]) for b in range(n - 1)], n, 0)
            worst = min(worst, closest_pair(M)[1])
    gate("2700 random J profiles at uniform gamma: never a degeneracy", worst > 1e-6,
         f"smallest gap seen = {worst:.4e}")
    for n in (3, 5, 7, 9):
        M = block_M([CANONICAL_G0] * n, chain_bonds(n, 1.0), n, 0)
        c = M @ M.conj().T - M.conj().T @ M
        # Exact route: at uniform gamma the two products cancel entry by entry, so this is
        # compared to literal zero. A nonzero residual here would be a finding, not noise.
        gate(f"N={n}: at uniform gamma [M, M^dag] is exactly zero", bool(np.all(c == 0)),
             f"max |entry| = {np.abs(c).max():.1e}")


def g10_star_control():
    """The control that isolates the right hypothesis. An earlier version used a star at
    UNIFORM gamma, where M is normal, so normality alone explained the semisimplicity and
    tridiagonality never entered: that fixture could not have failed. With a rate PROFILE the
    star is genuinely non-normal and the repetition is STILL semisimple, which is what
    separates 'not normal' from 'defective'."""
    print("\nG10 control: a PROFILED star is non-normal and still semisimple")
    # THREE equal leaf rates, not two. On a star the repetition comes from the leaf
    # DIFFERENCE modes, whose eigenvalue is -2*gamma_leaf because the hub coupling cancels;
    # two equal leaves span only one such mode and give a SIMPLE eigenvalue. A first version
    # used two and the gate said so, which is the fixture being checked rather than assumed.
    profiles = {5: [0.05, 0.3, 0.7, 0.7, 0.7],
                6: [0.05, 0.3, 0.9, 0.7, 0.7, 0.7],
                7: [0.05, 0.3, 0.9, 0.4, 0.7, 0.7, 0.7]}
    for n, gamma in profiles.items():
        M = block_M(gamma, [((0, k), 1.0) for k in range(1, n)], n, 0)
        c = M @ M.conj().T - M.conj().T @ M
        lam, sep = closest_pair(M)
        sv = np.linalg.svd(M - lam * np.eye(n), compute_uv=False)
        nullity = int(np.sum(sv < 1e-9 * sv[0]))
        gate(f"N={n}: the profiled star is NOT normal", np.abs(c).max() > 1e-3,
             f"||[M, M^dag]||_max = {np.abs(c).max():.3f}")
        gate(f"N={n}: it still has a repeated eigenvalue", sep < 1e-9, f"|split| = {sep:.1e}")
        gate(f"N={n}: and that eigenvalue is SEMISIMPLE", nullity >= 2,
             f"nullity = {nullity}, so not a Jordan block")


def g11_scale_law():
    """Exact route: both terms of the adjacency form are homogeneous of degree one in
    (gamma, J) jointly, so the EP condition depends on gamma/J alone."""
    print("\nG11 the XY block obeys M(gamma, J) = J * M(gamma/J, 1), exactly")
    rng = np.random.default_rng(17)
    for n in (3, 5, 7):
        gamma = rng.uniform(0.01, 2.0, n)
        j = CANONICAL_J
        lhs = block_M(gamma, chain_bonds(n, j), n, 0)
        rhs = j * block_M(gamma / j, chain_bonds(n, 1.0), n, 0)
        gate(f"N={n}: the two sides agree at the float noise",
             np.abs(lhs - rhs).max() < 10 * EPS * np.abs(lhs).max(),
             f"residual {np.abs(lhs - rhs).max():.1e}; the EP condition depends only on gamma/J")


def g12_what_the_contrast_costs():
    """The reachability reading, and it is a TRADE-OFF CURVE rather than the single number an
    earlier version printed. §(d) gives the N=2 law exactly, EP at |gamma_1 - gamma_0| = 2J,
    so with the hop 2J and the carrier ratio Q = 2J/gamma_0 the required contrast is
    1 + Q at N=2, and the N=4 witness family tracks it with a mild size correction. The
    canonical Q = 1.5 point is the one the repo actually operates at, and it is CHEAPER in
    contrast than the J = 0.075 line this file's witnesses use."""
    print("\nG12 what the contrast costs: it follows the coupling, and 1 + Q is the N=2 law")
    # Exact route at N=2: the rates and the coupling are exact decimals and the discriminant
    # is a rational, so the coalescence is asserted against literal zero. An earlier version
    # asked the eigensolver for the split and gated it at 1e-12, which is the wrong law: AT a
    # Jordan point the computed split is of order sqrt(eps), so that gate failed on healthy
    # arithmetic. G13 below is where the sqrt behaviour is the object rather than the noise.
    for q_carrier in (sp.Rational(12, 10), sp.Rational(15, 10), sp.Rational(2),
                      sp.Rational(3)):
        g0 = sp.Rational(5, 100)
        j_proof = q_carrier * g0 / 2                       # hop 2J = Q*gamma_0
        g1 = g0 + 2 * j_proof                              # the N=2 closed form
        M = sp.Matrix([[-2 * g0, 2 * sp.I * j_proof], [2 * sp.I * j_proof, -2 * g1]])
        d = exact_disc(M)
        # The contrast identity is true BY CONSTRUCTION (g1 was built as g0 + 2J = g0*(1+Q)),
        # so it is printed, not asserted. What the gate asserts is the part with content:
        # that this constructed point really is the coalescence, exactly.
        gate(f"Q = {q_carrier}: the N=2 point at contrast {float(g1 / g0):.2f} = 1 + Q IS the EP",
             d == 0,
             f"exact disc == {d}; contrast {g1 / g0} = 1 + Q holds by construction")
    # and the same trade-off on the N=4 witness family, where the size correction lives
    for q_carrier, expect in ((1.5, 2.30), (3.0, 4.42)):
        j_proof = q_carrier * 0.05 / 2.0
        build = lambda x: [0.05, x, 0.07, 0.09]
        f = lambda x: discriminant(block_M(build(x), chain_bonds(4, j_proof), 4, 0)).real
        lo, hi = 0.051, 1.5
        xs = np.linspace(lo, hi, 3000)
        root = None
        prev = f(xs[0])
        for x in xs[1:]:
            cur = f(x)
            if prev * cur < 0:
                root = brentq(f, x - (xs[1] - xs[0]), x, xtol=1e-15)
                break
            prev = cur
        contrast = max(build(root)) / min(build(root))
        gate(f"N=4 at Q = {q_carrier}: the contrast is {contrast:.2f}",
             abs(contrast - expect) < 0.005,
             f"EP at gamma_1 = {root:.6f}, contrast {contrast:.4f} against {expect}")


def g13_the_binding_constraint_is_the_tuning():
    """What actually gates a device, and it is NOT the contrast. At an EP the split opens like
    the square root of the detuning (G8 measures slope 1/2), so a relative error d in one rate
    opens a split of order sqrt(d) times the scale. Measured against the mode's own decay
    rate, one percent of rate-setting accuracy already splits the N=4 witness by about a fifth
    of |Re lambda|. The gate asserts the sqrt LAW rather than any one number."""
    print("\nG13 the binding constraint: the split opens like sqrt(detuning), so tuning gates it")
    label, n, build, bracket = WITNESSES[0]
    x0 = ep_root(build, n, bracket)
    M0 = block_M(build(x0), chain_bonds(n, CANONICAL_J), n, 0)
    lam0, _ = closest_pair(M0)
    rate = abs(lam0.real)
    # The percentages the proof's reachability section quotes, as LITERALS. An earlier
    # version interpolated them into the gate NAME and asserted only the slope, so the
    # document could have gone stale against a green gate.
    gate(f"{label}: the mode rate the percentages are measured against is 0.2124",
         abs(rate - 0.2124) < 5e-5, f"|Re lambda| = {rate:.6f}")
    prev_ratio = None
    ratios_seen = []
    for rel, pct in ((1e-6, 0.2), (1e-4, 1.9), (1e-2, 18.9)):
        M = block_M(build(x0 * (1 + rel)), chain_bonds(n, CANONICAL_J), n, 0)
        _, sep = closest_pair(M)
        ratio = sep / np.sqrt(rel)
        gate(f"{label}: a relative rate error {rel:.0e} splits by {sep:.1e}, "
             f"i.e. {100 * sep / rate:.1f}% of |Re lambda|, and the doc says {pct}%",
             abs(100 * sep / rate - pct) < 0.05
             and (prev_ratio is None or abs(ratio / prev_ratio - 1) < 0.05),
             f"split/sqrt(rel) = {ratio:.4f}"
             + ("" if prev_ratio is None else f", flat against {prev_ratio:.4f}"))
        prev_ratio = ratio
        ratios_seen.append(ratio)
    # The constant the proof CONVERTS through, asserted rather than printed. An earlier
    # version left 0.394 in a detail string only, so the document could have gone stale
    # against a green gate: the very defect this gate's literals were added to close.
    gate("the constant the proof converts through is 0.394, measured at the near-EP end",
         abs(ratios_seen[0] - 0.394) < 0.0005,
         f"split/sqrt(rel) = {ratios_seen[0]:.4f} at rel = 1e-6, drifting to "
         f"{ratios_seen[-1]:.4f} at rel = 1e-2, the far end of the sqrt branch")


def main():
    print("The edge block carries a defective EP under a per-site rate profile")
    print("=" * 90)
    g1_anchor()
    g2_non_derogatory_symbolic()
    g3_nullity_where_it_bites()
    g4_dimer_closed_form()
    g5_bipartite_reality()
    g6_existence_by_sign_change()
    g6b_negative_control()
    g6c_the_reality_floor_holds_across_scales()
    g7_window_and_positivity()
    g8_sqrt_branch()
    g9_the_rate_profile_is_what_does_it()
    g10_star_control()
    g11_scale_law()
    g12_what_the_contrast_costs()
    g13_the_binding_constraint_is_the_tuning()
    print("=" * 90)
    print(f"{PASSED}/{PASSED + FAILED} gates passed")
    if FAILED:
        raise AssertionError(f"{FAILED} gate(s) failed")


if __name__ == "__main__":
    main()
