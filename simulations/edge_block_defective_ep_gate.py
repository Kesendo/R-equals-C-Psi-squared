"""The edge block carries a defective EP under a per-site rate profile, and on the XY chain
the EP set has codimension one.

WHAT THIS SETTLES. The FIRST item of the (b2) inventory in the arc `site_resolved_vacuum_block`
(compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs:4118) says, verbatim, "'no defective EP
can live on an edge block' is no longer forbidden a priori ... Settle it before fencing, do not
assume either way." Settled in the positive: the EP exists, at real positive rates, at the
repo's canonical anchors, and at the sizes where the lemma that forbids it is used.

TWO OPERATORS, AND THE FIRST VERSION OF THIS GATE MEASURED THE WRONG ONE. On the (0,1)
coherence block, spanned by |0><j|:

    Delta = 1 (XXX)  M = -2i*(D_J - A_J) - 2*diag(gamma)     the Laplacian form, F152's
    Delta = 0 (XY)   M = +2i*A_J        - 2*diag(gamma)      the adjacency form

`PROOF_CODIM1_BY_ADDITIVITY.md` and F125, which carry the Edge lemma, are scoped to the **XY
chain, Delta = 0**. F152 is fenced to |Delta| = 1 and says so. A first version of this gate ran
every witness on the Laplacian form and was 81/81 green while measuring an operator the refuted
lemma is not about; on the adjacency form those same profiles have splits 1.54, 1.67 and 0.107,
which are not EPs at all. Both forms are gated here, and every WITNESS is on Delta = 0.

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
     rate space rather than two, so its solution set is a CURVE and turning a single rate
     crosses it. Existence follows from a sign change, by the intermediate value theorem on
     a real continuous function; no topological degree is needed and none is used.
     At Delta = 1 the degree diagonal D breaks that antisymmetry (S*D*S^-1 = D, not -D), the
     discriminant is genuinely complex, and the EPs are isolated points instead. G5 gates
     both halves of that contrast.

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
does not force defectiveness; tridiagonality does.

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
CANONICAL_J = 0.075        # docs/Q_REGIME_ANCHORS.md, the Q = 1.5 anchor
CANONICAL_G0 = 0.05        # the canonical base rate

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
    number, and the law had to be MEASURED rather than guessed. Against eps*||block|| alone
    the ratio grows from about 0.7 at N=3 to about 1.9 at N=7, because the block dimension
    is N and more summands enter; against eps*||block||*N it is flat. So the model is

        residual  <=  C * eps * ||block||_max * N,

    and the claim is that C does not grow. A first version gated the SPREAD max/min of the
    ratio, which is the wrong statistic for an upper bound: a configuration whose entries
    happen to be exactly representable makes the residual smaller than the model, and that
    is not a violation of anything. What is asserted instead is a bound on C and the absence
    of a trend in C across the size range."""
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
                c = np.abs(ref - mine).max() / (EPS * scale * n)
                ratios.append((n, c))
                gate(f"{name} n={n}, J={j}: closure leak is exactly zero", leak == 0.0,
                     f"leak = {leak:.1e}")
                gate(f"{name} n={n}, J={j}: residual within the noise model", c < 0.5,
                     f"C = residual/(eps*||block||*N) = {c:.3f}")
    small = max(c for n, c in ratios if n == 3)
    large = max(c for n, c in ratios if n == 7)
    gate("C does not grow across the size range (the law, not the bound)",
         large < 2.0 * small,
         f"max C at N=3 is {small:.3f}, at N=7 is {large:.3f}, growth {large / small:.2f}x")


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
    M_xxx = M_xy - 2 * sp.I * J * sp.eye(2)
    gate("at N=2 the two books differ by the scalar -2iJ, exactly",
         sp.simplify(M_xxx - (M_xy - 2 * sp.I * J * sp.eye(2))) == sp.zeros(2, 2),
         "so the EP condition is the same in both")
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
             f"max |Im coeff| = {np.abs(c0.imag).max():.1e} against |coeff| <= {np.abs(c0).max():.1e}")
        M1 = block_M(gamma, bonds, n, 1)
        gate(f"N={n}, Delta=1: the same similarity FAILS, so that book is different",
             np.abs(S @ M1.conj() @ S - M1).max() > 1e-3,
             f"residual {np.abs(S @ M1.conj() @ S - M1).max():.3f}, the degree diagonal breaks it")


def g6_existence_by_sign_change():
    """The certificate. A real continuous function with opposite signs at the ends of an
    interval has a zero inside; that is all that is used, and it needs no degree, no radius
    and no point count. The gate asserts the two signs, that the imaginary part is at the
    float noise (so the function really is the real object the argument is about), and that
    the bracketed root is a coalescence."""
    print("\nG6  existence: the REAL discriminant changes sign across each bracket")
    for label, n, build, (a, b) in WITNESSES:
        da = discriminant(block_M(build(a), chain_bonds(n, CANONICAL_J), n, 0))
        db = discriminant(block_M(build(b), chain_bonds(n, CANONICAL_J), n, 0))
        scale = max(abs(da), abs(db))
        gate(f"{label}: the discriminant is real at both ends",
             max(abs(da.imag), abs(db.imag)) < 1e-6 * scale,
             f"|Im|/|disc| <= {max(abs(da.imag), abs(db.imag)) / scale:.1e}")
        gate(f"{label}: and it changes sign across [{a:.6f}, {b:.6f}]", da.real * db.real < 0,
             f"disc {da.real:+.3e} -> {db.real:+.3e}")
        x = ep_root(build, n, (a, b))
        _, sep = closest_pair(block_M(build(x), chain_bonds(n, CANONICAL_J), n, 0))
        gate(f"{label}: the bracketed root is a coalescence", sep < 1e-6,
             f"gamma_1 = {x:.12f}, |split| = {sep:.2e}")


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


def main():
    print("The edge block carries a defective EP under a per-site rate profile")
    print("=" * 90)
    g1_anchor()
    g2_non_derogatory_symbolic()
    g3_nullity_where_it_bites()
    g4_dimer_closed_form()
    g5_bipartite_reality()
    g6_existence_by_sign_change()
    g7_window_and_positivity()
    g8_sqrt_branch()
    g9_the_rate_profile_is_what_does_it()
    g10_star_control()
    g11_scale_law()
    print("=" * 90)
    print(f"{PASSED}/{PASSED + FAILED} gates passed")
    if FAILED:
        raise AssertionError(f"{FAILED} gate(s) failed")


if __name__ == "__main__":
    main()
