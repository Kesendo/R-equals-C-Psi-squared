#!/usr/bin/env python3
"""The band-edge sector under a rate PROFILE: what the "whole sector at Re = -2*gamma" family may say.

A large family of statements in this repository says that the single-excitation band-edge sector,
the |vac><psi_k| coherences of the (0,1) block, sits entirely at Re = -2*gamma. Every one of them
is a UNIFORM-gamma statement, and the common ancestor F50 (docs/proofs/PROOF_WEIGHT1_DEGENERACY.md)
says so in its own first line; the consumers dropped the word at the citation. Fencing them needs an
instrument that can express the violation, which is what simulations/framework/weight_coherence_block.py
gained on 2026-08-21, and a fence with no such instrument is unfalsifiable rather than true.

WHAT THIS SCRIPT SETTLES, and the point is that the naive fence is wrong in BOTH directions. The
tempting repair, "under a profile the value 2*gamma becomes a window", is not what happens:

  (1) the DIAGONAL always comes apart. Under a profile the cell (0,j) carries -2*gamma_j, so the
      Absorption-Theorem diagonal is the site-resolved sum and never one number. That part is exact
      and needs no eigensolver.

  (2) the SPECTRUM need not follow it. On the R90 locus (an arithmetic ramp) at delta = 0 and J = 1
      the block is FLAT to machine zero, and flat at -2*gamma_bar: not at its FLOOR -2*gamma_min. A
      statement that survives as a value can still be false about WHICH value, which is the shape
      F153 records as "the block spreads is wrong in two directions".

  (3) the SET the statement quantifies over can be EMPTY. With one site detuned, no mode sits at
      -2*gamma_bar at all, so a claim of the form "the modes at Re = -2*gamma are ..." loses its
      subject rather than its value. PROOF_RING_GAP_DOMINANCE.md:126 states exactly this, measured,
      and this script reproduces it from the other side.

  (4) the MECHANISM fails outright. The gamma-protection argument is not that gamma cancels but that
      |vac><psi_k| is a simultaneous eigenoperator of L_D and L_H. Under any non-uniform profile it
      is an eigenoperator of neither: the residual is O(1), not a small correction.

So a class-B fence must say which of (1)-(4) it means. "It becomes a window" says none of them and
is false as often as it is true.

WHAT THIS SCRIPT DOES NOT DO. It does not touch the Edge lemma or the defective-EP question, which
are class A and live in docs/proofs/PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md. It does not settle
whether the band edge still sets the GAP of the full Liouvillian under a profile: that is D06's
separate question (Q > Q*_gap(N)) and needs the full 4^N spectrum, not one block. And it does not
claim the flat cases below are the only ones; they are the ones a fence must not contradict.

Tolerances follow the repo rule. Where an exact route exists it is taken: the rates are dyadic, so
the diagonal sums are exact in binary floating point and are compared with ==. Where an eigensolver
is unavoidable the test is a SEPARATION, not a threshold: the flat cases sit at 1e-15 and the spread
cases at 7.6e-2 or larger, thirteen decades apart, and each side is asserted against the other rather
than against a number chosen to pass.
"""

import sys

import numpy as np

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402

TOL_FLAT = 1e-12      # an eigensolver's noise floor on these sizes; the flat cases measure ~1e-15
MIN_SPREAD = 1e-3     # the spread cases measure 7.6e-2 .. 1.8e0, so this is a wall, not a fit


# --------------------------------------------------------------------------- topology + the block

def chain(n):
    return [(i, i + 1) for i in range(n - 1)]


def ring(n):
    return [(i, (i + 1) % n) for i in range(n)]


def star(n):
    return [(0, i) for i in range(1, n)]


def complete(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def band_edge_block(n, bonds, gammas, J=1.0, delta=0.0):
    """The (0,1) block on an arbitrary graph, built from the bond list.

    On |vac><j| the coherent part is the single-particle hopping (adjacency for XY; the ZZ term of
    an XXZ coupling adds the degree diagonal, which is why delta turns the adjacency into a
    Laplacian), and the dissipator is entry-wise -2*gamma_j because every such coherence disagrees
    in exactly the one bit j.
    """
    adj = np.zeros((n, n))
    for bond in bonds:
        a, b = bond[0], bond[1]
        w = bond[2] if len(bond) > 2 else 1.0
        adj[a, b] = adj[b, a] = w
    deg = np.diag(adj.sum(axis=1))
    hop = adj - delta * deg
    return -2j * J * hop - 2.0 * np.diag(np.asarray(gammas, float))


def re_spread(m):
    e = np.linalg.eigvals(m).real
    return float(np.ptp(e)), e


def r90_defect(g):
    """Distance from the R90 fixed set x -> 2*avg - reverse(x) (ParameterKlein's anti-palindromic class).

    A control profile has to PROVE it is generic: the first version of this script used
    [0.375, 0.875, 0.25, 0.75] as a "generic" profile at N=4 and it is R90-FIXED, so the control
    measured the same locus it was meant to contrast with, and the gate failed. That is the control
    working; a control that cannot fail proves nothing.
    """
    g = np.asarray(g, float)
    return float(np.abs(g - (2 * g.mean() - g[::-1])).max())


# --------------------------------------------------------------------------------------- the gates

def gate_diagonal_is_exact(report):
    """(1) The diagonal is the site-resolved sum, exactly, and only uniform gamma collapses it."""
    ok = True
    for n in (4, 5, 6, 7):
        prof = np.array([0.5, 0.25, 1.0, 0.125, 0.75, 0.375, 2.0][:n])
        A, _ = fw.weight_block_pencil(n, 0, 1, gamma=prof)
        # the (0,1) cell j disagrees in exactly bit j, so its rate is -2*gamma_j, exactly
        want = -2.0 * prof
        exact = bool(np.all(A == want))
        distinct = len(set(A.tolist())) == n
        good = exact and distinct
        ok &= good
        report(good, f"N={n} profile diagonal is -2*gamma_j EXACTLY (== , dyadic rates) "
                     f"and carries {len(set(A.tolist()))} distinct values, not one")
    for n in (4, 5, 6):
        A, _ = fw.weight_block_pencil(n, 0, 1, gamma=0.5)
        good = bool(np.all(A == -1.0))
        ok &= good
        report(good, f"N={n} uniform gamma=0.5 collapses the same diagonal to the single value -1.0 (==)")
    return ok


def gate_mechanism_fails(report):
    """(4) |vac><psi_k| is a simultaneous eigenoperator only at uniform gamma."""
    ok = True
    for n in (4, 5, 6):
        for label, g, expect_eig in (
            ("uniform", np.full(n, 0.5), True),
            ("ramp", np.linspace(0.125, 1.0, n), False),
            ("one site detuned", np.concatenate(([1.0], np.full(n - 1, 0.5))), False),
        ):
            A, K = fw.weight_block_pencil(n, 0, 1, gamma=g)
            _, V = np.linalg.eigh(K)
            D = np.diag(A)
            worst = 0.0
            for k in range(n):
                v = V[:, k]
                Dv = D @ v
                worst = max(worst, float(np.linalg.norm(Dv - (v @ Dv) * v)))
            is_eig = worst <= TOL_FLAT
            good = is_eig == expect_eig
            ok &= good
            report(good, f"N={n} {label:<16} H-eigenvector residual under D = {worst:.2e} -> "
                         f"{'eigenoperator' if is_eig else 'NOT an eigenoperator'}")
    return ok


def gate_spectrum_three_valued(report):
    """(2) and (3): flat at the wrong constant, spread, or an empty set, depending on the profile.

    The flatness is not a property of the locus alone; it also needs the coupling. Do NOT read the
    rows below as a coupling THRESHOLD: the pinning group further down shows the mechanism is
    per-eigenvalue, so what looks like an onset is the last Im-degeneracy lifting. Every row here
    fixes J against the rate scale, and the third row of each N shows the same on-locus profile
    losing its flatness where degeneracies are still present.
    """
    ok = True
    for n in (4, 5, 6):
        bonds = chain(n)
        ramp = np.linspace(0.125, 1.0, n)
        gbar, gmin = float(ramp.mean()), float(ramp.min())
        d_ramp = r90_defect(ramp)

        # the ramp is on the locus, and says so rather than being assumed to be
        good = d_ramp <= 1e-12
        ok &= good
        report(good, f"N={n} the ramp is R90-FIXED (defect {d_ramp:.1e}), which is the hypothesis "
                     f"the next row reads")

        # (2) flat, at -2*gbar, and NOT at the floor
        spread, e = re_spread(band_edge_block(n, bonds, ramp, J=1.0, delta=0.0))
        at_mean = abs(float(e.mean()) - (-2 * gbar))
        off_floor = abs(-2 * gbar - (-2 * gmin))
        good = spread <= TOL_FLAT and at_mean <= TOL_FLAT and off_floor > 0.5
        ok &= good
        report(good, f"N={n} on-locus, XY, J=1 (gamma <= 1): FLAT ({spread:.1e}) but at "
                     f"-2*gbar={-2*gbar:.4f}, a distance {off_floor:.3f} off its floor "
                     f"-2*gmin={-2*gmin:.4f}")

        # negative control A: OFF the locus at the same coupling, the block spreads
        generic = ramp + np.array([0.0, 0.31, 0.0, 0.17, 0.0, 0.23][:n])
        d_gen = r90_defect(generic)
        spread_g, _ = re_spread(band_edge_block(n, bonds, generic, J=1.0, delta=0.0))
        good = d_gen > 1e-9 and spread_g > MIN_SPREAD
        ok &= good
        report(good, f"N={n} OFF-locus (R90 defect {d_gen:.3f}), same J and delta: spread "
                     f"{spread_g:.4f} > {MIN_SPREAD}")

        # negative control B: the SAME on-locus profile at weak coupling is not flat either,
        # so the flat row is a REGIME and not a property of the locus
        spread_w, _ = re_spread(band_edge_block(n, bonds, ramp, J=0.05, delta=0.0))
        good = spread_w > MIN_SPREAD
        ok &= good
        report(good, f"N={n} the same on-locus profile at J=0.05: spread {spread_w:.4f}, so the flat "
                     f"row needs more than the locus; per the pinning group that more is the lifting "
                     f"of the Im-degeneracies, not a coupling scale as such")

        # (3) the set at -2*gbar is EMPTY once a single site is detuned
        det = np.concatenate(([1.0], np.full(n - 1, 0.5)))
        _, e_det = re_spread(band_edge_block(n, bonds, det, J=1.0, delta=0.0))
        tgt_mean, tgt_min = -2 * float(det.mean()), -2 * float(det.min())
        n_mean = int(np.sum(np.abs(e_det - tgt_mean) < 1e-9))
        n_min = int(np.sum(np.abs(e_det - tgt_min) < 1e-9))
        good = n_mean == 0 and n_min == 0
        ok &= good
        report(good, f"N={n} one site detuned: {n_mean} of {n} modes at -2*gbar and {n_min} at "
                     f"-2*gmin, so the SET the claim quantifies over is empty")
    return ok


def gate_matches_committed_numbers(report):
    """The Heisenberg spreads F153 publishes, reproduced from the other instrument."""
    ok = True
    for n, want in ((4, 0.4903), (5, 0.9477)):
        g = np.linspace(0.1, 1.0, n)
        L = fw.weight_block_build(n, 0, 1, 1.0, gamma=g, delta=1.0)
        spread = float(np.ptp(np.linalg.eigvals(L).real))
        good = abs(spread - want) < 5e-5
        ok &= good
        report(good, f"N={n} gamma=linspace(0.1,1), Delta=1, J=1: spread {spread:.4f} against "
                     f"F153's published {want}")
    for n, want in ((4, 1.77),):
        g = np.linspace(0.1, 1.0, n)
        L = fw.weight_block_build(n, 0, 1, 0.05, gamma=g, delta=0.0)
        spread = float(np.ptp(np.linalg.eigvals(L).real))
        good = abs(spread - want) < 5e-3
        ok &= good
        report(good, f"N={n} same profile, Delta=0, J=0.05: spread {spread:.4f} against F153's {want}, "
                     f"the flat XY case opening as the coupling falls")
    return ok


def gate_topology_general(report):
    """The class-B sites say "any graph", so the three-valued reading is checked off the chain too."""
    ok = True
    n = 5
    for name, bonds in (("chain", chain(n)), ("ring", ring(n)), ("star", star(n)),
                        ("complete", complete(n))):
        uni = np.full(n, 0.5)
        spread_u, e_u = re_spread(band_edge_block(n, bonds, uni, J=1.0, delta=0.0))
        flat_at = abs(float(e_u.mean()) - (-1.0))
        generic = np.array([0.125, 0.75, 0.3125, 0.9375, 0.5])
        d_gen = r90_defect(generic)
        spread_p, _ = re_spread(band_edge_block(n, bonds, generic, J=1.0, delta=0.0))
        good = (spread_u <= TOL_FLAT and flat_at <= TOL_FLAT and spread_p > MIN_SPREAD
                and d_gen > 1e-9)
        ok &= good
        report(good, f"N={n} {name:<9}: uniform flat at -2*gamma ({spread_u:.1e}), profile spreads "
                     f"({spread_p:.4f}); this group is a UNIFORM-gamma reading plus one generic "
                     f"profile, and says nothing about the R90 locus off the chain, reversal not "
                     f"being an automorphism of every graph here")
    return ok


def gate_cross_instrument(report):
    """The bond-list builder here and the framework's chain builder must agree, or neither is evidence."""
    ok = True
    for n in (4, 5, 6):
        g = np.array([0.5, 0.25, 1.0, 0.125, 0.75, 0.375][:n])
        mine = band_edge_block(n, chain(n), g, J=1.0, delta=0.0)
        theirs = fw.weight_block_build(n, 0, 1, 1.0, gamma=g, delta=0.0)
        d = float(np.abs(np.sort(np.linalg.eigvals(mine)) - np.sort(np.linalg.eigvals(theirs))).max())
        good = d <= 1e-12
        ok &= good
        report(good, f"N={n} this script's bond-list block agrees with weight_block_build to {d:.1e}")
    return ok


def gate_axis_membership(report):
    """WHICH modes sit on the palindrome axis, and the four fences that took five rounds to state.

    NOTHING HERE IS A NEW LAW, and it belongs to a CLAUSE of F153 rather than a number of its own.
    Committed already: the identity (*) `P^T L P = -conj(L) - c*Id` of
    experiments/WHAT_THE_R90_LOCUS_BUYS.md, with algebraic multiplicities and no normality needed,
    and the reading that for real c the pairing reflects Re and leaves Im UNCHANGED; F153's own
    statement that the symmetry holds at every REAL coupling while flatness sets in only above one;
    and the generator itself (F152). The self-pairing step is owned too, as a typed claim:
    F89BranchLocusPalindromeClaim says of its own antiunitary pairing that "every EP lies on the line
    or in a mirror pair across it, no orphan". This group reads that step on this pairing:

        a mode whose Im is shared with no other DISTINCT eigenvalue is its own image, hence on the
        axis exactly, threshold or no threshold; flatness is the case where every mode is lone.

    Four fences, each of which a draft here got wrong, each with a row below:
      - the coupling must be REAL (F153's fence, not a new one);
      - the identity is (*)'s conj form, NOT R M R = -M^dagger - c, which agrees only because a real
        hopping makes M symmetric;
      - the sharp hypothesis on H is `P A P = A^T`, WEAKER than source 2's "H reflection-symmetric":
        a flux ring satisfies the former and violates the latter, and (*) holds for it;
      - the hypotheses are SUFFICIENT, NOT NECESSARY: a ring whose RATES violate the locus fails (*)
        and can still sit wholly on the axis, so a failed identity forbids nothing;
      - the loneness test needs a tolerance at or above eps*|lambda|, and the axis residual is
        O(1)*eps*|lambda|, a LAW rather than a fixed number, so it is gated as a ratio.
    """
    ok = True

    def axis_ratio(e, axis):
        scale = np.finfo(float).eps * max(1.0, float(np.abs(e).max()))
        return float(np.abs(np.asarray(e).real - axis).max()) / scale

    def lone_indices(e, tol=None):
        n = len(e)
        if tol is None:
            tol = 1e3 * np.finfo(float).eps * max(1.0, float(np.abs(e).max()))
        return [k for k in range(n)
                if all(abs(e[k].imag - e[j].imag) > tol for j in range(n) if j != k)]

    # (*) in its own form, against the dagger form that coincides only for a real hopping
    n = 6
    g = np.linspace(0.125, 1.0, n)
    gbar = float(g.mean())
    P = np.eye(n)[::-1]
    for label, phase in (("real hopping", None), ("Peierls phase 0.4", 0.4)):
        A = np.zeros((n, n), complex)
        for i in range(n - 1):
            w = np.exp(1j * phase) if phase is not None else 1.0
            A[i, i + 1] = w
            A[i + 1, i] = np.conj(w)
        M = -2j * 1.0 * A - 2.0 * np.diag(g)
        star = float(np.abs(P @ M @ P - (-M.conj() - 4 * gbar * np.eye(n))).max())
        dagger = float(np.abs(P @ M @ P - (-M.conj().T - 4 * gbar * np.eye(n))).max())
        good = star <= 1e-14 and (dagger <= 1e-14) == (phase is None)
        ok &= good
        report(good, f"N={n} {label:<18}: (*) conj form {star:.1e}, dagger form {dagger:.1e} "
                     f"-> the two agree only for a real hopping")

    # axis membership, per mode, at Delta = 0 and Delta = 1, across the span
    for delta in (0.0, 1.0):
        for n in (5, 7):
            g = np.linspace(0.125, 1.0, n)
            axis = -2 * float(g.mean())
            for J in (0.3, 0.6, 1.0, 3.0):
                e = np.linalg.eigvals(band_edge_block(n, chain(n), g, J=J, delta=delta))
                lone = lone_indices(e)
                ratio = axis_ratio(e[lone], axis) if lone else 0.0
                good = ratio <= 50.0
                ok &= good
                report(good, f"N={n} Delta={delta} J={J:<4}: {len(lone)} of {n} modes Im-lone, on the "
                             f"axis to {ratio:.1f}*eps*|lambda| (a ratio, not a threshold), block span "
                             f"{float(np.ptp(e.real)):.4f}")

    # fence: the coupling must be REAL
    n = 7
    g = np.linspace(0.125, 1.0, n)
    gbar = float(g.mean())
    P = np.eye(n)[::-1]
    A = np.zeros((n, n))
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1.0
    for J in (1.0 + 0.3j, 1.0 + 0.0j):
        M = -2j * J * A - 2.0 * np.diag(g)
        resid = float(np.abs(P @ M @ P - (-M.conj() - 4 * gbar * np.eye(n))).max())
        e = np.linalg.eigvals(M)
        lone = lone_indices(e)
        off = max((abs(e[k].real - (-2 * gbar)) for k in lone), default=0.0)
        real_J = abs(J.imag) < 1e-15
        good = (resid <= 1e-14) == real_J and (off <= 1e-9) == real_J
        ok &= good
        report(good, f"N={n} J={J}: identity {resid:.1e}, worst lone mode off axis {off:.2f} "
                     f"-> the REAL-coupling fence is F153's and is load-bearing")

    # the hypothesis is SHARPER than source 2 states it. Source 2 asks for H reflection-symmetric;
    # what (*) actually needs is P A P = A^T, which a real symmetric hopping satisfies via R H R = H
    # but which a flux ring ALSO satisfies while R H R != H. So the flux ring is inside the theorem,
    # not a counterexample to it, and an earlier row here mis-read it as one by testing the dagger
    # form instead of (*).
    n = 6
    g = np.linspace(0.2, 1.0, n)
    g = 0.5 * (g + (2 * g.mean() - g[::-1]))
    gbar = float(g.mean())
    A = np.zeros((n, n), complex)
    ph = np.exp(0.7j)
    for i in range(n):
        A[i, (i + 1) % n] = ph
        A[(i + 1) % n, i] = np.conj(ph)
    M = -2j * 1.0 * A - 2.0 * np.diag(g)
    P = np.eye(n)[::-1]
    h_sym = float(np.abs(P @ A @ P - A).max())
    transpose_cond = float(np.abs(P @ A @ P - A.T).max())
    resid = float(np.abs(P @ M @ P - (-M.conj() - 4 * gbar * np.eye(n))).max())
    e = np.linalg.eigvals(M)
    lone = lone_indices(e)
    ratio = axis_ratio(e[lone], -2 * gbar) if lone else float("inf")
    good = h_sym > 1e-2 and transpose_cond <= 1e-14 and resid <= 1e-14 and ratio <= 50.0
    ok &= good
    report(good, f"N={n} flux ring: ||RHR-H|| = {h_sym:.3f} yet ||PAP-A^T|| = {transpose_cond:.1e}, "
                 f"so (*) HOLDS ({resid:.1e}) and all modes sit on the axis to "
                 f"{ratio:.1f}*eps*|lambda| -> the sharp hypothesis is P A P = A^T, weaker than "
                 f"source 2's 'H reflection-symmetric'")

    # fence: the hypotheses are SUFFICIENT, NOT NECESSARY. Break the RATE hypothesis instead: the
    # identity then genuinely fails and the block can still sit wholly on the axis.
    n = 4
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = 1.0
    g = np.array([0.75, 0.5, 0.25, 0.5])
    gbar = float(g.mean())
    P = np.eye(n)[::-1]
    M = -2j * 0.5 * A - 2.0 * np.diag(g)
    r90 = float(np.abs(g - (2 * g.mean() - g[::-1])).max())
    resid = float(np.abs(P @ M @ P - (-M.conj() - 4 * gbar * np.eye(n))).max())
    e = np.linalg.eigvals(M)
    off = float(np.abs(e.real - (-2 * gbar)).max())
    good = r90 > 1e-2 and resid > 1e-2 and off <= 1e-13
    ok &= good
    report(good, f"N={n} ring, R90 defect {r90:.2f} (rate hypothesis VIOLATED): (*) fails at "
                 f"{resid:.2f} and yet every mode sits on the axis to {off:.1e} -> the hypotheses "
                 f"are SUFFICIENT, not necessary; a failed identity forbids nothing")

    # fence: the loneness predicate is not decidable at exact equality
    n = 4
    g = np.array([0.9, 0.65, 0.35, 0.1])
    e = np.linalg.eigvals(band_edge_block(n, chain(n), g, J=0.3, delta=0.0))
    strict = [k for k in range(n) if all(e[k].imag != e[j].imag for j in range(n) if j != k)]
    tolerant = lone_indices(e)
    off_strict = max((abs(e[k].real - (-2 * float(g.mean()))) for k in strict), default=0.0)
    good = len(strict) == n and len(tolerant) == 0 and off_strict > 0.1
    ok &= good
    report(good, f"N={n} two Im values differ by 7e-16: a STRICT test calls {len(strict)} of {n} lone "
                 f"and the axis claim appears to fail by {off_strict:.3f}; with an eps-scaled "
                 f"tolerance {len(tolerant)} are lone and nothing is claimed")

    # the open question, executed rather than asserted
    rng = np.random.default_rng(20260822)
    draws = 3000
    hits = 0
    for _ in range(draws):
        m = int(rng.integers(4, 8))
        gg = rng.uniform(0.05, 1.5, size=m)
        gg = 0.5 * (gg + (2 * gg.mean() - gg[::-1]))
        if float(np.ptp(gg)) < 1e-6:
            continue
        J = float(rng.choice([0.5, 1.0, 2.0, 5.0]))
        e = np.linalg.eigvals(band_edge_block(m, chain(m), gg, J=J, delta=0.0))
        if len(lone_indices(e)) < m and axis_ratio(e, -2 * float(gg.mean())) <= 50.0:
            hits += 1
    good = hits == 0
    ok &= good
    report(good, f"the open question, searched here rather than asserted: {draws} non-uniform "
                 f"on-locus draws produced {hits} degenerate-Im blocks wholly on the axis")
    return ok


def main():
    passed = failed = 0

    def report(good, msg):
        nonlocal passed, failed
        if good:
            passed += 1
        else:
            failed += 1
        print(f"  [{'PASS' if good else 'FAIL'}] {msg}")

    ok = True
    for title, fn in (
        ("cross-instrument: the two builders agree", gate_cross_instrument),
        ("(1) the diagonal is the site-resolved sum, exactly", gate_diagonal_is_exact),
        ("(4) the eigenoperator mechanism fails under any non-uniform profile", gate_mechanism_fails),
        ("(2)+(3) the spectrum: flat at the wrong constant, spread, or an empty set",
         gate_spectrum_three_valued),
        ("the committed numbers F153 publishes, from the other instrument",
         gate_matches_committed_numbers),
        ("off the chain: ring, star and complete read the same way", gate_topology_general),
        ("which modes sit on the axis (a CLAUSE of F153, not a law of its own)",
         gate_axis_membership),
    ):
        print(f"\n{title}")
        ok &= fn(report)

    print("\n" + "=" * 90)
    print(f"{passed}/{passed + failed} gates passed" if ok else
          f"{failed} FAILURES ({passed} passed)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
