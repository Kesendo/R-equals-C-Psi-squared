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


def gate_locus_pinning_criterion(report):
    """The exact law behind the flat rows above, on the (0,1) block, PER EIGENVALUE.

    Not an "iff": one direction is proved, the converse is false, and F153 uses "criterion" for an iff
    so the word is avoided here.

    The IDENTITY is the repo's, not this script's: it is source 2 of
    experiments/WHAT_THE_R90_LOCUS_BUYS.md, and it carries TWO hypotheses, not one: gamma
    anti-palindromic under a site reflection AND H itself reflection-symmetric under the same
    reflection. The second is load-bearing and is exercised below, because dropping a hypothesis at
    the citation is the very failure this whole script exists to fence. Source 2 is gated at its home
    by simulations/two_sources_gate.py; the collapse it produces is that note's measured table.

    What is added here is why the collapse is EXACT, and the argument is per-eigenvalue:

        the pairing lambda -> -conj(lambda) - 4*gbar preserves Im lambda and reflects Re lambda
        about -2*gbar, so ANY eigenvalue whose Im is shared with no other is its own image and has
        Re lambda = -2*gbar exactly.

    Two consequences the first draft of this group got wrong and which the gate now measures. The
    statement is about ONE eigenvalue at a time, so pinning arrives MODE BY MODE as Im-degeneracies
    lift; there is no single J onset, and a block can carry pinned modes while its overall Re-spread
    is large (measured below: three of seven pinned to 1e-15 at a spread of 0.87). And it is ONE
    DIRECTION: the converse is false, a star at uniform gamma being pinned while N-2 of its modes
    share Im = 0, exhibited in the last rows, and false only trivially, the diagonal being scalar
    there for an elementary reason.

    This is not the route that note rejects, which argued from the palindrome plus the one-sided
    Absorption bound and merged FLAT with ON THE FLOOR; nothing here mentions the floor, and the axis
    sits a long way off it. It is not THE_SPREAD_IS_A_RESONANCE's X^N cancellation either, which
    needs an index at N/2 that (0,1) never has for N > 2.
    """
    ok = True

    # source 2's SECOND hypothesis: H reflection-symmetric. Dropping it breaks the identity.
    n = 7
    g = np.linspace(0.125, 1.0, n)
    gbar = float(g.mean())
    R = np.eye(n)[::-1]
    for label, w, want_identity in (
        ("palindromic bonds", [1.0, 1.25, 1.5, 1.5, 1.25, 1.0], True),
        ("ASYMMETRIC bonds", [1.0, 1.25, 1.5, 0.8, 1.1, 0.6], False),
    ):
        bonds = [(i, i + 1, w[i]) for i in range(n - 1)]
        M = band_edge_block(n, bonds, g, J=1.0, delta=0.0)
        resid = float(np.abs(R @ M @ R - (-M.conj().T - 4 * gbar * np.eye(n))).max())
        spread = float(np.ptp(np.linalg.eigvals(M).real))
        holds = resid <= 1e-14
        good = holds == want_identity
        ok &= good
        report(good, f"N={n} on-locus gamma, {label:<18}: identity residual {resid:.2e}, Re spread "
                     f"{spread:.4f} -> source 2 {'holds' if holds else 'FAILS'}, so H "
                     f"reflection-symmetric is a hypothesis and not decoration")

    # the criterion, per eigenvalue, across the region where the block is NOT flat
    for n in (5, 6, 7):
        g = np.linspace(0.125, 1.0, n)
        gbar = float(g.mean())
        bonds = chain(n)
        for J in (0.3, 0.45, 0.6, 0.78, 1.0, 3.0):
            e = np.linalg.eigvals(band_edge_block(n, bonds, g, J=J, delta=0.0))
            lone = [k for k in range(n)
                    if all(abs(e[k].imag - e[j].imag) > 1e-9 for j in range(n) if j != k)]
            worst = max((abs(e[k].real + 2 * gbar) for k in lone), default=0.0)
            spread = float(np.ptp(e.real))
            good = worst <= TOL_FLAT
            ok &= good
            report(good, f"N={n} J={J:<5}: {len(lone)} of {n} eigenvalues Im-lone, all pinned to "
                         f"{worst:.1e}, while the BLOCK spread is {spread:.4f}")

    # The converse is false, but only TRIVIALLY, and saying so is the point. The witness is uniform
    # gamma, where the diagonal is already the scalar -2*gamma*Id, so Re = -2*gamma holds for an
    # elementary reason and not through the pairing: this is an IDENTITY stated, not a measurement,
    # in the shape two_sources_gate.py uses for the same situation. A star is convenient because its
    # hopping spectrum is degenerate (N-2 modes at Im = 0; the spectrum is {-2J*sqrt(N-1), 0,
    # +2J*sqrt(N-1)}, three distinct values, NOT "fully degenerate").
    # OPEN, and searched: whether a NON-uniform on-locus profile can be pinned with a degenerate Im
    # spectrum. 48000 draws over chain, ring and complete at N = 4..7 and four couplings found none.
    # So the converse fails only where the conclusion has an independent reason.
    for n in (5, 6):
        g = np.full(n, 0.5)
        e = np.linalg.eigvals(band_edge_block(n, star(n), g, J=1.0, delta=0.0))
        im = np.sort(e.imag)
        min_gap = float(np.min(np.diff(im)))
        n_distinct = len(set(np.round(im, 9).tolist()))
        max_dev = float(np.abs(e.real - (-2 * float(g.mean()))).max())
        good = min_gap <= 1e-9 and max_dev <= TOL_FLAT and n_distinct == 3
        ok &= good
        report(good, f"N={n} star at UNIFORM gamma: {n_distinct} distinct Im values, min gap "
                     f"{min_gap:.1e}, pinned to {max_dev:.1e} -> the converse is false, trivially, "
                     f"the diagonal being scalar here (identity, not a measurement)")

    # and the star is not a source-2 graph at all: reversal is not one of its automorphisms
    n = 5
    gs = np.array([0.2, 0.4, 0.5, 0.6, 0.8])
    gs = 0.5 * (gs + (2 * gs.mean() - gs[::-1]))
    M = band_edge_block(n, star(n), gs, J=1.0, delta=0.0)
    resid = float(np.abs(np.eye(n)[::-1] @ M @ np.eye(n)[::-1]
                         - (-M.conj().T - 4 * float(gs.mean()) * np.eye(n))).max())
    good = resid > 1e-2
    ok &= good
    report(good, f"N={n} star with reversal-anti-palindromic gamma: identity residual {resid:.3f}, "
                 f"because reversal is not an automorphism of a star; r90_defect is the locus test "
                 f"only for graphs the reflection preserves")
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
        ("the locus pinning law: Im-lone => pinned to -2*gbar (one direction; converse false)",
         gate_locus_pinning_criterion),
    ):
        print(f"\n{title}")
        ok &= fn(report)

    print("\n" + "=" * 90)
    print(f"{passed}/{passed + failed} gates passed" if ok else
          f"{failed} FAILURES ({passed} passed)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
