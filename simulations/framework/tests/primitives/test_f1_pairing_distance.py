"""Tests for the F1 pairing distance primitive.

The point of this file is not that the metric works. It is that the metric's
BLIND SPOTS are pinned, because the greedy palindrome scorers this primitive
replaced were wrong in exactly the way an unpinned metric goes wrong: they
measured their own matcher, and every table they fed read as evidence.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def _palindromic_spectrum(sigma, n=8, seed=0):
    """A spectrum exactly closed under lambda -> -2*sigma - lambda."""
    rng = np.random.default_rng(seed)
    half = rng.normal(size=n) + 1j * rng.normal(size=n)
    return np.concatenate([half, -2.0 * sigma - half])


def test_exact_palindrome_at_zero_sigma_is_exactly_zero():
    """At sigma = 0 the reflection is exactly involutory in floating point."""
    spec = _palindromic_spectrum(0.0)
    assert fw.max_f1_pairing_distance(spec, 0.0) == 0.0


def test_residual_at_nonzero_sigma_is_the_double_negation_rounding():
    """A perfectly palindromic spectrum does NOT give 0.0 once sigma != 0, and
    the leftover is readable rather than something to gate.

    The metric forms `reflected = -2*sigma - spectrum`. A partner constructed as
    y = -2*sigma - x therefore reflects back to -2*sigma - y, which is x only up
    to one rounding of the subtraction. So the residual is a deterministic
    function of sigma and of the spectrum, computable exactly here and compared
    exactly. It is exactly 0.0 at sigma = 0 (the test above), which is what
    distinguishes an informative rounding route from a broken one.
    """
    for sigma in (0.25, 3.0):
        spec = _palindromic_spectrum(sigma)
        reflected = -2.0 * sigma - spec
        expected = max(abs(x - (-2.0 * sigma - y))
                       for x, y in zip(spec, reflected))
        assert fw.max_f1_pairing_distance(spec, sigma) == expected
        # and it is genuinely at the floor, not a real asymmetry
        assert expected < 8.0 * np.finfo(float).eps * max(abs(spec))


def test_wrong_sigma_is_detected():
    sigma = 0.25
    spec = _palindromic_spectrum(sigma)
    assert fw.max_f1_pairing_distance(spec, sigma + 1e-6) > 1e-7


def test_removal_is_what_catches_a_contested_partner():
    """The one test that fails if `taken[]` is dropped.

    Written after a mutation pass showed that the two tests below both pass
    against a matcher WITHOUT removal, i.e. neither pinned the property the
    docstring calls the metric's reason to exist. Removal is what makes it
    multiset-aware: here two copies of x compete for a single mirror point,
    so a matcher that lets both claim it reports ~0 and a matcher that
    consumes it reports the full mirror separation.
    """
    sigma = 0.25
    x = 1.0 + 0.5j
    spec = np.array([x, x, -2.0 * sigma - x])
    d = fw.max_f1_pairing_distance(spec, sigma)
    # without removal every point finds a partner at distance ~0
    reflected = -2.0 * sigma - spec
    no_removal = max(float(np.min(np.abs(p - reflected))) for p in spec)
    assert no_removal < 1e-12
    assert d > 1.0, "removal must leave the second copy of x unmatched"


def test_asymmetric_duplication_is_caught():
    sigma = 0.25
    spec = _palindromic_spectrum(sigma)
    bad = spec.copy()
    bad[3] = spec[0]          # duplicate one eigenvalue onto an unrelated slot
    assert fw.max_f1_pairing_distance(bad, sigma) > 0.1


def test_deletion_is_caught():
    sigma = 0.25
    spec = _palindromic_spectrum(sigma)
    assert fw.max_f1_pairing_distance(np.delete(spec, 3), sigma) > 0.1


def test_blind_to_a_reflection_respecting_multiplicity_error():
    """The documented blind spot, pinned so it cannot be forgotten.

    Dropping one mirror pair and duplicating another keeps the multiset closed
    under the reflection, so the metric cannot see it. This is what a BLOCK
    solver would produce, which is why this number does not certify a block
    spectrum against a dense one.
    """
    sigma = 0.25
    n = 8
    spec = _palindromic_spectrum(sigma, n=n)
    bad = spec.copy()
    bad[1], bad[n + 1] = spec[0], spec[n]     # pair 1 becomes a copy of pair 0
    assert fw.max_f1_pairing_distance(bad, sigma) == fw.max_f1_pairing_distance(spec, sigma)


def test_non_finite_raises_rather_than_understating():
    sigma = 0.25
    spec = _palindromic_spectrum(sigma)
    bad = spec.copy()
    bad[2] = np.nan
    with pytest.raises(ValueError):
        fw.max_f1_pairing_distance(bad, sigma)


def test_greedy_matcher_is_order_dependent_on_a_real_spectrum():
    """Not a defect of this port, a property of greedy matching.

    Pinned because the whole campaign that produced this primitive started from
    a score that was measuring its matcher. Permuting only the array order
    leaves the spectrum identical, so any spread here is bookkeeping, and two
    systems must not be ranked by differences of this size.
    """
    N = 4
    d = 2 ** N
    rng = np.random.default_rng(1)
    # A HEISENBERG chain, not a random symmetric matrix. F1 is a theorem about
    # the Heisenberg/XXZ generator under Z-dephasing; a generic H gives a
    # spectrum that is not palindromic at all, and an earlier version of this
    # test measured the order-spread of that meaningless number instead.
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Zp = np.array([[1, 0], [0, -1]], dtype=complex)

    def at(op, site):
        m = np.eye(1, dtype=complex)
        for s in range(N):
            m = np.kron(m, op if s == site else np.eye(2, dtype=complex))
        return m

    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y, Zp):
            H = H + at(P, b) @ at(P, b + 1)
    # A FIXTURE, not a calibration. These are the raw 1/T2* readings of the
    # Q85-Q94 line; a Z-dephasing model would take gamma = 1/(2*T2*) instead
    # (docs/GLOSSARY.md, the T2 -> gamma table). Both assertions below are
    # scale-free, so the book does not enter: one compares the spread to
    # itself, the other to a corruption of the same spectrum. Do not cite
    # these numbers as the profile, and do not halve them expecting a change.
    gammas = np.array([0.2681, 0.0163, 0.0103, 0.0147])[:N]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    for k in range(N):
        op = np.eye(1, dtype=complex)
        for m in range(N):
            op = np.kron(op, Z if m == k else np.eye(2, dtype=complex))
        Lk = np.sqrt(gammas[k]) * op
        LdL = Lk.conj().T @ Lk
        L = L + np.kron(Lk.conj(), Lk) - 0.5 * np.kron(Id, LdL) - 0.5 * np.kron(LdL.T, Id)
    evals = np.linalg.eigvals(L)
    sigma = float(sum(gammas))
    vals = [fw.f1_distance_in_eps(rng.permutation(evals), sigma)[0] for _ in range(40)]
    assert max(vals) > min(vals), "order dependence is a property of the matcher"
    # The bound that means something: the whole order-induced spread must stay
    # far below what an actual corruption produces, or the metric could not
    # separate bookkeeping from a real violation. A fixed numeric ceiling here
    # would be a threshold that cannot fail; this one is measured against the
    # signal it has to beat.
    corrupted = evals.copy()
    corrupted[3] = evals[0]
    signal = fw.f1_distance_in_eps(corrupted, sigma)[0]
    assert max(vals) < signal / 1000.0


def test_zero_spectrum_at_nonzero_sigma_is_a_maximal_violation():
    """The trap: a zero spectrum is the WORST case, not a clean one.

    Every eigenvalue at 0 puts the whole reflection at -2*sigma, so the
    distance is 2|sigma|. An earlier version of the guard returned 0.0 here,
    which would have reported a perfect score for a maximal violation.
    """
    # A synthetic centre, not the Q85-Q94 profile: on that line 0.3199 is
    # 2*Sigma gamma, not Sigma gamma. Nothing here reads the hardware.
    sigma = 0.3199
    dist, radius = fw.f1_distance_in_eps(np.zeros(4, dtype=complex), sigma)
    assert radius == 0.0
    assert dist == pytest.approx(2 * sigma, rel=0, abs=0)


def test_zero_spectrum_at_zero_sigma_is_zero():
    dist, radius = fw.f1_distance_in_eps(np.zeros(4, dtype=complex), 0.0)
    assert radius == 0.0
    assert dist == 0.0
