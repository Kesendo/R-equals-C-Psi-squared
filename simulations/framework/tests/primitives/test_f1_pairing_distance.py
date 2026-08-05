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
    H = rng.normal(size=(d, d))
    H = H + H.T
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
    assert max(vals) / min(vals) < 10.0, "but it stays inside the solver floor"


def test_zero_spectrum_at_nonzero_sigma_is_a_maximal_violation():
    """The trap: a zero spectrum is the WORST case, not a clean one.

    Every eigenvalue at 0 puts the whole reflection at -2*sigma, so the
    distance is 2|sigma|. An earlier version of the guard returned 0.0 here,
    which would have reported a perfect score for a maximal violation.
    """
    sigma = 0.3199
    dist, radius = fw.f1_distance_in_eps(np.zeros(4, dtype=complex), sigma)
    assert radius == 0.0
    assert dist == pytest.approx(2 * sigma, rel=0, abs=0)


def test_zero_spectrum_at_zero_sigma_is_zero():
    dist, radius = fw.f1_distance_in_eps(np.zeros(4, dtype=complex), 0.0)
    assert radius == 0.0
    assert dist == 0.0
