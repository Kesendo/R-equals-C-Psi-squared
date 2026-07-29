"""Tests for F77 trichotomy classifier (truly / soft / hard) at k=2 and k≥3."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F77_classify_pauli_pair_k2_known_cases():
    """k=2 trichotomy on known V-Effect cases."""
    chain = fw.ChainSystem(N=4)

    # Truly: M = 0
    assert fw.classify_pauli_pair(chain, [('X', 'X'), ('Y', 'Y')]) == 'truly'
    assert fw.classify_pauli_pair(chain, [('X', 'X'), ('Z', 'Z')]) == 'truly'

    # Soft: M ≠ 0, eigenvalue pairing intact
    assert fw.classify_pauli_pair(chain, [('X', 'Y'), ('Y', 'X')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('X', 'Z'), ('Z', 'X')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('Y', 'Z'), ('Z', 'Y')]) == 'soft'

    # Hard: eigenvalue pairing broken
    assert fw.classify_pauli_pair(chain, [('X', 'Y'), ('X', 'Z')]) == 'hard'
    assert fw.classify_pauli_pair(chain, [('X', 'X'), ('X', 'Y')]) == 'hard'


def test_F77_classify_pauli_pair_k3_kbody_dispatch():
    """k=3 chain Hamiltonians use the chain sliding-window k-body builder.

    The I-padded cases below are NOT an equivalence check. An identity in a
    length-3 term does not reduce it to the 2-body bilinear: on an N=4 chain the
    k=2 form covers three bonds while the I-padded form covers only the two
    windows l=0,1. For the first pair below that is ‖H_kbody − H_bilinear‖_F =
    5.657 against ‖H_bilinear‖_F = 9.798.
    What they pin is weaker and still worth pinning: the class survives dropping
    the boundary bond.
    """
    chain = fw.ChainSystem(N=4)

    # Same class as 2-body XY+YX, on two windows instead of three: soft
    assert fw.classify_pauli_pair(chain, [('X', 'Y', 'I'), ('Y', 'X', 'I')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('I', 'X', 'Y'), ('I', 'Y', 'X')]) == 'soft'

    # Same class as 2-body XX+YY, on two windows instead of three: truly
    assert fw.classify_pauli_pair(chain, [('X', 'X', 'I'), ('Y', 'Y', 'I')]) == 'truly'
    assert fw.classify_pauli_pair(chain, [('I', 'X', 'X'), ('I', 'Y', 'Y')]) == 'truly'

    # Genuine 3-body XYZ pair: empirically soft (Klein-homogeneous)
    assert fw.classify_pauli_pair(chain, [('X', 'Y', 'Z'), ('Y', 'Z', 'X')]) == 'soft'

    # All-truly 3-body: M = 0 → truly
    assert fw.classify_pauli_pair(chain, [('X', 'X', 'X'), ('Y', 'I', 'Y')]) == 'truly'


def test_F77_classify_pauli_pair_k3_klein_homogeneity_not_strict():
    """Klein-homogeneity does not settle the class. Some Klein-homogeneous k=3
    Hamiltonians are F87-hard.

    Empirical verification: the 294-pair Z₂³-homogeneous sweep at k=3 N=4
    (PROOF_F103_F87_Z2_CUBED_REFINEMENT §3.2) finds 50 hard in the Klein-(0,1)
    cell under Z-dephasing. The two-letter case is not a clean contrast: the 6
    identity-free Klein-homogeneous pairs are 0/6 hard, but the 15 carrying an
    identity are 5/15 hard.

    The assertion below is the counterexample itself: a Klein-homogeneous k=3
    pair that classifies hard, contradicting the blanket
    "Klein-homogeneous ⟹ never hard" reading.
    """
    chain = fw.ChainSystem(N=4)

    # Both strings sit in the same Klein cell, and the pair is F77-hard.
    hard_cases = [
        [('I', 'I', 'Z'), ('I', 'X', 'Y')],
        [('I', 'I', 'Z'), ('X', 'I', 'Y')],
    ]
    for terms in hard_cases:
        klein = {fw.klein_index(t) for t in terms}
        assert len(klein) == 1, f"{terms}: not Klein-homogeneous, cells {klein}"
        assert fw.classify_pauli_pair(chain, terms) == 'hard', \
            f"k=3 {terms} is Klein-homogeneous and must classify hard"

    # Klein-homogeneity does not force hard either: these are soft.
    for terms in [[('X', 'Y', 'Z'), ('X', 'Z', 'Y')], [('X', 'X', 'Y'), ('Y', 'Y', 'Y')]]:
        assert fw.classify_pauli_pair(chain, terms) == 'soft'


def test_F77_classify_mixed_body_count():
    """Hamiltonian with mixed body counts (one k=2 term + one k=3 term) builds H
    by summing both contributions and classifies the combined L.

    Scope: the k=2 part is built on `chain.bonds` (topology-aware) while the k≥3
    part goes through the sliding-window chain builder, which ignores topology.
    On a non-chain topology the two halves therefore sit on different graphs, so
    this test pins the chain case only.
    """
    chain = fw.ChainSystem(N=4)

    # k=2 XY and k=3 IXY are each soft alone; the summed Hamiltonian stays soft.
    assert fw.classify_pauli_pair(chain, [('X', 'Y')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('I', 'X', 'Y')]) == 'soft'
    assert fw.classify_pauli_pair(chain, [('X', 'Y'), ('I', 'X', 'Y')]) == 'soft'


def test_F77_classify_validation_errors():
    """Classifier raises ValueError on invalid term body counts."""
    chain = fw.ChainSystem(N=3)

    with pytest.raises(ValueError):
        fw.classify_pauli_pair(chain, [('X',)])  # k=1 < 2

    with pytest.raises(ValueError):
        fw.classify_pauli_pair(chain, [('X', 'Y', 'Z', 'I')])  # k=4 > N=3


def test_F77_classify_empty_terms():
    """Empty term list is trivially truly (zero Hamiltonian)."""
    chain = fw.ChainSystem(N=3)
    assert fw.classify_pauli_pair(chain, []) == 'truly'


# ----------------------------------------------------------------------
# The spectrum-pairing error: exact bottleneck, not the old greedy pass.
# ----------------------------------------------------------------------

def test_pairing_error_is_zero_on_an_exactly_mirrored_spectrum():
    """A spectrum built to be closed under λ ↦ −λ − 2Σγ must score 0."""
    from framework.diagnostics.f77_trichotomy import spectrum_pairing_error
    rng = np.random.default_rng(2)
    sigma = 0.17
    half = rng.normal(size=9) + 1j * rng.normal(size=9)
    evals = np.concatenate([half, -half - 2 * sigma])
    assert spectrum_pairing_error(evals, sigma) < 1e-12


def test_pairing_error_never_exceeds_the_greedy_pass():
    """The exact bottleneck is a minimum over bijections; greedy realises one of
    them and additionally forces an involution, so greedy can only be larger.

    This is what makes a greedy classifier conservative (soft may be reported as
    hard, never the reverse), and it is the relation the search relies on to use
    the greedy value as its upper bound.
    """
    from framework.diagnostics.f77_trichotomy import (
        spectrum_pairing_error, _greedy_pairing_error,
    )
    rng = np.random.default_rng(11)
    seen_strict = False
    for _ in range(60):
        n = int(rng.integers(4, 12))
        sigma = float(rng.uniform(0.05, 0.4))
        evals = rng.normal(size=n) + 1j * rng.normal(size=n)
        exact = spectrum_pairing_error(evals, sigma)
        greedy = _greedy_pairing_error(evals, sigma)
        assert exact <= greedy + 1e-12, f"exact {exact} exceeded greedy {greedy}"
        if exact < greedy - 1e-9:
            seen_strict = True
    assert seen_strict, "greedy never lost; the comparison would be vacuous"


def test_pairing_error_detects_a_broken_mirror():
    """Displace one eigenvalue and the bottleneck must see it, at its own size."""
    from framework.diagnostics.f77_trichotomy import spectrum_pairing_error
    sigma = 0.1
    half = np.array([0.3 + 0.2j, -1.1 + 0.7j, 0.05 - 0.9j])
    evals = np.concatenate([half, -half - 2 * sigma])
    evals[0] += 0.25
    err = spectrum_pairing_error(evals, sigma)
    assert 0.1 < err < 0.6, f"expected an error of order the 0.25 displacement, got {err}"


def test_pairing_error_return_always_admits_a_perfect_matching():
    """The returned ε must be a threshold at which a bijection actually exists.

    Regression: D was built as |λ_i − (−λ_j − 2σ)| while the greedy bound used
    |λ_j − (−λ_i − 2σ)|. Algebraically equal, ~1 ulp apart in IEEE, so the prune
    `D <= greedy` could drop the very edge greedy had used; no candidate was then
    feasible and the binary search returned its last candidate untested. Both are
    written as sums now, which is exactly symmetric. A 4000-trial sweep of this
    shape found 27 failures before the fix; the 400 below are the same shape.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching

    from framework.diagnostics.f77_trichotomy import spectrum_pairing_error

    rng = np.random.default_rng(20260729)
    for _ in range(400):
        n = int(rng.integers(4, 10))
        # integer-valued spectra force exact ties, where the ulp gap bites
        evals = (rng.integers(-4, 5, size=n) + 1j * rng.integers(-4, 5, size=n)).astype(complex)
        sigma = float(rng.choice([0.1, 0.15, 0.3]))
        eps = spectrum_pairing_error(evals, sigma)
        D = np.abs(evals[:, None] + evals[None, :] + 2 * sigma)
        matching = maximum_bipartite_matching(csr_matrix(D <= eps), perm_type='column')
        assert np.all(matching >= 0), \
            f"returned eps={eps!r} admits no perfect matching (n={n}, sigma={sigma})"
