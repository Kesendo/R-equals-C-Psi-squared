"""Tests for F80 Bloch sign-walk and predict_M_spectrum_pi2_odd / k-body identity."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F80_bloch_signwalk_chain_pi2_odd():
    """F80: chain Π²-odd 2-body M-cluster values follow the open-chain
    free-fermion Bloch sign-walk formula (γ-independent by Master Lemma):

        cluster(N) = 2|c|·|Σ_{k=1..⌊N/2⌋} σ_k · 2cos(πk/(N+1))|

    Equivalent direct identity (discovered 2026-04-29 via data sweep):

        Spec(M)_{nontrivial} = ±2i · Spec(H)_{nontrivial, many-body}

    M's spectrum is directly 2i times the chain Hamiltonian's many-body
    eigenvalues. The Bloch sign-walk form is just the free-fermion many-body
    spectrum written out using Bogoliubov mode energies E_k = 4|c|·cos(πk/(N+1)).

    Verified at N=4, 5 (small enough for fast pytest); N=6, 7 verified
    in scripts (see _pi2_odd_universality_data_sweep.py and
    _n7_bloch_signwalk_verification.txt).
    """
    from itertools import product as iproduct
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_bilinear

    def predict(N, c=1.0):
        eps = [2.0 * np.cos(np.pi * k / (N + 1)) for k in range(1, N // 2 + 1)]
        sign_combos = list(iproduct([1, -1], repeat=len(eps)))
        sums = [abs(sum(s * e for s, e in zip(sigs, eps))) for sigs in sign_combos]
        # Distinct values
        distinct = []
        for v in sums:
            if not any(abs(v - d) < 1e-9 for d in distinct):
                distinct.append(v)
        return sorted([2 * c * v for v in distinct], reverse=True)

    for N in [4, 5]:
        bonds = [(i, i + 1) for i in range(N - 1)]
        # Test all 4 Π²-odd Pauli pairs (universality)
        for (P, Q) in [('X', 'Y'), ('X', 'Z'), ('Y', 'X'), ('Z', 'X')]:
            H = _build_bilinear(N, bonds, [(P, Q, 1.0)])
            L = lindbladian_z_dephasing(H, [1.0] * N)
            M = palindrome_residual(L, N * 1.0, N)
            svs = np.linalg.svd(M, compute_uv=False)

            # Distinct cluster values (above zero)
            observed = []
            for s in svs:
                if s > 1e-6 and not any(abs(s - o) < 1e-5 for o in observed):
                    observed.append(s)
            observed = sorted(observed, reverse=True)

            predicted = predict(N)
            assert len(observed) == len(predicted), \
                f"N={N} ({P},{Q}): {len(observed)} observed clusters vs {len(predicted)} predicted"
            for o, p in zip(observed, predicted):
                assert abs(o - p) < 1e-6, \
                    f"N={N} ({P},{Q}): observed cluster {o} vs predicted {p}"

            # Verify multiplicity 4^N / num_distinct (excluding zero)
            n_distinct = len(predicted)
            expected_mult = (4 ** N) // n_distinct
            for pred in predicted:
                actual_mult = int(np.sum(np.abs(svs - pred) < 1e-6))
                assert actual_mult == expected_mult, \
                    f"N={N} ({P},{Q}) cluster {pred}: mult {actual_mult} vs {expected_mult}"


def test_F80_predict_M_spectrum_pi2_odd_method():
    """ChainSystem.predict_M_spectrum_pi2_odd reproduces actual M's spectrum
    bit-exact for chain Π²-odd 2-body bilinears.

    Verifies the F80 structural identity Spec(M) = ±2i · Spec(H_non-truly)
    with multiplicity ×2^N: prediction (computed from H eigenvalues only)
    matches numerical M-eigenvalues from palindrome_residual.

    Also covers the trichotomy edge cases: truly-only returns {0: 4^N},
    identity letters and Π²-even non-truly bilinears raise ValueError.
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_bilinear

    def actual_M_spectrum(N, bonds, terms_with_c):
        H = _build_bilinear(N, bonds, terms_with_c)
        L = lindbladian_z_dephasing(H, [0.0] * N)  # γ=0 isolates structural M
        M = palindrome_residual(L, 0.0, N)
        evs = np.linalg.eigvals(M)
        out = {}
        for ev in evs:
            assert abs(ev.real) < 1e-7, f"M eigenvalue must be purely imaginary, got {ev}"
            key = round(ev.imag, 6)
            out[key] = out.get(key, 0) + 1
        return out

    def normalize_pred(pred):
        return {round(k.imag, 6): v for k, v in pred.items()}

    # Test 1-4: chain N=3 various Π²-odd cases at γ=0
    chain3 = fw.ChainSystem(N=3)
    bonds3 = [(0, 1), (1, 2)]

    for label, terms, c in [
        ('XY+YX', [('X', 'Y'), ('Y', 'X')], 1.0),
        ('XY', [('X', 'Y')], 1.0),
        ('ZX', [('Z', 'X')], 1.0),
        ('XY c=0.5', [('X', 'Y')], 0.5),
    ]:
        pred = normalize_pred(chain3.predict_M_spectrum_pi2_odd(terms, c=c))
        actual = actual_M_spectrum(3, bonds3, [(a, b, c) for (a, b) in terms])
        assert pred == actual, f"N=3 {label}: pred {pred} vs actual {actual}"

    # Test 5-6: chain N=4
    chain4 = fw.ChainSystem(N=4)
    bonds4 = [(0, 1), (1, 2), (2, 3)]

    for label, terms in [
        ('XY+YX', [('X', 'Y'), ('Y', 'X')]),
        ('XY', [('X', 'Y')]),
    ]:
        pred = normalize_pred(chain4.predict_M_spectrum_pi2_odd(terms, c=1.0))
        actual = actual_M_spectrum(4, bonds4, [(a, b, 1.0) for (a, b) in terms])
        assert pred == actual, f"N=4 {label}: pred {pred} vs actual {actual}"

    # Test 7: truly-only returns {0: 4^N}
    pred_truly = chain3.predict_M_spectrum_pi2_odd([('X', 'X'), ('Y', 'Y')], c=1.0)
    assert pred_truly == {0 + 0j: 4 ** 3}, f"truly-only: expected {{0: 64}}, got {pred_truly}"

    # Test 8: mixed truly + Π²-odd (XX + XY) drops the truly contribution
    pred_mixed = chain3.predict_M_spectrum_pi2_odd([('X', 'X'), ('X', 'Y')], c=1.0)
    pred_xy = chain3.predict_M_spectrum_pi2_odd([('X', 'Y')], c=1.0)
    assert pred_mixed == pred_xy, f"truly XX should be dropped: mixed {pred_mixed} vs XY {pred_xy}"

    # Test 9-10: identity raises (single-body falls under F78, not F80)
    with pytest.raises(ValueError, match="single-body"):
        chain3.predict_M_spectrum_pi2_odd([('I', 'Y')], c=1.0)
    with pytest.raises(ValueError, match="single-body"):
        chain3.predict_M_spectrum_pi2_odd([('Z', 'I')], c=1.0)

    # Test 11-12: Π²-even non-truly raises (YZ, ZY out of F80 scope)
    with pytest.raises(ValueError, match="Π²-even"):
        chain3.predict_M_spectrum_pi2_odd([('Y', 'Z')], c=1.0)
    with pytest.raises(ValueError, match="Π²-even"):
        chain3.predict_M_spectrum_pi2_odd([('Z', 'Y')], c=1.0)


def test_F80_kbody_spectrum_identity():
    """F80 (Spec(M) = ±2i·Spec(H_non-truly), mult ×2^N) generalizes to k-body
    chain Π²-odd Hamiltonians.

    PROOF_F80 had bit-exact verification at 2-body chain N=3..7. The proof
    structure (JW + Bogoliubov + per-mode tensor sum) carries verbatim to
    k-body via k-fold Majorana products. This test empirically verifies the
    generalization at k=3 (N=4) and k=4 (N=5) for representative Π²-odd
    Hamiltonians.

    Counterpart to test_F80_bloch_signwalk_chain_pi2_odd (which tests the
    cluster-value sign-walk formula at 2-body via SVD); this test verifies
    the underlying spectral identity at higher k.
    """
    from collections import Counter
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_kbody_chain

    def verify(N, letters):
        H = _build_kbody_chain(N, [letters + (1.0,)])
        L = lindbladian_z_dephasing(H, [0.0] * N)
        M = palindrome_residual(L, 0.0, N)
        M_evs = np.linalg.eigvals(M)
        H_evs = np.linalg.eigvalsh(H)
        # Verify M is purely imaginary (anti-Hermitian)
        max_real = max(abs(ev.real) for ev in M_evs)
        assert max_real < 1e-6, f'k={len(letters)} {letters} N={N}: M has real part {max_real}'
        # F80: Spec(M) = 2i · Spec(H), mult ×2^N
        h_counts = Counter(round(float(e), 8) for e in H_evs)
        m_im_counts = Counter(round(ev.imag, 6) for ev in M_evs)
        bra_factor = 2 ** N
        predicted = Counter()
        for lam, mult in h_counts.items():
            predicted[round(2 * lam, 6)] += mult * bra_factor
        assert predicted == m_im_counts, \
            f'k={len(letters)} {letters} N={N}: F80 fails. ' \
            f'Predicted {dict(predicted)}, measured {dict(m_im_counts)}'

    # k=3 chain Π²-odd at N=4 (lightweight)
    for letters in [('X', 'X', 'Y'), ('Y', 'Y', 'Y'), ('X', 'X', 'Z'), ('X', 'Y', 'X')]:
        verify(4, letters)

    # k=4 chain Π²-odd at N=5
    verify(5, ('X', 'X', 'X', 'Y'))
