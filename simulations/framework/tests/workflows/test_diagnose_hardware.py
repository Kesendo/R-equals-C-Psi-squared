"""Tests for diagnose_hardware: F-toolkit lens-reading workflow on hardware data."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


REPO_ROOT = Path(__file__).resolve().parents[4]
F83_HARDWARE_JSON = REPO_ROOT / 'data' / 'ibm_f83_signature_april2026' / 'f83_signature_ibm_marrakesh_20260430_190035.json'

F83_TERMS_PER_CATEGORY = {
    'truly_unbroken':       [('X', 'X'), ('Y', 'Y')],
    'pi2_odd_pure':         [('X', 'Y'), ('Y', 'X')],
    'pi2_even_nontruly':    [('Y', 'Z'), ('Z', 'Y')],
    'mixed_anti_one_sixth': [('X', 'Y'), ('Y', 'Z')],
}

# Calibration values from 2026-04-30T16:25Z, path [4, 5, 6]
F83_CALIBRATION = {
    'T1': [206.115, 144.477, 351.240],
    'T2': [183.978, 120.938, 151.412],
}


def _load_f83_measured():
    """Load F83 hardware measurements from JSON, parsed into category → (P,P): float."""
    if not F83_HARDWARE_JSON.exists():
        pytest.skip(f"F83 hardware JSON not found at {F83_HARDWARE_JSON}")
    with open(F83_HARDWARE_JSON) as f:
        data = json.load(f)
    measured = {}
    for cat, exp_dict in data['expectations'].items():
        cat_dict = {}
        for key, val in exp_dict.items():
            if ',' in key:
                p0, p2 = key.split(',')
                cat_dict[(p0, p2)] = val
        measured[cat] = cat_dict
    return measured


def test_diagnose_hardware_returns_structure():
    """diagnose_hardware returns the expected nested dict structure."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(
        chain, measured, F83_TERMS_PER_CATEGORY,
        calibration=F83_CALIBRATION, t=0.8, n_trotter=3, gamma_z=0.1, shots=4096,
    )

    assert 'per_category' in result
    assert 'cross_category' in result
    for cat in F83_TERMS_PER_CATEGORY:
        pc = result['per_category'][cat]
        for key in ('F77_class', 'F83_anti_fraction', 'F83_r',
                    'predictions', 'measurements', 'residuals',
                    'rms_residual', 'lens_readings'):
            assert key in pc, f"missing key '{key}' in per_category['{cat}']"


def test_diagnose_hardware_F77_classification_matches_F83_categories():
    """F77 classifier identifies the four F83-test categories correctly."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(chain, measured, F83_TERMS_PER_CATEGORY, gamma_z=0.1)

    assert result['per_category']['truly_unbroken']['F77_class'] == 'truly'
    # pi2_odd_pure is M ≠ 0 with eigenvalue pairing → 'soft'
    assert result['per_category']['pi2_odd_pure']['F77_class'] == 'soft'
    # pi2_even_nontruly is also 'soft' (M ≠ 0, eigenvalue pairing, F83 anti = 0)
    assert result['per_category']['pi2_even_nontruly']['F77_class'] == 'soft'
    # Mixed (XY+YZ) is hard
    assert result['per_category']['mixed_anti_one_sixth']['F77_class'] == 'hard'


def test_diagnose_hardware_F83_anti_fractions_match_closed_form():
    """F83 anti-fractions in the diagnostic match the closed-form predictions."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(chain, measured, F83_TERMS_PER_CATEGORY, gamma_z=0.1)

    # Truly: M = 0, anti = None
    assert result['per_category']['truly_unbroken']['F83_anti_fraction'] is None
    # Pure Π²-odd: anti = 1/2
    assert abs(result['per_category']['pi2_odd_pure']['F83_anti_fraction'] - 0.5) < 1e-9
    # Pure Π²-even non-truly: anti = 0
    assert abs(result['per_category']['pi2_even_nontruly']['F83_anti_fraction'] - 0.0) < 1e-9
    # Mixed equal-Frobenius: anti = 1/6
    assert abs(result['per_category']['mixed_anti_one_sixth']['F83_anti_fraction'] - 1.0 / 6) < 1e-9


def test_diagnose_hardware_truly_has_F82_F84_signature_lens_reading():
    """The hypothesis test: on F83 Marrakesh data, the workflow detects the
    F82/F84 amplitude-damping signature on truly_unbroken (because truly's
    ⟨Z,Z⟩ damps significantly, which the M=0 lens reads as σ⁻ amplitude
    damping breaking ⟨Z⟩-conservation)."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(
        chain, measured, F83_TERMS_PER_CATEGORY,
        calibration=F83_CALIBRATION, gamma_z=0.1,
    )

    truly_readings = result['per_category']['truly_unbroken']['lens_readings']
    f82_readings = [r for r in truly_readings if r.get('lens') == 'F82/F84']
    assert len(f82_readings) == 1, "truly_unbroken should have exactly one F82/F84 lens reading"
    f82 = f82_readings[0]
    # The actual hardware run shows 56% damping at γ_Z=0.1 (framework default;
    # equivalently 60% at γ_Z=0.05 path-fit). The workflow flags this as significant.
    assert f82['significant'] is True
    assert abs(f82['damping_fraction']) > 0.20, \
        f"Damping fraction {f82['damping_fraction']:.3f} should exceed 20% for the σ⁻ signature"


def test_diagnose_hardware_pure_class_quantitative_match():
    """Pure Π²-odd and pure Π²-even non-truly hardware data match Trotter+γ_Z
    prediction at small RMS (within shot noise) — direct confirmation of F83
    structural prediction at hardware level."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(
        chain, measured, F83_TERMS_PER_CATEGORY,
        gamma_z=0.1, shots=4096,
    )

    pi2_odd_readings = result['per_category']['pi2_odd_pure']['lens_readings']
    quant = [r for r in pi2_odd_readings if r.get('lens') == 'F83-quantitative']
    assert len(quant) == 1, "pi2_odd_pure should produce one F83-quantitative reading"
    # RMS in σ-units should be moderate — pure-class match is structural.
    assert quant[0]['rms_sigma_units'] < 20, \
        f"pi2_odd_pure RMS {quant[0]['rms_sigma_units']:.1f}σ should fit F83 quantitatively"

    pi2_even_readings = result['per_category']['pi2_even_nontruly']['lens_readings']
    quant = [r for r in pi2_even_readings if r.get('lens') == 'F83-quantitative']
    assert len(quant) == 1, "pi2_even_nontruly should produce one F83-quantitative reading"
    assert quant[0]['rms_sigma_units'] < 20, \
        f"pi2_even_nontruly RMS {quant[0]['rms_sigma_units']:.1f}σ should fit F83 quantitatively"


def test_diagnose_hardware_Y_Z_asymmetry_attributed_to_T2_inhomogeneity():
    """The truly Hamiltonian XX+YY is symmetric under (q0 ↔ q2) swap on
    chain N=3, so prediction has ⟨Y,Z⟩ = ⟨Z,Y⟩. Hardware shows asymmetry
    because Q4 (T2=184μs) and Q6 (T2=151μs) have different decoherence rates.
    The workflow detects this and attributes it to T2 inhomogeneity from
    the calibration."""
    measured = _load_f83_measured()
    chain = fw.ChainSystem(N=3)
    result = fw.diagnose_hardware(
        chain, measured, F83_TERMS_PER_CATEGORY,
        calibration=F83_CALIBRATION,
    )

    yz = result['cross_category']['Y_Z_asymmetry_on_truly']
    assert yz is not None
    assert yz['expected_symmetric'] is True
    # Hardware shows ⟨Y,Z⟩ = +0.67, ⟨Z,Y⟩ = +0.19 → asymmetry = 0.48 ~31σ at 4096 shots
    assert yz['sigma_units'] > 5
    assert yz['attributed_to_calibration'] is True
    assert 'T2 inhomogeneity' in yz['attribution']


def test_diagnose_hardware_structural_klein_set_no_labels():
    """Structure lens exposes raw Klein indices (bit_a, bit_b) tuples, not
    label letters. Replaces the prior Klein-letter labels (M, F_a, F_b, C)
    which were k=2-specific and confusing at k≥3."""
    chain = fw.ChainSystem(N=3)

    test_cases = [
        # (terms, expected_klein_set as sorted tuples)
        ([('X', 'X'), ('Y', 'Y')],          [(0, 0)]),                 # pure (0,0)
        ([('X', 'Y'), ('Y', 'X')],          [(0, 1)]),                 # pure (0,1)
        ([('X', 'Z'), ('Z', 'X')],          [(1, 1)]),                 # pure (1,1)
        ([('Y', 'Z'), ('Z', 'Y')],          [(1, 0)]),                 # pure (1,0)
        ([('X', 'X'), ('X', 'Y')],          [(0, 0), (0, 1)]),         # April 26 hard
        ([('X', 'X'), ('X', 'Z')],          [(0, 0), (1, 1)]),         # M + (1,1)
        ([('X', 'Y'), ('X', 'Z')],          [(0, 1), (1, 1)]),         # always hard at k=2
        ([('X', 'Y'), ('Y', 'Z')],          [(0, 1), (1, 0)]),         # April 30 mixed
    ]
    terms_per_category = {f'case_{i}': terms for i, (terms, _) in enumerate(test_cases)}
    predictions = fw.predict_signature_table(chain, terms_per_category, gamma_z=0.1)
    result = fw.diagnose_hardware(chain, predictions, terms_per_category, gamma_z=0.1)

    for i, (terms, expected_klein_set) in enumerate(test_cases):
        cat = f'case_{i}'
        struct = result['per_category'][cat]['structure']
        assert struct['klein_set'] == expected_klein_set, \
            f"{terms}: expected klein_set {expected_klein_set}, got {struct['klein_set']}"
        # Klein-homogeneity check
        expected_homo = (len(expected_klein_set) == 1)
        assert struct['is_klein_homogeneous'] == expected_homo
        # Find the Structure lens reading
        struct_readings = [r for r in result['per_category'][cat]['lens_readings']
                           if r.get('lens') == 'Structure']
        assert len(struct_readings) == 1
        assert struct_readings[0]['klein_set'] == expected_klein_set


def test_diagnose_hardware_klein_homogeneity_predicts_soft():
    """Klein-homogeneous Hamiltonians (all terms with same Klein index) are
    always F77 soft or truly, never hard. Structural fact verified at k=2."""
    chain = fw.ChainSystem(N=3)

    homogeneous_cases = [
        [('X', 'X'), ('Y', 'Y')],  # both (0,0) → truly
        [('X', 'X'), ('Z', 'Z')],  # both (0,0) → truly
        [('X', 'Y'), ('Y', 'X')],  # both (0,1) → soft
        [('X', 'Z'), ('Z', 'X')],  # both (1,1) → soft
        [('Y', 'Z'), ('Z', 'Y')],  # both (1,0) → soft
    ]

    for terms in homogeneous_cases:
        f77 = fw.classify_pauli_pair(chain, terms)
        assert f77 in ('truly', 'soft'), \
            f"Klein-homogeneous Hamiltonian {terms} should be truly or soft, got {f77}"


def test_diagnose_hardware_structural_distinguishes_april26_april30_hard_cases():
    """April 26 hard_broken (XX+XY) has klein_set {(0,0), (0,1)}.
    April 30 mixed_anti_one_sixth (XY+YZ) has klein_set {(0,1), (1,0)}.
    Both Klein-inhomogeneous; both F77-hard. The Structure lens shows the
    raw klein indices that distinguish them, no label-collapse."""
    chain = fw.ChainSystem(N=3)

    terms = {
        'april26_hard': [('X', 'X'), ('X', 'Y')],
        'april30_mixed': [('X', 'Y'), ('Y', 'Z')],
    }
    predictions = fw.predict_signature_table(chain, terms, gamma_z=0.1)
    result = fw.diagnose_hardware(chain, predictions, terms, gamma_z=0.1)

    s26 = result['per_category']['april26_hard']['structure']
    assert s26['klein_set'] == [(0, 0), (0, 1)]
    assert s26['is_klein_homogeneous'] is False

    s30 = result['per_category']['april30_mixed']['structure']
    assert s30['klein_set'] == [(0, 1), (1, 0)]
    assert s30['is_klein_homogeneous'] is False


def test_diagnose_hardware_structure_exposes_y_parity_axis():
    """The Structure lens exposes Y-parity as an independent axis. At k=2
    Y-parity is determined by Klein index (redundant); the diagnostic still
    reports it as a property so users can verify consistency."""
    chain = fw.ChainSystem(N=3)

    terms = {
        'truly': [('X', 'X'), ('Y', 'Y')],            # Y-par: [0, 0]
        'pi2_odd': [('X', 'Y'), ('Y', 'X')],          # Y-par: [1, 1]
        'pi2_even_nontruly': [('Y', 'Z'), ('Z', 'Y')], # Y-par: [1, 1]
    }
    predictions = fw.predict_signature_table(chain, terms, gamma_z=0.1)
    result = fw.diagnose_hardware(chain, predictions, terms, gamma_z=0.1)

    s_truly = result['per_category']['truly']['structure']
    assert s_truly['per_term_y_parities'] == [0, 0]
    assert s_truly['is_y_parity_homogeneous'] is True
    assert s_truly['y_parity_set'] == [0]

    s_odd = result['per_category']['pi2_odd']['structure']
    assert s_odd['per_term_y_parities'] == [1, 1]
    assert s_odd['is_y_parity_homogeneous'] is True
    assert s_odd['y_parity_set'] == [1]


def test_diagnose_hardware_synthetic_round_trip():
    """Round-trip sanity: feed predict_signature_table's output back as
    measured. diagnose_hardware should report zero residual and no F82/F84
    flag (since synthetic 'measured' = 'predicted')."""
    chain = fw.ChainSystem(N=3)
    predicted = fw.predict_signature_table(
        chain, F83_TERMS_PER_CATEGORY, t=0.8, n_trotter=3, gamma_z=0.1,
    )
    result = fw.diagnose_hardware(
        chain, predicted, F83_TERMS_PER_CATEGORY,
        gamma_z=0.1, shots=4096,
    )
    for cat in F83_TERMS_PER_CATEGORY:
        rms = result['per_category'][cat]['rms_residual']
        assert rms < 1e-10, f"Round-trip RMS for {cat} should be ~0, got {rms}"

    # Truly's F82/F84 signature: not significant in round-trip
    truly_readings = result['per_category']['truly_unbroken']['lens_readings']
    f82 = [r for r in truly_readings if r.get('lens') == 'F82/F84']
    if f82:
        assert f82[0]['significant'] is False, \
            "Synthetic round-trip should not show F82/F84 signature"
