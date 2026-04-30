"""Tests for F82 closed-form T1 dissipator violation and the predict/invert primitives."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F82_closed_form_T1_dissipator():
    """F82: f81_violation = ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1).

    Verifies the closed-form scaling of the F81 violation under T1 amplitude
    damping for various N, uniform and non-uniform per-site rates, and across
    multiple Hamiltonians (the violation is H-independent).
    """
    soft = [('X', 'Y'), ('Y', 'X')]

    # N-scaling: ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N−1) for uniform γ_T1
    expected = {
        2: 0.10 * (2 ** 0.5) * (2 ** 1),  # 0.10 · √2 · 2 = 0.2828
        3: 0.10 * (3 ** 0.5) * (2 ** 2),  # 0.10 · √3 · 4 = 0.6928
        4: 0.10 * (4 ** 0.5) * (2 ** 3),  # 0.10 · √4 · 8 = 1.6000
        5: 0.10 * (5 ** 0.5) * (2 ** 4),  # 0.10 · √5 · 16 = 3.5777
    }
    import warnings
    for N, exp_val in expected.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # silence N=2 degeneracy warning
            chain = fw.ChainSystem(N=N)
        d = fw.pi_decompose_M(chain,[('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=0.1)
        assert abs(d['f81_violation'] - exp_val) < 1e-9, \
            f"N={N}: F82 closed form predicts {exp_val:.6f}, got {d['f81_violation']:.6f}"

    # Non-uniform per-site γ_T1 at N=3
    chain3 = fw.ChainSystem(N=3)
    test_cases = [
        ([0.10, 0.0, 0.0], 0.10 * (2 ** 2)),                       # single-site
        ([0.10, 0.10, 0.0], (2 * 0.10 ** 2) ** 0.5 * (2 ** 2)),    # two-site
        ([0.05, 0.10, 0.15], (0.05**2 + 0.10**2 + 0.15**2)**0.5 * (2 ** 2)),  # all different
    ]
    for gt1_l, expected_violation in test_cases:
        d = fw.pi_decompose_M(chain3,[('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
        assert abs(d['f81_violation'] - expected_violation) < 1e-9, \
            f"γ_T1_l={gt1_l}: predicted {expected_violation:.6f}, got {d['f81_violation']:.6f}"

    # H-independence at fixed γ_T1 (uniform 0.1 at N=3): violation = 0.6928 for any H
    expected_violation = 0.10 * (3 ** 0.5) * (2 ** 2)
    for label, terms in [('truly XX+YY', [('X', 'X'), ('Y', 'Y')]),
                         ('soft XY+YX', soft),
                         ('hard XX+XY', [('X', 'X'), ('X', 'Y')]),
                         ('YZ+ZY (Π²-even non-truly)', [('Y', 'Z'), ('Z', 'Y')])]:
        d = fw.pi_decompose_M(chain3,terms, gamma_z=0.1, gamma_t1=0.1)
        assert abs(d['f81_violation'] - expected_violation) < 1e-9, \
            f"{label}: violation should be H-independent, got {d['f81_violation']:.6f}"

    # γ_z-independence at fixed γ_T1: violation = 0.6928 for γ_z ∈ {0, 0.1, 1.0}
    for gz in [0.0, 0.1, 1.0]:
        d = fw.pi_decompose_M(chain3,soft, gamma_z=gz, gamma_t1=0.1)
        assert abs(d['f81_violation'] - expected_violation) < 1e-9, \
            f"γ_z={gz}: violation should be γ_z-independent, got {d['f81_violation']:.6f}"


def test_F82_predict_and_invert_primitives():
    """F82 framework primitives: predict_T1_dissipator_violation and
    estimate_T1_from_violation are forward/inverse closed-form pairs that
    match what pi_decompose_M(gamma_t1=...) returns.

    These primitives wrap the F82 closed form ‖D_{T1, odd}‖_F = √(Σ γ²) ·
    2^(N−1) for direct prediction and the inverse γ_T1, RMS = violation /
    (√N · 2^(N−1)) for hardware-T1 readout.
    """
    # Forward direction: predict matches numerical pi_decompose_M output
    for N in [3, 4]:
        chain = fw.ChainSystem(N=N)
        for gt1 in [0.05, 0.1, 0.5]:
            predicted = fw.predict_T1_dissipator_violation(chain,gt1)
            d = fw.pi_decompose_M(chain,[('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1)
            assert abs(predicted - d['f81_violation']) < 1e-9, \
                f"N={N} γ_T1={gt1}: predict={predicted}, measured={d['f81_violation']}"

        # Non-uniform rates
        gt1_l = [0.05, 0.1, 0.15][:N] + [0.0] * max(0, N - 3)
        predicted_nu = fw.predict_T1_dissipator_violation(chain,gt1_l)
        d_nu = fw.pi_decompose_M(chain,[('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
        assert abs(predicted_nu - d_nu['f81_violation']) < 1e-9, \
            f"N={N} non-uniform: predict={predicted_nu}, measured={d_nu['f81_violation']}"

    # Inverse direction: estimate_T1_from_violation recovers uniform γ_T1
    chain3 = fw.ChainSystem(N=3)
    for gt1_true in [0.0, 0.05, 0.1, 0.5]:
        d = fw.pi_decompose_M(chain3,[('X', 'X'), ('Y', 'Y')], gamma_z=0.1, gamma_t1=gt1_true)
        gt1_rms = fw.estimate_T1_from_violation(chain3,d['f81_violation'])
        assert abs(gt1_rms - gt1_true) < 1e-9, \
            f"γ_T1 inversion: true={gt1_true}, recovered={gt1_rms}"

    # Inverse on non-uniform recovers RMS, not individual rates
    gt1_l = [0.05, 0.1, 0.15]
    rms_expected = (sum(g * g for g in gt1_l) / 3) ** 0.5
    d_nu = fw.pi_decompose_M(chain3,[('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
    rms_recovered = fw.estimate_T1_from_violation(chain3,d_nu['f81_violation'])
    assert abs(rms_recovered - rms_expected) < 1e-9, \
        f"non-uniform RMS: expected={rms_expected}, recovered={rms_recovered}"

    # Edge cases
    assert fw.predict_T1_dissipator_violation(chain3,0.0) == 0.0
    assert fw.estimate_T1_from_violation(chain3,0.0) == 0.0

    # Length validation
    with pytest.raises(ValueError, match="must have length"):
        fw.predict_T1_dissipator_violation(chain3,[0.1, 0.2])  # length 2 ≠ 3

    # Negative violation rejected
    with pytest.raises(ValueError, match="non-negative"):
        fw.estimate_T1_from_violation(chain3,-0.1)
