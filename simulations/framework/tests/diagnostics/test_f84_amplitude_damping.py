"""Tests for F84 amplitude-damping (cooling+heating, Pauli-channel cancellation)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_F84_amplitude_damping_thermal_bath():
    """F84: f81_violation = √(Σ(γ_↓_l − γ_↑_l)²)·2^(N-1) for combined cooling
    and heating. Special cases: cooling only (= F82), heating only (= F82
    by symmetry), detailed balance γ_↓ = γ_↑ (violation = 0), arbitrary mix.
    Plus: Pauli-channel dissipators (D[Z], D[X], D[Y]) are Π²-symmetric and
    contribute zero to f81_violation.
    """
    chain = fw.ChainSystem(N=3)
    soft = [('X', 'Y'), ('Y', 'X')]

    # Cooling only (= F82)
    d_cool = chain.pi_decompose_M(soft, gamma_z=0.1, gamma_t1=0.1)
    expected = 0.10 * np.sqrt(3) * (2 ** 2)
    assert abs(d_cool['f81_violation'] - expected) < 1e-9

    # Heating only (symmetric to cooling)
    d_heat = chain.pi_decompose_M(soft, gamma_z=0.1, gamma_pump=0.1)
    assert abs(d_heat['f81_violation'] - expected) < 1e-9

    # Detailed balance: γ_↓ = γ_↑, violation must be zero
    d_balance = chain.pi_decompose_M(soft, gamma_z=0.1, gamma_t1=0.1, gamma_pump=0.1)
    assert d_balance['f81_violation'] < 1e-10, \
        f"Detailed balance: violation should be 0, got {d_balance['f81_violation']}"

    # Net cooling at intermediate (γ_↓=0.10, γ_↑=0.05)
    d_net = chain.pi_decompose_M(soft, gamma_z=0.1, gamma_t1=0.1, gamma_pump=0.05)
    expected_net = 0.05 * np.sqrt(3) * (2 ** 2)
    assert abs(d_net['f81_violation'] - expected_net) < 1e-9

    # Net HEATING at intermediate (γ_↓=0.05, γ_↑=0.10), violation symmetric to net cooling
    d_net_heat = chain.pi_decompose_M(soft, gamma_z=0.1, gamma_t1=0.05, gamma_pump=0.10)
    assert abs(d_net_heat['f81_violation'] - expected_net) < 1e-9, \
        f"Net heating must give same violation as net cooling at same |Δγ|, got {d_net_heat['f81_violation']}"

    # Closed-form forward primitive matches numerical
    pred = chain.predict_amplitude_damping_violation(0.1, 0.05)
    assert abs(pred - expected_net) < 1e-12

    # Inverse closed form: estimate_net_cooling_from_violation
    delta_recovered = chain.estimate_net_cooling_from_violation(d_net['f81_violation'])
    assert abs(delta_recovered - 0.05) < 1e-9, \
        f"net cooling rate inversion: expected 0.05, got {delta_recovered}"

    # Non-uniform cooling and heating
    gt1_l = [0.10, 0.05, 0.15]
    gp_l = [0.0, 0.02, 0.08]
    delta_l = [d - u for d, u in zip(gt1_l, gp_l)]
    expected_nu = np.sqrt(sum(d * d for d in delta_l)) * (2 ** 2)
    d_nu = chain.pi_decompose_M(soft, gamma_z=0.0, gamma_t1=gt1_l, gamma_pump=gp_l)
    assert abs(d_nu['f81_violation'] - expected_nu) < 1e-9
    pred_nu = chain.predict_amplitude_damping_violation(gt1_l, gp_l)
    assert abs(pred_nu - expected_nu) < 1e-12

    # Heating-only is recovered by predict_amplitude_damping_violation when
    # pump is given but t1 is zero
    pred_heat = chain.predict_amplitude_damping_violation([0.0]*3, [0.1]*3)
    assert abs(pred_heat - expected) < 1e-12

    # Backward compatibility: predict_T1_dissipator_violation = predict_amplitude_damping_violation(γ_T1, None)
    pred_t1 = chain.predict_T1_dissipator_violation(0.1)
    pred_amp_only = chain.predict_amplitude_damping_violation(0.1)
    assert abs(pred_t1 - pred_amp_only) < 1e-12


def test_F84_pauli_channels_pi2_symmetric():
    """F84 corollary: pure D[Z], D[X], D[Y] dissipators are Π²-symmetric and
    give zero F81 violation. Verifies the Pauli-Channel Cancellation Lemma
    by constructing each dissipator and checking the violation is at machine
    precision regardless of the rate.

    Note: D[Z] is built into Z-dephasing (already tested via F81); this test
    verifies D[X] and D[Y] explicitly via lindbladian_general.
    """
    from framework.lindblad import lindbladian_general, palindrome_residual
    from framework.pauli import _vec_to_pauli_basis_transform, site_op, _build_bilinear
    from framework.symmetry import build_pi_full

    N = 3
    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)
    T = _vec_to_pauli_basis_transform(N)
    d = 2 ** N

    H_soft = _build_bilinear(N, bonds, [('X', 'Y', 1.0), ('Y', 'X', 1.0)])
    L_H_vec = -1j * (np.kron(H_soft, np.eye(d, dtype=complex)) -
                     np.kron(np.eye(d, dtype=complex), H_soft.T))
    L_H_p = (T.conj().T @ L_H_vec @ T) / d

    for letter in ['X', 'Y']:
        c_ops = [np.sqrt(0.1) * site_op(N, l, letter) for l in range(N)]
        L_full = lindbladian_general(H_soft, c_ops)
        M = palindrome_residual(L_full, 0.0, N)
        PiMPi = Pi @ M @ Pi_inv
        M_anti = (M - PiMPi) / 2
        violation = float(np.linalg.norm(M_anti - L_H_p))
        assert violation < 1e-10, \
            f"D[{letter}] should give zero F81 violation (Pauli-Channel Cancellation), got {violation}"
