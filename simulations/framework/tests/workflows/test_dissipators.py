"""Tests for dissipator c1/c2 closed forms and the HARDWARE_DISSIPATORS table."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_dissipator_c1_c2_pure_paulis():
    """Pure single-Pauli dissipators."""
    # X: c1 = c2 = 16·|α|⁴ = 16
    c1, c2 = fw.dissipator_c1_c2_from_pauli(1, 0, 0)
    assert (c1, c2) == (16.0, 16.0)
    # Y, Z: c1 = 0, c2 = 16
    assert fw.dissipator_c1_c2_from_pauli(0, 1, 0) == (0.0, 16.0)
    assert fw.dissipator_c1_c2_from_pauli(0, 0, 1) == (0.0, 16.0)
    # I-only: both zero
    assert fw.dissipator_c1_c2_from_pauli(0, 0, 0) == (0.0, 0.0)


def test_dissipator_c1_c2_scaling():
    """Operator scaling: c (α=2) gives (c1, c2) scaled by 16 = 2⁴."""
    c1_x, c2_x = fw.dissipator_c1_c2_from_pauli(1, 0, 0)
    c1_2x, c2_2x = fw.dissipator_c1_c2_from_pauli(2, 0, 0)
    assert (c1_2x, c2_2x) == (16 * c1_x, 16 * c2_x)
    c1_half, c2_half = fw.dissipator_c1_c2_from_pauli(0.5, 0, 0)
    assert abs(c1_half - c1_x / 16) < 1e-12
    assert abs(c2_half - c2_x / 16) < 1e-12


def test_dissipator_c1_c2_sigma_minus_plus():
    """σ⁻ = (X-iY)/2 and σ⁺ = (X+iY)/2 give identical (c1, c2) = (3, 4)."""
    c1_sm, c2_sm = fw.dissipator_c1_c2_from_pauli(0.5, -0.5j, 0)
    c1_sp, c2_sp = fw.dissipator_c1_c2_from_pauli(0.5, 0.5j, 0)
    assert abs(c1_sm - 3.0) < 1e-12 and abs(c2_sm - 4.0) < 1e-12
    assert abs(c1_sp - 3.0) < 1e-12 and abs(c2_sp - 4.0) < 1e-12


def test_dissipator_c1_c2_phase_sensitivity():
    """X+Y real and X+iY differ in c1 (32 vs 48), but share c2 = 64."""
    c1_real, c2_real = fw.dissipator_c1_c2_from_pauli(1, 1, 0)
    c1_imag, c2_imag = fw.dissipator_c1_c2_from_pauli(1, 1j, 0)
    assert (c1_real, c2_real) == (32.0, 64.0)
    assert (c1_imag, c2_imag) == (48.0, 64.0)


def test_dissipator_c1_c2_matches_numerical_M():
    """Closed form matches numerical ‖M(L)‖² for H=0, single-class dissipator."""
    from framework import lindbladian_general, palindrome_residual
    sigma_x_2 = np.array([[0,1],[1,0]], dtype=complex)
    sigma_y_2 = np.array([[0,-1j],[1j,0]], dtype=complex)
    sigma_z_2 = np.array([[1,0],[0,-1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    def site_op_general(N, l, op):
        ops = [I2]*N
        ops[l] = op
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out
    N = 3
    test_ops = [
        (1, 0, 0),       # X
        (1, 1j, 0),      # X+iY
        (0.5, -0.5j, 0), # σ⁻
        (1, 1, 1),       # X+Y+Z
        (2, 1, 1),       # 2X+Y+Z
    ]
    for alpha, beta, delta in test_ops:
        c1_pred, c2_pred = fw.dissipator_c1_c2_from_pauli(alpha, beta, delta)
        c_op = alpha*sigma_x_2 + beta*sigma_y_2 + delta*sigma_z_2
        # Two probes: single-site γ=0.1 vs uniform γ=0.1
        for gamma_l, sum_g_sq, sum_g_then_sq in [
            ([0.1, 0, 0],     0.01,  0.01),
            ([0.1, 0.1, 0.1], 0.03,  0.09),
        ]:
            c_ops = [np.sqrt(g) * site_op_general(N, l, c_op)
                     for l, g in enumerate(gamma_l) if g != 0]
            H = np.zeros((2**N, 2**N), dtype=complex)
            L = lindbladian_general(H, c_ops)
            M = palindrome_residual(L, 0.0, N)
            actual = float(np.linalg.norm(M)**2)
            predicted = (4**(N-1)) * (c1_pred * sum_g_sq + c2_pred * sum_g_then_sq)
            assert abs(actual - predicted) < 1e-6, \
                f"({alpha},{beta},{delta}) γ={gamma_l}: actual={actual} pred={predicted}"


def test_hardware_dissipators_table_consistency():
    """HARDWARE_DISSIPATORS table c1/c2 match dissipator_c1_c2_from_pauli."""
    for name, spec in fw.HARDWARE_DISSIPATORS.items():
        a, b, d = spec['pauli']
        c1_pred, c2_pred = fw.dissipator_c1_c2_from_pauli(a, b, d)
        assert c1_pred == spec['c1'], f"{name}: c1 table {spec['c1']} vs computed {c1_pred}"
        assert c2_pred == spec['c2'], f"{name}: c2 table {spec['c2']} vs computed {c2_pred}"
