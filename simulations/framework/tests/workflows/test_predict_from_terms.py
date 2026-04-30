"""Tests for ChainSystem.predict_residual_norm_squared_from_terms and predict_residual_with_hardware_noise."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_predict_residual_with_hardware_noise_t1_only():
    """T1 only: matches the existing T1 closed form."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    T1_l = [0.005] * 3
    result = chain.predict_residual_with_hardware_noise(T1_l=T1_l)
    # T1 alone: 4^(N-1) · [3·Σγ² + 4·(Σγ)²] = 16·[3·3·0.005² + 4·0.015²]
    expected = 16 * (3 * 3 * 0.005**2 + 4 * 0.015**2)
    assert abs(result['per_class']['T1'] - expected) < 1e-12
    assert abs(result['total'] - expected) < 1e-12
    assert result['cross'] == {}


def test_predict_residual_with_hardware_noise_t1_plus_tphi():
    """T1 + Tphi simultaneous: matches numerical ‖M‖² with σ_offset=0."""
    from framework import lindbladian_general, palindrome_residual

    sigma_x_2 = np.array([[0,1],[1,0]], dtype=complex)
    sigma_y_2 = np.array([[0,-1j],[1j,0]], dtype=complex)
    sigma_z_2 = np.array([[1,0],[0,-1]], dtype=complex)
    sigma_minus = (sigma_x_2 - 1j*sigma_y_2)/2

    def site_op_local(N, l, op):
        I2 = np.eye(2, dtype=complex)
        ops = [I2]*N
        ops[l] = op
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    N = 3
    chain = fw.ChainSystem(N=N, gamma_0=0.05)
    gT1, gTphi = 0.005, 0.005

    # Numerical
    c_ops = []
    for l in range(N):
        c_ops.append(np.sqrt(gT1) * site_op_local(N, l, sigma_minus))
        c_ops.append(np.sqrt(gTphi) * site_op_local(N, l, sigma_z_2))
    H = np.zeros((2**N, 2**N), dtype=complex)
    L = lindbladian_general(H, c_ops)
    M = palindrome_residual(L, 0.0, N)
    actual = float(np.linalg.norm(M)**2)

    # Predicted via cockpit method
    result = chain.predict_residual_with_hardware_noise(
        T1_l=[gT1]*N, Tphi_l=[gTphi]*N)
    assert abs(result['total'] - actual) < 1e-9, \
        f"actual={actual} predicted={result['total']}"
    # Check the cross-term is non-zero (cross_T1_Tphi has d1=0, d2=16, only d2 contributes)
    assert ('T1', 'Tphi') in result['cross']


def test_predict_residual_with_hardware_noise_full_stack():
    """Full hardware stack: T1 + Tphi + Xnoise; closed form matches numerical."""
    from framework import lindbladian_general, palindrome_residual

    sigma_x_2 = np.array([[0,1],[1,0]], dtype=complex)
    sigma_y_2 = np.array([[0,-1j],[1j,0]], dtype=complex)
    sigma_z_2 = np.array([[1,0],[0,-1]], dtype=complex)
    sigma_minus = (sigma_x_2 - 1j*sigma_y_2)/2
    I2 = np.eye(2, dtype=complex)

    def site_op_local(N, l, op):
        ops = [I2]*N
        ops[l] = op
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    N = 3
    chain = fw.ChainSystem(N=N, gamma_0=0.05)
    rates = {'T1': [0.001, 0.002, 0.003],
             'Tphi': [0.0005, 0.001, 0.0008],
             'Xnoise': [0.0001, 0.0002, 0.0001]}

    c_ops = []
    op_map = {'T1': sigma_minus, 'Tphi': sigma_z_2, 'Xnoise': sigma_x_2}
    for cls, gl in rates.items():
        op = op_map[cls]
        for l, g in enumerate(gl):
            c_ops.append(np.sqrt(g) * site_op_local(N, l, op))
    H = np.zeros((2**N, 2**N), dtype=complex)
    L = lindbladian_general(H, c_ops)
    M = palindrome_residual(L, 0.0, N)
    actual = float(np.linalg.norm(M)**2)

    result = chain.predict_residual_with_hardware_noise(
        T1_l=rates['T1'], Tphi_l=rates['Tphi'], Xnoise_l=rates['Xnoise'])
    assert abs(result['total'] - actual) < 1e-9, \
        f"actual={actual}, predicted={result['total']}"


def test_predict_from_terms_matches_numerical_chain():
    """Frobenius prediction matches numerical ||M||^2 across topologies and N."""
    cases = [
        (3, 'chain', [('I','Y')], 512.0),         # 2^5·1·16  (||H||²_F=2B·d=2·2·8=16)
        (4, 'chain', [('I','Y')], 3072.0),        # 2^6·1·48
        (4, 'chain', [('Y','Z')], 6144.0),        # 2^6·2·48
        (4, 'chain', [('Y','Z'),('Z','Y')], 12288.0),
        (4, 'chain', [('I','Y'),('Y','I')], 10240.0),
        (4, 'ring',  [('Y','Z'),('Z','Y')], 16384.0),
        (4, 'star',  [('Y','Z'),('Z','Y')], 12288.0),
        (4, 'complete', [('I','Y'),('Y','I')], 36864.0),
        (5, 'chain', [('Y','Z'),('Z','Y')], 65536.0),
    ]
    for N, topo, terms, expected in cases:
        chain = fw.ChainSystem(N=N, topology=topo)
        pred = chain.predict_residual_norm_squared_from_terms(terms)
        assert abs(pred - expected) < 1e-6, \
            f"N={N} topo={topo} terms={terms}: predicted {pred}, expected {expected}"


def test_predict_from_terms_truly_returns_zero():
    chain = fw.ChainSystem(N=4)
    assert chain.predict_residual_norm_squared_from_terms(
        [('X','X'),('Y','Y'),('Z','Z')]) == 0.0
    assert chain.predict_residual_norm_squared_from_terms([('I','X'),('X','I')]) == 0.0


def test_predict_from_terms_handles_mixed_n_yz_classes():
    """Per-term Frobenius sum handles mixed n_YZ classes automatically (no split needed)."""
    chain = fw.ChainSystem(N=4)
    # Mixed n_YZ=(2,2) for YZ+ZY plus n_YZ=(1,1) for IY+YI.
    # Each term contributes 2^(N+2)·n_YZ_k·||H_k||²_F separately.
    pred = chain.predict_residual_norm_squared_from_terms(
        [('Y','Z'),('Z','Y'),('I','Y'),('Y','I')])
    actual = chain.residual_norm_squared(
        [('Y','Z'),('Z','Y'),('I','Y'),('Y','I')])
    assert abs(pred - actual) < 1e-6


def test_predict_from_terms_decomposes_additively():
    """For mixed-class H, sum of per-class predictions equals total numerical ||M||^2."""
    chain = fw.ChainSystem(N=4)
    p_2yz = chain.predict_residual_norm_squared_from_terms([('Y','Z'),('Z','Y')])
    p_1yz = chain.predict_residual_norm_squared_from_terms([('I','Y'),('Y','I')])
    actual = chain.residual_norm_squared(
        [('Y','Z'),('Z','Y'),('I','Y'),('Y','I')])
    assert abs(p_2yz + p_1yz - actual) < 1e-6


def test_predict_from_terms_t1_truly_hamiltonian():
    """For truly Hamiltonians, the T1 contribution is the *only* term."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    pred = chain.predict_residual_norm_squared_from_terms(
        [('X','X'),('Y','Y'),('Z','Z')], gamma_t1=0.005)
    actual = chain.residual_norm_squared(
        [('X','X'),('Y','Y'),('Z','Z')], gamma_t1=0.005)
    # 4^(3-1) * [3·3·0.005² + 4·(3·0.005)²] = 16·(0.000225 + 0.0009) = 16·0.001125 = 0.018
    assert abs(pred - 0.018) < 1e-9
    assert abs(pred - actual) < 1e-9


def test_predict_from_terms_t1_soft_hamiltonian_additive():
    """T1 additively extends the Frobenius result for soft Hamiltonians."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    z_only = chain.predict_residual_norm_squared_from_terms(
        [('Y','Z'),('Z','Y')])  # 2048
    with_t1 = chain.predict_residual_norm_squared_from_terms(
        [('Y','Z'),('Z','Y')], gamma_t1=0.005)
    # T1 part = 16·[3·3·0.005² + 4·(0.015)²] = 16·0.001125 = 0.018
    assert abs(with_t1 - z_only - 0.018) < 1e-9
    actual = chain.residual_norm_squared(
        [('Y','Z'),('Z','Y')], gamma_t1=0.005)
    assert abs(with_t1 - actual) < 1e-9


def test_predict_from_terms_t1_nonuniform_distribution():
    """T1 formula handles arbitrary per-site distributions."""
    chain = fw.ChainSystem(N=4, gamma_0=0.05)
    gamma_t1 = [0.001, 0.005, 0.01, 0.002]
    pred = chain.predict_residual_norm_squared_from_terms(
        [('I','Y'),('Y','I')], gamma_t1=gamma_t1)
    actual = chain.residual_norm_squared(
        [('I','Y'),('Y','I')], gamma_t1=gamma_t1)
    assert abs(pred - actual) < 1e-9


def test_predict_from_terms_t1_no_hamiltonian():
    """Pure-T1 (no Hamiltonian) gives the analytic dissipator formula."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    # H = 0 → only T1 part survives.
    pred_empty = chain.predict_residual_norm_squared_from_terms(
        [], gamma_t1=0.005)
    # 4^(3-1) · [3·3·0.005² + 4·(0.015)²] = 16 · 0.001125 = 0.018
    assert abs(pred_empty - 0.018) < 1e-9


def test_predict_from_terms_v_effect_mixed_truly_nontruly():
    """V-Effect 36-combos exposed: per-term truly handling matters.

    YY+YZ has n_YZ=(2,2) homogeneously, but YY is truly (M=0) so only YZ
    contributes. Old gross-list logic predicted 2048; correct is 1024.
    """
    chain = fw.ChainSystem(N=3, gamma_0=0.1)
    # YY+YZ: YY truly (M=0), YZ contributes 32·2·16 = 1024
    val = chain.predict_residual_norm_squared_from_terms(
        [('Y','Y'),('Y','Z')])
    actual = chain.residual_norm_squared([('Y','Y'),('Y','Z')])
    assert abs(val - actual) < 1e-6
    assert abs(val - 1024.0) < 1e-6
    # XX+YZ: XX truly, YZ contributes
    val = chain.predict_residual_norm_squared_from_terms(
        [('X','X'),('Y','Z')])
    actual = chain.residual_norm_squared([('X','X'),('Y','Z')])
    assert abs(val - actual) < 1e-6


def test_predict_from_terms_v_effect_full_36_combos():
    """All 36 V-Effect combos at N=3: predict matches numerical exactly."""
    from itertools import combinations
    SINGLE = ['XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ']
    chain = fw.ChainSystem(N=3, gamma_0=0.1)
    for t1, t2 in combinations(SINGLE, 2):
        terms = [(t1[0], t1[1]), (t2[0], t2[1])]
        pred = chain.predict_residual_norm_squared_from_terms(terms)
        num = chain.residual_norm_squared(terms)
        assert abs(pred - num) < 1e-6, \
            f"{t1}+{t2}: pred={pred}, num={num}"
