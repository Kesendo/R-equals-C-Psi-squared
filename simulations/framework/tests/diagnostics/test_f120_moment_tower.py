"""Tests for F120 moment-tower pump channel: the deg-1 girth ladder read linearly."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw
from framework.lindblad import lindbladian_general
from framework.pauli import _build_kbody_chain, site_op


SM = np.array([[0, 1], [0, 0]], dtype=complex)  # σ⁻ = (X+iY)/2 (framework convention)
SP = SM.conj().T                                 # σ⁺


def _site_2x2(N, l, op2):
    ops = [np.eye(2, dtype=complex)] * N
    ops[l] = op2
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def _build_L(H, N, g_deph, g_down, g_up):
    """L = −i[H,·] + Σ γ^deph D[Z_l] + Σ γ↓ D[σ⁻_l] + Σ γ↑ D[σ⁺_l] (dense vec form)."""
    c_ops = []
    for l in range(N):
        if g_deph[l]:
            c_ops.append(np.sqrt(g_deph[l]) * site_op(N, l, 'Z'))
    for l in range(N):
        if g_down[l]:
            c_ops.append(np.sqrt(g_down[l]) * _site_2x2(N, l, SM))
    for l in range(N):
        if g_up[l]:
            c_ops.append(np.sqrt(g_up[l]) * _site_2x2(N, l, SP))
    return lindbladian_general(H, c_ops)


def _slope_dense(A, L, d):
    """d/dt Tr(A ρ)|_{ρ = I/d} from the dense Liouvillian."""
    vecI = np.eye(d, dtype=complex).flatten() / d
    return complex(np.vdot(A.conj().T.flatten(), L @ vecI))


def _girth2_witness(N=3):
    """H = X₀ + X₀Z₁ + 0.7·X₁X₂: t₁ ≡ 0, t₂ fires at site 1 with value 16."""
    return site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')


def test_F120_girth2_witness():
    """moment_tower on the girth-2 witness: t₁ ≡ 0, t₂ = [0, 16, 0], girth = 2,
    verdict 'hard, m* = 5' per the F87 girth dichotomy (m* = 2ℓ+1)."""
    N = 3
    H = _girth2_witness(N)
    tower = fw.moment_tower(H, N, j_max=4)
    assert all(t == 0 for t in tower['t'][1]), f"t_1 = {tower['t'][1]} != 0"
    assert tower['t'][2][0] == 0 and tower['t'][2][2] == 0, "t_2 off-site entries nonzero"
    assert tower['t'][2][1] == 16.0, f"t_2(1) = {tower['t'][2][1]} != 16"
    assert tower['girth'] == 2
    assert tower['deg1_verdict'] == 'hard, m* = 5'

    # girth-1 control: a Z-drive fires at j = 1 with m* = 3
    H1 = 0.5 * site_op(N, 0, 'Z')
    tower1 = fw.moment_tower(H1, N, j_max=2)
    assert tower1['girth'] == 1
    assert tower1['deg1_verdict'] == 'hard, m* = 3'


def test_F120_slope_law_vs_dense_L():
    """predict_pump_slope == the slope read off the dense Liouvillian at N = 3,
    site-dependent rates, all three channel types together, j = 1..4."""
    N, d = 3, 8
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    dg = [a - b for a, b in zip(g_dn, g_up)]
    rng = np.random.default_rng(11)
    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    for H in (_girth2_witness(N), (M + M.conj().T) / 2):
        L = _build_L(H, N, g_deph, g_dn, g_up)
        for j in range(1, 5):
            measured = _slope_dense(np.linalg.matrix_power(H, j), L, d).real
            predicted = fw.predict_pump_slope(H, j, dg)
            assert abs(measured - predicted) <= 1e-14 * max(1.0, abs(measured)), \
                f"j={j}: measured {measured}, predicted {predicted}"

    # scalar Δγ broadcast (uniform) follows the F82 convention
    H = _girth2_witness(N)
    assert fw.predict_pump_slope(H, 2, 0.1) == pytest.approx(
        fw.predict_pump_slope(H, 2, [0.1] * N), abs=1e-15)

    # length validation
    with pytest.raises(ValueError, match="must have length"):
        fw.predict_pump_slope(H, 1, [0.1, 0.2])


def test_F120_f113_bridge_vs_polarity_coordinates():
    """f113_bridge_asymmetry_from_slope == the F113 closed form == the actual Frobenius
    asymmetry from polarity_coordinates_from_hc, for the F113-scope generator
    H = Σ (ω_l/2)·Z_l with per-site σ⁻ (γ_T1) and σ⁺ (γ_pump) channels."""
    for N in (2, 3):
        omega = [0.13, 0.29, 0.07][:N]
        gt1 = [0.10, 0.04, 0.16][:N]
        gpu = [0.03, 0.11, 0.02][:N]
        H = sum((omega[l] / 2.0) * site_op(N, l, 'Z') for l in range(N))
        dg = [a - b for a, b in zip(gt1, gpu)]

        closed = (4 ** N / 2.0) * sum(omega[l] * (gpu[l] - gt1[l]) for l in range(N))
        bridged = fw.f113_bridge_asymmetry_from_slope(H, dg)
        c_ops = [_site_2x2(N, l, SM) for l in range(N)] + \
                [_site_2x2(N, l, SP) for l in range(N)]
        actual = fw.polarity_coordinates_from_hc(H, c_ops, list(gt1) + list(gpu), N)['asymmetry']

        assert abs(bridged - closed) <= 1e-15, f"N={N}: bridge {bridged} != closed {closed}"
        assert abs(actual - closed) <= 1e-12, f"N={N}: Frobenius {actual} != closed {closed}"


def test_F120_honest_negative_control_k4_silent():
    """The k = 4 pair IIXY+ZXZY (N = 5): the deg-1 tower is silent through j = 5 (t_j = 0
    exactly, slope = 0), yet the pair is HARD at m* = 11 with p₁₁ = 86507520·γ⁵ (deg-5,
    pinned in f87_girth_dichotomy). The verdict must carry the honest one-sided line."""
    N = 5
    H = _build_kbody_chain(N, [('I', 'I', 'X', 'Y', 1.0), ('Z', 'X', 'Z', 'Y', 1.0)])
    tower = fw.moment_tower(H, N, j_max=5)
    for j in range(1, 6):
        assert max(abs(t) for t in tower['t'][j]) == 0.0, f"t_{j} != 0"
        assert fw.predict_pump_slope(H, j, [0.1, 0.27, 0.05, 0.18, 0.07]) == 0.0
    assert tower['girth'] is None
    assert tower['deg1_verdict'] == \
        'deg-1 silent through j_max = 5 (NOT a softness certificate)'
