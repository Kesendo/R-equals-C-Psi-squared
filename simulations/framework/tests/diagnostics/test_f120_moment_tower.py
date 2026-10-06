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


def _rung2_witness(N=3):
    """H = X₀ + X₀Z₁ + 0.7·X₁X₂: t₁ ≡ 0, t₂ fires at site 1 with value 16."""
    return site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')


def test_F120_rung2_witness():
    """moment_tower on the rung-2 witness: t₁ ≡ 0, t₂ = [0, 16, 0], first nonzero rung 2,
    verdict 'hard, m* <= 5' (tight, as at every k = 2: PROOF_MOMENT_TOWER_PUMP_CHANNEL §4)."""
    N = 3
    H = _rung2_witness(N)
    tower = fw.moment_tower(H, N, j_max=4)
    assert all(t == 0 for t in tower['t'][1]), f"t_1 = {tower['t'][1]} != 0"
    assert tower['t'][2][0] == 0 and tower['t'][2][2] == 0, "t_2 off-site entries nonzero"
    assert tower['t'][2][1] == 16.0, f"t_2(1) = {tower['t'][2][1]} != 16"
    assert tower['first_firing_rung'] == 2
    assert tower['deg1_verdict'] == 'hard, m* <= 5'

    # girth-1 control: a Z-drive fires at j = 1 with m* = 3
    H1 = 0.5 * site_op(N, 0, 'Z')
    tower1 = fw.moment_tower(H1, N, j_max=2)
    assert tower1['first_firing_rung'] == 1
    assert tower1['deg1_verdict'] == 'hard, m* <= 3'


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
    for H in (_rung2_witness(N), (M + M.conj().T) / 2):
        # Pin the operands. Every comparison below is slope-vs-law, and both read 0.0 for
        # a zero H or at detailed balance, so the test would otherwise certify nothing.
        assert np.linalg.norm(H) > 1e-6, "the witness H is zero"
        assert any(abs(a - b) > 1e-6 for a, b in zip(g_dn, g_up)), (
            "the rates are at detailed balance, so every comparison below reads 0 == 0")
        L = _build_L(H, N, g_deph, g_dn, g_up)
        for j in range(1, 5):
            measured = _slope_dense(np.linalg.matrix_power(H, j), L, d).real
            predicted = fw.predict_pump_slope(H, j, dg)
            assert abs(measured - predicted) <= 1e-14 * max(1.0, abs(measured)), \
                f"j={j}: measured {measured}, predicted {predicted}"

    # scalar Δγ broadcast (uniform) follows the F82 convention
    H = _rung2_witness(N)
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
        # a three-route CONSISTENCY check agrees at zero too, in two ways: pin both
        assert any(abs(w) > 1e-6 for w in omega), "the drive is zero"
        assert any(abs(x - y) > 1e-6 for x, y in zip(gt1, gpu)), (
            "detailed balance: all three routes read 0.0 and the check is vacuous")
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
    # A silence certificate has to show the operator is there before it shows the tower
    # is not: everything below asserts a vanishing, so the zero Hamiltonian passed it.
    assert np.trace(H @ H).real == 128.0, f"not the witness ({np.trace(H @ H).real})"
    tower = fw.moment_tower(H, N, j_max=5)
    for j in range(1, 6):
        assert max(abs(t) for t in tower['t'][j]) == 0.0, f"t_{j} != 0"
        assert fw.predict_pump_slope(H, j, [0.1, 0.27, 0.05, 0.18, 0.07]) == 0.0
    assert tower['first_firing_rung'] is None
    assert tower['deg1_verdict'] == \
        'deg-1 silent through j_max = 5 (NOT a softness certificate)'


def test_F120_first_firing_rung_is_not_the_f87_girth():
    """The field this diagnostic returns used to be called 'girth'. It is not one.

    f87's ℓ (effective_ell) is 1 for a nonzero diagonal and otherwise the shortest odd
    cycle, so it is never even; the first firing rung routinely is. H = Z₀Z₁ + Z₀Z₁Z₂ is
    the cleanest separation: it is diagonal, so ℓ = 1, but H = z₀z₁(1 + z₂) squares to
    2I + 2Z₂, giving t₁ ≡ 0 and t₂ = (0, 0, 16). The certificate is m* = 2·2+1 = 5, not
    the 2ℓ+1 = 3 the girth reading would report.
    """
    N = 3
    Z = [site_op(N, l, 'Z') for l in range(N)]
    H = Z[0] @ Z[1] + Z[0] @ Z[1] @ Z[2]

    tower = fw.moment_tower(H, N, j_max=3)
    assert max(abs(t) for t in tower['t'][1]) == 0.0, "the girth rung must be silent here"
    assert [complex(t).real for t in tower['t'][2]] == [0.0, 0.0, 16.0]
    assert tower['first_firing_rung'] == 2
    assert tower['deg1_verdict'] == 'hard, m* <= 5'

    # H is diagonal with nonzero entries, which is exactly f87's ℓ = 1 branch, so a
    # reading that identified the two would have emitted m* = 3 instead.
    assert np.count_nonzero(np.diag(H)) > 0
    assert np.allclose(H, np.diag(np.diag(H))), "H must be diagonal for the ℓ = 1 branch"


def test_F120_the_rung_bound_is_not_an_equality():
    """2k+1 BOUNDS m*, it does not give it, and the docstring must not promise equality.

    H = Y₂ + Z₀Z₁Y₂ + Z₀Z₁Z₂ at N = 3 has t₁ = t₂ = 0 and t₃ = (0, 0, 16), so the first
    nonzero rung is k = 3 and the bound reads 2k+1 = 7. But H carries the weight-3 Z-string
    Z₀Z₁Z₂, so the deg-3 class already fires at 5: p₅ = 120·4³·γ³ = 7680·γ³
    (PROOF_MOMENT_TOWER_PUMP_CHANNEL §4). The true hardness moment is 5.

    p₇ is not a monomial and changes sign, so the moment the bound names is not the one
    that certifies hardness throughout.
    """
    N = 3
    Y2 = site_op(N, 2, 'Y')
    Z0, Z1, Z2 = (site_op(N, l, 'Z') for l in range(3))
    H = Y2 + Z0 @ Z1 @ Y2 + Z0 @ Z1 @ Z2
    assert np.allclose(H, H.conj().T)

    tower = fw.moment_tower(H, N, j_max=3)
    assert max(abs(t) for t in tower['t'][1]) == 0.0
    assert max(abs(t) for t in tower['t'][2]) == 0.0
    assert [complex(t).real for t in tower['t'][3]] == [0.0, 0.0, 16.0]
    assert tower['first_firing_rung'] == 3
    assert tower['deg1_verdict'] == 'hard, m* <= 7'

    # p_5 fires below the bound, so m* = 5 < 7. M = A + γQ, the f87 recentred generator.
    d = 2 ** N
    Q = sum(np.kron(site_op(N, l, 'Z'), site_op(N, l, 'Z').conj()) for l in range(N))
    A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))

    def p(m, g):
        return float(np.real(np.trace(np.linalg.matrix_power(A + g * Q, m))))

    for g in (0.1, 0.25, 0.4):
        assert abs(p(1, g)) < 1e-7 and abs(p(3, g)) < 1e-7, f"p1/p3 should be silent at {g}"
        assert p(5, g) == pytest.approx(7680.0 * g ** 3, rel=1e-9)   # fires, monomial γ³
    # and the bound's own moment is not sign-definite
    assert p(7, 0.25) > 0 and p(7, 0.4) < 0


def test_F120_p53_closed_form_and_k3_tightness():
    """PROOF_MOMENT_TOWER_PUMP_CHANNEL §4: P_{5,3} = 20·4^N·[6·Σ_{|S|=3} h_S² + (3N−2)·Σ_l c_l²] for every
    Hermitian H, read exactly (integer Tr(M^5) at γ = 0..5, interpolated); and at k = 3 the bound is slack
    exactly when H carries a weight-3 Z-string (the Y₂ row carries Z₀Z₁Z₂, p₅ = 7680·γ³; IXX + IZZ + ZYY,
    k = 3 with none, has p₅ ≡ 0, so m* = 7 = 2k+1)."""
    from fractions import Fraction
    from itertools import product

    def string(s):
        out = np.array([[1]], dtype=complex)
        for c in s:
            out = np.kron(out, {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}[c])
        return out

    def p5_coeffs(H, N):
        d = 2 ** N
        A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
        Q = sum(np.kron(string('I' * l + 'Z' + 'I' * (N - l - 1)), string('I' * l + 'Z' + 'I' * (N - l - 1)))
                for l in range(N)).real
        Ar, Ai, Qi = (np.round(A.real).astype(object), np.round(A.imag).astype(object), np.round(Q).astype(object))
        vals = []
        for g in range(6):
            R, I = Ar + g * Qi, Ai
            PR, PI = R.copy(), I.copy()
            for _ in range(4):
                PR, PI = PR.dot(R) - PI.dot(I), PR.dot(I) + PI.dot(R)
            assert sum(PI[i, i] for i in range(d * d)) == 0
            vals.append(Fraction(int(sum(PR[i, i] for i in range(d * d)))))
        coef = vals[:]                                # Newton divided differences on γ = 0..5
        for j in range(1, 6):
            for i in range(5, j - 1, -1):
                coef[i] = (coef[i] - coef[i - 1]) / j
        poly, basis = [Fraction(0)] * 6, [Fraction(1)]
        for i in range(6):
            for k_, b in enumerate(basis):
                poly[k_] += coef[i] * b
            nb = [Fraction(0)] * (len(basis) + 1)
            for k_, b in enumerate(basis):
                nb[k_ + 1] += b
                nb[k_] -= i * b
            basis = nb
        return [int(c) for c in poly]

    N = 3
    rng = np.random.default_rng(20261006)
    labels = [''.join(t) for t in product('IXYZ', repeat=N) if set(t) != {'I'}]
    for _ in range(3):
        coef = {s: int(rng.integers(-3, 4)) for s in rng.choice(labels, 6, replace=False)}
        coef.update({'ZZZ': 2, 'IZI': -1})
        H = sum(c * string(s) for s, c in coef.items())
        w3 = sum(c * c for s, c in coef.items() if set(s) <= {'I', 'Z'} and s.count('Z') == 3)
        w1 = sum(c * c for s, c in coef.items() if set(s) <= {'I', 'Z'} and s.count('Z') == 1)
        assert p5_coeffs(H, N)[3] == 20 * 4 ** N * (6 * w3 + (3 * N - 2) * w1)
    with_zzz = string('IIY') + string('ZZY') + string('ZZZ')
    without = string('IXX') + string('IZZ') + string('ZYY')
    assert p5_coeffs(with_zzz, N) == [0, 0, 0, 7680, 0, 0]
    assert p5_coeffs(without, N) == [0] * 6
    tower = fw.moment_tower(without, N, j_max=3)
    assert tower['first_firing_rung'] == 3


def test_F120_tightness_examples_exact():
    """The rows PROOF_MOMENT_TOWER_PUMP_CHANNEL §4 and F120 (c) quote, read exactly at N = 3: (k, m*, p_m*) for a
    slack and a tight k = 4 row, a slack k = 5 row under F H F = -H with no weight-3 Z-string, and the girth
    dichotomy failing outside F H F = -H (ell = 3 with t_3 != 0, yet k = 2 and m* = 5)."""
    from fractions import Fraction
    from itertools import product

    P1 = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}

    def string(s):
        out = np.array([[1]], dtype=complex)
        for c in s:
            out = np.kron(out, P1[c])
        return out

    def moments(H, N, mmax):
        d = 2 ** N
        A = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
        Q = sum(np.kron(string('I' * l + 'Z' + 'I' * (N - l - 1)), string('I' * l + 'Z' + 'I' * (N - l - 1)))
                for l in range(N)).real
        Ar, Ai, Qi = (np.round(A.real).astype(object), np.round(A.imag).astype(object), np.round(Q).astype(object))
        tr = []
        for g in range(mmax + 1):
            R = Ar + g * Qi
            PR, PI = R.copy(), Ai.copy()
            row = [sum(PR[i, i] for i in range(d * d))]
            for _ in range(mmax - 1):
                PR, PI = PR.dot(R) - PI.dot(Ai), PR.dot(Ai) + PI.dot(R)
                row.append(sum(PR[i, i] for i in range(d * d)))
            tr.append(row)
        out = {}
        for m in range(1, mmax + 1, 2):
            vals = [Fraction(int(tr[g][m - 1])) for g in range(m + 1)]
            coef = vals[:]
            for j in range(1, m + 1):
                for i in range(m, j - 1, -1):
                    coef[i] = (coef[i] - coef[i - 1]) / j
            poly, basis = [Fraction(0)] * (m + 1), [Fraction(1)]
            for i in range(m + 1):
                for k_, b in enumerate(basis):
                    poly[k_] += coef[i] * b
                nb = [Fraction(0)] * (len(basis) + 1)
                for k_, b in enumerate(basis):
                    nb[k_ + 1] += b
                    nb[k_] -= i * b
                basis = nb
            out[m] = [int(c) for c in poly]
        return out

    def read(coef, N=3):
        H = sum(c * string(s) for s, c in coef.items())
        k = fw.moment_tower(H, N, j_max=8)['first_firing_rung']
        pm = moments(H, N, 2 * k + 1)
        ms = next(m for m in range(1, 2 * k + 2, 2) if any(pm[m]))
        return k, ms, {j: c for j, c in enumerate(pm[ms]) if c}

    assert read({'XYZ': -3, 'IYY': -3, 'IXX': 2, 'YXI': -1}) == (4, 7, {3: 580608})
    assert read({'ZIY': 1, 'XYI': -3, 'ZXX': 2, 'XIZ': 1}) == (4, 9, {1: 92897280})
    k5 = {'IXY': -2, 'XIY': 3, 'IYX': 2, 'YYY': 1, 'IXZ': 1}
    H5 = sum(c * string(s) for s, c in k5.items())
    F = string('XXX')
    assert np.array_equal(F @ H5 @ F, -H5)
    assert read(k5) == (5, 9, {3: 4644864})
    dich = {'XXZ': 1, 'IXZ': 1, 'XIY': 1, 'XYX': 1, 'IXI': 1}
    Hd = sum(c * string(s) for s, c in dich.items())
    t3 = [int(round(np.trace(string('I' * l + 'Z' + 'I' * (2 - l)) @ np.linalg.matrix_power(Hd, 3)).real)) for l in range(3)]
    assert t3 == [0, -16, 0] and not np.array_equal(F @ Hd @ F, -Hd)
    assert read(dich) == (2, 5, {1: 7680})
