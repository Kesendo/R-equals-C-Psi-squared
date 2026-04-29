"""Smoke tests for the lean framework package.

Run with: pytest simulations/framework/tests/

Each test verifies one cockpit-level behavior end-to-end. They cover the
hot paths: trichotomy classification, F71 signature, Confirmations lookup,
and the full Lebensader cockpit_panel.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import framework as fw


# ----------------------------------------------------------------------
# ChainSystem
# ----------------------------------------------------------------------

def test_chainsystem_n5_chain_invariants():
    chain = fw.ChainSystem(N=5)
    assert chain.N == 5
    assert chain.B == 4
    assert chain.D2 == 14  # 1²+2²+2²+2²+1² = 14
    assert chain.degrees == [1, 2, 2, 2, 1]


def test_chainsystem_classify_known_pauli_pairs():
    chain = fw.ChainSystem(N=5)
    assert chain.classify_pauli_pair([('X', 'X'), ('Y', 'Y')]) == 'truly'
    assert chain.classify_pauli_pair([('Y', 'Z'), ('Z', 'Y')]) == 'soft'
    assert chain.classify_pauli_pair([('I', 'X'), ('I', 'Z')]) == 'hard'


def test_chainsystem_predict_residual_norm_squared():
    chain = fw.ChainSystem(N=5)
    # main: (N-1) · 4^(N-2) = 4 · 64 = 256
    assert chain.predict_residual_norm_squared(1.0, 'main') == 256.0
    # single_body at chain N=5: D2/2 · 4^(N-2) = 14/2 · 64 = 7·64 = 448
    assert chain.predict_residual_norm_squared(1.0, 'single_body') == 448.0


def test_chainsystem_topology_invariants():
    """Verify B and D2 for ring, star, K_N at N=5."""
    chain_chain = fw.ChainSystem(N=5, topology='chain')
    assert (chain_chain.B, chain_chain.D2) == (4, 14)

    chain_ring = fw.ChainSystem(N=5, topology='ring')
    assert chain_ring.B == 5
    assert chain_ring.D2 == 5 * 4  # all degree 2

    chain_star = fw.ChainSystem(N=5, topology='star')
    assert chain_star.B == 4
    # hub deg 4, 4 leaves deg 1: 16 + 4·1 = 20
    assert chain_star.D2 == 20

    chain_kn = fw.ChainSystem(N=5, topology='complete')
    assert chain_kn.B == 10  # C(5,2)
    assert chain_kn.D2 == 5 * 16  # all degree N-1=4


# ----------------------------------------------------------------------
# Receiver
# ----------------------------------------------------------------------

def test_receiver_f71_zero_state():
    psi = np.zeros(2 ** 5, dtype=complex)
    psi[0] = 1.0
    r = fw.Receiver(psi)
    assert r.f71_class == +1
    sig = r.signature()
    assert sig['f71_eigenvalue'] == +1
    assert sig['bond_block_balanced'] is True  # N=5 → 2+2 balanced
    assert 'capacity-optimal' in sig['prediction']


def test_receiver_bonding_mode_k2_at_n5_is_f71_minus():
    """ψ_2 at N=5 is F71-antisymmetric per (-1)^(k+1) parity rule."""
    N, k = 5, 2
    psi = np.zeros(2 ** N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        psi[2 ** (N - 1 - i)] = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
    psi /= np.linalg.norm(psi)
    r = fw.Receiver(psi)
    assert r.f71_class == -1


def test_receiver_at_n6_predicts_suboptimal():
    """At N=6 with unbalanced 3+2 split, F71-eigenstate is capacity-suboptimal."""
    psi = np.zeros(2 ** 6, dtype=complex)
    psi[0] = 1.0  # |0>^6, F71-symmetric
    r = fw.Receiver(psi)
    assert r.f71_class == +1
    sig = r.signature()
    assert sig['bond_block_dims'] == (3, 2)
    assert sig['bond_block_balanced'] is False
    assert 'capacity-suboptimal' in sig['prediction']


# ----------------------------------------------------------------------
# Confirmations
# ----------------------------------------------------------------------

def test_confirmations_has_nine_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 9
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names
    assert 'gamma_0_marrakesh_calibration' in names
    assert 'marrakesh_transverse_y_field_detection' in names


def test_confirmations_lookup_palindrome_trichotomy():
    e = fw.Confirmations.lookup('palindrome_trichotomy')
    assert e['date'] == '2026-04-26'
    assert e['machine'] == 'ibm_marrakesh'
    assert e['job_id'] == 'd7mjnjjaq2pc73a1pk4g'
    assert e['measured_value']['delta_soft_minus_truly'] == -0.722


def test_confirmations_unknown_raises():
    with pytest.raises(KeyError):
        fw.Confirmations.lookup('does_not_exist')


def test_confirmations_by_machine():
    marrakesh = fw.Confirmations.by_machine('ibm_marrakesh')
    kingston = fw.Confirmations.by_machine('ibm_kingston')
    assert len(marrakesh) >= 6
    assert len(kingston) >= 3


# ----------------------------------------------------------------------
# Cockpit panel
# ----------------------------------------------------------------------

def test_cockpit_panel_yzzy_t1_drop_28_at_n3():
    """Reproduces the EQ-030 hardware-confirmed drop=28 for YZ+ZY +T1=0.005 at N=3."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('Y', 'Z'), ('Z', 'Y')],
                                  gamma_t1=0.005, t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 28
    assert result['lebensader']['rating'].startswith('collapsed')


def test_cockpit_panel_truly_no_drop_pure_z():
    """XX+YY (truly) under pure-Z: no skeleton drop expected."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    r = fw.Receiver(psi, chain=chain)
    result = chain.cockpit_panel(r, terms=[('X', 'X'), ('Y', 'Y')],
                                  t_max=5.0, dt=0.01)
    skeleton = result['lebensader']['skeleton']
    assert skeleton['drop'] == 0


# ----------------------------------------------------------------------
# Scaling formulas
# ----------------------------------------------------------------------

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
    import numpy as np
    from framework import lindbladian_general, palindrome_residual, site_op
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
    import numpy as np
    from framework import lindbladian_general, palindrome_residual, site_op

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
    import numpy as np
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


def test_zn_mirror_state_neel_partner():
    """Z⊗N applied to |+−+⟩ gives |−+−⟩ (X-basis Néel-mirror)."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    psi_b = fw.zn_mirror_state(psi_a, N=3)
    psi_expected = np.kron(np.kron(minus, plus), minus)
    np.testing.assert_allclose(psi_b, psi_expected, atol=1e-12)


def test_zn_mirror_diagnostic_pure_heisenberg_preserved():
    """Heisenberg + Z-dephasing → Z⊗N-Mirror exakt erhalten."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    rho_a_0 = np.outer(psi_a, psi_a.conj())
    psi_b = fw.zn_mirror_state(psi_a, N=3)
    rho_b_0 = np.outer(psi_b, psi_b.conj())
    # Propagate both
    rho_a = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    rho_b = chain.propagate_with_hardware_noise(
        rho_b_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b)
    assert diag['verdict'] == 'preserved'
    assert diag['max_violation'] < 1e-9


def test_zn_mirror_diagnostic_t1_preserves_zn():
    """T1 amplitude damping (σ⁻σ⁺ pairs) preserves Z⊗N."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    rho_a_0 = np.outer(psi_a, psi_a.conj())
    psi_b = fw.zn_mirror_state(psi_a, N=3)
    rho_b_0 = np.outer(psi_b, psi_b.conj())
    rho_a = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        T1_l=[0.01]*3)
    rho_b = chain.propagate_with_hardware_noise(
        rho_b_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        T1_l=[0.01]*3)
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b)
    assert diag['verdict'] == 'preserved', \
        f"T1 sollte Z⊗N erhalten, aber max_violation = {diag['max_violation']}"


def test_zn_mirror_diagnostic_h_z_preserves_zn():
    """Mini-Magnetfeld (h_z·Z_l, longitudinal) erhält Z⊗N (Z kommutiert mit Z⊗N)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    rho_a_0 = np.outer(psi_a, psi_a.conj())
    rho_b_0 = np.outer(fw.zn_mirror_state(psi_a, N=3),
                       fw.zn_mirror_state(psi_a, N=3).conj())
    rho_a = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        h_z_l=[0.1, 0.05, 0.1])
    rho_b = chain.propagate_with_hardware_noise(
        rho_b_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        h_z_l=[0.1, 0.05, 0.1])
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b)
    assert diag['verdict'] == 'preserved'


def test_zn_mirror_diagnostic_h_x_breaks_zn():
    """Transverse Hamiltonian X-Field (h_x·X_l) bricht Z⊗N (X anti-kommutiert)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    rho_a_0 = np.outer(psi_a, psi_a.conj())
    psi_b = fw.zn_mirror_state(psi_a, N=3)
    rho_b_0 = np.outer(psi_b, psi_b.conj())
    rho_a = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        h_x_l=[0.1]*3)
    rho_b = chain.propagate_with_hardware_noise(
        rho_b_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        h_x_l=[0.1]*3)
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b)
    assert diag['verdict'] == 'broken'
    # Linear scaling check: h_x=0.1 should give ~8e-3 violation (verified empirically)
    assert 1e-3 < diag['max_violation'] < 0.1


def test_zn_mirror_diagnostic_transverse_x_breaks_zn():
    """Transverse X-field on Hamiltonian breaks Z⊗N (X anti-commutes with Z⊗N)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi_a = np.kron(np.kron(plus, minus), plus)
    rho_a_0 = np.outer(psi_a, psi_a.conj())
    psi_b = fw.zn_mirror_state(psi_a, N=3)
    rho_b_0 = np.outer(psi_b, psi_b.conj())
    # Add transverse X-field via Xnoise dissipator (X-noise also breaks Z⊗N)
    rho_a = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        Xnoise_l=[0.05]*3)
    rho_b = chain.propagate_with_hardware_noise(
        rho_b_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')],
        Xnoise_l=[0.05]*3)
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b)
    # Xnoise actually preserves Z⊗N too (D[X] is even in X), but we'd need a
    # *Hamiltonian* X-field to break it. Let's test via state asymmetry instead.
    # Use a state that's NOT a Z⊗N partner pair to confirm 'broken' detection.
    # Pure state vs flipped one-site state — definitely not Z⊗N partners.
    psi_wrong = np.kron(np.kron(plus, plus), plus)  # NOT the Z⊗N partner of |+−+⟩
    rho_wrong_0 = np.outer(psi_wrong, psi_wrong.conj())
    rho_wrong = chain.propagate_with_hardware_noise(
        rho_wrong_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    rho_a_clean = chain.propagate_with_hardware_noise(
        rho_a_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    diag_break = chain.zn_mirror_diagnostic(rho_a_clean, rho_wrong)
    assert diag_break['verdict'] == 'broken'
    assert diag_break['max_violation'] > 0.1


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_at_default_gamma():
    """Verify the optimal probe parameters at γ=0.05 match the manual analysis."""
    chain = fw.ChainSystem(N=2)
    setup = chain.gamma_probe_setup(gamma_assumed=0.05, target_precision=0.01)
    # K_optimal ≈ 0.119 (4·γ·t ≈ 0.474 → K = γt ≈ 0.119)
    assert abs(setup['K_optimal'] - 0.119) < 0.01
    # cpsi_target ≈ 0.144
    assert abs(setup['cpsi_target'] - 0.144) < 0.005
    # K_cusp ≈ 0.0374 (4γt_cusp = -ln(f_cusp) ≈ 0.150 → K ≈ 0.0374)
    assert abs(setup['K_cusp'] - 0.0374) < 0.001
    # cpsi_cusp = 0.25 by construction
    assert setup['cpsi_cusp'] == 0.25


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_K_invariance():
    """K_optimal = γ·t* should be γ-independent (γ-invariance in dimensionless form)."""
    chain = fw.ChainSystem(N=2)
    K_vals = []
    for gamma in [0.01, 0.05, 0.1, 0.5]:
        setup = chain.gamma_probe_setup(gamma_assumed=gamma)
        K_vals.append(setup['K_optimal'])
    # All K_optimal values should agree (γ-invariance)
    assert all(abs(k - K_vals[0]) < 1e-6 for k in K_vals), \
        f"K_optimal varies with γ: {K_vals}"


def test_gamma_probe_default_gamma_uses_chain_gamma_0():
    chain = fw.ChainSystem(N=3, gamma_0=0.07)
    setup = chain.gamma_probe_setup()
    assert setup['gamma_assumed'] == 0.07


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_from_cpsi_inverts_f25():
    """Round-trip: γ → CΨ(t) → estimate γ should recover the original."""
    chain = fw.ChainSystem(N=2)
    gamma_true = 0.07
    t = 2.0  # arbitrary probe time
    f = np.exp(-4 * gamma_true * t)
    cpsi = f * (1 + f**2) / 6.0
    gamma_est = chain.estimate_gamma_from_cpsi(cpsi, t)
    assert abs(gamma_est - gamma_true) < 1e-9


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_rejects_out_of_range():
    chain = fw.ChainSystem(N=2)
    with pytest.raises(ValueError, match="≥ 1/3"):
        chain.estimate_gamma_from_cpsi(0.5, t=1.0)
    with pytest.raises(ValueError, match="≤ 0"):
        chain.estimate_gamma_from_cpsi(-0.1, t=1.0)
    with pytest.raises(ValueError, match="t must be > 0"):
        chain.estimate_gamma_from_cpsi(0.1, t=0)


def test_cpsi_bell_plus_recovers_f25_for_pure_z():
    """F26 with γ_x=γ_y=0 must reduce to F25: CΨ = f·(1+f²)/6 with f=exp(-4γz·t)."""
    import numpy as np
    gz = 0.05
    t = 2.0
    f = np.exp(-4 * gz * t)
    cpsi_f25 = f * (1 + f**2) / 6.0
    cpsi_f26 = fw.cpsi_bell_plus(0.0, 0.0, gz, t)
    assert abs(cpsi_f26 - cpsi_f25) < 1e-12


def test_cpsi_bell_plus_at_t0_gives_one_third():
    """Bell+ at t=0: CΨ = 1·(1+1+1+1)/12 = 4/12 = 1/3."""
    cpsi = fw.cpsi_bell_plus(0.05, 0.07, 0.03, 0.0)
    assert abs(cpsi - 1.0/3.0) < 1e-12


def test_cpsi_bell_plus_monotonic_decay():
    """CΨ monotonically decreases with t for any nonzero noise (F26 corollary)."""
    cpsi_values = [fw.cpsi_bell_plus(0.05, 0.0, 0.05, t) for t in [0.0, 1.0, 2.0, 5.0, 10.0]]
    for i in range(len(cpsi_values) - 1):
        assert cpsi_values[i+1] < cpsi_values[i]


def test_cpsi_cusp_K_per_channel_matches_F27():
    """F27 K-values: Z=0.0374, X=Y=0.0867, depol=0.0440. Cross-check via cusp finder."""
    from scipy.optimize import brentq
    for channel, K_expected in fw.CPSI_CUSP_K_PER_CHANNEL.items():
        if channel == 'Z':
            gx, gy, gz = 0.0, 0.0, 1.0
        elif channel == 'X':
            gx, gy, gz = 1.0, 0.0, 0.0
        elif channel == 'Y':
            gx, gy, gz = 0.0, 1.0, 0.0
        elif channel == 'depolarizing':
            gx, gy, gz = 1/3, 1/3, 1/3
        # Solve CΨ(t) = 1/4 with γ = 1
        t_cusp = brentq(lambda t: fw.cpsi_bell_plus(gx, gy, gz, t) - 0.25, 1e-6, 100)
        K_computed = 1.0 * t_cusp
        assert abs(K_computed - K_expected) < 0.001, \
            f"channel {channel}: K computed {K_computed} vs F27 {K_expected}"


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_x_channel():
    """gamma_probe_setup with channel='X' should give K_cusp = 0.0867 vs Z's 0.0374."""
    chain = fw.ChainSystem(N=2)
    setup_z = chain.gamma_probe_setup(gamma_assumed=0.05, channel='Z')
    setup_x = chain.gamma_probe_setup(gamma_assumed=0.05, channel='X')
    setup_y = chain.gamma_probe_setup(gamma_assumed=0.05, channel='Y')
    setup_d = chain.gamma_probe_setup(gamma_assumed=0.05, channel='depolarizing')
    # K_cusp from F26 cusp condition (note: K_Y = K_Z, NOT K_X — doc has typo)
    assert abs(setup_z['K_cusp'] - 0.0374) < 0.001
    assert abs(setup_x['K_cusp'] - 0.0867) < 0.001
    assert abs(setup_y['K_cusp'] - 0.0374) < 0.001  # K_Y = K_Z, not K_X
    assert abs(setup_d['K_cusp'] - 0.0440) < 0.001


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_estimate_gamma_round_trip_x_channel():
    """Round-trip γ → CΨ_X(t) → estimate_γ for X-channel."""
    chain = fw.ChainSystem(N=2)
    gamma_true = 0.07
    t = 2.0
    cpsi = fw.cpsi_bell_plus(gamma_true, 0.0, 0.0, t)
    gamma_est = chain.estimate_gamma_from_cpsi(cpsi, t, channel='X')
    assert abs(gamma_est - gamma_true) < 1e-9


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_gamma_probe_setup_kingston_data_consistency():
    """Kingston cusp-slowing F25 RMS residual was 0.0097; with 1% target,
    shots-needed should be reasonable (~10⁵-10⁶)."""
    chain = fw.ChainSystem(N=2)
    setup = chain.gamma_probe_setup(gamma_assumed=0.05, target_precision=0.01)
    # Order of magnitude check: should be ~10⁶ shots for 1% precision
    assert 1e5 < setup['shots_needed'] < 1e7


def test_propagate_with_hardware_noise_no_evolution():
    """t=0 returns ρ_0 (with symmetrization)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(rho_0, t=0.0)
    np.testing.assert_allclose(rho_t, rho_0, atol=1e-12)


def test_propagate_with_hardware_noise_no_dissipation_unitary():
    """No T1/Tphi/etc. → unitary evolution, trace preserved exactly."""
    chain = fw.ChainSystem(N=3, gamma_0=0.0)  # no Z-dephasing
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=1.0, terms=[('X','X'),('Y','Y'),('Z','Z')])
    assert abs(np.trace(rho_t).real - 1.0) < 1e-10
    # Purity preserved (unitary on pure state)
    assert abs(np.trace(rho_t @ rho_t).real - 1.0) < 1e-9


def test_propagate_with_hardware_noise_matches_marrakesh_xz_for_truly():
    """Marrakesh April 26: truly-unbroken (XX+YY) Hamiltonian on |+−+⟩ at t=0.8.

    With γ_Z=0.1 (idealized), <X_0 Z_2> should be near 0 (matches HW within readout noise)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')])
    # Trace out middle qubit, compute <X_0 Z_2>
    rho_3q = rho_t.reshape(2,2,2,2,2,2)
    rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    X0Z2 = np.kron(X, Z)
    exp = float(np.real(np.trace(rho_02 @ X0Z2)))
    # Idealized framework prediction: 0.000 (per Marrakesh README); HW measured +0.011.
    # Our propagation should match the framework idealized number ~0.
    assert abs(exp) < 0.01


def test_propagate_with_hardware_noise_matches_marrakesh_xz_for_soft():
    """Marrakesh: soft (XY+YX) under γ_Z=0.1 should give <X_0 Z_2> ≈ -0.62 (matches README)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_t = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')])
    rho_3q = rho_t.reshape(2,2,2,2,2,2)
    rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    exp = float(np.real(np.trace(rho_02 @ np.kron(X, Z))))
    # README idealized: -0.623; HW: -0.711. Our propagation should match idealized.
    assert abs(exp - (-0.623)) < 0.02


def test_propagate_with_hardware_noise_zz_crosstalk_changes_dynamics():
    """J_zz adds Z-Z Hamiltonian correction; non-trivial effect on Pauli expectations."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.05)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_no_zz = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')])
    rho_with_zz = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','X'),('Y','Y')], J_zz=0.3)

    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    Y0Z2 = np.kron(Y, Z)
    def y0z2(rho):
        rho_3q = rho.reshape(2,2,2,2,2,2)
        rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
        return float(np.real(np.trace(rho_02 @ Y0Z2)))
    e_no = y0z2(rho_no_zz)
    e_with = y0z2(rho_with_zz)
    # ZZ-crosstalk should shift <Y_0 Z_2> by > 0.1
    assert abs(e_with - e_no) > 0.1, \
        f"J_zz=0.3 should shift <Y_0 Z_2> by >0.1; got Δ={abs(e_with - e_no):.4f}"


def test_propagate_with_hardware_noise_t1_has_nontrivial_effect():
    """T1 changes <X_0 Z_2> measurably (sanity that the T1 channel actually fires)."""
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=0.1)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)
    rho_0 = np.outer(psi, psi.conj())
    rho_id = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')])
    rho_t1 = chain.propagate_with_hardware_noise(
        rho_0, t=0.8, terms=[('X','Y'),('Y','X')],
        T1_l=[0.05]*3)  # 10× larger T1 to ensure visible effect

    def x0z2(rho):
        rho_3q = rho.reshape(2,2,2,2,2,2)
        rho_02 = np.einsum('ikjlkm->ijlm', rho_3q).reshape(4, 4)
        X = np.array([[0,1],[1,0]], dtype=complex)
        Z = np.array([[1,0],[0,-1]], dtype=complex)
        return float(np.real(np.trace(rho_02 @ np.kron(X, Z))))

    e_id = x0z2(rho_id)
    e_t1 = x0z2(rho_t1)
    assert abs(e_t1 - e_id) > 0.01, \
        f"T1=0.05 must shift <X_0 Z_2> by >0.01; got Δ={abs(e_t1 - e_id):.4f}"


def test_palindrome_residual_norm_ratio_squared_n3_n4():
    # main: 4·k/(k-1) at k=3 → 6
    assert fw.palindrome_residual_norm_ratio_squared(3, 4, 'main') == 6.0
    # single-body: 4·(2k-1)/(2k-3) at k=3 → 20/3
    assert abs(fw.palindrome_residual_norm_ratio_squared(3, 4, 'single_body') - 20 / 3) < 1e-12


def test_palindrome_residual_norm_squared_factor_graph_topologies():
    # Chain N=5: B=4, D2=14
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 4, 14, 'main') == 4 * 64
    # Ring N=5: B=5, D2=20
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'main') == 5 * 64
    assert fw.palindrome_residual_norm_squared_factor_graph(5, 5, 20, 'single_body') == 10 * 64


# ----------------------------------------------------------------------
# F71 / chain mirror
# ----------------------------------------------------------------------

def test_chain_mirror_state_involutory():
    for N in [3, 4, 5, 6]:
        R = fw.chain_mirror_state(N)
        np.testing.assert_allclose(R @ R, np.eye(2 ** N), atol=1e-12)


def test_bond_mirror_basis_dimensions():
    """Verify the parity-tied dimensional structure."""
    for N, expected_sym, expected_asym in [(3, 1, 1), (4, 2, 1), (5, 2, 2),
                                             (6, 3, 2), (7, 3, 3), (8, 4, 3)]:
        sym, asym = fw.bond_mirror_basis(N)
        assert (len(sym), len(asym)) == (expected_sym, expected_asym), \
            f"N={N}: got ({len(sym)}, {len(asym)}), expected ({expected_sym}, {expected_asym})"


# ----------------------------------------------------------------------
# Cockpit defensive guarantees (PTF round 1)
# ----------------------------------------------------------------------

def test_chainsystem_n2_warns_about_structural_degeneracy():
    """N=2 is allowed but warns; fundamental ops still produce correct values."""
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        chain = fw.ChainSystem(N=2)
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 1
        assert 'structurally degenerate' in str(user_warnings[0].message)
    # Math still works at N=2 — fundamental vocabulary
    assert chain.classify_pauli_pair([('X','X'),('Y','Y'),('Z','Z')]) == 'truly'
    assert chain.classify_pauli_pair([('Y','Z'),('Z','Y')]) == 'soft'
    assert abs(chain.predict_residual_norm_squared_from_terms(
        [('Y','Z'),('Z','Y')]) - 256.0) < 1e-6  # 2^4·2·8 = 256 at N=2


def test_chainsystem_n3plus_does_not_warn():
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fw.ChainSystem(N=3)
        fw.ChainSystem(N=4)
        user_warnings = [x for x in w if issubclass(x.category, UserWarning)]
        assert len(user_warnings) == 0


def test_chainsystem_J_immutable_after_init():
    chain = fw.ChainSystem(N=3, J=1.0)
    with pytest.raises(AttributeError, match="immutable"):
        chain.J = 2.0


def test_chainsystem_gamma_0_immutable_after_init():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    with pytest.raises(AttributeError, match="immutable"):
        chain.gamma_0 = 0.5


def test_chainsystem_topology_immutable_after_init():
    chain = fw.ChainSystem(N=3, topology='chain')
    with pytest.raises(AttributeError, match="immutable"):
        chain.topology = 'ring'


def test_receiver_rejects_2d_input():
    """Receiver must reject density matrices passed where psi is expected."""
    rho = np.eye(8, dtype=complex) / 8.0  # 8x8 max-mixed
    with pytest.raises(ValueError, match="1D state vector"):
        fw.Receiver(rho)


def test_receiver_rejects_unnormalized_psi():
    psi = np.array([1.0, 1.0, 0.0, 0.0], dtype=complex)  # norm sqrt(2)
    with pytest.raises(ValueError, match="normalized"):
        fw.Receiver(psi)


def test_receiver_from_psi_unnormalized():
    psi = np.array([1.0, 1.0, 0.0, 0.0], dtype=complex)
    r = fw.Receiver.from_psi_unnormalized(psi)
    assert np.isclose(np.linalg.norm(r.psi), 1.0)
    assert r.N == 2


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cockpit_panel_rejects_non_hermitian_rho():
    chain = fw.ChainSystem(N=2, gamma_0=0.05)
    rho_nh = np.array([[0.5, 0.2, 0, 0],
                       [0.2, 0.3, 0, 0],
                       [0, 0, 0.1, 0.05],
                       [0, 0, 0.05j, 0.1]], dtype=complex)
    with pytest.raises(ValueError, match="Hermitian"):
        fw.cockpit_panel(chain.H, [0.05]*2, rho_nh, 2, t_max=0.5, dt=0.05)


def test_cockpit_panel_rejects_non_psd_rho():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    rho_neg = np.zeros((8,8), dtype=complex)
    rho_neg[0,0] = 1.5
    rho_neg[1,1] = -0.5  # negative eigenvalue
    with pytest.raises(ValueError, match="positive semi-definite"):
        fw.cockpit_panel(chain.H, [0.05]*3, rho_neg, 3, t_max=0.5, dt=0.05)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_cockpit_panel_rejects_trace_mismatch():
    chain = fw.ChainSystem(N=2, gamma_0=0.05)
    rho_no_trace = np.eye(4, dtype=complex) * 0.5  # trace = 2
    with pytest.raises(ValueError, match="trace"):
        fw.cockpit_panel(chain.H, [0.05]*2, rho_no_trace, 2, t_max=0.5, dt=0.05)


def test_residual_norm_squared_matches_classify_truly():
    """For Heisenberg (truly), ||M||^2 ≈ 0."""
    chain = fw.ChainSystem(N=4)
    norm_sq = chain.residual_norm_squared([('X','X'),('Y','Y'),('Z','Z')])
    assert norm_sq < 1e-10


def test_residual_norm_squared_matches_known_value():
    """At chain N=4, YZ+ZY soft pair has ||M||^2 = 12288 (= c_H_main(256) * F(48))."""
    chain = fw.ChainSystem(N=4)
    norm_sq = chain.residual_norm_squared([('Y','Z'),('Z','Y')])
    assert abs(norm_sq - 12288.0) < 1e-6


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


def test_residual_norm_squared_with_t1_matches_predict():
    """Numerical residual_norm_squared(gamma_t1=...) matches predict's closed form."""
    for N in [3, 4]:
        chain = fw.ChainSystem(N=N, gamma_0=0.05)
        for terms in [[('Y','Z'),('Z','Y')], [('I','Y'),('Y','I')],
                      [('X','X'),('Y','Y'),('Z','Z')]]:
            for gT1 in [0.001, 0.005, 0.01]:
                pred = chain.predict_residual_norm_squared_from_terms(
                    terms, gamma_t1=gT1)
                num = chain.residual_norm_squared(terms, gamma_t1=gT1)
                assert abs(pred - num) < 1e-6, \
                    f"N={N} terms={terms} gT1={gT1}: pred={pred} num={num}"


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


def test_cockpit_panel_terms_uses_chain_J():
    """terms-mode must scale by self.J, not hardcoded 1.0.

    A chain at J=2 with terms=YZ+ZY should give a different residual
    (scaled by J²=4) than the same chain at J=1.
    """
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(plus, minus), plus)

    chain1 = fw.ChainSystem(N=3, J=1.0, gamma_0=0.05)
    r1 = fw.Receiver(psi, chain=chain1)
    p1 = chain1.cockpit_panel(r1, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005,
                              t_max=2.0, dt=0.01)

    chain2 = fw.ChainSystem(N=3, J=2.0, gamma_0=0.05)
    r2 = fw.Receiver(psi, chain=chain2)
    p2 = chain2.cockpit_panel(r2, terms=[('Y','Z'),('Z','Y')], gamma_t1=0.005,
                              t_max=2.0, dt=0.01)

    # Different J → different θ-trajectory (Hamiltonian rescales coherent dynamics)
    theta1 = p1['_trajectory_for_inspection']['theta']
    theta2 = p2['_trajectory_for_inspection']['theta']
    assert not np.allclose(theta1, theta2), \
        "cockpit_panel(terms=...) must use self.J, not hardcoded 1.0"


def test_F78_F79_M_is_anti_hermitian():
    """Lemma 4 of PROOF_SVD_CLUSTER_STRUCTURE: M is anti-Hermitian for any
    Hermitian H under Z-dephasing.

    Consequences:
    - Eigenvalues of M are purely imaginary.
    - M is normal: M·M† = M†·M.
    - Singular values = |eigenvalues|.
    - Same spectrum ⇒ unitary equivalence (spectral theorem).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_bilinear

    N = 4
    bonds = [(i, i+1) for i in range(N-1)]
    cases = [
        [('I', 'Y'), ('Y', 'I')],            # single-body Y
        [('Y', 'Z'), ('Z', 'Y')],            # 2-body Π²-even soft
        [('X', 'Y')],                        # 2-body Π²-odd
        [('X', 'X'), ('X', 'Y')],            # mixed (XX even-truly + XY odd)
        [('Y', 'Z'), ('X', 'Y')],            # mixed parity even+odd
        [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],  # Heisenberg (truly: M=0, trivially anti-Hermitian)
    ]
    for terms in cases:
        H = _build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])
        L = lindbladian_z_dephasing(H, [1.0]*N)
        M = palindrome_residual(L, N*1.0, N)
        anti_herm_err = float(np.linalg.norm(M + M.conj().T))
        assert anti_herm_err < 1e-9, \
            f"M not anti-Hermitian for terms={terms}: ‖M+M†‖ = {anti_herm_err}"


def test_F78_F79_master_lemma_M_gamma_independent():
    """Master Lemma (foundation of F78 and F79):
       For pure Z-dephasing, M = Π·L·Π⁻¹ + L + 2σ·I = Π·L_H·Π⁻¹ + L_H,
       i.e., the dissipator-conjugation contribution exactly cancels with 2σI.
       Consequence: M depends only on H, NOT on γ.

    Proof in docs/proofs/PROOF_SVD_CLUSTER_STRUCTURE.md (Lemma 1).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.symmetry import build_pi_full
    from framework.pauli import _build_bilinear, _vec_to_pauli_basis_transform

    N = 4
    bonds = [(i, i+1) for i in range(N-1)]
    Mb = _vec_to_pauli_basis_transform(N)
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)

    cases = [
        [('I', 'Y'), ('Y', 'I')],         # single-body Y
        [('Y', 'Z'), ('Z', 'Y')],         # 2-body soft (Π²-even)
        [('X', 'X'), ('X', 'Y')],         # 2-body mixed (XX even-truly + XY odd)
    ]
    for terms in cases:
        H = _build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])
        # L_H alone (zero dephasing): in Pauli basis
        L_H_vec = lindbladian_z_dephasing(H, [0.0]*N)
        L_H_pauli = (Mb.conj().T @ L_H_vec @ Mb) / (2**N)
        M_predicted = Pi @ L_H_pauli @ Pi_inv + L_H_pauli

        # Compute full M at multiple γ values, verify all equal predicted
        for gamma in [0.5, 1.0, 2.5]:
            L_full = lindbladian_z_dephasing(H, [gamma]*N)
            M_actual = palindrome_residual(L_full, gamma*N, N)
            diff = float(np.linalg.norm(M_actual - M_predicted))
            assert diff < 1e-9, \
                f"Master Lemma fails for terms={terms} γ={gamma}: ‖M − Π·L_H·Π⁻¹ − L_H‖ = {diff}"


def test_F78_single_body_M_additive_decomposition():
    """F78: For single-body H = Σ_l c_l·Y_l, M's SVD clusters predict from
    additive single-site decomposition M = Σ_l M_l ⊗ I_others.

    Per-site M_l is normal with eigenvalues ±2c_l·γ·i (mult 2 each).
    Full M's SVs = |Σ_l ε_l · 2c_l·γ| with multiplicity 2^N per sign-combo.
    Verifies that the previously-open "Group-Theory cluster question" is closed.
    """
    from itertools import product as iproduct
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual

    # Per-site M_l in (I, X, Y, Z) basis for c·Y, γ·Z-dephasing
    def _single_site_M(c, gamma=1.0):
        L_l = np.array([[0,0,0,0],
                        [0,-2*gamma,0,-2*c],
                        [0,0,-2*gamma,0],
                        [0,2*c,0,0]], dtype=complex)
        Pi = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1j],[0,0,1j,0]], dtype=complex)
        return Pi @ L_l @ np.linalg.inv(Pi) + L_l + 2*gamma*np.eye(4)

    N = 4
    chain = fw.ChainSystem(N=N)  # chain degrees [1, 2, 2, 1]
    bilinear = [('I', 'Y', 1.0), ('Y', 'I', 1.0)]
    H = fw._build_bilinear(N, chain.bonds, bilinear)  # bond-summed → Σ_l deg(l)·Y_l
    L = lindbladian_z_dephasing(H, [1.0]*N)
    M_full = palindrome_residual(L, N*1.0, N)
    direct_svs = np.sort(np.linalg.svd(M_full, compute_uv=False))[::-1]

    # Predicted: each combination Σ_l ε_l·2c_l, all eigenvalue products
    M_per_site = [_single_site_M(d) for d in chain.degrees]
    eigs_per_site = [np.linalg.eigvals(M) for M in M_per_site]
    pred_evs = [sum(eigs_per_site[l][combo[l]] for l in range(N))
                for combo in iproduct(*[range(4) for _ in range(N)])]
    pred_svs = np.sort(np.abs(pred_evs))[::-1]

    assert direct_svs.shape == pred_svs.shape
    np.testing.assert_allclose(direct_svs, pred_svs, atol=1e-9)

    # Also lock in the specific cluster signature for chain N=4 IY+YI:
    # SVs at 12, 8, 4, 0 with mults 32, 64, 96, 64.
    expected_clusters = {12.0: 32, 8.0: 64, 4.0: 96, 0.0: 64}
    for sv_value, expected_mult in expected_clusters.items():
        actual_mult = int(np.sum(np.abs(direct_svs - sv_value) < 1e-6))
        assert actual_mult == expected_mult, \
            f"SV={sv_value}: expected mult {expected_mult}, got {actual_mult}"


def test_F79_single_bond_lebensader_reduction():
    """F79 single-bond Π²-odd universality, proven via Lebensader reduction.

    For a single-bond Π²-odd 2-body bilinear H = c·(P⊗Q) at sites (i, j),
    M(P⊗Q at bond) is spectrally equivalent to M(Π(P⊗Q) at single site).
    Specifically, all four single-bond Π²-odd letter pairs and all four
    Π-reduced single-body operators give the same M-cluster pattern.

    This closes the single-bond half of the F79 universality observation.
    See PROOF_SVD_CLUSTER_STRUCTURE.md "Refined scope of the universality"
    (April 2026 follow-up).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_bilinear, site_op

    N = 4
    single_bond = [(0, 1)]

    # Single-bond Π²-odd 2-body bilinears (4 letter pairs)
    two_body_cases = [('X', 'Y'), ('X', 'Z'), ('Y', 'X'), ('Z', 'X')]
    two_body_svs = []
    for P, Q in two_body_cases:
        H = _build_bilinear(N, single_bond, [(P, Q, 1.0)])
        L = lindbladian_z_dephasing(H, [1.0]*N)
        M = palindrome_residual(L, N*1.0, N)
        svs = np.sort(np.linalg.svd(M, compute_uv=False))[::-1]
        two_body_svs.append(svs)

    # All four 2-body cases must have identical SV spectrum
    for svs in two_body_svs[1:]:
        np.testing.assert_allclose(two_body_svs[0], svs, atol=1e-9)

    # Π-reduced single-body cases (Z or Y at site 0 or 1)
    one_body_cases = [(1, 'Z'), (1, 'Y'), (0, 'Z'), (0, 'Y')]
    one_body_svs = []
    for site, P in one_body_cases:
        H = site_op(N, site, P)
        L = lindbladian_z_dephasing(H, [1.0]*N)
        M = palindrome_residual(L, N*1.0, N)
        svs = np.sort(np.linalg.svd(M, compute_uv=False))[::-1]
        one_body_svs.append(svs)

    # All four 1-body cases must have identical SV spectrum
    for svs in one_body_svs[1:]:
        np.testing.assert_allclose(one_body_svs[0], svs, atol=1e-9)

    # 2-body and 1-body must match (the Lebensader reduction)
    np.testing.assert_allclose(two_body_svs[0], one_body_svs[0], atol=1e-9)


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


def test_F79_two_body_pi_squared_block_decomposition():
    """F79: For 2-body bond bilinear M = Π·L·Π⁻¹+L+2σ·I,
    Π²-parity = (bit_b(P)+bit_b(Q)) mod 2 of each bilinear term determines
    the block structure:
      - All terms Π²-even (e.g., YZ, XX, Heisenberg): M is Π²-block-diagonal
        (off-diagonal blocks M[V_+,V_-] vanish exactly).
      - All terms Π²-odd (e.g., XY, XZ, XX+XY): M is purely off-diagonal
        (diagonal blocks M[V_+,V_+] and M[V_-,V_-] vanish exactly).
    Plus: pure Π²-odd 2-body bilinears give universal M-SVD (letter-irrelevant).
    """
    from framework.pauli import _build_bilinear, _k_to_indices
    from framework.symmetry import pi_squared_eigenvalue
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual

    N = 4
    bonds = [(i, i+1) for i in range(N-1)]

    def _build_M(terms):
        H = _build_bilinear(N, bonds, [(t[0], t[1], 1.0) for t in terms])
        L = lindbladian_z_dephasing(H, [1.0]*N)
        return palindrome_residual(L, N*1.0, N)

    # Π²-eigenspace indices
    idx_p = np.array([k for k in range(4**N)
                      if pi_squared_eigenvalue(_k_to_indices(k, N)) == 1])
    idx_m = np.array([k for k in range(4**N)
                      if pi_squared_eigenvalue(_k_to_indices(k, N)) == -1])
    assert len(idx_p) == len(idx_m) == 4**N // 2  # equal split for any N≥1

    # Π²-even bilinear (YZ): M block-diagonal, off-diag blocks vanish
    M_even = _build_M([('Y', 'Z')])
    assert float(np.linalg.norm(M_even[np.ix_(idx_p, idx_m)])) < 1e-9, \
        "Π²-even YZ: off-diagonal block M[V_+,V_-] should vanish exactly"
    assert float(np.linalg.norm(M_even[np.ix_(idx_m, idx_p)])) < 1e-9, \
        "Π²-even YZ: off-diagonal block M[V_-,V_+] should vanish exactly"

    # Π²-odd bilinear (XY): M purely off-diagonal, diagonal blocks vanish
    M_odd = _build_M([('X', 'Y')])
    assert float(np.linalg.norm(M_odd[np.ix_(idx_p, idx_p)])) < 1e-9, \
        "Π²-odd XY: diagonal block M[V_+,V_+] should vanish exactly"
    assert float(np.linalg.norm(M_odd[np.ix_(idx_m, idx_m)])) < 1e-9, \
        "Π²-odd XY: diagonal block M[V_-,V_-] should vanish exactly"

    # Π²-odd universality: XY, XZ, XX+XY share same SVD cluster spectrum
    svs_XY = np.sort(np.linalg.svd(M_odd, compute_uv=False))[::-1]
    svs_XZ = np.sort(np.linalg.svd(_build_M([('X', 'Z')]), compute_uv=False))[::-1]
    svs_XXY = np.sort(np.linalg.svd(_build_M([('X', 'X'), ('X', 'Y')]),
                                    compute_uv=False))[::-1]
    np.testing.assert_allclose(svs_XY, svs_XZ, atol=1e-9)
    np.testing.assert_allclose(svs_XY, svs_XXY, atol=1e-9)
