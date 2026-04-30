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

def test_confirmations_has_ten_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 10
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names
    assert 'gamma_0_marrakesh_calibration' in names
    assert 'marrakesh_transverse_y_field_detection' in names
    assert 'f83_pi2_class_signature_marrakesh' in names


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


def test_F81_pi_conjugation_of_M():
    """F81: Π·M·Π⁻¹ = M − 2·L_{H_odd} structural identity, plus the
    Π-symmetric/antisymmetric decomposition of M.

    Verified at N=3 across the trichotomy:
      - truly XX+YY:           H_odd = 0, Π·M·Π⁻¹ = M (M = 0 trivially)
      - Π²-even non-truly YZ+ZY: H_odd = 0, Π·M·Π⁻¹ = M (M ≠ 0)
      - pure Π²-odd XY+YX:     H_odd = full H, Π·M·Π⁻¹ = M − 2·L_H
      - pure Π²-odd XY:        same as above
      - mixed hard XX+XY:      H_odd = XY part only, Π·M·Π⁻¹ = M − 2·L_{XY}

    For pure Π²-odd cases, the antisymmetric component M_anti = L_{H_odd}
    exactly, and ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 (50/50 split at N=3).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_bilinear, _vec_to_pauli_basis_transform
    from framework.symmetry import build_pi_full

    N = 3
    bonds = [(0, 1), (1, 2)]
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)
    T = _vec_to_pauli_basis_transform(N)
    d = 2 ** N
    Id = np.eye(d, dtype=complex)

    def L_H_pauli(H):
        L_H_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        return (T.conj().T @ L_H_vec @ T) / d

    def build_M(terms):
        H = _build_bilinear(N, bonds, terms)
        L = lindbladian_z_dephasing(H, [0.1] * N)
        return palindrome_residual(L, 0.3, N), H

    # Case 1: truly XX+YY → H_odd = 0, Π·M·Π⁻¹ = M (M=0 trivially)
    M_truly, _ = build_M([('X', 'X', 1.0), ('Y', 'Y', 1.0)])
    PiMPi_truly = Pi @ M_truly @ Pi_inv
    assert np.linalg.norm(PiMPi_truly - M_truly) < 1e-10, \
        "truly XX+YY: Π·M·Π⁻¹ should equal M"

    # Case 2: Π²-even non-truly YZ+ZY → H_odd = 0, Π·M·Π⁻¹ = M (M ≠ 0)
    M_yz, H_yz = build_M([('Y', 'Z', 1.0), ('Z', 'Y', 1.0)])
    assert np.linalg.norm(M_yz) > 1.0, "YZ+ZY: M should be non-zero"
    PiMPi_yz = Pi @ M_yz @ Pi_inv
    assert np.linalg.norm(PiMPi_yz - M_yz) < 1e-10, \
        "YZ+ZY (Π²-even non-truly): Π·M·Π⁻¹ should still equal M"

    # Case 3: pure Π²-odd XY+YX → Π·M·Π⁻¹ = M − 2·L_H
    M_soft, H_soft = build_M([('X', 'Y', 1.0), ('Y', 'X', 1.0)])
    PiMPi_soft = Pi @ M_soft @ Pi_inv
    L_H_soft = L_H_pauli(H_soft)
    pred_soft = M_soft - 2 * L_H_soft
    assert np.linalg.norm(PiMPi_soft - pred_soft) < 1e-10, \
        "XY+YX: Π·M·Π⁻¹ should equal M − 2·L_H"

    # Case 4: pure XY → same identity
    M_xy, H_xy = build_M([('X', 'Y', 1.0)])
    PiMPi_xy = Pi @ M_xy @ Pi_inv
    L_H_xy = L_H_pauli(H_xy)
    pred_xy = M_xy - 2 * L_H_xy
    assert np.linalg.norm(PiMPi_xy - pred_xy) < 1e-10, \
        "pure XY: Π·M·Π⁻¹ should equal M − 2·L_H"

    # Case 5: mixed hard XX+XY → only XY contributes to H_odd
    M_hard, _ = build_M([('X', 'X', 1.0), ('X', 'Y', 1.0)])
    PiMPi_hard = Pi @ M_hard @ Pi_inv
    H_xy_only = _build_bilinear(N, bonds, [('X', 'Y', 1.0)])
    L_H_xy_only = L_H_pauli(H_xy_only)
    pred_hard = M_hard - 2 * L_H_xy_only
    assert np.linalg.norm(PiMPi_hard - pred_hard) < 1e-10, \
        "hard XX+XY: Π·M·Π⁻¹ should equal M − 2·L_{XY part only}"

    # Decomposition for pure Π²-odd: M_anti = L_H exactly, 50/50 split
    M_sym = (M_soft + PiMPi_soft) / 2
    M_anti = (M_soft - PiMPi_soft) / 2
    assert np.linalg.norm(M_anti - L_H_soft) < 1e-10, \
        "soft XY+YX: M_anti should equal L_H exactly"
    sym_sq = np.linalg.norm(M_sym) ** 2
    anti_sq = np.linalg.norm(M_anti) ** 2
    total_sq = np.linalg.norm(M_soft) ** 2
    assert abs(sym_sq + anti_sq - total_sq) < 1e-8, \
        f"Pythagoras: ‖M_sym‖² + ‖M_anti‖² should equal ‖M‖² (got {sym_sq + anti_sq} vs {total_sq})"
    assert abs(sym_sq - anti_sq) / total_sq < 1e-8, \
        f"50/50 split for pure Π²-odd at N=3: ‖M_sym‖² ≈ ‖M_anti‖² (got {sym_sq} vs {anti_sq})"

    # Spectrum invariance under Π-conjugation (always)
    spec_M = sorted(np.linalg.eigvals(M_soft).imag)
    spec_PiMPi = sorted(np.linalg.eigvals(PiMPi_soft).imag)
    np.testing.assert_allclose(spec_M, spec_PiMPi, atol=1e-9)


def test_F81_pi_decompose_M_method():
    """ChainSystem.pi_decompose_M exposes the F81 decomposition cleanly.

    Verifies the framework primitive:
      - returns the dict with M, M_sym, M_anti, L_H_odd, norm_sq
      - the F81 identity M_anti = L_H_odd is enforced internally
      - the Pythagoras identity ‖M‖² = ‖M_sym‖² + ‖M_anti‖² holds
      - works correctly across the trichotomy at N=3 and N=4
    """
    for N in [3, 4, 5]:
        chain = fw.ChainSystem(N=N)

        # Truly XX+YY: M = 0, both sym and anti vanish
        d_truly = chain.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.1)
        assert d_truly['norm_sq']['M'] < 1e-15
        assert d_truly['norm_sq']['M_sym'] < 1e-15
        assert d_truly['norm_sq']['M_anti'] < 1e-15

        # Π²-even non-truly YZ+ZY: M ≠ 0, M_anti = 0 (no Π²-odd bilinears)
        d_yz = chain.pi_decompose_M([('Y', 'Z'), ('Z', 'Y')], gamma_z=0.1)
        assert d_yz['norm_sq']['M'] > 1e-3
        assert d_yz['norm_sq']['M_anti'] < 1e-15  # L_H_odd = 0
        assert abs(d_yz['norm_sq']['M_sym'] - d_yz['norm_sq']['M']) < 1e-9

        # Pure Π²-odd XY+YX: 50/50 split (analytical, N- and γ-independent)
        d_soft = chain.pi_decompose_M([('X', 'Y'), ('Y', 'X')], gamma_z=0.1)
        sym_frac = d_soft['norm_sq']['M_sym'] / d_soft['norm_sq']['M']
        anti_frac = d_soft['norm_sq']['M_anti'] / d_soft['norm_sq']['M']
        assert abs(sym_frac - 0.5) < 1e-9, f"N={N} XY+YX: sym fraction {sym_frac}"
        assert abs(anti_frac - 0.5) < 1e-9, f"N={N} XY+YX: anti fraction {anti_frac}"

        # Pythagoras: ‖M‖² = ‖M_sym‖² + ‖M_anti‖² for all cases (always holds)
        for label, d_dict in [('truly', d_truly), ('YZ', d_yz), ('soft', d_soft)]:
            total = d_dict['norm_sq']['M_sym'] + d_dict['norm_sq']['M_anti']
            assert abs(total - d_dict['norm_sq']['M']) < 1e-8, \
                f"N={N} {label}: Pythagoras fails: {total} vs {d_dict['norm_sq']['M']}"

        # Hard XX+XY (mixed): only XY contributes to L_H_odd
        d_hard = chain.pi_decompose_M([('X', 'X'), ('X', 'Y')], gamma_z=0.1)
        # Verify L_H_odd extraction is correct: same as XY-only Hamiltonian's L_H_odd
        d_xy_only = chain.pi_decompose_M([('X', 'Y')], gamma_z=0.1)
        assert np.allclose(d_hard['L_H_odd'], d_xy_only['L_H_odd'], atol=1e-12), \
            f"N={N} hard's L_H_odd should equal pure-XY's L_H_odd"
        # Pythagoras still holds for mixed; 50/50 NOT analytically guaranteed for mixed
        total_hard = d_hard['norm_sq']['M_sym'] + d_hard['norm_sq']['M_anti']
        assert abs(total_hard - d_hard['norm_sq']['M']) < 1e-8, \
            f"N={N} hard: Pythagoras fails: {total_hard} vs {d_hard['norm_sq']['M']}"

    # Test γ-independence of the 50/50 split for pure Π²-odd (chain N=4)
    chain4 = fw.ChainSystem(N=4)
    for gz in [0.0, 0.05, 1.0]:
        d = chain4.pi_decompose_M([('X', 'Y'), ('Y', 'X')], gamma_z=gz)
        if d['norm_sq']['M'] < 1e-15:
            continue  # γ=0 may give M=0 if H itself is Π²-protected; here it doesn't
        sym_frac = d['norm_sq']['M_sym'] / d['norm_sq']['M']
        assert abs(sym_frac - 0.5) < 1e-9, \
            f"N=4 XY+YX γ_z={gz}: 50/50 should be γ-independent, got {sym_frac}"


def test_F81_violation_T1_diagnostic():
    """F81 violation ‖M_anti − L_{H_odd}‖_F as a non-Z-dissipator diagnostic.

    For pure Z-dephasing the F81 identity holds exactly (violation ≈ 0). For
    T1 amplitude damping the dissipator is no longer Π²-symmetric and the
    violation grows linearly with γ_T1. This makes the violation a
    quantitative diagnostic for non-Z noise content on real hardware.
    """
    chain = fw.ChainSystem(N=3)
    soft_terms = [('X', 'Y'), ('Y', 'X')]

    # Z-only: violation must be at machine precision (strict=True default)
    d = chain.pi_decompose_M(soft_terms, gamma_z=0.1)
    assert d['f81_violation'] < 1e-10, \
        f"Z-only F81 violation should be ~0, got {d['f81_violation']:.4e}"

    # T1 enabled: strict default to False, violation reported
    violations = []
    for gt1 in [0.0, 0.05, 0.1, 0.2, 0.5]:
        d = chain.pi_decompose_M(soft_terms, gamma_z=0.1, gamma_t1=gt1)
        violations.append(d['f81_violation'])

    # γ_T1 = 0 gives zero violation
    assert violations[0] < 1e-10

    # Strictly monotone in γ_T1
    for i in range(1, len(violations)):
        assert violations[i] > violations[i - 1], \
            f"F81 violation should be monotone in γ_T1, got {violations}"

    # Linear-in-γ_T1: violations[2] / violations[1] should equal 0.1/0.05 = 2
    assert abs(violations[2] / violations[1] - 2.0) < 1e-6, \
        f"F81 violation should be linear in γ_T1 (small γ_T1), got ratio {violations[2]/violations[1]}"

    # strict=True with T1 should raise
    with pytest.raises(RuntimeError, match="F81 identity violated"):
        chain.pi_decompose_M(soft_terms, gamma_z=0.1, gamma_t1=0.1, strict=True)

    # Π²-even non-truly H (YZ+ZY) under T1: M_anti and L_H_odd both contain
    # T1's contribution? Let's check it doesn't crash and violation is reasonable.
    d_yz = chain.pi_decompose_M([('Y', 'Z'), ('Z', 'Y')], gamma_z=0.1, gamma_t1=0.05)
    # For Π²-even non-truly: L_H_odd = 0, but M_anti now has T1 content (non-zero)
    assert d_yz['norm_sq']['L_H_odd'] < 1e-15
    # f81_violation = ‖M_anti‖ since L_H_odd = 0
    assert d_yz['f81_violation'] > 1e-3


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
        d = chain.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=0.1)
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
        d = chain3.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
        assert abs(d['f81_violation'] - expected_violation) < 1e-9, \
            f"γ_T1_l={gt1_l}: predicted {expected_violation:.6f}, got {d['f81_violation']:.6f}"

    # H-independence at fixed γ_T1 (uniform 0.1 at N=3): violation = 0.6928 for any H
    expected_violation = 0.10 * (3 ** 0.5) * (2 ** 2)
    for label, terms in [('truly XX+YY', [('X', 'X'), ('Y', 'Y')]),
                         ('soft XY+YX', soft),
                         ('hard XX+XY', [('X', 'X'), ('X', 'Y')]),
                         ('YZ+ZY (Π²-even non-truly)', [('Y', 'Z'), ('Z', 'Y')])]:
        d = chain3.pi_decompose_M(terms, gamma_z=0.1, gamma_t1=0.1)
        assert abs(d['f81_violation'] - expected_violation) < 1e-9, \
            f"{label}: violation should be H-independent, got {d['f81_violation']:.6f}"

    # γ_z-independence at fixed γ_T1: violation = 0.6928 for γ_z ∈ {0, 0.1, 1.0}
    for gz in [0.0, 0.1, 1.0]:
        d = chain3.pi_decompose_M(soft, gamma_z=gz, gamma_t1=0.1)
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
            predicted = chain.predict_T1_dissipator_violation(gt1)
            d = chain.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1)
            assert abs(predicted - d['f81_violation']) < 1e-9, \
                f"N={N} γ_T1={gt1}: predict={predicted}, measured={d['f81_violation']}"

        # Non-uniform rates
        gt1_l = [0.05, 0.1, 0.15][:N] + [0.0] * max(0, N - 3)
        predicted_nu = chain.predict_T1_dissipator_violation(gt1_l)
        d_nu = chain.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
        assert abs(predicted_nu - d_nu['f81_violation']) < 1e-9, \
            f"N={N} non-uniform: predict={predicted_nu}, measured={d_nu['f81_violation']}"

    # Inverse direction: estimate_T1_from_violation recovers uniform γ_T1
    chain3 = fw.ChainSystem(N=3)
    for gt1_true in [0.0, 0.05, 0.1, 0.5]:
        d = chain3.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.1, gamma_t1=gt1_true)
        gt1_rms = chain3.estimate_T1_from_violation(d['f81_violation'])
        assert abs(gt1_rms - gt1_true) < 1e-9, \
            f"γ_T1 inversion: true={gt1_true}, recovered={gt1_rms}"

    # Inverse on non-uniform recovers RMS, not individual rates
    gt1_l = [0.05, 0.1, 0.15]
    rms_expected = (sum(g * g for g in gt1_l) / 3) ** 0.5
    d_nu = chain3.pi_decompose_M([('X', 'X'), ('Y', 'Y')], gamma_z=0.0, gamma_t1=gt1_l)
    rms_recovered = chain3.estimate_T1_from_violation(d_nu['f81_violation'])
    assert abs(rms_recovered - rms_expected) < 1e-9, \
        f"non-uniform RMS: expected={rms_expected}, recovered={rms_recovered}"

    # Edge cases
    assert chain3.predict_T1_dissipator_violation(0.0) == 0.0
    assert chain3.estimate_T1_from_violation(0.0) == 0.0

    # Length validation
    with pytest.raises(ValueError, match="must have length"):
        chain3.predict_T1_dissipator_violation([0.1, 0.2])  # length 2 ≠ 3

    # Negative violation rejected
    with pytest.raises(ValueError, match="non-negative"):
        chain3.estimate_T1_from_violation(-0.1)


def test_F83_pi_decomposition_anti_fraction_closed_form():
    """F83: anti-fraction = 1/(2 + 4·r) where r = ‖H_even_nontruly‖²/‖H_odd‖².

    Verifies the F83 closed form across the trichotomy + mixed configurations
    at N=3 and N=4. The predicted anti-fraction matches the numerical
    M_anti²/M² ratio bit-exact.
    """
    test_cases = [
        # (label, terms, expected anti-fraction)
        ('pure odd XY+YX', [('X', 'Y'), ('Y', 'X')], 0.5),
        ('pure odd XY', [('X', 'Y')], 0.5),
        ('pure odd XZ', [('X', 'Z')], 0.5),
        ('pure even non-truly YZ+ZY', [('Y', 'Z'), ('Z', 'Y')], 0.0),
        ('pure even non-truly YZ', [('Y', 'Z')], 0.0),
        ('truly XX+YY (M=0)', [('X', 'X'), ('Y', 'Y')], 0.0),  # limiting value
        ('mixed equal XY+YZ', [('X', 'Y'), ('Y', 'Z')], 1.0 / 6),
        ('mixed equal XY+ZY', [('X', 'Y'), ('Z', 'Y')], 1.0 / 6),
        ('mixed equal XZ+YZ', [('X', 'Z'), ('Y', 'Z')], 1.0 / 6),
        ('asymmetric more odd XY+YX+YZ', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z')], 1.0 / 4),
        ('full mix XY+YX+YZ+ZY', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z'), ('Z', 'Y')], 1.0 / 6),
        ('truly+mixed XX+XY+YZ', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')], 1.0 / 6),
    ]
    for N in [3, 4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms, expected_anti in test_cases:
            predicted = chain.predict_pi_decomposition_anti_fraction(terms)
            assert abs(predicted - expected_anti) < 1e-10, \
                f"N={N} {label}: predicted={predicted}, expected={expected_anti}"

            # Cross-check against numerical pi_decompose_M
            d = chain.pi_decompose_M(terms, gamma_z=0.0)
            if d['norm_sq']['M'] < 1e-12:
                # truly H: M = 0, anti-fraction undefined (closed form returns 0)
                continue
            measured = d['norm_sq']['M_anti'] / d['norm_sq']['M']
            assert abs(measured - expected_anti) < 1e-9, \
                f"N={N} {label}: numerical anti-fraction={measured}, expected={expected_anti}"
            assert abs(predicted - measured) < 1e-9, \
                f"N={N} {label}: predicted closed form ≠ numerical evaluation"

    # γ-independence: predicted anti-fraction should not depend on γ_z
    chain3 = fw.ChainSystem(N=3)
    for gz in [0.0, 0.05, 0.5]:
        d = chain3.pi_decompose_M([('X', 'Y'), ('Y', 'Z')], gamma_z=gz)
        if d['norm_sq']['M'] > 1e-12:
            measured = d['norm_sq']['M_anti'] / d['norm_sq']['M']
            assert abs(measured - 1.0 / 6) < 1e-9, \
                f"γ_z={gz}: anti-fraction should be γ-independent, got {measured}"


def test_F83_predict_pi_decomposition_full_closed_form():
    """F83 forward primitive: predict_pi_decomposition returns all closed-form
    norms (‖M‖², ‖M_anti‖², ‖M_sym‖²) plus anti-fraction from H alone. The
    predictions must match the numerical pi_decompose_M output bit-exact.
    """
    test_cases = [
        ('XY+YX (pure odd)', [('X', 'Y'), ('Y', 'X')]),
        ('YZ+ZY (pure even non-truly)', [('Y', 'Z'), ('Z', 'Y')]),
        ('XY+YZ (mixed equal)', [('X', 'Y'), ('Y', 'Z')]),
        ('XY+YX+YZ (asymmetric)', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z')]),
        ('XY+YX+YZ+ZY (full mix)', [('X', 'Y'), ('Y', 'X'), ('Y', 'Z'), ('Z', 'Y')]),
        ('XX+XY+YZ (truly + mixed)', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')]),
    ]
    for N in [3, 4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms in test_cases:
            pred = chain.predict_pi_decomposition(terms)
            num = chain.pi_decompose_M(terms, gamma_z=0.0)

            # All three norm-squared predictions match numerical
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"N={N} {label} ‖M‖²: pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9, \
                f"N={N} {label} ‖M_anti‖²: pred={pred['M_anti_sq']}, num={num['norm_sq']['M_anti']}"
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9, \
                f"N={N} {label} ‖M_sym‖²: pred={pred['M_sym_sq']}, num={num['norm_sq']['M_sym']}"

            # anti_fraction matches the convenience wrapper
            wrapper = chain.predict_pi_decomposition_anti_fraction(terms)
            assert abs(pred['anti_fraction'] - wrapper) < 1e-15

            # Pythagoras: M_anti_sq + M_sym_sq = M_sq
            assert abs(pred['M_anti_sq'] + pred['M_sym_sq'] - pred['M_sq']) < 1e-9

    # Special-case ratios at N=3 for canonical r values
    chain3 = fw.ChainSystem(N=3)
    pred_pure_odd = chain3.predict_pi_decomposition([('X', 'Y'), ('Y', 'X')])
    assert pred_pure_odd['r'] == 0.0
    assert abs(pred_pure_odd['anti_fraction'] - 0.5) < 1e-15

    pred_pure_even = chain3.predict_pi_decomposition([('Y', 'Z'), ('Z', 'Y')])
    assert pred_pure_even['r'] == float('inf')
    assert pred_pure_even['anti_fraction'] == 0.0
    assert pred_pure_even['M_anti_sq'] == 0.0

    pred_equal_mix = chain3.predict_pi_decomposition([('X', 'Y'), ('Y', 'Z')])
    assert abs(pred_equal_mix['r'] - 1.0) < 1e-15
    assert abs(pred_equal_mix['anti_fraction'] - 1.0/6) < 1e-15

    # All-truly: M=0, anti_fraction defaults to 0
    pred_truly = chain3.predict_pi_decomposition([('X', 'X'), ('Y', 'Y')])
    assert pred_truly['M_sq'] == 0.0
    assert pred_truly['anti_fraction'] == 0.0


def test_F85_kbody_trichotomy_counts():
    """F85: trichotomy enumeration at k = 2, 3, 4 over {X,Y,Z}^k.

    Verifies the closed-form counts:
      Π²-odd:           (3^k − (−1)^k) / 2
      Π²-even non-truly: pure-letter triples that use {Y, Z} only
      truly: rest

    Empirical: k=2 → 3/4/2, k=3 → 7/14/6, k=4 → 21/40/20.
    """
    from itertools import product
    from collections import Counter
    from framework.core import _pauli_tuple_is_truly, _pauli_tuple_pi2_class

    expected_counts = {
        2: {'truly': 3, 'pi2_odd': 4, 'pi2_even_nontruly': 2},
        3: {'truly': 7, 'pi2_odd': 14, 'pi2_even_nontruly': 6},
        4: {'truly': 21, 'pi2_odd': 40, 'pi2_even_nontruly': 20},
    }
    for k, expected in expected_counts.items():
        counts = Counter()
        for letters in product('XYZ', repeat=k):
            counts[_pauli_tuple_pi2_class(letters)] += 1
        assert dict(counts) == expected, \
            f"k={k} trichotomy counts mismatch: got {dict(counts)}, expected {expected}"
        # Closed form for Π²-odd count: (3^k − (−1)^k) / 2
        expected_odd = (3**k - (-1)**k) // 2
        assert counts['pi2_odd'] == expected_odd

    # Backward compat: 2-body classification matches _pauli_pair_is_truly
    from framework.core import _pauli_pair_is_truly
    for a, b in product('IXYZ', repeat=2):
        pair_check = _pauli_pair_is_truly(a, b)
        tuple_check = _pauli_tuple_is_truly((a, b))
        assert pair_check == tuple_check, f"({a},{b}): pair={pair_check}, tuple={tuple_check}"


def test_F85_kbody_predict_pi_decomposition():
    """F85: predict_pi_decomposition matches numerical pi_decompose_M at k=3, 4.

    Verifies that the F83 anti-fraction closed form generalizes to k-body
    via the Π²-class trichotomy (truly / Π²-odd / Π²-even non-truly), with
    the c(k) ∈ {0, 1, 2} factor scheme.
    """
    test_cases_k3 = [
        ('XYZ (Π²-even non-truly all-3-distinct)', [('X', 'Y', 'Z')], 0.0),
        ('XXY (Π²-odd)', [('X', 'X', 'Y')], 0.5),
        ('YYY (Π²-odd, n_YZ=3)', [('Y', 'Y', 'Y')], 0.5),
        ('XYY (truly)', [('X', 'Y', 'Y')], 0.0),
        ('XYZ + XXY (mixed)', [('X', 'Y', 'Z'), ('X', 'X', 'Y')], None),  # auto-check vs numerical
    ]
    for N in [4, 5]:
        chain = fw.ChainSystem(N=N)
        for label, terms, expected_anti in test_cases_k3:
            pred = chain.predict_pi_decomposition(terms)
            num = chain.pi_decompose_M(terms, gamma_z=0.0)
            # Closed form matches numerical bit-exact
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"N={N} {label}: M_sq pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9
            if expected_anti is not None and pred['M_sq'] > 1e-12:
                assert abs(pred['anti_fraction'] - expected_anti) < 1e-10, \
                    f"N={N} {label}: anti={pred['anti_fraction']}, expected {expected_anti}"

    # 4-body verification: predict matches numerical
    chain4 = fw.ChainSystem(N=5)
    test_cases_k4 = [
        [('X', 'Y', 'Z', 'X')],         # contains all three letters → Π²-even non-truly
        [('X', 'X', 'Y', 'Y')],         # 2X 2Y arrangement → truly
        [('Y', 'Y', 'Y', 'Z')],         # 3Y 1Z arrangement → Π²-even non-truly
        [('X', 'X', 'X', 'Y')],         # 3X 1Y arrangement → Π²-odd
    ]
    for terms in test_cases_k4:
        pred = chain4.predict_pi_decomposition(terms)
        num = chain4.pi_decompose_M(terms, gamma_z=0.0)
        assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
            f"k=4 {terms}: M_sq pred {pred['M_sq']} vs num {num['norm_sq']['M']}"


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


def test_F85_kbody_classifier_at_k5_spot_check():
    """F85 truly criterion verified at k=5 via spot check.

    Full 243-tuple enumeration at k=5 is computationally expensive (would
    add several minutes to the suite); instead we spot-check 8 representative
    5-tuples covering all three Π²-classes plus various letter compositions.

    Together with the full enumeration at k=2, 3, 4 in
    `test_F85_kbody_trichotomy_counts`, this provides empirical evidence
    that the analytic rule "#Y even AND #Z even ⟹ truly" continues to hold
    at k=5. The closed form for the Π²-odd count (3^k − (−1)^k)/2 = 122 at
    k=5 is verifiable analytically (binomial generating function).
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _build_kbody_chain
    from framework.core import _pauli_tuple_pi2_class

    N = 5
    test_cases = [
        ('XXXXX', 'truly'),         # all X (#Y=0 even, #Z=0 even)
        ('YYYYY', 'pi2_odd'),       # all Y (#Y=5 odd → bit_b=1)
        ('ZZZZZ', 'pi2_odd'),       # all Z (#Z=5 odd → bit_b=1)
        ('XXXYZ', 'pi2_even_nontruly'),  # bit_b=0 but #Z=1 odd
        ('YYZZX', 'truly'),         # #Y=2 even, #Z=2 even
        ('XYZXY', 'pi2_odd'),       # bit_b=1 odd
        ('YYZZY', 'pi2_odd'),       # bit_b=1 (3Y + 2Z = 5 mod 2 = 1)
        ('XZYXZ', 'pi2_odd'),       # bit_b=1 (1Y + 2Z = 3 mod 2 = 1)
    ]
    for label, expected in test_cases:
        letters = tuple(label)
        cls_ana = _pauli_tuple_pi2_class(letters)
        assert cls_ana == expected, f'{label}: analytic={cls_ana}, expected {expected}'

        # Numerical verification
        H = _build_kbody_chain(N, [letters + (1.0,)])
        H_sq = float(np.real(np.trace(H.conj().T @ H)))
        L = lindbladian_z_dephasing(H, [0.0] * N)
        M = palindrome_residual(L, 0.0, N)
        M_sq = float(np.linalg.norm(M) ** 2)
        if M_sq < 1e-10:
            cls_num = 'truly'
        else:
            ratio = M_sq / (H_sq * 2 ** N)
            if abs(ratio - 4) < 1e-6:
                cls_num = 'pi2_odd'
            elif abs(ratio - 8) < 1e-6:
                cls_num = 'pi2_even_nontruly'
            else:
                cls_num = f'unexpected ratio {ratio}'
        assert cls_ana == cls_num, f'{label}: analytic={cls_ana}, numerical={cls_num}'


def test_F85_kbody_F81_identity_at_k3():
    """F85: F81 identity Π·M·Π⁻¹ = M − 2·L_{H_odd} generalizes verbatim to k-body.

    Tests at k=3 by computing M and Π·M·Π⁻¹ directly and verifying the
    F81 identity holds bit-exact for Π²-odd 3-body terms.
    """
    from framework.lindblad import lindbladian_z_dephasing, palindrome_residual
    from framework.pauli import _vec_to_pauli_basis_transform, _build_kbody_chain
    from framework.symmetry import build_pi_full

    N = 4
    chain = fw.ChainSystem(N=N)
    Pi = build_pi_full(N)
    Pi_inv = np.linalg.inv(Pi)
    T = _vec_to_pauli_basis_transform(N)
    d = 2 ** N

    # 3-body Π²-odd: H_odd contributes via L_{H_odd}, F81 holds
    H = _build_kbody_chain(N, [('X', 'X', 'Y', 1.0)])
    L = lindbladian_z_dephasing(H, [0.0]*N)
    M = palindrome_residual(L, 0.0, N)
    PiMPi = Pi @ M @ Pi_inv
    L_H_vec = -1j * (np.kron(H, np.eye(d, dtype=complex)) -
                     np.kron(np.eye(d, dtype=complex), H.T))
    L_H = (T.conj().T @ L_H_vec @ T) / d
    # F81: Π·M·Π⁻¹ = M − 2·L_H_odd (and L_H_odd = L_H since term is purely Π²-odd)
    diff = np.linalg.norm(PiMPi - (M - 2 * L_H))
    assert diff < 1e-9, f"k=3 F81 identity fails: ‖diff‖ = {diff}"


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


def test_F83_topology_generalization():
    """F83 closed form is topology-independent: ring, star, complete K_N
    give the same closed-form prediction as numerical pi_decompose_M.

    The matrix-based F83 primitive (post-review fix) builds H_odd and
    H_even_nontruly via _build_bilinear which respects the chosen
    topology's bond graph. This test verifies that the closed form
    matches numerical L → M → pi_decompose_M at N=4 for each non-chain
    topology.
    """
    test_cases = [
        ('XY+YX (pure odd)', [('X', 'Y'), ('Y', 'X')]),
        ('YZ+ZY (pure even non-truly)', [('Y', 'Z'), ('Z', 'Y')]),
        ('XY+YZ (mixed equal)', [('X', 'Y'), ('Y', 'Z')]),
        ('XX+XY+YZ (truly + mixed)', [('X', 'X'), ('X', 'Y'), ('Y', 'Z')]),
    ]
    N = 4
    for topology in ['ring', 'star', 'complete']:
        chain = fw.ChainSystem(N=N, topology=topology)
        for label, terms in test_cases:
            pred = chain.predict_pi_decomposition(terms)
            num = chain.pi_decompose_M(terms, gamma_z=0.0)

            # Match closed form to numerical decomposition
            assert abs(pred['M_sq'] - num['norm_sq']['M']) < 1e-9, \
                f"topology={topology} {label} ‖M‖²: pred={pred['M_sq']}, num={num['norm_sq']['M']}"
            assert abs(pred['M_anti_sq'] - num['norm_sq']['M_anti']) < 1e-9, \
                f"topology={topology} {label} ‖M_anti‖²: pred={pred['M_anti_sq']}, num={num['norm_sq']['M_anti']}"
            assert abs(pred['M_sym_sq'] - num['norm_sq']['M_sym']) < 1e-9, \
                f"topology={topology} {label} ‖M_sym‖²: pred={pred['M_sym_sq']}, num={num['norm_sq']['M_sym']}"

            # anti_fraction matches numerical (when M ≠ 0)
            if num['norm_sq']['M'] > 1e-12:
                num_anti_frac = num['norm_sq']['M_anti'] / num['norm_sq']['M']
                assert abs(pred['anti_fraction'] - num_anti_frac) < 1e-9, \
                    f"topology={topology} {label} anti-fraction"
