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

def test_confirmations_has_seven_entries():
    names = fw.Confirmations.list_names()
    assert len(names) == 7
    assert 'palindrome_trichotomy' in names
    assert 'lebensader_skeleton_trace_decoupling' in names


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
    assert len(marrakesh) >= 4
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
