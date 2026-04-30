"""Tests for F78/F79 M-structure (anti-Hermitian, master lemma, single-body and 2-body decomposition)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


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
