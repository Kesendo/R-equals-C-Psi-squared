"""Tests for F81 Π-decomposition (Π·M·Π⁻¹ = M − 2·L_{H_odd}, sym/anti split, T1 violation diagnostic)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


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
