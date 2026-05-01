"""Tests for the PTF (Perspectival Time Field) workflow + diagnostics.

Validates:
  - bond_perturbation primitive structure
  - pt_matrix_elements on slow modes (Π-invariance protection at first order)
  - ptf_alpha_fit returns trivial α=1 under no-defect (J_mod=1) and sensible
    α-pattern under a real defect

References: hypotheses/PERSPECTIVAL_TIME_FIELD.md (Tier 2). The N=7 reference
α-pattern (1.095, 1.182, 1.051, 0.991, 0.845, 0.923, 0.997) at J_mod=1.1, ψ_1
bonding mode, defect bond (0,1) is the canonical anchor — too slow for the
default test suite, validated separately.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def _bonding_mode_state(N):
    """ψ = (|vac⟩ + |ψ_1⟩) / √2 — the canonical PTF initial state."""
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * (i + 1) / (N + 1))
        psi[2 ** (N - 1 - i)] += amp
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def test_bond_perturbation_shape_and_antihermitian():
    """V_L is anti-Hermitian (it generates unitary evolution: V_L = -i[H, ·])."""
    N = 4
    V_L = fw.bond_perturbation(N, (0, 1), kind='XY')
    d2 = 4 ** N
    assert V_L.shape == (d2, d2)
    # Anti-Hermitian: V_L† = -V_L
    assert np.allclose(V_L.conj().T, -V_L, atol=1e-12)


def test_bond_perturbation_kinds():
    """All supported kinds build without error and produce non-zero V_L."""
    N = 3
    for kind in ('XY', 'heisenberg', 'XX', 'YY', 'ZZ'):
        V_L = fw.bond_perturbation(N, (0, 1), kind=kind)
        assert np.linalg.norm(V_L) > 0.1, f"{kind} V_L is zero"


def test_bond_perturbation_invalid_kind():
    """Invalid kind raises ValueError."""
    with pytest.raises(ValueError):
        fw.bond_perturbation(3, (0, 1), kind='ZX')


def test_pt_matrix_elements_shape():
    """Matrix elements have shape (n_slow, n_slow)."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    sm = fw.slow_modes(chain, n_slow=10)
    V_L = fw.bond_perturbation(4, (0, 1), kind='XY')
    ME = fw.pt_matrix_elements(sm, V_L)
    assert ME.shape == (10, 10)


def test_pt_matrix_elements_diagonal_protection_for_stationary():
    """For modes at λ ≈ 0 (stationary), ⟨W_s | V_L | M_s⟩ = 0 by U(1) excitation
    conservation: H_pert (XX+YY)/2 conserves total excitation, and the stationary
    modes are excitation-sector projectors. V_L · projector = 0 identically.
    """
    chain = fw.ChainSystem(N=4, H_type='xy')
    # Don't exclude stationary so we can test their protection
    sm = fw.slow_modes(chain, n_slow=20, exclude_stationary=False)
    V_L = fw.bond_perturbation(4, (0, 1), kind='XY')
    ME = fw.pt_matrix_elements(sm, V_L)
    # Find rows where the eigenvalue is near zero (stationary)
    stationary_mask = np.abs(sm['eigenvalues']) < 1e-9
    if stationary_mask.any():
        diag_at_stationary = np.diag(ME)[stationary_mask]
        assert np.all(np.abs(diag_at_stationary) < 1e-9), \
            f"stationary diagonal shifts {diag_at_stationary} not protected"


def test_pt_eigvec_shift_orthogonal_to_self():
    """First-order eigenvector shift δM_s is orthogonal to M_s (no self-mixing)."""
    chain = fw.ChainSystem(N=4, H_type='xy')
    sm = fw.slow_modes(chain, n_slow=8)
    V_L = fw.bond_perturbation(4, (0, 1), kind='XY')
    M_s = sm['right_eigvecs'][:, 0]
    W_s = sm['left_covecs'][0]
    delta = fw.pt_eigvec_shift(sm, V_L, mode_idx=0)
    # ⟨W_s | δM_s⟩ should be 0 by construction (the s'=s term is excluded)
    self_overlap = float(np.abs(W_s @ delta))
    assert self_overlap < 1e-10, \
        f"self-overlap {self_overlap}, expected 0 (s'=s excluded)"


def test_ptf_alpha_fit_no_defect_returns_unity():
    """When J_mod = 1.0 (no defect), α_i should all be 1.0 (no rescaling)."""
    N = 4
    chain = fw.ChainSystem(N=N, H_type='xy')
    rho_0 = _bonding_mode_state(N)
    result = fw.ptf_alpha_fit(chain, rho_0, defect_bond=0, J_mod=1.0,
                               t_max=10.0, n_t=100)
    assert np.allclose(result['alphas'], 1.0, atol=1e-3), \
        f"no-defect α should be 1; got {result['alphas']}"
    assert abs(result['sigma_log_alpha']) < 1e-3
    assert not result['on_boundary'].any()


def test_ptf_alpha_fit_defect_n4_sensible():
    """N=4 ψ_1 bonding mode + bond-0 defect at J_mod=1.1 gives α-pattern with:
      - α_0 > 1 (defect side speeds up)
      - all α in [0.5, 2.0] (perturbative range)
      - small RMSE (<0.01)
      - no boundary hits
    """
    N = 4
    chain = fw.ChainSystem(N=N, H_type='xy')
    rho_0 = _bonding_mode_state(N)
    result = fw.ptf_alpha_fit(chain, rho_0, defect_bond=0, J_mod=1.1,
                               t_max=20.0, n_t=200)
    alphas = result['alphas']
    assert alphas[0] > 1.0, f"α_0 = {alphas[0]}, expected > 1 at defect side"
    assert np.all(alphas > 0.5), f"some α < 0.5: {alphas}"
    assert np.all(alphas < 2.0), f"some α > 2: {alphas}"
    assert np.all(result['rmses'] < 0.01), \
        f"RMSEs too large: {result['rmses']}"
    assert not result['on_boundary'].any()


def test_ptf_painter_panel_with_matrix_elements():
    """Full painter panel includes V_L and slow-mode matrix elements."""
    N = 4
    chain = fw.ChainSystem(N=N, H_type='xy')
    rho_0 = _bonding_mode_state(N)
    panel = fw.ptf_painter_panel(chain, rho_0, defect_bond=0, J_mod=1.1,
                                  t_max=10.0, n_t=100,
                                  include_matrix_elements=True, n_slow=10)
    assert 'V_L' in panel
    assert 'slow_modes' in panel
    assert 'matrix_elements' in panel
    assert 'eigenvalue_shifts' in panel
    assert panel['V_L'].shape == (4 ** N, 4 ** N)
    assert panel['matrix_elements'].shape == (10, 10)
    # First-order shifts on slow modes: Re(δλ_s) ≈ 0 (Π-invariance + U(1)
    # excitation conservation protect the decay-rate component). Im(δλ_s)
    # is the energy shift, which is not generically protected.
    slow_re_shift = np.real(panel['eigenvalue_shifts'][:5])
    assert np.all(np.abs(slow_re_shift) < 1e-9), \
        f"slow-mode Re(δλ_s) not protected: {slow_re_shift}"


def test_ptf_alpha_fit_invalid_bond():
    """Out-of-range defect_bond raises."""
    N = 4
    chain = fw.ChainSystem(N=N, H_type='xy')
    rho_0 = _bonding_mode_state(N)
    with pytest.raises(ValueError):
        fw.ptf_alpha_fit(chain, rho_0, defect_bond=10, J_mod=1.1)
