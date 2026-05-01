"""Tests for d=0 substrate diagnostics: stationary modes (kernel of L) +
d=0/d=2 state decomposition.

The d²−2d=0 condition forces d=0 or d=2. d=2 is the qubit. d=0 is the
substrate axis. The kernel of L gives direct access to that axis: for
uniform XY/Heisenberg + Z-dephasing it is exactly the F4 set of N+1
sector projectors P_n in the {I,Z}^N Pauli sublattice (n_xy = 0).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def _n_xy_of_pauli_index(k, N):
    """Count of X/Y letters in the Pauli string with framework index k."""
    from framework.pauli import _pauli_label
    label = _pauli_label(k, N)
    return sum(1 for ch in label if ch in ('X', 'Y'))


def test_stationary_modes_n3_xy_chain_kernel_dim_matches_f4():
    """N=3 XY chain + Z-dephasing has kernel dimension N+1 = 4 (F4 prediction:
    4 excitation-sector projectors P_0, P_1, P_2, P_3)."""
    chain = fw.ChainSystem(N=3, H_type='xy', gamma_0=0.05)
    sm = fw.stationary_modes(chain)
    assert sm['kernel_dimension'] == 4, \
        f"expected kernel_dim = 4 (F4), got {sm['kernel_dimension']}"


def test_stationary_modes_kernel_eigenvalues_are_zero():
    """All kernel-mode eigenvalues satisfy |λ| < tol."""
    chain = fw.ChainSystem(N=3, H_type='xy', gamma_0=0.05)
    sm = fw.stationary_modes(chain)
    assert np.all(np.abs(sm['eigenvalues']) < 1e-9)


def test_stationary_modes_kernel_lives_in_iz_sublattice():
    """For uniform XY/Heisenberg + Z-dephasing, every kernel mode has Pauli
    decomposition supported only on n_xy = 0 strings ({I,Z}^N sublattice).

    This is the F4 / lens-immune property: stationary modes are XOR-trivial
    under the polarity layer.
    """
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    sm = fw.stationary_modes(chain)
    decomp = sm['pauli_decomposition']  # n_kernel × 4^N
    for ii in range(decomp.shape[0]):
        for k in range(4 ** N):
            if _n_xy_of_pauli_index(k, N) > 0:
                assert abs(decomp[ii, k]) < 1e-9, \
                    f"kernel mode {ii} has nonzero coeff on n_xy>0 string k={k}"


def test_stationary_modes_n4_heisenberg_kernel_dim():
    """N=4 Heisenberg chain: F4 still gives N+1 = 5 sector projectors."""
    chain = fw.ChainSystem(N=4, H_type='heisenberg', gamma_0=0.05)
    sm = fw.stationary_modes(chain)
    assert sm['kernel_dimension'] == 5


def test_d_zero_decomposition_reconstructs_rho():
    """ρ = ρ_d0 + ρ_d2 exactly (up to numerical roundoff)."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    psi = fw.polarity_state(N, +1)
    rho = np.outer(psi, psi.conj())
    d = fw.d_zero_decomposition(rho, chain)
    rho_reconstructed = d['rho_d0'] + d['rho_d2']
    assert np.allclose(rho_reconstructed, rho, atol=1e-10)


def test_d_zero_decomposition_preserves_trace():
    """Tr(ρ_d0) = Tr(ρ) since L preserves trace and ρ_d0 = lim e^{Lt} ρ."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    psi = fw.polarity_state(N, +1)
    rho = np.outer(psi, psi.conj())
    d = fw.d_zero_decomposition(rho, chain)
    assert abs(d['d0_weight'] - 1.0) < 1e-9


def test_d_zero_decomposition_zero_excitation_state_fully_in_d0():
    """|0⟩^N⟨0|^N = P_0 (zero-excitation sector projector) is itself a kernel
    mode → fully d=0, d2_norm ≈ 0."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    rho = np.zeros((2 ** N, 2 ** N), dtype=complex)
    rho[0, 0] = 1.0  # |0...0⟩⟨0...0| = P_0
    d = fw.d_zero_decomposition(rho, chain)
    assert d['d2_norm'] < 1e-9, \
        f"|0⟩^N is P_0, should be fully d=0, but d2_norm = {d['d2_norm']:.2e}"
    assert np.allclose(d['rho_d0'], rho, atol=1e-9)


def test_d_zero_decomposition_full_excitation_state_fully_in_d0():
    """|1⟩^N⟨1|^N = P_N (fully-excited sector projector) is also a kernel mode."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    rho = np.zeros((2 ** N, 2 ** N), dtype=complex)
    rho[-1, -1] = 1.0  # |1...1⟩⟨1...1| = P_N
    d = fw.d_zero_decomposition(rho, chain)
    assert d['d2_norm'] < 1e-9


def test_d_zero_decomposition_polarity_state_has_nonzero_d2():
    """|+⟩^N has equal amplitude on all computational basis states, so it
    has off-diagonal coherences that decay under Z-dephasing — must have
    nonzero d=2 content."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    psi = fw.polarity_state(N, +1)
    rho = np.outer(psi, psi.conj())
    d = fw.d_zero_decomposition(rho, chain)
    assert d['d2_norm'] > 1e-3, \
        f"|+⟩^N has coherences, expected nonzero d2_norm; got {d['d2_norm']:.2e}"
    # Trace conservation still holds
    assert abs(d['d0_weight'] - 1.0) < 1e-9


def test_d_zero_decomposition_single_basis_state_partly_in_d0():
    """|k⟩⟨k| for k in the n=1 sector (e.g., |001⟩) is NOT P_1 by itself —
    XY mixes within the 1-excitation sector. Steady state is (1/3)·P_1.
    So d_zero_decomposition gives nonzero d=2 part (the difference)."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    rho = np.zeros((2 ** N, 2 ** N), dtype=complex)
    rho[1, 1] = 1.0  # |001⟩⟨001| (one excitation, but not the whole P_1)
    d = fw.d_zero_decomposition(rho, chain)
    # Trace preserved
    assert abs(d['d0_weight'] - 1.0) < 1e-9
    # But the full state is not in the kernel — Z-dephasing leaves it
    # invariant but XY mixes the n=1 manifold; only the symmetric mixture
    # (1/3)·P_1 is stationary.
    assert d['d2_norm'] > 1e-6


def test_d_zero_decomposition_steady_state_of_basis_state_is_sector_projector():
    """The d=0 part of |001⟩⟨001| is (1/3)·P_1 = (1/3)·(|001⟩⟨001| +
    |010⟩⟨010| + |100⟩⟨100|). Verify the diagonal is uniform 1/3 across
    the n=1 sector."""
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    rho = np.zeros((2 ** N, 2 ** N), dtype=complex)
    rho[1, 1] = 1.0  # |001⟩⟨001|
    d = fw.d_zero_decomposition(rho, chain)
    rho_d0 = d['rho_d0']
    # n=1 basis indices are 1 (=001), 2 (=010), 4 (=100)
    diag = np.real(np.diag(rho_d0))
    assert abs(diag[1] - 1.0 / 3.0) < 1e-6
    assert abs(diag[2] - 1.0 / 3.0) < 1e-6
    assert abs(diag[4] - 1.0 / 3.0) < 1e-6
    # Other sectors should be empty
    for idx in (0, 3, 5, 6, 7):
        assert abs(diag[idx]) < 1e-6


def test_sector_populations_basis_states():
    """|0⟩^N has p_0 = 1; |1⟩^N has p_N = 1; mean and variance trivial."""
    N = 3
    d_phys = 2 ** N

    psi_zero = np.eye(d_phys)[0]   # |000⟩
    sp = fw.sector_populations(psi_zero, N=N)
    assert sp['p'][0] == pytest.approx(1.0)
    assert sp['p'][1:].sum() == pytest.approx(0.0)
    assert sp['mean_n'] == pytest.approx(0.0)
    assert sp['var_n'] == pytest.approx(0.0)
    assert sp['entropy'] == pytest.approx(0.0)

    psi_one = np.eye(d_phys)[-1]   # |111⟩
    sp = fw.sector_populations(psi_one, N=N)
    assert sp['p'][N] == pytest.approx(1.0)
    assert sp['mean_n'] == pytest.approx(N)


def test_sector_populations_plus_state_is_binomial():
    """|+⟩^N has p_n = C(N,n) / 2^N (binomial), ⟨n⟩ = N/2, var = N/4."""
    from math import comb
    N = 4
    psi = fw.polarity_state(N, +1)
    sp = fw.sector_populations(psi, N=N)
    expected = np.array([comb(N, n) / (2 ** N) for n in range(N + 1)])
    assert np.allclose(sp['p'], expected, atol=1e-12)
    assert sp['mean_n'] == pytest.approx(N / 2)
    assert sp['var_n'] == pytest.approx(N / 4)


def test_sector_populations_within_sector_superposition():
    """(|001⟩+|010⟩)/√2 sits entirely in n=1; p_1 = 1, ⟨n⟩ = 1, var = 0."""
    N = 3
    d_phys = 2 ** N
    psi = (np.eye(d_phys)[1] + np.eye(d_phys)[2]) / np.sqrt(2)
    sp = fw.sector_populations(psi, N=N)
    assert sp['p'][1] == pytest.approx(1.0)
    assert sp['mean_n'] == pytest.approx(1.0)
    assert sp['var_n'] == pytest.approx(0.0)


def test_sector_populations_density_matrix_input():
    """Both ψ and ρ inputs supported, give the same answer."""
    N = 3
    psi = fw.polarity_state(N, +1)
    rho = np.outer(psi, psi.conj())
    sp_psi = fw.sector_populations(psi, N=N)
    sp_rho = fw.sector_populations(rho, N=N)
    assert np.allclose(sp_psi['p'], sp_rho['p'], atol=1e-12)


def test_sector_populations_conserved_under_xy_z_dephasing():
    """The d=0 axis observable: p_n is invariant under uniform XY+Z-dephasing
    Liouvillian. This is the operational meaning of "sector populations live
    in the kernel of L".

    Drift in p_n(t) would indicate non-{Z, XY/Heisenberg} noise (T1, σ±,
    transverse fields), which is the basis of the d=0 hardware diagnostic.
    """
    from scipy.linalg import expm
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    psi_0 = fw.polarity_state(N, +1)
    rho_0 = np.outer(psi_0, psi_0.conj())
    p_0 = fw.sector_populations(rho_0, N=N)['p']

    # Propagate to several timepoints and verify p stays fixed
    L = chain.L
    rho_vec_0 = rho_0.flatten('F')
    d_phys = 2 ** N
    for t in (0.1, 1.0, 5.0):
        rho_vec_t = expm(L * t) @ rho_vec_0
        rho_t = rho_vec_t.reshape(d_phys, d_phys, order='F')
        p_t = fw.sector_populations(rho_t, N=N)['p']
        assert np.allclose(p_t, p_0, atol=1e-9), \
            f"sector populations drifted at t={t}: max |Δp| = {np.max(np.abs(p_t - p_0)):.2e}"


def test_psi_vanishes_on_d_zero_substrate():
    """Algebraic bridge: Ψ(ρ_d0) = 0 for uniform XY/Heisenberg + Z-dephasing.

    The kernel of L lives in {I,Z}^N (sector projectors, no off-diagonal
    coherence in the computational basis), so ℓ₁(ρ_d0) = 0 by construction
    ⇒ Ψ(ρ_d0) = 0 exactly. All ℓ₁-coherence sits in ρ_d2.

    This closes the algebraic bridge between R = CΨ² (the Ψ-bearing
    recursion R = C(Ψ+R)² with bifurcation at CΨ = 1/4) and d²−2d = 0
    (the Ψ=0 dimension equation). Setting Ψ=0 and C=1/2 in the recursion
    gives R(R−2)=0 — the same quadratic, with d=0 and d=2 as its roots.
    The d=0 substrate IS the Ψ=0 fixed-point of the family.
    """
    N = 3
    chain = fw.ChainSystem(N=N, H_type='xy', gamma_0=0.05)
    d_phys = 2 ** N

    def psi_norm(rho_):
        # Ψ = ℓ₁/(d²-1), where ℓ₁ = sum of |off-diagonal entries|
        off = rho_ - np.diag(np.diag(rho_))
        return float(np.sum(np.abs(off))) / (d_phys ** 2 - 1)

    # Two states with different Ψ-content profiles
    plus_N = fw.polarity_state(N, +1)                       # uniform full coherence
    superpos_n1 = (np.eye(d_phys)[1] + np.eye(d_phys)[2]) / np.sqrt(2)  # |001⟩+|010⟩, within-sector

    for name, psi in (('|+⟩^N', plus_N), ('(|001⟩+|010⟩)/√2', superpos_n1)):
        rho = np.outer(psi, psi.conj())
        Psi_total = psi_norm(rho)
        decomp = fw.d_zero_decomposition(rho, chain)
        Psi_d0 = psi_norm(decomp['rho_d0'])
        Psi_d2 = psi_norm(decomp['rho_d2'])

        # ρ_d0 is in {I,Z}^N → no off-diagonals → Ψ = 0 exactly
        assert Psi_d0 < 1e-9, \
            f"[{name}] Ψ(ρ_d0) should vanish; got {Psi_d0:.2e}"
        # All ℓ₁-coherence survives intact in ρ_d2
        assert abs(Psi_d2 - Psi_total) < 1e-9, \
            f"[{name}] Ψ(ρ_d2)={Psi_d2:.6f} should match Ψ(ρ)={Psi_total:.6f}"
        # Sanity: state actually has nonzero coherence
        assert Psi_total > 1e-6, \
            f"[{name}] test only meaningful for Ψ(ρ)>0; got {Psi_total:.2e}"
