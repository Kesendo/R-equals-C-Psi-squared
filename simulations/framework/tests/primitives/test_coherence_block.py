"""Tests for coherence_block primitives (block-restricted L for chromaticity / Q-scale analysis)."""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_chromaticity_F74_formula():
    # F74: c(n, N) = min(n, N-1-n) + 1
    assert fw.chromaticity(7, 0) == 1
    assert fw.chromaticity(7, 1) == 2
    assert fw.chromaticity(7, 2) == 3
    assert fw.chromaticity(7, 3) == 4   # central, c_max for odd N
    assert fw.chromaticity(7, 4) == 3
    assert fw.chromaticity(7, 6) == 1
    # even N: two adjacent c_max blocks
    assert fw.chromaticity(8, 3) == 4
    assert fw.chromaticity(8, 4) == 4


def test_block_basis_dimensions():
    P_n, P_np1 = fw.block_basis(5, 2)
    assert len(P_n) == comb(5, 2) == 10
    assert len(P_np1) == comb(5, 3) == 10
    # all states have correct popcount
    assert all(bin(p).count("1") == 2 for p in P_n)
    assert all(bin(q).count("1") == 3 for q in P_np1)


def test_block_L_split_xy_dimensions_and_dephasing():
    """D matrix has correct diagonal -2γ₀·HD entries; structure consistent."""
    N, n, gamma_0 = 4, 1, 0.05
    D, M_H_per_bond, P_n, P_np1 = fw.block_L_split_xy(N, n, gamma_0)
    Mn, Mnp1 = len(P_n), len(P_np1)
    Mtot = Mn * Mnp1
    assert D.shape == (Mtot, Mtot)
    assert len(M_H_per_bond) == N - 1
    for Mb in M_H_per_bond:
        assert Mb.shape == (Mtot, Mtot)
    # D is diagonal
    assert np.allclose(D - np.diag(np.diag(D)), 0)
    # diagonal entries are -2γ₀·HD with HD ∈ {1, 3} (chromaticity 2 here)
    diag = np.real(np.diag(D))
    for v in diag:
        assert abs(v + 0.10) < 1e-10 or abs(v + 0.30) < 1e-10  # -2γ₀·1 or -2γ₀·3


def test_block_L_split_xy_matches_full_construction():
    """Block-L matches full L on the (n, n+1) sector for N=4 XY chain."""
    N, n, gamma_0, J = 4, 1, 0.05, 1.0
    chain = fw.ChainSystem(N=N, gamma_0=gamma_0, J=J, H_type='xy')
    L_full = chain.L

    D, M_H_per_bond, P_n, P_np1 = fw.block_L_split_xy(N, n, gamma_0)
    L_block = D + J * sum(M_H_per_bond)

    # Block eigenvalues should appear in the full spectrum
    block_evals = np.linalg.eigvals(L_block)
    full_evals = np.linalg.eigvals(L_full)
    for be in block_evals:
        dist = float(np.min(np.abs(full_evals - be)))
        assert dist < 1e-9, f"block eigenvalue {be} not in full spectrum (min dist {dist})"


def test_hd_channel_basis_orthonormal():
    """Channel-uniform projectors are orthonormal."""
    N, n = 7, 3
    P, HDs = fw.hd_channel_basis(N, n)
    c = fw.chromaticity(N, n)
    assert P.shape[1] == c
    assert HDs == [1, 3, 5, 7]
    gram = P.conj().T @ P
    assert np.allclose(gram, np.eye(c), atol=1e-12)


def test_M_H_eff_diagonal_in_hd_channel_basis():
    """EQ-022 (b1) finding: M_H_eff is purely diagonal in HD-channel basis
    (off-diagonals exactly zero). Generalises F73 to all c.
    """
    N, n, gamma_0 = 6, 2, 0.05  # c = 3
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P, _ = fw.hd_channel_basis(N, n)
    M_H_eff = P.conj().T @ M_H_total @ P
    # Diagonal is purely imaginary (anti-Hermitian shape from L_H = -i[H, ·])
    # Off-diagonal should be exactly zero
    off_diag = M_H_eff - np.diag(np.diag(M_H_eff))
    assert np.allclose(off_diag, 0, atol=1e-12), \
        f"M_H_eff has nonzero off-diagonals: max |off| = {np.max(np.abs(off_diag))}"


def test_dicke_block_probe_lives_in_channel_subspace():
    """Dicke probe is a specific linear combination of channel-uniform vectors;
    its projection onto channel-uniform subspace recovers full norm.
    """
    N, n = 7, 3  # c = 4
    rho = fw.dicke_block_probe(N, n)
    P, _ = fw.hd_channel_basis(N, n)
    chan_components = P.conj().T @ rho
    norm_in_chan = np.linalg.norm(chan_components)
    norm_full = np.linalg.norm(rho)
    assert abs(norm_in_chan - norm_full) < 1e-12


def test_spatial_sum_kernel_F73_at_c1():
    """F73: at c=1 ((0, 1) block), the spatial-sum coherence purity for the
    Dicke probe |vac⟩⟨S_1|/2 evolves as (1/2)·exp(-4γ₀·t), independent of J.
    Evaluate at t=0: S(0) = 1/2 exactly.
    """
    N, n = 4, 0  # vac-SE block, c=1
    assert fw.chromaticity(N, n) == 1
    rho0 = fw.dicke_block_probe(N, n)
    S_kernel = fw.spatial_sum_coherence_kernel(N, n)
    S0 = float(np.real(rho0.conj() @ S_kernel @ rho0))
    assert abs(S0 - 0.5) < 1e-12, f"F73 S(t=0) should be 1/2, got {S0}"


def test_spatial_sum_kernel_hermitian_psd():
    """spatial_sum_coherence_kernel is Hermitian PSD per its docstring claim."""
    S = fw.spatial_sum_coherence_kernel(5, 2)
    assert np.allclose(S, S.conj().T, atol=1e-12)
    eigvals = np.linalg.eigvalsh(S)
    assert eigvals.min() >= -1e-12, f"min eigenvalue {eigvals.min()} < 0; not PSD"


# F86 EP-derived constant
# ---------------------------------------------------------------------------

def test_t_peak_inverse_4_gamma():
    """t_peak(γ₀) = 1/(4γ₀) — universal EP timescale (Tier 1 from EP algebra)."""
    assert abs(fw.t_peak(0.05) - 5.0) < 1e-12
    assert abs(fw.t_peak(0.1) - 2.5) < 1e-12
    assert abs(fw.t_peak(1.0) - 0.25) < 1e-12
