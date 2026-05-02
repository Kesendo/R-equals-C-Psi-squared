#!/usr/bin/env python3
"""_eq022_b1_step_a_verify_blockL.py — EQ-022 (b1) Step (a).

Verify the Python block-L construction in `_eq022_b1_channel_projection.py`
against the C#-precomputed full N=7 XY chain Liouvillian eigendecomposition
saved in `eq014_*.bin` (J=1, γ=0.05, ChainXY).

Steps:
  1. Load eigvals (262 KB) + metadata.
  2. Build Python block-L for (N=7, n=3, J=1, γ=0.05). Block dim = 35*35
     = 1225. Diagonalize.
  3. Check every block-L eigenvalue appears in the saved full-L spectrum
     (within tolerance). This confirms Python construction.
  4. Optional: load right eigenvectors (4.3 GB), find which full-L
     eigenvectors live in the (n=3, n=4) popcount block, compare those
     eigenvalues to Python block-L spectrum.
  5. Compute Dicke probe |S_3⟩⟨S_4|/2 overlap with each block-L eigenvector.

Goal: prepare for step (b) (extended projection) by understanding the
block-L spectrum and where the Dicke probe lives in the eigenmode basis
at Q=20 (plateau).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "simulations" / "results"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def load_eigvals(path):
    return np.fromfile(path, dtype=np.complex128)


def main():
    # ---- 1. Load metadata + eigvals ----
    meta_path = RESULTS_DIR / "eq014_metadata.json"
    meta = json.loads(meta_path.read_text())
    print(f"# Metadata: N={meta['N']}, J={meta['J']}, γ={meta['gamma']}")
    print(f"# Topology: {meta['topology']}, Q = J/γ = {meta['J']/meta['gamma']:.1f}")
    print()

    N, J, gamma_0 = meta["N"], meta["J"], meta["gamma"]
    d = 2 ** N
    d2 = d * d

    eigvals_path = RESULTS_DIR / "eq014_eigvals_n7.bin"
    eigvals_full = load_eigvals(eigvals_path)
    print(f"Loaded {len(eigvals_full)} full-L eigenvalues from {eigvals_path.name}")
    print(f"  Re range: [{eigvals_full.real.min():.4f}, {eigvals_full.real.max():.4f}]")
    print()

    # ---- 2. Python block-L for (n=3, n=4) at N=7 ----
    n = 3  # gives c=4 at N=7
    c = fw.chromaticity(N, n)
    print(f"Building Python block-L for (N={N}, n={n}, J={J}, γ={gamma_0}), c={c}")
    t0 = time.time()
    D_block, M_H_per_bond, P_n_states, P_np1_states = fw.block_L_split_xy(N, n, gamma_0)
    L_block = D_block + J * sum(M_H_per_bond)
    Mn, Mnp1 = len(P_n_states), len(P_np1_states)
    print(f"  block dim = {Mn} x {Mnp1} = {Mn * Mnp1}")
    print(f"  built in {time.time() - t0:.2f} s")

    t0 = time.time()
    block_eigvals = np.linalg.eigvals(L_block)
    print(f"  diagonalized in {time.time() - t0:.2f} s")
    print(f"  Re range: [{block_eigvals.real.min():.4f}, {block_eigvals.real.max():.4f}]")
    print()

    # ---- 3. Verify each block eigenvalue appears in full-L spectrum ----
    print("## Sanity: each block-L eigenvalue should match a full-L eigenvalue")
    matched = 0
    max_dist = 0.0
    for be in block_eigvals:
        d_arr = np.abs(eigvals_full - be)
        idx = int(np.argmin(d_arr))
        dist = float(d_arr[idx])
        if dist < 1e-8:
            matched += 1
        max_dist = max(max_dist, dist)
    print(f"  matched within 1e-8: {matched} / {len(block_eigvals)}")
    print(f"  max nearest-neighbor distance: {max_dist:.2e}")
    if matched == len(block_eigvals):
        print("  → Python block-L construction VERIFIED against C# full-L.")
    else:
        print("  → MISMATCH — need to investigate basis convention.")
    print()

    # ---- 4. Block-L spectrum structure ----
    print("## Block-L spectrum at Q=20 (plateau region):")
    print("  Slowest 10 modes (largest Re(λ)):")
    sorted_idx = np.argsort(-block_eigvals.real)
    for i, j in enumerate(sorted_idx[:10]):
        bl = block_eigvals[j]
        print(f"    λ_{i:02d} = {bl.real:+.5f} {bl.imag:+.5f}j")
    print()
    print("  Fastest 5 modes (most negative Re(λ)):")
    for i, j in enumerate(sorted_idx[-5:]):
        bl = block_eigvals[j]
        print(f"    λ = {bl.real:+.5f} {bl.imag:+.5f}j")
    print()

    # ---- 5. Pure-rate ladder check (J=0 limit reference) ----
    print("## Reference: at J=0 the spectrum would collapse to pure-rate ladder")
    print(f"   c={c} → rates 2γ·HD for HD in [1, 3, 5, 7] = "
          f"{[2*gamma_0*hd for hd in [1, 3, 5, 7]]}")
    print()

    # ---- 6. Dicke probe overlap — eigenvalue identification ----
    # The probe |S_n><S_{n+1}|/2 is a uniform superposition; its overlap
    # with right eigenvectors v_j (in the block basis) tells us which modes
    # carry the Dicke weight at this Q.
    print("## Dicke probe |S_n><S_{n+1}|/2 in block-L eigenbasis (Q=20)")
    rho_probe = fw.dicke_block_probe(N, n)
    # Diagonalize block-L with eigenvectors
    t0 = time.time()
    block_eigvals_full, block_R = np.linalg.eig(L_block)
    print(f"  diagonalized with eigvecs in {time.time() - t0:.2f} s")

    # Inverse for left covectors (biorthogonal)
    block_W_T = np.linalg.inv(block_R)  # rows = left covectors

    # Probe overlap: c_j = w_j^T · rho_probe (where w_j is row j of W^T)
    overlap = block_W_T @ rho_probe
    overlap_abs = np.abs(overlap)

    # Sort by overlap magnitude
    order = np.argsort(-overlap_abs)
    print("  Top 10 modes by |c_j| = |⟨w_j | ρ_probe⟩|:")
    print(f"  {'idx':>5} {'|c_j|':>10} {'Re(λ)':>10} {'Im(λ)':>10}")
    cumul = 0.0
    total_w = float(np.sum(overlap_abs ** 2))
    for i in range(min(10, len(order))):
        j = order[i]
        lj = block_eigvals_full[j]
        c = overlap_abs[j]
        cumul += float(c ** 2)
        print(f"  {j:>5d} {c:>10.5f} {lj.real:>+10.5f} {lj.imag:>+10.5f}")
    print(f"  Cumulative weight in top 10: {cumul/total_w:.4f}")
    print()

    # ---- 7. Dicke probe in HD-channel-uniform basis (link to step b1 channel projection) ----
    print("## Dicke probe in HD-channel-uniform basis")
    P_chan, HDs = fw.hd_channel_basis(N, n)
    chan_amplitudes = P_chan.conj().T @ rho_probe
    print(f"  c×1 amplitude vector (channel components):")
    for k, hd in enumerate(HDs):
        print(f"    HD={hd}: amplitude = {chan_amplitudes[k]:+.5f}")
    print(f"  Norm in channel-uniform basis: {np.linalg.norm(chan_amplitudes):.5f}")
    print(f"  Norm of full probe: {np.linalg.norm(rho_probe):.5f}")
    print(f"  → fraction of probe in channel-uniform subspace: "
          f"{(np.linalg.norm(chan_amplitudes) / np.linalg.norm(rho_probe))**2:.4f}")
    print()
    print("  If ≈ 1: probe lives entirely in c-dim channel subspace.")
    print("  If < 1: probe has weight outside, where mixing happens.")


if __name__ == "__main__":
    main()
