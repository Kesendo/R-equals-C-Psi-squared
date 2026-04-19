#!/usr/bin/env python3
"""eq014_debug.py: quick sanity check on P_A values."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps

sys.stdout.reconfigure(encoding='utf-8')

sys.path.insert(0, str(Path(__file__).parent))
from eq014_step4567_closure import (
    load_bin, apply_biorth_fix, mode_to_matrix,
    partial_trace_all_sites, build_initial_states,
    compute_unperturbed_purity, build_V_L,
    TIME_WINDOW, DELTA_J,
)

RESULTS_DIR = Path(__file__).parent / "results"

with open(RESULTS_DIR / "eq014_biorth_metadata.json") as f:
    meta = json.load(f)
N = meta["N"]; d = meta["d"]; d2 = meta["d2"]

fix = np.load(RESULTS_DIR / "eq014_biorth_fix.npz")
c_diag = fix["c_diag"]
cluster_flat = fix["cluster_flat"]
cluster_offsets = fix["cluster_offsets"]
slow_indices = fix["slow_indices"]

vals = np.fromfile(RESULTS_DIR / "eq014_eigvals_n7.bin", dtype=np.complex128)
R = load_bin(RESULTS_DIR / "eq014_right_eigvecs_n7.bin", (d2, d2))
W = load_bin(RESULTS_DIR / "eq014_left_eigvecs_n7.bin", (d2, d2))
apply_biorth_fix(R, W, c_diag, cluster_flat, cluster_offsets)

R_slow = R[:, slow_indices]
W_slow = W[:, slow_indices]
slow_vals = vals[slow_indices]
n_slow = len(slow_indices)

print(f"Slow modes: {n_slow}")
print(f"Slow λ range: Re[{slow_vals.real.min():.4f}, {slow_vals.real.max():.4f}], "
      f"|Im|max={np.abs(slow_vals.imag).max():.4f}")

# Unperturbed purity at t=0 and t=1
slow_marginals = np.zeros((n_slow, N, 2, 2), dtype=complex)
for idx in range(n_slow):
    M = mode_to_matrix(R_slow[:, idx], d)
    slow_marginals[idx] = partial_trace_all_sites(M, d, N)

states = build_initial_states(N, d)

for name, rho_flat in states.items():
    c_s = W_slow.conj().T @ rho_flat
    print(f"\nState {name}:")
    print(f"  |c_s| top 5: {sorted(np.abs(c_s), reverse=True)[:5]}")
    print(f"  Total |c|^2 on slow subspace: {np.sum(np.abs(c_s)**2):.4f}")

    # P_A at a few times
    times_probe = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
    P_A = compute_unperturbed_purity(c_s, slow_vals, slow_marginals, times_probe)
    print(f"  P_A per site at t=0: {P_A[:, 0]}")
    print(f"  P_A per site at t=1: {P_A[:, 2]}")
    print(f"  P_A per site at t=4: {P_A[:, 4]}")

# Also check exact purity via direct evolution for comparison (t=1)
print("\n--- Exact vs slow-mode purity (psi_1_bonding, t=1) ---")
rho0 = states["psi_1_bonding"]
c_all = W.conj().T @ rho0
rho_t = R @ (c_all * np.exp(vals * 1.0))
rho_mat = rho_t.reshape((d, d))

margs = partial_trace_all_sites(rho_mat, d, N)
P_exact = np.array([float(np.real(np.trace(m @ m.conj().T))) for m in margs])
print(f"P_A exact (via all modes) at t=1:   {P_exact}")
P_A_psi1 = compute_unperturbed_purity(
    W_slow.conj().T @ rho0, slow_vals, slow_marginals, np.array([1.0]))
print(f"P_A slow-only at t=1:               {P_A_psi1[:, 0]}")

# Now check the perturbed side. Build V_L and compute mixing.
from eq014_step4567_closure import compute_perturbed_purity

V_L = build_V_L((0, 1), N)
print(f"\nV_L built, nnz={V_L.nnz}")

VR = V_L @ R_slow
A = W.conj().T @ VR
print(f"A shape: {A.shape}, |A| max: {np.abs(A).max():.3e}, mean: {np.abs(A).mean():.3e}")

delta_R_slow = np.zeros((d2, n_slow), dtype=complex)
for s_pos, s_glob in enumerate(slow_indices):
    denom = slow_vals[s_pos] - vals
    denom[s_glob] = 1.0
    coeffs = A[:, s_pos] / denom
    coeffs[s_glob] = 0.0
    delta_R_slow[:, s_pos] = R @ coeffs
print(f"delta_R_slow shape: {delta_R_slow.shape}, |max|: {np.abs(delta_R_slow).max():.3e}")

delta_slow_marginals = np.zeros((n_slow, N, 2, 2), dtype=complex)
for idx in range(n_slow):
    dM = mode_to_matrix(delta_R_slow[:, idx], d)
    delta_slow_marginals[idx] = partial_trace_all_sites(dM, d, N)
print(f"delta_slow_marginals |max|: {np.abs(delta_slow_marginals).max():.3e}")

WV = W_slow.conj().T @ V_L
B = WV @ R
print(f"B shape: {B.shape}, |B| max: {np.abs(B).max():.3e}")

rho0 = states["psi_1_bonding"]
c_s = W_slow.conj().T @ rho0
c_all = W.conj().T @ rho0

dc = np.zeros(n_slow, dtype=complex)
for s_pos, s_glob in enumerate(slow_indices):
    denom = slow_vals[s_pos] - vals
    denom[s_glob] = 1.0
    term = B[s_pos, :] / denom * c_all
    term[s_glob] = 0.0
    dc[s_pos] = np.sum(term)
print(f"dc |max|: {np.abs(dc).max():.3e}")
print(f"dc contains NaN: {np.any(np.isnan(dc))}")
print(f"Does A contain NaN? {np.any(np.isnan(A))}")
print(f"Does delta_R_slow contain NaN? {np.any(np.isnan(delta_R_slow))}")

# Check for denom near zero (degenerate eigenvalues)
for s_pos, s_glob in enumerate(slow_indices[:5]):
    denom = slow_vals[s_pos] - vals
    denom[s_glob] = 1.0
    # count how many denom values are tiny
    tiny = np.sum(np.abs(denom) < 1e-10)
    print(f"  slow[{s_pos}]=idx {s_glob}, λ={slow_vals[s_pos]:.4e}, "
          f"tiny denoms={tiny}, min|denom|={np.abs(denom).min():.2e}")

# Compute dP
dP_psi1 = compute_perturbed_purity(
    c_s, dc, slow_vals, slow_marginals, delta_slow_marginals, TIME_WINDOW)
print(f"\ndP shape: {dP_psi1.shape}, dP contains NaN: {np.any(np.isnan(dP_psi1))}")
print(f"dP |max|: {np.abs(dP_psi1).max():.3e}")
print(f"dP at t=0: {dP_psi1[:, 0]}")
print(f"dP at t=1: {dP_psi1[:, 10]}")
