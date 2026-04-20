#!/usr/bin/env python3
"""eq014_validate_groundtruth.py

Validate first-order PT prediction against exact time evolution of ρ_0
under L_B = L_A + δJ V_L, both for ψ_1 bonding mode.

Method:
 1. Load L_A eigendecomp data (XY chain).
 2. Compute P_A(t) exactly via Σ_s c_s e^{λ_s t} v_s + marginal + |·|².
 3. Build L_B = L_A + δJ V_L as sparse/dense d² × d² matrix.
 4. Evolve ρ_0 under L_B via RK4 to get P_B(t).
 5. Fit α_i from (P_A, P_B) with same algorithm as first-order.
 6. Compare to first-order PT prediction and to PTF empirical values.
"""
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
    partial_trace_all_sites, build_initial_states, build_V_L,
    fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results"
DELTA_J = 0.1
# Match PTF: integrate to T=80, fit over [0, 20], DT=0.2 sample.
TIME_WINDOW = np.arange(401) * 0.2
T_FIT = 20.0


def rebuild_L_A(N, J, gamma):
    """Build L_A sparse for XY chain (PTF convention) to cross-check and
    to serve as the base for L_B = L_A + δJ V_L.

    H_A = Σ_i (J/2)(X_i X_{i+1} + Y_i Y_{i+1}),  big-endian site indexing."""
    d = 2 ** N
    dd = d * d
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)

    def site_op(op, k):
        full = sps.eye(1, dtype=complex, format='csr')
        I2 = sps.eye(2, dtype=complex, format='csr')
        for j in range(N):
            full = sps.kron(full, sps.csr_matrix(op) if j == k else I2, format='csr')
        return full

    H = sps.csr_matrix((d, d), dtype=complex)
    for i in range(N - 1):
        H = H + (J / 2) * (site_op(X, i) @ site_op(X, i + 1)
                           + site_op(Y, i) @ site_op(Y, i + 1))
    Id = sps.eye(d, dtype=complex, format='csr')
    L_h = -1j * (sps.kron(H, Id, format='csr') - sps.kron(Id, H.T, format='csr'))
    a = np.arange(dd) // d
    b = np.arange(dd) % d
    h_vec = np.array([bin(int(av ^ bv)).count('1') for av, bv in zip(a, b)])
    L_d = sps.diags(-2.0 * gamma * h_vec.astype(complex), format='csr')
    return (L_h + L_d).tocsr()


def rk4_step(L, rho, dt):
    """RK4 step for drho/dt = L rho."""
    k1 = L @ rho
    k2 = L @ (rho + 0.5 * dt * k1)
    k3 = L @ (rho + 0.5 * dt * k2)
    k4 = L @ (rho + dt * k3)
    return rho + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def evolve_and_sample(L, rho0_flat, d, N, sample_times, dt_small=0.005):
    """RK4-evolve ρ_0 under L, sampling P_i(t) at sample_times.
    Use adaptive subdivision of time intervals so dt ≤ dt_small."""
    sample_times = np.asarray(sample_times)
    n_samples = len(sample_times)
    out = np.zeros((N, n_samples), dtype=float)
    current_rho = rho0_flat.astype(complex).copy()
    t_current = 0.0
    for si, t_target in enumerate(sample_times):
        dt_total = t_target - t_current
        if dt_total > 0:
            n_steps = max(1, int(np.ceil(dt_total / dt_small)))
            dt = dt_total / n_steps
            for _ in range(n_steps):
                current_rho = rk4_step(L, current_rho, dt)
            t_current = t_target
        M = current_rho.reshape((d, d))
        margs = partial_trace_all_sites(M, d, N)
        for i in range(N):
            out[i, si] = float(np.real(np.trace(margs[i] @ margs[i].conj().T)))
    return out


def main():
    with open(RESULTS_DIR / "eq014_biorth_metadata.json") as f:
        meta = json.load(f)
    N = meta["N"]; d = meta["d"]; d2 = meta["d2"]
    gamma_0 = meta["gamma"]; J = meta["J"]

    # Load the eigendecomp + biorth
    fix = np.load(RESULTS_DIR / "eq014_biorth_fix.npz")
    c_diag = fix["c_diag"]
    cluster_flat = fix["cluster_flat"]
    cluster_offsets = fix["cluster_offsets"]

    vals = np.fromfile(RESULTS_DIR / "eq014_eigvals_n7.bin", dtype=np.complex128)
    print("Loading R...")
    R = load_bin(RESULTS_DIR / "eq014_right_eigvecs_n7.bin", (d2, d2))
    print("Loading W...")
    W = load_bin(RESULTS_DIR / "eq014_left_eigvecs_n7.bin", (d2, d2))
    apply_biorth_fix(R, W, c_diag, cluster_flat, cluster_offsets)
    print("Biorth applied.")

    # Build L_A (sparse) to cross-check and to enable L_B construction
    print("Rebuilding L_A (sparse, cross-check)...")
    L_A = rebuild_L_A(N, J, gamma_0)
    print(f"L_A nnz: {L_A.nnz}")

    # Initial states
    states = build_initial_states(N, d)

    summary = []
    # Loop over all 5 initial states x 2 bonds
    for state_name, rho0 in states.items():
        print(f"\n##### State: {state_name} #####")
        # P_A exact via eigendecomp
        c_all = W.conj().T @ rho0
        P_A_exact = np.zeros((N, len(TIME_WINDOW)))
        for ti, t in enumerate(TIME_WINDOW):
            rho_t_flat = R @ (c_all * np.exp(vals * t))
            M = rho_t_flat.reshape((d, d))
            margs = partial_trace_all_sites(M, d, N)
            for i in range(N):
                P_A_exact[i, ti] = float(np.real(np.trace(margs[i] @ margs[i].conj().T)))

        for bond in [(0, 1), (3, 4)]:
            V_L = build_V_L(bond, N)
            L_B = (L_A + DELTA_J * V_L).tocsr()
            P_B_exact = evolve_and_sample(L_B, rho0, d, N, TIME_WINDOW, dt_small=0.01)
            alpha_exact = fit_alpha(P_A_exact, P_B_exact, TIME_WINDOW, t_max=T_FIT)
            valid = np.all(np.isfinite(alpha_exact) & (alpha_exact > 0))
            sum_ln_exact = float(np.sum(np.log(alpha_exact))) if valid else float('nan')
            print(f"  Bond {bond}: α = {alpha_exact}")
            print(f"            Σ ln α = {sum_ln_exact:+.4e}")
            summary.append((state_name, bond, alpha_exact.tolist(), sum_ln_exact))

    # Compare to PTF empirical and to my first-order PT (from closure_test.txt)
    print("\n\n=== SUMMARY TABLE: RK4 EXACT Σ ln α_i (this run) ===")
    print(f"{'state':<20} {'bond':<10} Σ ln(α_i) (RK4 exact)")
    for name, bond, _, s in summary:
        print(f"  {name:<20} {str(bond):<10} {s:+.4e}")
    return

    # (legacy single-state path kept for backward compat, but main loop runs all 5)
    pass


if __name__ == "__main__":
    main()
