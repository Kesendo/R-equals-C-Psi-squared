#!/usr/bin/env python3
"""EQ-014 chiral mirror law: verify at N ∈ {5, 8, 9}.

Established at N=7 (commit referenced in EQ-014 closure):
  Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k})    (chiral mirror law)

Mechanism: K_1 = ⊗_i (-1)^i is the chiral / sublattice symmetry of the open
XY chain. Sine-mode eigenvectors satisfy K_1 ψ_k = (−1)^? ψ_{N+1−k} (up to a
staggered phase). Per-site purity P_i = |ψ(i)|² is K_1-invariant under this
mapping, so α_i and Σ f_i are equal for K_1-paired modes.

This script extends to N=5 (3 mirror pairs: (1,5), (2,4), 3 fixed) and
N=8 (4 mirror pairs (1,8), (2,7), (3,6), (4,5), no fixed point), and N=9
(4 mirror pairs + k=5 fixed point).

Method: standalone RK4 evolution + α-fit, no biorth/dense decomposition
needed. Self-contained — does not depend on eq014 step files.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# Pauli matrices
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)


def site_op_sparse(op, k, N):
    full = sps.eye(1, dtype=complex, format='csr')
    I2 = sps.eye(2, dtype=complex, format='csr')
    for j in range(N):
        full = sps.kron(full, sps.csr_matrix(op) if j == k else I2,
                          format='csr')
    return full


def build_xy_chain_H(N, J):
    """H_A = Σ_i (J/2)·(X_i X_{i+1} + Y_i Y_{i+1}), big-endian site indexing."""
    d = 2 ** N
    H = sps.csr_matrix((d, d), dtype=complex)
    for i in range(N - 1):
        H = H + (J / 2) * (
            site_op_sparse(X, i, N) @ site_op_sparse(X, i + 1, N)
            + site_op_sparse(Y, i, N) @ site_op_sparse(Y, i + 1, N)
        )
    return H


def build_L(H, N, gamma):
    """L = -i(H ⊗ I − I ⊗ H^T) − Σ_i 2γ · (XOR-popcount of basis indices)."""
    d = 2 ** N
    dd = d * d
    Id = sps.eye(d, dtype=complex, format='csr')
    L_h = -1j * (sps.kron(H, Id, format='csr') - sps.kron(Id, H.T, format='csr'))
    a = np.arange(dd) // d
    b = np.arange(dd) % d
    h_vec = np.array([bin(int(av ^ bv)).count('1') for av, bv in zip(a, b)])
    L_d = sps.diags(-2.0 * gamma * h_vec.astype(complex), format='csr')
    return (L_h + L_d).tocsr()


def build_V_L_bond(bond, N):
    """δJ-coefficient of L: V_L = -i [bond_H, .] for bond (i, i+1) of XX+YY/2."""
    d = 2 ** N
    b0, b1 = bond
    Hb = 0.5 * (site_op_sparse(X, b0, N) @ site_op_sparse(X, b1, N)
                + site_op_sparse(Y, b0, N) @ site_op_sparse(Y, b1, N))
    Id = sps.eye(d, dtype=complex, format='csr')
    return (-1j * (sps.kron(Hb, Id, format='csr')
                   - sps.kron(Id, Hb.T, format='csr'))).tocsr()


def sine_mode_state(N, k, d):
    """Single-excitation sine mode psi_k as 2^N-dim ket; site i ↔ bit (N-1-i)."""
    psi = np.zeros(d, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        idx = 1 << (N - 1 - i)
        psi[idx] = amp
    return psi


def bonding_mode_density(N, k, d):
    """ρ_0 = |φ⟩⟨φ| for |φ⟩ = (|vac⟩ + |ψ_k⟩)/√2, returned flat."""
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = sine_mode_state(N, k, d)
    phi = (vac + psi) / np.sqrt(2)
    rho = np.outer(phi, phi.conj())
    return rho.flatten()


def rk4_step(L, rho, dt):
    k1 = L @ rho
    k2 = L @ (rho + 0.5 * dt * k1)
    k3 = L @ (rho + 0.5 * dt * k2)
    k4 = L @ (rho + dt * k3)
    return rho + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def partial_trace_all_sites(M, d, N):
    """Per-site marginals as (N, 2, 2) array."""
    out = np.zeros((N, 2, 2), dtype=complex)
    shape = [2] * (2 * N)
    T = M.reshape(shape)
    letters = "abcdefghijklmnop"
    for i in range(N):
        row_labels = list(letters[:N])
        col_labels = list(letters[N:2 * N])
        for j in range(N):
            if j != i:
                col_labels[j] = row_labels[j]
        in_spec = "".join(row_labels) + "".join(col_labels)
        out_spec = row_labels[i] + col_labels[i]
        out[i] = np.einsum(f"{in_spec}->{out_spec}", T)
    return out


def evolve_and_sample(L, rho0_flat, d, N, sample_times, dt_small=0.01):
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


def fit_alpha(P_A, P_B, times, t_max=None, zero_tol=1e-8,
              alpha_bounds=(0.1, 10.0)):
    if t_max is None:
        t_max = times[-1]
    mask = times <= t_max
    N, T = P_A.shape
    alpha = np.zeros(N)
    for i in range(N):
        pa_range = np.max(P_A[i]) - np.min(P_A[i])
        pab_range = np.max(np.abs(P_B[i] - P_A[i]))
        if pa_range < zero_tol and pab_range < zero_tol:
            alpha[i] = 1.0
            continue
        interp = interp1d(times, P_A[i], kind='cubic', bounds_error=False,
                          fill_value=(float(P_A[i, 0]), float(P_A[i, -1])))
        t_eval = times[mask]
        b = P_B[i][mask]

        def mse(a):
            d = interp(a * t_eval) - b
            return float(np.mean(d * d))

        try:
            res = minimize_scalar(mse, bounds=alpha_bounds, method='bounded',
                                    options={'xatol': 1e-6})
            alpha[i] = float(res.x)
        except Exception:
            alpha[i] = np.nan
    return alpha


def fourier_overlap(N, k):
    return float(np.sin(np.pi * k / (N + 1)) * np.sin(2 * np.pi * k / (N + 1)))


def run_N(N, J=1.0, gamma=0.05, delta_J=0.01, t_total=None, t_fit=20.0,
            n_samples=None, bond=(0, 1), dt_small=0.02):
    # Default: t_total=30 (covers fit window t≤20) for N≥8 to keep runtime
    # manageable; t_total=80 at smaller N to match PTF reference exactly.
    if t_total is None:
        t_total = 30.0 if N >= 8 else 80.0
    if n_samples is None:
        n_samples = 151 if N >= 8 else 401
    print()
    print(f"  Building L_A, V_L for N={N}...")
    d = 2 ** N
    H_A = build_xy_chain_H(N, J)
    L_A = build_L(H_A, N, gamma)
    V_L = build_V_L_bond(bond, N)
    L_B = (L_A + delta_J * V_L).tocsr()
    print(f"  L_A shape: {L_A.shape}, nnz: {L_A.nnz}")

    times = np.linspace(0, t_total, n_samples)

    print(f"  {'k':>3}  {'k_mirror':>9}  {'Σ f_i':>10}  {'M_k':>10}  "
          f"{'|M_k|':>8}  {'time':>6}")
    print("  " + "-" * 55)

    results = []
    for k in range(1, N + 1):
        rho0 = bonding_mode_density(N, k, d)
        m_k = fourier_overlap(N, k)
        t0 = time.time()
        P_A = evolve_and_sample(L_A, rho0, d, N, times, dt_small=dt_small)
        P_B = evolve_and_sample(L_B, rho0, d, N, times, dt_small=dt_small)
        alpha = fit_alpha(P_A, P_B, times, t_max=t_fit)
        valid = np.all(np.isfinite(alpha) & (alpha > 0))
        sum_ln = float(np.sum(np.log(alpha))) if valid else float('nan')
        f_coeff = sum_ln / delta_J if np.isfinite(sum_ln) else float('nan')
        elapsed = time.time() - t0
        k_mirror = N + 1 - k
        results.append({
            "N": N, "k": k, "k_mirror": k_mirror,
            "sum_f": f_coeff, "M_k": m_k,
            "elapsed": elapsed,
        })
        print(f"  {k:>3d}  {k_mirror:>9d}  {f_coeff:>+10.4f}  {m_k:>+10.4f}  "
              f"{abs(m_k):>8.4f}  {elapsed:>5.1f}s")

    # Verify chiral mirror law
    print()
    print(f"  Chiral mirror law verification: Σ f_i(k) =? Σ f_i(N+1−k)")
    for k in range(1, (N + 1) // 2 + 1):
        k_mirror = N + 1 - k
        if k > k_mirror:
            continue
        f_k = next(r["sum_f"] for r in results if r["k"] == k)
        f_m = next(r["sum_f"] for r in results if r["k"] == k_mirror)
        diff = abs(f_k - f_m)
        verdict = "OK" if diff < 1e-3 else f"MIRROR-BREAK ({diff:.2e})"
        if k == k_mirror:
            print(f"    k={k}: chiral fixed point (k = N+1-k); Σ f_i = {f_k:+.4f}")
        else:
            print(f"    k={k} ↔ k={k_mirror}: |Σ f_i diff| = {diff:.2e}  {verdict}")
    return results


def main():
    print("=" * 80)
    print("EQ-014 chiral mirror law: extension to N ∈ {5, 7-prior, 8, 9}")
    print("=" * 80)
    print()
    print("Hypothesis: Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k}) at every N (chiral K_1).")
    print("Verified at N=7: peak at chiral fixed-point k=4 (ψ_4 = K_1-eigenstate).")
    print()

    Ns = [int(x) for x in sys.argv[1:] if x.isdigit()] or [5]
    all_results = {}
    for N in Ns:
        print(f"##### N = {N} #####")
        all_results[N] = run_N(N)
        print()


if __name__ == "__main__":
    main()
