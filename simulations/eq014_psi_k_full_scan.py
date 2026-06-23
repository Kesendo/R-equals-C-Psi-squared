#!/usr/bin/env python3
"""EQ-014 surviving sub-question: full ψ_k scan for k = 1..N at N=7.

Existing data (eq014_first_order_from_rk4.py output, 2026-04-19):
  ψ_1 (k=1 bonding)  Σ f_i ≈ +0.97
  ψ_2 (k=2 bonding)  Σ f_i ≈ +0.05
  ψ_3 (k=3 bonding)  Σ f_i ≈ +0.36
  |+⟩^7              Σ f_i ≈ +1.29

Three out of four are O(1); ψ_2 is an order of magnitude smaller. The
candidate explanations were:
  (a) Fourier-content selection rule on V_L matrix element
  (b) Node structure of ψ_k at the defect bond
  (c) Excitation-sector weight (slow subspace dominance)

This script extends to k = 4, 5, 6, 7 to look for the full k-dependence
pattern. For each ψ_k bonding state at N=7, bond (0,1), δJ = 0.01:
  - run RK4 evolution for L_A (unperturbed) and L_B = L_A + δJ V_L
  - fit α_i per site
  - extract Σ f_i = Σ ln(α_i) / δJ

Compare against the Fourier-overlap predictor:
  M_k = sin(π k / (N+1)) · sin(2π k / (N+1))
which vanishes at k = (N+1)/2 (k=4 for N=7) — the "node-at-defect" case.

If Σ f_i scales with |M_k|, that's selection rule (a) confirmed.
If not, the pattern is in the bilinear purity expansion's eigenvector mixing.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from eq014_validate_groundtruth import rebuild_L_A, evolve_and_sample
from eq014_step4567_closure import build_V_L, fit_alpha


N = 7
D = 128
J = 1.0
GAMMA = 0.05
DELTA_J = 0.01
TIMES = np.arange(401) * 0.2
T_FIT = 20.0
BOND = (0, 1)


def build_psi_k_bonding(k, N, d):
    """ψ_k = single-excitation momentum k state, then (vac + ψ_k)/√2."""
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = np.zeros(d, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        idx = 1 << (N - 1 - i)
        psi[idx] = amp
    phi = (vac + psi) / np.sqrt(2)
    rho = np.outer(phi, phi.conj())
    return rho.flatten()


def fourier_overlap_predictor(k, N):
    """V_L matrix element on bond (0,1) for ψ_k: depends on
    sin(π k / (N+1)) · sin(2π k / (N+1)) (the bond-amplitude product)."""
    return float(
        np.sin(np.pi * k / (N + 1)) * np.sin(2 * np.pi * k / (N + 1))
    )


def main():
    print(f"=== EQ-014 ψ_k full scan, N={N}, bond={BOND}, δJ={DELTA_J} ===")
    print(f"  Existing data (2026-04-19, eq014_first_order_from_rk4.py):")
    print(f"    ψ_1: +0.97   ψ_2: +0.05   ψ_3: +0.36   |+⟩⁷: +1.29")
    print()

    print("Building L_A...")
    L_A = rebuild_L_A(N, J, GAMMA).tocsr()
    V_L = build_V_L(BOND, N)
    L_B = (L_A + DELTA_J * V_L).tocsr()
    print(f"  L_A nnz: {L_A.nnz}")
    print()

    # We'll re-evolve P_A from scratch for each ψ_k since each state has
    # a different P_A. Reuse-pattern: cache P_A per state across runs.
    # But for THIS script we only need one δJ, so just do A and B per state.

    print(f"  {'state':<14s}  {'Σ f_i':>10s}  {'M_k Fourier':>12s}  "
          f"{'|M_k|':>8s}  {'Σ ln α (δJ)':>13s}  {'time':>7s}")
    print('-' * 80)

    results = []
    for k in range(1, N + 1):
        rho0 = build_psi_k_bonding(k, N, D)
        m_k = fourier_overlap_predictor(k, N)

        t0 = time.time()
        P_A = evolve_and_sample(L_A, rho0, D, N, TIMES, dt_small=0.01)
        P_B = evolve_and_sample(L_B, rho0, D, N, TIMES, dt_small=0.01)
        alpha = fit_alpha(P_A, P_B, TIMES, t_max=T_FIT)
        valid = np.all(np.isfinite(alpha) & (alpha > 0))
        sum_ln = float(np.sum(np.log(alpha))) if valid else float('nan')
        f_coeff = sum_ln / DELTA_J if np.isfinite(sum_ln) else float('nan')
        elapsed = time.time() - t0

        print(f"  ψ_{k}_bonding    {f_coeff:>+10.4f}  {m_k:>+12.4f}  "
              f"{abs(m_k):>8.4f}  {sum_ln:>+13.6e}  {elapsed:>5.1f}s")

        results.append({
            'k': k, 'sum_f': f_coeff, 'M_k': m_k, 'sum_ln': sum_ln,
            'alpha': alpha.tolist(), 'time': elapsed,
        })

    print()
    print("Summary table:")
    print(f"  {'k':>2s}  {'Σ f_i':>9s}  {'M_k':>9s}  {'|M_k|':>7s}  {'ratio f_i/M_k':>14s}")
    for r in results:
        ratio = r['sum_f'] / r['M_k'] if abs(r['M_k']) > 1e-9 else float('nan')
        print(f"  {r['k']:>2d}  {r['sum_f']:>+9.4f}  {r['M_k']:>+9.4f}  "
              f"{abs(r['M_k']):>7.4f}  {ratio:>+14.3f}")

    # Pearson correlation
    sum_fs = np.array([r['sum_f'] for r in results])
    abs_mks = np.array([abs(r['M_k']) for r in results])
    if np.all(np.isfinite(sum_fs)) and abs_mks.std() > 0 and sum_fs.std() > 0:
        rho_abs = float(np.corrcoef(sum_fs, abs_mks)[0, 1])
        print()
        print(f"  Pearson(Σ f_i, |M_k|) = {rho_abs:+.4f}")
        rho_signed = float(np.corrcoef(sum_fs, np.array([r['M_k'] for r in results]))[0, 1])
        print(f"  Pearson(Σ f_i, M_k signed) = {rho_signed:+.4f}")


if __name__ == "__main__":
    main()
