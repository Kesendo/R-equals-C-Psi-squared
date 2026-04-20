#!/usr/bin/env python3
"""c1_bond_scan.py

Approach A: does the Dicke-state c_1 signature depend on which bond
the perturbation sits on?

Setup: N=5, Dicke states |S_n> for n=0..5, four bonds tested b in
{0, 1, 2, 3} where bond b is (b, b+1). Bonds come in Pi-pairs under
site reflection: (0,1) <-> (3,4) and (1,2) <-> (2,3).

Test:
- If c_1(|S_n>, bond b) is the same across all b, K_diag is bond-invariant
  (state-intrinsic signature).
- If it varies across b, K_diag is bond-local (depends on perturbation
  position in the chain).
- By bond-mirror symmetry, bonds (0,1) and (3,4) should give identical
  c_1 for each Dicke state (Dicke states are bit-flip invariant;
  combined with chain reflection they are site-reflection invariant
  up to sign). Similarly (1,2) and (2,3).

12 eigendecompositions at d^2 = 1024, ~30 seconds total.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    per_site_purity, fit_alpha,
)
from c1_bilinearity_test import dicke_state

RESULTS_DIR = Path(__file__).parent / "results" / "c1_bond_scan"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 5
DJ_EXTRACT = 0.01


def eig_and_inv(L):
    eigvals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return eigvals, V_R, V_Linv


def propagate(eigvals, V_R, V_Linv, rho_0, times):
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_Linv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def measure_c1(rho_0, decA, decBp, decBm, N, times):
    rho_A = propagate(*decA, rho_0, times)
    rho_Bp = propagate(*decBp, rho_0, times)
    rho_Bm = propagate(*decBm, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    for i in range(N):
        a, _ = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        alpha_p[i] = a
        a, _ = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        alpha_m[i] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    return (closure_p - closure_m) / (2 * DJ_EXTRACT), alpha_p.tolist(), alpha_m.tolist()


def main():
    print("=" * 70)
    print(f"Bond-scan of Dicke c_1 at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}\n")

    J_A = [J_UNIFORM] * (N - 1)
    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    # L_A is the same regardless of which bond is perturbed
    print(f"  Building L_A and eigendecomposing...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    decA = eig_and_inv(L_A)
    del L_A
    print(f"    done in {time.time()-t0:.1f} s")

    # Precompute Dicke states
    S = {n: dicke_state(N, n) for n in range(N + 1)}
    state_densities = {n: density_matrix(S[n]) for n in range(N + 1)}

    # For each bond b in 0..N-2, build L_B+ and L_B- and measure c_1
    # for each Dicke state
    bond_results = {}
    for b in range(N - 1):
        print(f"\n  Bond b = {b} ({b}, {b+1}):")
        J_Bp = list(J_A); J_Bp[b] = J_UNIFORM + DJ_EXTRACT
        J_Bm = list(J_A); J_Bm[b] = J_UNIFORM - DJ_EXTRACT
        t0 = time.time()
        L_Bp = build_liouvillian_matrix(build_H_XY(J_Bp, N), GAMMA_0, N)
        decBp = eig_and_inv(L_Bp); del L_Bp
        L_Bm = build_liouvillian_matrix(build_H_XY(J_Bm, N), GAMMA_0, N)
        decBm = eig_and_inv(L_Bm); del L_Bm
        print(f"    decomps in {time.time()-t0:.1f} s")

        c1_by_n = {}
        for n in range(N + 1):
            c1, _ap, _am = measure_c1(state_densities[n], decA, decBp, decBm,
                                        N, times_arr)
            c1_by_n[n] = c1
        bond_results[b] = c1_by_n

        # Print
        labels = "  ".join(f"|S_{n}>" for n in range(N + 1))
        values = "  ".join(f"{c1_by_n[n]:+.4f}" for n in range(N + 1))
        print(f"    c_1 per Dicke state:")
        print(f"    state:     {labels}")
        print(f"    c_1:       {values}")

    # Bond-mirror check. At N=5: bonds (0,1) <-> (3,4) are mirror pairs
    # under i <-> N-1-i, and (1,2) <-> (2,3) likewise.
    print(f"\n  Bond-mirror symmetry check (bond b <-> N-2-b):")
    for b in range((N - 1) // 2 + 1):
        b_mirror = (N - 2) - b
        if b > b_mirror:
            continue
        diffs = []
        for n in range(N + 1):
            diff = abs(bond_results[b][n] - bond_results[b_mirror][n])
            diffs.append(diff)
        if b == b_mirror:
            print(f"    bond {b} is self-mirror (center)")
        else:
            print(f"    bond {b} <-> bond {b_mirror}: "
                  f"max |c_1 diff over Dicke states| = {max(diffs):.2e}")

    # Variation check: does c_1 vary with bond for each Dicke state?
    print(f"\n  c_1 variation across bonds per Dicke state:")
    print(f"  {'state':>6}  " + "  ".join(f"b={b}" for b in range(N - 1))
          + f"  {'range':>10}  {'rel var':>10}")
    for n in range(N + 1):
        vals = [bond_results[b][n] for b in range(N - 1)]
        rng = max(vals) - min(vals)
        mean = sum(vals) / len(vals)
        rel_var = rng / abs(mean) if abs(mean) > 1e-10 else float('inf')
        vals_str = "  ".join(f"{v:+.4f}" for v in vals)
        print(f"  {'|S_'+str(n)+'>':>6}  {vals_str}  {rng:>10.4f}  "
              f"{rel_var:>10.2f}")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "bonds_tested": list(range(N - 1)),
        "bond_results": {
            f"bond_{b}_{b+1}": {f"S_{n}": bond_results[b][n]
                                 for n in range(N + 1)}
            for b in range(N - 1)
        },
    }
    path = RESULTS_DIR / "bond_scan.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
