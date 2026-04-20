#!/usr/bin/env python3
"""c1_veffect_scaling.py

Scaling of the closure-breaking coefficient c_1 with chain length N, and
comparison against the V-Effect V(N) = 1 + cos(π/N).

From ClaudeTasks/TASK_C1_VEFFECT_SCALING.md.

Hypothesis: c_1 at the endpoint bond (0,1) might scale with V(N) because
endpoints are the cavity windows where gamma_0 couples to the chain.

Procedure per N in {3, 4, 5, 6, 7}:
  * Build XY chain L_A (uniform J=1) and L_B (J_bond(0,1) = 1 + dJ).
  * Eigendecompose L_A ONCE; reuse for propagation.
  * Propagate rho(t) under L_A and L_B for initial state (|vac> + |psi_1>)/sqrt(2).
  * Fit per-site alpha_i and compute Sum_i ln(alpha_i).
  * Extract c_1 = [closure(+dJ) - closure(-dJ)] / (2 * dJ) (symmetric difference
    cancels the c_2 contribution to leading order).

Also compute full c_1 vector (one per bond) at N in {3, 4, 5, 6} for spatial profile.
Skip full vector at N=7 (too expensive).

Also repeat the N-scaling for initial state psi_2+vac (state universality test).

Rules:
  * XY Hamiltonian H = (J/2)(XX+YY).
  * Numbers from script output only.
  * UTF-8 stdout.
  * No em-dashes.
  * For N=7: dense eig but use eigendecomp-based propagation, not expm.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX, ALPHA_BOUNDS,
    X, Y, Z, I2, site_op, build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode,
    density_matrix, per_site_purity, fit_alpha,
)

RESULTS_DIR = Path(__file__).parent / "results" / "c1_veffect_scaling"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01   # probe amplitude for c_1 extraction
DJ_PROBE = [-DJ_EXTRACT, +DJ_EXTRACT]


def v_effect(N):
    return 1.0 + np.cos(np.pi / N)


def bonding_plus_vacuum(N, k=1):
    v = vacuum_ket(N)
    psi = single_excitation_mode(N, k=k)
    return (v + psi) / np.sqrt(2.0)


def propagate_via_eig(eigvals, V_R, V_L_inv, rho_0, times):
    """Propagate rho(t) given precomputed eigendecomposition of L_A."""
    d = rho_0.shape[0]
    rho0_vec = rho_0.flatten(order='F')
    c0 = V_L_inv @ rho0_vec
    out = np.empty((len(times), d, d), dtype=complex)
    for k, t in enumerate(times):
        rho_vec_t = V_R @ (np.exp(eigvals * t) * c0)
        out[k] = rho_vec_t.reshape(d, d, order='F')
    return out


def build_and_decompose(N, J_list):
    """Build L, eigendecompose. Returns (eigvals, V_R, V_L_inv)."""
    H = build_H_XY(J_list, N)
    L = build_liouvillian_matrix(H, GAMMA_0, N)
    eigvals, V_R = eig(L)
    V_L_inv = np.linalg.inv(V_R)
    return eigvals, V_R, V_L_inv


def compute_c1_for_bond(N, bond_idx, state, times, dJ_list):
    """Compute c_1 for a specified bond. Reuses L_A eigendecomp implicitly
    by building L_B fresh each time (cheaper than re-eig L_A).

    Returns dict with closures per dJ and extracted c_1."""
    J_A = [J_UNIFORM] * (N - 1)
    t0 = time.time()
    ev_A, VR_A, VLinv_A = build_and_decompose(N, J_A)
    t_decomp_A = time.time() - t0
    rho_0 = density_matrix(state)
    rho_A = propagate_via_eig(ev_A, VR_A, VLinv_A, rho_0, times)
    P_A = per_site_purity(rho_A, N)

    entries = []
    for dJ in dJ_list:
        J_B = list(J_A); J_B[bond_idx] += dJ
        t1 = time.time()
        ev_B, VR_B, VLinv_B = build_and_decompose(N, J_B)
        t_decomp_B = time.time() - t1
        rho_B = propagate_via_eig(ev_B, VR_B, VLinv_B, rho_0, times)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
        closure = float(np.sum(np.log(alpha)))
        entries.append({"dJ": float(dJ), "closure": closure,
                        "alpha": alpha.tolist(),
                        "decomp_time_s": t_decomp_B})

    # Extract c_1 from symmetric difference
    plus = next(e for e in entries if e["dJ"] > 0)
    minus = next(e for e in entries if e["dJ"] < 0)
    c1 = (plus["closure"] - minus["closure"]) / (plus["dJ"] - minus["dJ"])
    c2 = (plus["closure"] + minus["closure"]) / (plus["dJ"]**2 + minus["dJ"]**2)

    return {"entries": entries, "c_1": float(c1), "c_2_approx": float(c2),
            "L_A_decomp_time_s": t_decomp_A}


def main():
    print("="*70)
    print("c_1 vs N scaling + V-Effect comparison")
    print("="*70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}")
    print(f"  dJ for c_1 extraction = +/- {DJ_EXTRACT}")
    print(f"  Initial state: (|vac> + |psi_k>) / sqrt(2)")
    print(f"  Results -> {RESULTS_DIR}")

    all_results = {}
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    # -----------------------------------------------------------
    # Step 1 + 2: c_1(endpoint) vs N, state = psi_1 + vac
    # -----------------------------------------------------------
    print(f"\n{'-'*70}\nStep 1+2: c_1 at bond (0,1), psi_1+vac, N=3..7\n{'-'*70}")
    print(f"  {'N':>3} {'V(N)':>8} {'c_1 endpoint':>15} {'c_2 approx':>12}"
          f" {'closure@+dJ':>14} {'closure@-dJ':>14} {'decomp_A':>10} {'decomp_B':>10}")
    step12 = {}
    for N in [3, 4, 5, 6, 7]:
        state = bonding_plus_vacuum(N, k=1)
        t0 = time.time()
        r = compute_c1_for_bond(N, bond_idx=0, state=state,
                                times=times, dJ_list=DJ_PROBE)
        total = time.time() - t0
        plus_cl = next(e["closure"] for e in r["entries"] if e["dJ"] > 0)
        minus_cl = next(e["closure"] for e in r["entries"] if e["dJ"] < 0)
        decomp_B = np.mean([e["decomp_time_s"] for e in r["entries"]])
        print(f"  {N:>3d} {v_effect(N):>8.4f} {r['c_1']:>+15.5f} "
              f"{r['c_2_approx']:>+12.4f} "
              f"{plus_cl:>+14.6f} {minus_cl:>+14.6f} "
              f"{r['L_A_decomp_time_s']:>10.2f} {decomp_B:>10.2f}")
        step12[N] = {"V_N": v_effect(N), "c_1": r["c_1"],
                     "c_2_approx": r["c_2_approx"],
                     "closure_plus": plus_cl, "closure_minus": minus_cl,
                     "alpha_plus": r["entries"][1]["alpha"],
                     "alpha_minus": r["entries"][0]["alpha"],
                     "total_time_s": total}
    all_results["step12_psi1_endpoint_scan"] = step12

    # Scaling fits
    Ns = sorted(step12.keys())
    c1s = np.array([step12[N]["c_1"] for N in Ns])
    Vs = np.array([v_effect(N) for N in Ns])
    # Power fit c_1 ~ A * N^p
    log_c1 = np.log(np.abs(c1s))
    log_N = np.log(Ns)
    p_N, logA_N = np.polyfit(log_N, log_c1, 1)
    # Power fit c_1 ~ B * V(N)^q
    log_V = np.log(Vs)
    p_V, logB_V = np.polyfit(log_V, log_c1, 1)
    # Test c_1 ~ A * N^2 * V(N) ?
    ansatz_N2V = np.array([N**2 * v_effect(N) for N in Ns])
    log_N2V = np.log(ansatz_N2V)
    p_N2V, logA_N2V = np.polyfit(log_N2V, log_c1, 1)

    print(f"\n  Power-law fits (log-log):")
    print(f"    c_1 ~ A * N^p:           p = {p_N:.3f}, "
          f"A = {np.exp(logA_N):.4f}, residuals:")
    for N, c in zip(Ns, c1s):
        pred = np.exp(logA_N) * N**p_N
        print(f"      N={N}: actual={c:+.4f}, pred={pred:+.4f}, "
              f"ratio={c/pred:.3f}")
    print(f"    c_1 ~ B * V(N)^q:        q = {p_V:.3f}, "
          f"B = {np.exp(logB_V):.4f}, residuals:")
    for N, c in zip(Ns, c1s):
        pred = np.exp(logB_V) * v_effect(N)**p_V
        print(f"      N={N}: actual={c:+.4f}, pred={pred:+.4f}, "
              f"ratio={c/pred:.3f}")
    print(f"    c_1 ~ A * (N^2 V(N))^r:  r = {p_N2V:.3f}, "
          f"A = {np.exp(logA_N2V):.4f}, residuals:")
    for N, c in zip(Ns, c1s):
        pred = np.exp(logA_N2V) * (N**2 * v_effect(N))**p_N2V
        print(f"      N={N}: actual={c:+.4f}, pred={pred:+.4f}, "
              f"ratio={c/pred:.3f}")
    all_results["step12_fits"] = {
        "c_1_vs_N_power": {"p": float(p_N), "A": float(np.exp(logA_N))},
        "c_1_vs_V_power": {"q": float(p_V), "B": float(np.exp(logB_V))},
        "c_1_vs_N2V_power": {"r": float(p_N2V), "A": float(np.exp(logA_N2V))},
    }

    # -----------------------------------------------------------
    # Step 3: full c_1 vector at N = 3, 4, 5, 6
    # -----------------------------------------------------------
    print(f"\n{'-'*70}\nStep 3: full c_1 vector over all bonds, psi_1+vac, N=3..6\n{'-'*70}")
    step3 = {}
    for N in [3, 4, 5, 6]:
        state = bonding_plus_vacuum(N, k=1)
        c1_vec = np.zeros(N - 1)
        print(f"\n  N = {N}, bonds = {list(range(N-1))}")
        for b in range(N - 1):
            r = compute_c1_for_bond(N, bond_idx=b, state=state,
                                    times=times, dJ_list=DJ_PROBE)
            c1_vec[b] = r["c_1"]
            print(f"    bond ({b},{b+1}): c_1 = {r['c_1']:+.5f}")
        # Mirror symmetry check
        mirror_diff = np.linalg.norm(c1_vec - c1_vec[::-1])
        endpoint_c1 = c1_vec[0]
        interior_c1 = c1_vec[(N-1)//2]
        ratio = endpoint_c1 / interior_c1 if abs(interior_c1) > 1e-10 else float('inf')
        print(f"    Full c_1 vector: {c1_vec}")
        print(f"    Mirror residual (||c_1 - reversed||): {mirror_diff:.2e}")
        print(f"    endpoint/interior ratio: {ratio:+.3f}")
        step3[N] = {"c_1_vector": c1_vec.tolist(),
                    "mirror_residual": float(mirror_diff),
                    "endpoint_over_interior": float(ratio)}
    all_results["step3_full_vector_psi1"] = step3

    # -----------------------------------------------------------
    # Step 4: N-scaling for psi_2 + vac (endpoint only)
    # -----------------------------------------------------------
    print(f"\n{'-'*70}\nStep 4: c_1 at bond (0,1), psi_2+vac, N=3..6\n{'-'*70}")
    step4 = {}
    print(f"  {'N':>3} {'V(N)':>8} {'c_1 endpoint':>15}")
    for N in [3, 4, 5, 6]:
        # psi_2 requires N >= 2
        state = (vacuum_ket(N) + single_excitation_mode(N, k=2)) / np.sqrt(2.0)
        r = compute_c1_for_bond(N, bond_idx=0, state=state,
                                times=times, dJ_list=DJ_PROBE)
        print(f"  {N:>3d} {v_effect(N):>8.4f} {r['c_1']:>+15.5f}")
        step4[N] = {"V_N": v_effect(N), "c_1": r["c_1"]}
    all_results["step4_psi2_endpoint_scan"] = step4

    # Save
    out = RESULTS_DIR / "c1_vs_N.json"
    with open(out, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved: {out}")

    # Summary TXT
    summary_path = RESULTS_DIR / "c1_scaling_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"c_1 scaling with N and V-Effect comparison\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Setup: XY chain, gamma_0 = {GAMMA_0}, J = {J_UNIFORM}\n")
        f.write(f"dJ probe = +/- {DJ_EXTRACT}\n\n")
        f.write(f"Step 1+2: psi_1+vac, endpoint bond (0,1)\n")
        f.write(f"{'N':>3} {'V(N)':>8} {'c_1':>12}\n")
        for N in Ns:
            f.write(f"{N:>3} {v_effect(N):>8.4f} "
                    f"{step12[N]['c_1']:>+12.5f}\n")
        f.write(f"\nFits:\n")
        f.write(f"  c_1 ~ N^p:       p = {p_N:.3f}\n")
        f.write(f"  c_1 ~ V(N)^q:    q = {p_V:.3f}\n")
        f.write(f"  c_1 ~ (N^2 V)^r: r = {p_N2V:.3f}\n\n")
        f.write(f"Step 3: full c_1 vector, psi_1+vac\n")
        for N in [3, 4, 5, 6]:
            f.write(f"  N={N}: c_1 vector = {step3[N]['c_1_vector']}\n")
        f.write(f"\nStep 4: psi_2+vac, endpoint\n")
        f.write(f"{'N':>3} {'c_1':>12}\n")
        for N in [3, 4, 5, 6]:
            f.write(f"{N:>3} {step4[N]['c_1']:>+12.5f}\n")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
