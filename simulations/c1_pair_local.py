#!/usr/bin/env python3
"""c1_pair_local.py

Approach B: pair-local observable analog of the per-site alpha_i.

Per-pair purity P_{ij}(t) = Tr(rho_{ij}^2) for the 2-site reduced
density matrix rho_{ij} = Tr_{not {i, j}}(rho). Fit alpha_{ij} such
that P_B(i, j, t) ~ P_A(i, j, alpha_{ij} * t). Closure:
Sum_{(i,j)} ln(alpha_{ij}), with c_1 = first derivative in dJ.

Kinematic prediction: single-site alpha_i captures only |Delta N| <= 1
sector blocks of rho_0. Pair-local alpha_{ij} captures |Delta N| <= 2.
Therefore coherence blocks previously invisible (|Delta N| = 2) should
now contribute, while |Delta N| >= 3 remains invisible.

Empirical test:
- (|vac> + |S_2>)/sqrt(2) at N=5: Delta N = 2 coherence. Site-local
  c_1 = -0.29 (all from diagonal-cross); coherence contribution was 0.
  Pair-local should now show a coherence contribution.
- (|vac> + |S_3>)/sqrt(2) at N=5: Delta N = 3. Pair-local still
  invisible. Coherence contribution should stay 0.

N=5, d=32, d^2 = 1024. 3 eigendecompositions ~20s. 10 pairs per state.
~2 minutes total.
"""
from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import eig

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX, ALPHA_BOUNDS,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, density_matrix,
    fit_alpha,
)
from c1_bilinearity_test import dicke_state

RESULTS_DIR = Path(__file__).parent / "results" / "c1_pair_local"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N = 5
DJ_EXTRACT = 0.01


# ---------------------------------------------------------------------------
# Two-site partial trace
# ---------------------------------------------------------------------------
def partial_trace_keep_sites(rho, sites, N):
    """Trace out all sites not in `sites`. Returns 2^k x 2^k where k = |sites|."""
    sites = tuple(sorted(sites))
    shape_2N = [2] * (2 * N)
    out = rho.reshape(shape_2N)
    ket_axes = list(range(N))
    bra_axes = list(range(N, 2 * N))
    for j in range(N - 1, -1, -1):
        if j in sites:
            continue
        a_k = ket_axes[j]
        a_b = bra_axes[j]
        out = np.trace(out, axis1=a_k, axis2=a_b)
        lo, hi = sorted((a_k, a_b))
        for k in range(N):
            if k == j:
                continue
            if ket_axes[k] > hi:
                ket_axes[k] -= 2
            elif ket_axes[k] > lo:
                ket_axes[k] -= 1
            if bra_axes[k] > hi:
                bra_axes[k] -= 2
            elif bra_axes[k] > lo:
                bra_axes[k] -= 1
    k = len(sites)
    d_sub = 2**k
    # After tracing, the remaining axes are the ket and bra axes of the kept
    # sites. With k=2, shape is (2, 2, 2, 2). We need to collapse to (4, 4)
    # with kets as rows and bras as columns. The order of kets in the final
    # shape follows the remaining ket_axes[sites] values; same for bras.
    # For k=2 sites (i, j), after all other traces, remaining axes ordering
    # follows the remaining ket_axes values (ascending i, j) and bra values.
    # Simpler: reorder axes so kets are first, bras are last, then reshape.
    remaining_axes_ket = [ket_axes[i] for i in sites]
    remaining_axes_bra = [bra_axes[i] for i in sites]
    current_order = remaining_axes_ket + remaining_axes_bra
    # Map current_order to new positions
    perm = np.argsort(current_order)  # move these axes to positions 0..2k-1
    out = np.transpose(out, perm)
    return out.reshape(d_sub, d_sub)


def per_pair_purity(rho_traj, N):
    """Return array of shape (T, C(N,2)) of pair purities."""
    T = rho_traj.shape[0]
    pairs = list(combinations(range(N), 2))
    out = np.zeros((T, len(pairs)))
    for t_idx in range(T):
        rho_t = rho_traj[t_idx]
        for p_idx, (i, j) in enumerate(pairs):
            rho_ij = partial_trace_keep_sites(rho_t, (i, j), N)
            out[t_idx, p_idx] = float(np.trace(rho_ij @ rho_ij).real)
    return out, pairs


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


def measure_c1_pair(label, ket, decomps, N, times):
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_0 = density_matrix(ket)
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_A, pairs = per_pair_purity(rho_A, N)
    P_Bp, _ = per_pair_purity(rho_Bp, N)
    P_Bm, _ = per_pair_purity(rho_Bm, N)
    n_pairs = len(pairs)
    alpha_p = np.zeros(n_pairs); alpha_m = np.zeros(n_pairs)
    for k in range(n_pairs):
        a, _ = fit_alpha(times, P_A[:, k], P_Bp[:, k])
        alpha_p[k] = a
        a, _ = fit_alpha(times, P_A[:, k], P_Bm[:, k])
        alpha_m[k] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    c_1 = (closure_p - closure_m) / (2 * DJ_EXTRACT)
    return {
        "label": label, "c_1_pair_local": c_1,
        "closure_plus": closure_p, "closure_minus": closure_m,
        "pairs": [list(p) for p in pairs],
        "alpha_plus": alpha_p.tolist(),
        "alpha_minus": alpha_m.tolist(),
    }


# ---------------------------------------------------------------------------
# Comparison: single-site c_1 computed from the same trajectories
# ---------------------------------------------------------------------------
def per_site_purity_new(rho_traj, N):
    from pi_pair_closure_investigation import partial_trace_keep_site_fast
    T = rho_traj.shape[0]
    P = np.zeros((T, N))
    for t_idx in range(T):
        for i in range(N):
            rho_i = partial_trace_keep_site_fast(rho_traj[t_idx], i, N)
            P[t_idx, i] = float(np.trace(rho_i @ rho_i).real)
    return P


def measure_c1_site(ket, decomps, N, times):
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_0 = density_matrix(ket)
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_A = per_site_purity_new(rho_A, N)
    P_Bp = per_site_purity_new(rho_Bp, N)
    P_Bm = per_site_purity_new(rho_Bm, N)
    alpha_p = np.zeros(N); alpha_m = np.zeros(N)
    for i in range(N):
        a, _ = fit_alpha(times, P_A[:, i], P_Bp[:, i])
        alpha_p[i] = a
        a, _ = fit_alpha(times, P_A[:, i], P_Bm[:, i])
        alpha_m[i] = a
    closure_p = float(np.sum(np.log(alpha_p)))
    closure_m = float(np.sum(np.log(alpha_m)))
    return (closure_p - closure_m) / (2 * DJ_EXTRACT)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print(f"Pair-local c_1 at N = {N}")
    print("=" * 70)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}\n")

    J_A = [J_UNIFORM] * (N - 1)
    J_B_plus = list(J_A); J_B_plus[0] = J_UNIFORM + DJ_EXTRACT
    J_B_minus = list(J_A); J_B_minus[0] = J_UNIFORM - DJ_EXTRACT

    print("  Building Liouvillians...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    decA = eig_and_inv(L_A); del L_A
    L_Bp = build_liouvillian_matrix(build_H_XY(J_B_plus, N), GAMMA_0, N)
    decBp = eig_and_inv(L_Bp); del L_Bp
    L_Bm = build_liouvillian_matrix(build_H_XY(J_B_minus, N), GAMMA_0, N)
    decBm = eig_and_inv(L_Bm); del L_Bm
    print(f"  Built in {time.time()-t0:.1f} s")

    decomps = (decA, decBp, decBm)
    times_arr = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    # Test states
    S = {n: dicke_state(N, n) for n in range(N + 1)}
    vac = vacuum_ket(N)
    psi_1 = single_excitation_mode(N, k=1)

    def norm(v):
        return v / np.linalg.norm(v)

    test_states = [
        # (label, ket, expected_delta_N_coherence)
        ("(|vac>+|psi_1>)/sqrt(2) [PTF]",          norm(vac + psi_1),   1),
        ("(|vac>+|S_1>)/sqrt(2) [W+vac]",          norm(vac + S[1]),    1),
        ("(|vac>+|S_2>)/sqrt(2) [DelN=2, new!]",   norm(vac + S[2]),    2),
        ("(|vac>+|S_3>)/sqrt(2) [DelN=3, still invisible]", norm(vac + S[3]), 3),
        ("(|vac>+|S_4>)/sqrt(2) [DelN=4]",         norm(vac + S[4]),    4),
        ("(|S_1>+|S_2>)/sqrt(2) [DelN=1]",         norm(S[1] + S[2]),   1),
        ("(|S_1>+|S_3>)/sqrt(2) [DelN=2, new!]",   norm(S[1] + S[3]),   2),
        ("(|S_2>+|S_4>)/sqrt(2) [DelN=2, new!]",   norm(S[2] + S[4]),   2),
    ]

    print(f"\n  Per-state c_1 comparison: site-local vs pair-local")
    print(f"  {'label':<50} {'c_1 site':>10} {'c_1 pair':>10} {'|ΔN|':>5}")
    results = []
    for label, ket, dN in test_states:
        t0 = time.time()
        c1_site = measure_c1_site(ket, decomps, N, times_arr)
        r_pair = measure_c1_pair(label, ket, decomps, N, times_arr)
        c1_pair = r_pair["c_1_pair_local"]
        elapsed = time.time() - t0
        print(f"  {label:<50} {c1_site:>+10.4f} {c1_pair:>+10.4f} {dN:>5d}"
              f"   [{elapsed:.1f}s]")
        results.append({
            "label": label, "delta_N_coherence": dN,
            "c_1_site": c1_site, "c_1_pair": c1_pair,
            "alpha_plus": r_pair["alpha_plus"],
            "alpha_minus": r_pair["alpha_minus"],
            "pairs": r_pair["pairs"],
        })

    # Analysis: did the ΔN = 2 coherences contribute differently?
    print(f"\n  Analysis:")
    for r in results:
        dN = r["delta_N_coherence"]
        diff = r["c_1_pair"] - r["c_1_site"]
        print(f"    {r['label']:<50} ΔN={dN}: c_pair - c_site = {diff:+.4f}")

    # Save
    out = {
        "N": N, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "defect_bond": [0, 1],
        "results": results,
    }
    path = RESULTS_DIR / "pair_local.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
