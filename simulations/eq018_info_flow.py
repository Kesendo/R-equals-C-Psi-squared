#!/usr/bin/env python3
"""eq018_info_flow.py

Step 2 from ORTHOGONALITY_SELECTION_FAMILY.md §6.2 plan: measure the
information-flow landscape between sites under bond perturbation.

Observable:
    C_ij(t) = Tr(rho_i(t) · rho_j(t))    site-to-site cross-correlation
    dC_ij/dJ = (C_ij^(B+) - C_ij^(B-)) / (2 dJ)    flow response

Interpretation:
  - C_ij(t) measures how similar the reduced states at sites i, j are
    at time t. For i=j, it reduces to per-site purity Tr(rho_i^2).
    For i!=j, it captures how information about one site's state is
    correlated with another's.
  - dC_ij/dJ is the dynamical response of this correlation to a
    bond-b perturbation. It shows how the perturbation's effect
    propagates across the chain.

Key things to extract:
  1. Lieb-Robinson front: dC_ij(t) is zero for |i-j| > v*t (v ~ 2J for XY).
  2. Pi-pair structure: dC_ij(t) should show structure paired under
     spatial reflection (i,j) <-> (N-1-i, N-1-j) for uniform chain.
  3. Flow range as a function of N: does the "reach" of perturbation
     scale linearly, saturate, etc.?

Initial state: PTF bonding (|vac> + |psi_1>)/sqrt(2) for clean signal.
Bond perturbation at b=0. Scan N in {4, 5, 6}.

Rules: UTF-8 stdout, no em-dashes, hyphens only.
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
    partial_trace_keep_site_fast,
)
from eq018_kernel_extract import build_decomps, propagate

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_info_flow"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01


def compute_C_ij_trajectory(rho_traj, N):
    """Compute C_ij(t) = Tr(rho_i(t) * rho_j(t)) for all t, i, j.

    rho_traj: (T, d, d)
    Returns: (T, N, N) array (real-valued)
    """
    T = rho_traj.shape[0]
    out = np.zeros((T, N, N))
    for t_idx in range(T):
        # Compute all rho_i at this time
        rhos_i = np.empty((N, 2, 2), dtype=complex)
        for i in range(N):
            rhos_i[i] = partial_trace_keep_site_fast(rho_traj[t_idx], i, N)
        # C_ij = Tr(rhos_i[i] @ rhos_i[j])
        for i in range(N):
            for j in range(N):
                out[t_idx, i, j] = np.trace(rhos_i[i] @ rhos_i[j]).real
    return out


def measure_dC_ij_dJ(rho_0, decomps, times, N, dj=DJ_EXTRACT):
    """Perturbative dC_ij/dJ via symmetric difference."""
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    C_A = compute_C_ij_trajectory(rho_A, N)
    C_Bp = compute_C_ij_trajectory(rho_Bp, N)
    C_Bm = compute_C_ij_trajectory(rho_Bm, N)
    dC_dJ = (C_Bp - C_Bm) / (2 * dj)  # (T, N, N)
    return C_A, dC_dJ


def run_N(N, bond=0, t_subsample=None):
    print(f"\n{'=' * 78}")
    print(f"N = {N}, bond = {bond}")
    print(f"{'=' * 78}")
    t0 = time.time()
    decomps = build_decomps(N, bond)
    print(f"  L_A, L_B+/- decomposed in {time.time() - t0:.1f} s")

    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    dt = times[1] - times[0]

    # Initial state: PTF bonding
    phi = (vacuum_ket(N) + single_excitation_mode(N, k=1)) / np.sqrt(2.0)
    rho_0 = density_matrix(phi)

    t0 = time.time()
    C_A, dC_dJ = measure_dC_ij_dJ(rho_0, decomps, times, N)
    print(f"  dC/dJ landscape computed in {time.time() - t0:.1f} s")
    print(f"  C_A shape: {C_A.shape},   dC_dJ shape: {dC_dJ.shape}")

    # Summary: find peak |dC_ij/dJ| vs distance |i-j| and time
    max_per_dist = np.zeros((len(times), N))  # max over pairs at each distance
    for t_idx in range(len(times)):
        for d in range(N):  # distance d = |i-j|
            vals = [abs(dC_dJ[t_idx, i, j]) for i in range(N) for j in range(N)
                    if abs(i - j) == d]
            max_per_dist[t_idx, d] = max(vals) if vals else 0.0

    # Find when each distance "lights up": first t where max_per_dist > threshold
    threshold = 1e-4
    first_light_t = {}
    for d in range(N):
        mask = max_per_dist[:, d] > threshold
        if mask.any():
            first_light_t[d] = float(times[mask][0])
        else:
            first_light_t[d] = None

    # Find peak value per distance (when does it max out?)
    peak_val_per_dist = {}
    peak_time_per_dist = {}
    for d in range(N):
        if max_per_dist[:, d].max() > threshold:
            peak_idx = int(np.argmax(max_per_dist[:, d]))
            peak_val_per_dist[d] = float(max_per_dist[peak_idx, d])
            peak_time_per_dist[d] = float(times[peak_idx])
        else:
            peak_val_per_dist[d] = None
            peak_time_per_dist[d] = None

    print(f"\n  Information-flow summary:")
    print(f"  {'distance':>10} {'first light (t)':>18} "
          f"{'peak value':>14} {'peak time':>12}")
    for d in range(N):
        fl = first_light_t[d]
        pv = peak_val_per_dist[d]
        pt = peak_time_per_dist[d]
        fl_s = f"{fl:.3f}" if fl is not None else "never"
        pv_s = f"{pv:.4e}" if pv is not None else "never"
        pt_s = f"{pt:.2f}" if pt is not None else "-"
        print(f"  {d:>10d} {fl_s:>18} {pv_s:>14} {pt_s:>12}")

    # Pi-pair symmetry check: dC_ij = dC_{N-1-i, N-1-j}?
    print(f"\n  Pi-pair (spatial reflection) symmetry check:")
    max_residual_pi = 0.0
    for t_idx in [len(times) // 4, len(times) // 2, 3 * len(times) // 4,
                  len(times) - 1]:
        diff = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                diff[i, j] = abs(dC_dJ[t_idx, i, j]
                                 - dC_dJ[t_idx, N - 1 - i, N - 1 - j])
        r = float(diff.max())
        max_residual_pi = max(max_residual_pi, r)
        print(f"    t = {times[t_idx]:>6.2f}: max |dC_ij - dC_{{N-1-i, N-1-j}}| "
              f"= {r:.2e}")
    # Note: pi symmetry is broken by the perturbation at bond 0
    # Unless we use the mirror bond N-2 perturbation.
    # For bond 0, the expected reflection is dC_ij(bond=0) = dC_{N-1-j, N-1-i}(bond=N-2)
    # Under the full L, [L_B(b), R_sup] = L_B(N-2-b)
    # So the "perfect" pi-pair check would be two simulations (b=0 and b=N-2-b=N-2)
    # and then dC_ij(b=0) = dC_{N-1-i, N-1-j}(b=N-2)

    # Save compact summary + sampled time slices
    if t_subsample is None:
        # Subsample at 10 time points across the window
        t_indices = list(np.linspace(0, len(times) - 1, 10).astype(int))
    else:
        t_indices = t_subsample

    out = {
        "N": N,
        "bond": bond,
        "gamma_0": GAMMA_0,
        "J": J_UNIFORM,
        "dJ": DJ_EXTRACT,
        "times_sampled": [float(times[i]) for i in t_indices],
        "C_A_sampled": [C_A[i].tolist() for i in t_indices],
        "dC_dJ_sampled": [dC_dJ[i].tolist() for i in t_indices],
        "max_per_dist_time_series": max_per_dist.tolist(),
        "first_light_t_by_dist": first_light_t,
        "peak_val_by_dist": peak_val_per_dist,
        "peak_time_by_dist": peak_time_per_dist,
        "max_pi_pair_residual": max_residual_pi,
    }
    path = RESULTS_DIR / f"info_flow_N{N}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n  Saved: {path}")

    return out


def main():
    start = time.time()
    print("=" * 78)
    print("EQ-018 info-flow landscape: C_ij(t) = Tr(rho_i(t) * rho_j(t))")
    print("=" * 78)
    print(f"  Initial state: (|vac> + |psi_1>)/sqrt(2), perturbation on bond 0")
    print(f"  Gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}")

    results = {}
    for N in [4, 5, 6]:
        results[f"N_{N}"] = run_N(N, bond=0)

    # Cross-N comparison
    print(f"\n{'=' * 78}")
    print(f"CROSS-N INFORMATION-FLOW REACH")
    print(f"{'=' * 78}")
    print(f"  N    distance   first light t   peak value   peak time")
    for N in [4, 5, 6]:
        r = results[f"N_{N}"]
        for d in range(N):
            fl = r["first_light_t_by_dist"].get(d) or r["first_light_t_by_dist"].get(str(d))
            pv = r["peak_val_by_dist"].get(d) or r["peak_val_by_dist"].get(str(d))
            pt = r["peak_time_by_dist"].get(d) or r["peak_time_by_dist"].get(str(d))
            fl_s = f"{fl:.3f}" if fl is not None else "never"
            pv_s = f"{pv:.3e}" if pv is not None else "never"
            pt_s = f"{pt:.2f}" if pt is not None else "-"
            print(f"  {N}   {d:>2d}         {fl_s:>10} {pv_s:>12} {pt_s:>8}")

    print(f"\n  Total walltime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
