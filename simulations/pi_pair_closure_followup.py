#!/usr/bin/env python3
"""Follow-up to pi_pair_closure_investigation: sharper characterization.

Addresses three non-tested items from the first pass:
  1. Defect location scan: is c_1 (linear coefficient of closure vs delta_J)
     invariant under defect position along the chain?
  2. Systematic c_1 extraction across initial states: for which states
     does c_1 -> 0 (strict protection)?
  3. Voice-mode-projected alpha: subtract stationary-subspace contribution
     from rho(t), refit alpha_i, check if closure is cleaner.

Shared setup: XY chain, uniform Z-dephasing gamma_0 = 0.05, N in {3, 5}.
Extract c_1 numerically via 3-point finite difference at small delta_J.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).parent))
from pi_pair_closure_investigation import (
    GAMMA_0, J_UNIFORM, T_FINAL, N_STEPS, T_FIT_MAX, ALPHA_BOUNDS,
    X, Y, Z, I2, site_op,
    build_H_XY, build_liouvillian_matrix,
    vacuum_ket, single_excitation_mode, bonding_plus_vacuum,
    density_matrix, propagate_vectorised,
    per_site_purity, partial_trace_keep_site_fast,
    fit_alpha,
)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results" / "pi_pair_closure_investigation"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Small delta_J for linear coefficient extraction (avoids quadratic contamination)
DJ_LIN = 0.01


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def compute_c1_and_c2(N, defect_bond, initial_state, times, L_A, P_A):
    """Fit c_1, c_2 of closure(dJ) = c_0 + c_1 * dJ + c_2 * dJ^2 via 5-point stencil."""
    J_A = [J_UNIFORM] * (N - 1)
    rho_0 = density_matrix(initial_state)
    dJ_values = [-2 * DJ_LIN, -DJ_LIN, 0, +DJ_LIN, +2 * DJ_LIN]
    closures = []
    for dJ in dJ_values:
        if abs(dJ) < 1e-12:
            closures.append(0.0)
            continue
        J_B = list(J_A); J_B[defect_bond[0]] += dJ
        L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
        rho_B = propagate_vectorised(L_B, rho_0, times)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
        closures.append(float(np.sum(np.log(alpha))))
    # Fit cubic: c_0 + c_1 dJ + c_2 dJ^2 + c_3 dJ^3
    coeffs = np.polyfit(dJ_values, closures, 3)
    # coeffs = [c_3, c_2, c_1, c_0]
    return {"dJ": dJ_values, "closure": closures,
            "c_0": float(coeffs[3]), "c_1": float(coeffs[2]),
            "c_2": float(coeffs[1]), "c_3": float(coeffs[0])}


def alpha_from_Pdiff(times, P_A, P_B, N):
    """Alpha fit on the subtracted profile P(t) - P(infty).
    P(infty) is taken as P at the last time point (approximation for T large)."""
    alpha = np.zeros(N)
    for i in range(N):
        P_A_inf = float(P_A[-1, i])
        P_B_inf = float(P_B[-1, i])
        pA_diff = P_A[:, i] - P_A_inf
        pB_diff = P_B[:, i] - P_B_inf
        # Avoid divide-by-zero: if pA_diff is uniformly small, flag NaN
        scale = float(np.max(np.abs(pA_diff)))
        if scale < 1e-10:
            alpha[i] = np.nan
            continue
        # Interpolate and fit
        interp = interp1d(times, pA_diff, bounds_error=False,
                          fill_value=(float(pA_diff[0]), float(pA_diff[-1])),
                          kind='cubic')
        mask = times <= T_FIT_MAX
        te = times[mask]
        b = pB_diff[mask]
        def mse(a):
            d = interp(a * te) - b
            return float(np.mean(d * d))
        res = minimize_scalar(mse, bounds=ALPHA_BOUNDS, method='bounded',
                              options={'xatol': 1e-6})
        alpha[i] = float(res.x)
    return alpha


# ---------------------------------------------------------------------------
# Part 1: Defect-location scan
# ---------------------------------------------------------------------------
def defect_location_scan(N):
    """For each bond in the chain, extract c_1 of closure(dJ) and the
    alpha_i pattern at dJ = +0.01 (linear regime). Check mirror symmetry."""
    print(f"\n{'='*60}")
    print(f"Defect location scan at N={N}")
    print(f"{'='*60}")
    J_A = [J_UNIFORM] * (N - 1)
    H_A = build_H_XY(J_A, N)
    L_A = build_liouvillian_matrix(H_A, GAMMA_0, N)
    phi = bonding_plus_vacuum(N, k=1)
    rho_0 = density_matrix(phi)
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    rho_A = propagate_vectorised(L_A, rho_0, times)
    P_A = per_site_purity(rho_A, N)

    results = []
    print(f"  Bond   c_0         c_1          c_2         c_3         alpha_i@dJ=+0.01")
    for b in range(N - 1):
        r = compute_c1_and_c2(N, (b, b+1), phi, times, L_A, P_A)
        # Also compute alpha at dJ=+0.01 specifically
        J_B = list(J_A); J_B[b] += DJ_LIN
        L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
        rho_B = propagate_vectorised(L_B, rho_0, times)
        P_B = per_site_purity(rho_B, N)
        alpha = np.zeros(N)
        for i in range(N):
            a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
            alpha[i] = a
        print(f"  ({b},{b+1}) {r['c_0']:+.4f}    {r['c_1']:+.4f}    "
              f"{r['c_2']:+.4f}    {r['c_3']:+.4f}    "
              + " ".join(f"{a:.4f}" for a in alpha))
        results.append({
            "bond": [b, b+1],
            "c_0": r["c_0"], "c_1": r["c_1"],
            "c_2": r["c_2"], "c_3": r["c_3"],
            "alpha_at_dJ_0.01": alpha.tolist(),
        })
    # Mirror symmetry check: bond (b, b+1) mirrors to (N-2-b, N-1-b)
    print(f"\n  Mirror symmetry (alpha_i at bond b should equal reversed alpha_i at mirror bond):")
    for b in range(N - 1):
        b_mirror = N - 2 - b
        if b >= b_mirror:
            continue
        alpha_b = np.array(results[b]["alpha_at_dJ_0.01"])
        alpha_m = np.array(results[b_mirror]["alpha_at_dJ_0.01"])[::-1]
        diff = np.linalg.norm(alpha_b - alpha_m)
        c1_diff = abs(results[b]["c_1"] - results[b_mirror]["c_1"])
        print(f"    bond ({b},{b+1}) <-> ({b_mirror},{b_mirror+1}): "
              f"||alpha - reversed(alpha_m)|| = {diff:.3e}, "
              f"|c_1 - c_1_mirror| = {c1_diff:.3e}")
    return results


# ---------------------------------------------------------------------------
# Part 2: Systematic c_1 extraction across initial states
# ---------------------------------------------------------------------------
def initial_state_c1_scan(N, defect_bond=(0, 1)):
    """Extract c_1 for a suite of initial states."""
    print(f"\n{'='*60}")
    print(f"c_1 scan over initial states at N={N}, defect bond = {defect_bond}")
    print(f"{'='*60}")
    d = 2**N
    J_A = [J_UNIFORM] * (N - 1)
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    states = {}
    for k in range(1, N + 1):
        states[f"psi_{k}+vac"] = (vacuum_ket(N) + single_excitation_mode(N, k)) / np.sqrt(2.0)
        states[f"psi_{k}_only"] = single_excitation_mode(N, k)
    # |+>^N
    plus_N = np.ones(d, dtype=complex) / np.sqrt(d)
    states["plus_N"] = plus_N
    # Haar-random single-excitation state
    rng = np.random.default_rng(42)
    coeffs = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    coeffs /= np.linalg.norm(coeffs)
    haar_1exc = np.zeros(d, dtype=complex)
    for i in range(N):
        haar_1exc[2**(N - 1 - i)] = coeffs[i]
    states["haar_1exc_only"] = haar_1exc
    states["haar_1exc+vac"] = (vacuum_ket(N) + haar_1exc) / np.sqrt(2.0)

    print(f"  State                   c_1         c_2         closure(+0.01)  closure(-0.01)")
    rows = []
    for name, ket in states.items():
        rho_0 = density_matrix(ket)
        rho_A = propagate_vectorised(L_A, rho_0, times)
        P_A = per_site_purity(rho_A, N)
        r = compute_c1_and_c2(N, defect_bond, ket, times, L_A, P_A)
        # Direct closure at ±0.01
        cl_pos = r["closure"][3]  # dJ=+DJ_LIN
        cl_neg = r["closure"][1]  # dJ=-DJ_LIN
        print(f"  {name:22s}  {r['c_1']:+.4f}    {r['c_2']:+.4f}    "
              f"{cl_pos:+.5f}    {cl_neg:+.5f}")
        rows.append({"state": name, "c_1": r["c_1"], "c_2": r["c_2"],
                     "closure_plus": cl_pos, "closure_neg": cl_neg})
    return rows


# ---------------------------------------------------------------------------
# Part 3: Voice-mode-projected alpha
# ---------------------------------------------------------------------------
def voice_projected_alpha(N, defect_bond=(0, 1), J_mod=1.1):
    """Subtract the asymptotic (t -> infty) site marginals, refit alpha.

    Rationale: the strict stationary subspace is J-invariant (verified above),
    so any 'J-signature' in the alpha_i fit comes from the transient (voice)
    part of rho(t). Fitting P_B(i,t) - P_B(i, infty) against P_A(i,t) - P_A(i, infty)
    removes the trivial memory contribution."""
    print(f"\n{'='*60}")
    print(f"Voice-mode projected alpha at N={N}, bond={defect_bond}, J_mod={J_mod}")
    print(f"{'='*60}")
    J_A = [J_UNIFORM] * (N - 1)
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    J_B = list(J_A); J_B[defect_bond[0]] = J_mod
    L_B = build_liouvillian_matrix(build_H_XY(J_B, N), GAMMA_0, N)
    phi = bonding_plus_vacuum(N, k=1)
    rho_0 = density_matrix(phi)
    # Use a long time horizon so P(infty) is reached
    times = np.linspace(0.0, 200.0, 801)
    rho_A = propagate_vectorised(L_A, rho_0, times)
    rho_B = propagate_vectorised(L_B, rho_0, times)
    P_A = per_site_purity(rho_A, N)
    P_B = per_site_purity(rho_B, N)

    # Asymptotic purity (at t = 200)
    P_A_inf = P_A[-1]
    P_B_inf = P_B[-1]
    diff_inf = np.abs(P_A_inf - P_B_inf)
    print(f"  P_A(infty) = " + "  ".join(f"{p:.5f}" for p in P_A_inf))
    print(f"  P_B(infty) = " + "  ".join(f"{p:.5f}" for p in P_B_inf))
    print(f"  |P_A(inf) - P_B(inf)| = " + "  ".join(f"{d:.1e}" for d in diff_inf)
          + f"   max = {diff_inf.max():.2e}")

    # Standard alpha fit (baseline)
    alpha_std = np.zeros(N)
    for i in range(N):
        a, _ = fit_alpha(times, P_A[:, i], P_B[:, i])
        alpha_std[i] = a
    closure_std = float(np.sum(np.log(alpha_std)))
    print(f"  Standard alpha:      " + "  ".join(f"{a:.4f}" for a in alpha_std)
          + f"     Sum ln = {closure_std:+.5f}")

    # Voice-projected alpha (subtract infty)
    alpha_voice = alpha_from_Pdiff(times, P_A, P_B, N)
    if np.all(np.isfinite(alpha_voice)):
        closure_voice = float(np.sum(np.log(alpha_voice)))
        print(f"  Voice alpha:         " + "  ".join(f"{a:.4f}" for a in alpha_voice)
              + f"     Sum ln = {closure_voice:+.5f}")
    else:
        print(f"  Voice alpha:         (some NaN, skipped closure)")
        closure_voice = None
    return {
        "alpha_std": alpha_std.tolist(), "closure_std": closure_std,
        "alpha_voice": alpha_voice.tolist(),
        "closure_voice": closure_voice,
        "P_A_inf": P_A_inf.tolist(), "P_B_inf": P_B_inf.tolist(),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("="*60)
    print("Pi-pair closure follow-up")
    print(f"  delta_J for linear fit = +/- {DJ_LIN}")
    print("="*60)

    out = {}

    for N in [3, 5]:
        print(f"\n\n{'#'*60}\n# N = {N}\n{'#'*60}")
        out[f"n{N}_defect_scan"] = defect_location_scan(N)
        out[f"n{N}_state_c1"] = initial_state_c1_scan(N)
        out[f"n{N}_voice"] = voice_projected_alpha(N)

    path = RESULTS_DIR / "followup.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")


if __name__ == "__main__":
    main()
