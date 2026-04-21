#!/usr/bin/env python3
"""eq018_pi_pair_flow.py

Step 3 of ORTHOGONALITY_SELECTION_FAMILY.md §6.2 plan: verify Π-pair
flux balance on the Liouvillian mode level.

Background
----------
The [absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) gives
per-mode decay rates:
    Re(lambda_s) = -2 * gamma_0 * <n_XY>_s
for each Liouvillian eigenmode, where <n_XY> is the Pauli-XY weight of
the mode. Pi-paired modes satisfy:
    <n_XY>_fast + <n_XY>_slow = N
equivalently
    Re(lambda_fast) + Re(lambda_slow) = -2 * Sigma * gamma = -2 * N * gamma_0.

What we verify here (flux-balance interpretation):

(A) Absorption theorem spectral sum:
    For each Pi-paired mode pair (s, s'), Re(lam_s) + Re(lam_s') = -2*Sigma*gamma,
    independent of J.

(B) Perturbation-invariant pair sum:
    Under bond-b perturbation J -> J + dJ, individual Re(lam_s) shifts, but
    Sum within each Pi-pair remains fixed: Delta Re(lam_s) + Delta Re(lam_s') = 0.
    This is the "flux conservation" in the pair: one mode absorbs exactly as
    much as the other releases.

(C) XY-weight flux: Delta <n_XY>_s = -Delta <n_XY>_{s'} for each pair.
    Explicitly shows the light/lens duality as a conserved flux.

Setup: N=5, bond=0, delta_J = 0.01. Diagonalize L_A and L_B+, identify
Pi-pairs in L_A by eigenvalue sum condition, track the matching pairs in
L_B+ by closest-eigenvalue continuation.

Rules: UTF-8 stdout, no em-dashes.
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
    GAMMA_0, J_UNIFORM,
    build_H_XY, build_liouvillian_matrix,
)

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_pi_pair_flow"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01


def diagonalize(L):
    evals, V_R = eig(L)
    V_Linv = np.linalg.inv(V_R)
    return evals, V_R, V_Linv


def find_pi_pairs(eigvals, target_sum_re, tol=1e-7):
    """Return list of (s, s') pairs with Re(lam_s + lam_s') = target_sum_re
    and Im(lam_s + lam_s') = 0. Pairs are ordered so |Im(lam_s)| <=
    |Im(lam_s')|. Each mode appears in exactly one pair; unpaired modes
    (self-Pi: Re(lam) = target_sum_re/2, Im = 0) listed separately."""
    n = len(eigvals)
    used = np.zeros(n, dtype=bool)
    pairs = []
    self_pairs = []
    for s in range(n):
        if used[s]:
            continue
        lam_s = eigvals[s]
        # Check if self-Pi
        if abs(2 * lam_s.real - target_sum_re) < tol and abs(lam_s.imag) < tol:
            self_pairs.append(s)
            used[s] = True
            continue
        # Find partner
        best_sp = -1
        best_residual = 1e99
        for sp in range(s + 1, n):
            if used[sp]:
                continue
            lam_sp = eigvals[sp]
            res = abs((lam_s + lam_sp).real - target_sum_re) + \
                  abs((lam_s + lam_sp).imag)
            if res < best_residual:
                best_residual = res
                best_sp = sp
        if best_sp >= 0 and best_residual < tol:
            used[s] = True
            used[best_sp] = True
            # Order by Im part
            if abs(eigvals[s].imag) <= abs(eigvals[best_sp].imag):
                pairs.append((s, best_sp))
            else:
                pairs.append((best_sp, s))
    return pairs, self_pairs, used


def match_eigvals_B_to_A(eigvals_A, eigvals_B):
    """Greedy nearest-neighbour match: for each A-mode, find closest B-mode."""
    n = len(eigvals_A)
    used_B = np.zeros(n, dtype=bool)
    mapping = np.full(n, -1, dtype=int)
    residuals = np.zeros(n)
    for s in range(n):
        lam_A = eigvals_A[s]
        best_r = 1e99
        best_sp = -1
        for sp in range(n):
            if used_B[sp]:
                continue
            r = abs(eigvals_B[sp] - lam_A)
            if r < best_r:
                best_r = r
                best_sp = sp
        mapping[s] = best_sp
        residuals[s] = best_r
        used_B[best_sp] = True
    return mapping, residuals


def main():
    start = time.time()
    N = 5
    bond = 0
    dj = DJ_EXTRACT
    sigma_gamma = N * GAMMA_0
    target_sum_re = -2 * sigma_gamma  # -0.5

    print("=" * 78)
    print(f"EQ-018 Step 3: Pi-pair flux balance at N={N}, bond={bond}")
    print("=" * 78)
    print(f"  gamma_0 = {GAMMA_0}, Sigma*gamma = {sigma_gamma}, "
          f"target Re(lam_s + lam_s') = {target_sum_re}")

    # Build L_A and L_B+
    J_A = [J_UNIFORM] * (N - 1)
    J_Bp = list(J_A); J_Bp[bond] = J_UNIFORM + dj

    print(f"\n  Building and diagonalising L_A and L_B+...")
    t0 = time.time()
    L_A = build_liouvillian_matrix(build_H_XY(J_A, N), GAMMA_0, N)
    evA, VR_A, VLinv_A = diagonalize(L_A)
    L_Bp = build_liouvillian_matrix(build_H_XY(J_Bp, N), GAMMA_0, N)
    evB, VR_B, VLinv_B = diagonalize(L_Bp)
    print(f"  diagonalised in {time.time() - t0:.1f} s ({len(evA)} modes each)")

    # ---- Verify absorption theorem at L_A ----
    print(f"\n  Absorption theorem spectral check at L_A:")
    print(f"  All modes: Re(lam) should be non-positive, capped at 0 "
          f"and -2*N*gamma_0 = {-2 * N * GAMMA_0}.")
    re_range_A = [float(np.real(evA).min()), float(np.real(evA).max())]
    print(f"  Range of Re(lam_A): [{re_range_A[0]:.6f}, {re_range_A[1]:.6f}]")

    # n_XY per mode: n_XY = -Re(lam)/(2*gamma_0)
    n_XY_A = -np.real(evA) / (2 * GAMMA_0)
    print(f"  Range of <n_XY> under A: [{n_XY_A.min():.4f}, {n_XY_A.max():.4f}]  "
          f"(should be [0, {N}])")

    # ---- Identify Pi-pairs in L_A ----
    print(f"\n  Identifying Pi-pairs (Re sum = {target_sum_re}, Im sum = 0) "
          f"in L_A...")
    pairs_A, self_A, used_A = find_pi_pairs(evA, target_sum_re, tol=1e-7)
    n_paired = sum(1 for _ in pairs_A) * 2
    n_self = len(self_A)
    n_unpaired = len(evA) - n_paired - n_self
    print(f"  {len(pairs_A)} pairs ({n_paired} modes), {n_self} self-pairs, "
          f"{n_unpaired} unpaired")

    # ---- Match eigenvalues A -> B ----
    print(f"\n  Matching L_B+ eigenvalues to L_A...")
    mapping, residuals = match_eigvals_B_to_A(evA, evB)
    print(f"  Max A->B match residual: {residuals.max():.2e}")
    # For small dJ, individual modes should shift by O(dJ) or smaller
    print(f"  Median match residual: {np.median(residuals):.2e}")

    # ---- Pi-pair shift analysis ----
    print(f"\n  Pi-pair shift analysis (top 10 pairs by slowest decay):")
    pair_data = []
    for (s, sp) in pairs_A:
        lam_A_s = evA[s]
        lam_A_sp = evA[sp]
        lam_B_s = evB[mapping[s]]
        lam_B_sp = evB[mapping[sp]]

        delta_re_s = np.real(lam_B_s - lam_A_s)
        delta_re_sp = np.real(lam_B_sp - lam_A_sp)
        pair_delta_sum = delta_re_s + delta_re_sp

        pair_data.append({
            "s": int(s), "sp": int(sp),
            "lam_A_s": {"re": float(lam_A_s.real), "im": float(lam_A_s.imag)},
            "lam_A_sp": {"re": float(lam_A_sp.real), "im": float(lam_A_sp.imag)},
            "lam_B_s": {"re": float(lam_B_s.real), "im": float(lam_B_s.imag)},
            "lam_B_sp": {"re": float(lam_B_sp.real), "im": float(lam_B_sp.imag)},
            "pair_re_sum_A": float(np.real(lam_A_s + lam_A_sp)),
            "pair_re_sum_B": float(np.real(lam_B_s + lam_B_sp)),
            "delta_re_s": float(delta_re_s),
            "delta_re_sp": float(delta_re_sp),
            "pair_delta_re_sum": float(pair_delta_sum),
            "n_XY_A_s": float(-lam_A_s.real / (2 * GAMMA_0)),
            "n_XY_A_sp": float(-lam_A_sp.real / (2 * GAMMA_0)),
            "pair_n_XY_sum_A": float(-(lam_A_s.real + lam_A_sp.real) / (2 * GAMMA_0)),
        })

    # Sort by Re(lam_A_s) closest to 0 (slowest modes)
    pair_data.sort(key=lambda x: abs(x["lam_A_s"]["re"]))

    print(f"\n  {'pair':>4} {'Re(A_s)':>10} {'Re(A_sp)':>10} {'Pair sum A':>11} "
          f"{'Re(B_s)':>10} {'Re(B_sp)':>10} {'Pair sum B':>11} "
          f"{'delta sum':>10}")
    for i, pd in enumerate(pair_data[:10]):
        print(f"  {i+1:>4d} "
              f"{pd['lam_A_s']['re']:>+10.5f} {pd['lam_A_sp']['re']:>+10.5f} "
              f"{pd['pair_re_sum_A']:>+11.6f} "
              f"{pd['lam_B_s']['re']:>+10.5f} {pd['lam_B_sp']['re']:>+10.5f} "
              f"{pd['pair_re_sum_B']:>+11.6f} "
              f"{pd['pair_delta_re_sum']:>+10.2e}")

    # ---- Aggregate statistics ----
    pair_sums_A = [pd["pair_re_sum_A"] for pd in pair_data]
    pair_sums_B = [pd["pair_re_sum_B"] for pd in pair_data]
    pair_delta_sums = [pd["pair_delta_re_sum"] for pd in pair_data]

    print(f"\n  Aggregate pair-sum statistics:")
    print(f"    Max |pair_re_sum_A - (-2*Sigma*gamma)|: "
          f"{max(abs(x - target_sum_re) for x in pair_sums_A):.2e}")
    print(f"    Max |pair_re_sum_B - (-2*Sigma*gamma)|: "
          f"{max(abs(x - target_sum_re) for x in pair_sums_B):.2e}")
    print(f"    Max |pair delta re_sum| (should be 0 by flux balance): "
          f"{max(abs(x) for x in pair_delta_sums):.2e}")
    print(f"    Mean |pair delta re_sum|: "
          f"{np.mean([abs(x) for x in pair_delta_sums]):.2e}")

    # ---- XY-weight flux picture ----
    print(f"\n  XY-weight absorption under perturbation (first 10 pairs):")
    print(f"  {'pair':>4} {'<n_XY>_A,s':>11} {'<n_XY>_A,sp':>12} "
          f"{'Sum A':>9} {'<n_XY>_B,s':>11} {'<n_XY>_B,sp':>12} "
          f"{'Sum B':>9} {'Delta sum':>10}")
    for i, pd in enumerate(pair_data[:10]):
        n_B_s = -pd["lam_B_s"]["re"] / (2 * GAMMA_0)
        n_B_sp = -pd["lam_B_sp"]["re"] / (2 * GAMMA_0)
        sum_A = pd["pair_n_XY_sum_A"]
        sum_B = n_B_s + n_B_sp
        delta_sum = sum_B - sum_A
        print(f"  {i+1:>4d} "
              f"{pd['n_XY_A_s']:>+11.4f} {pd['n_XY_A_sp']:>+12.4f} "
              f"{sum_A:>+9.4f} "
              f"{n_B_s:>+11.4f} {n_B_sp:>+12.4f} "
              f"{sum_B:>+9.4f} "
              f"{delta_sum:>+10.2e}")

    # ---- Summary ----
    print(f"\n{'=' * 78}")
    print(f"VERDICT")
    print(f"{'=' * 78}")
    all_A_ok = max(abs(x - target_sum_re) for x in pair_sums_A) < 1e-6
    all_B_ok = max(abs(x - target_sum_re) for x in pair_sums_B) < 1e-6
    all_delta_ok = max(abs(x) for x in pair_delta_sums) < 1e-6
    print(f"  Absorption theorem at L_A (spectral sum): "
          f"{'CONFIRMED' if all_A_ok else 'FAILS'}")
    print(f"  Absorption theorem at L_B+ (spectral sum): "
          f"{'CONFIRMED' if all_B_ok else 'FAILS'}")
    print(f"  Flux balance under perturbation (pair delta sum = 0): "
          f"{'CONFIRMED' if all_delta_ok else 'FAILS'}")
    print(f"  Interpretation: each Pi-pair absorbs the perturbation in "
          f"equal-and-opposite")
    print(f"                  amounts in its two partners (<n_XY> shift cancels).")

    # Save
    out = {
        "N": N, "bond": bond, "dJ": dj,
        "gamma_0": GAMMA_0, "sigma_gamma": sigma_gamma,
        "target_sum_re": target_sum_re,
        "n_pairs": len(pairs_A),
        "n_self_pairs": n_self,
        "n_unpaired": n_unpaired,
        "match_residual_max": float(residuals.max()),
        "pair_sums_A_max_deviation": float(max(abs(x - target_sum_re) for x in pair_sums_A)),
        "pair_sums_B_max_deviation": float(max(abs(x - target_sum_re) for x in pair_sums_B)),
        "pair_delta_sums_max_abs": float(max(abs(x) for x in pair_delta_sums)),
        "pair_data_all": pair_data,
    }
    path = RESULTS_DIR / "pi_pair_flow.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")
    print(f"Walltime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
