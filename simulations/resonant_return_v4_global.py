#!/usr/bin/env python3
"""
Resonant Return V4: Global Optimization via Differential Evolution
==================================================================
Uses scipy.optimize.differential_evolution with C# RK4 backend.
Tests whether the Nelder-Mead sacrifice-zone (SumMI=0.1439) is
the global optimum or just a local one.

Script:  simulations/resonant_return_v4_global.py
Output:  simulations/results/resonant_return_v4_global.txt
"""

import numpy as np
from scipy.optimize import differential_evolution
import subprocess
import json
import os
import sys
import time as _time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
CSHARP_DIR = os.path.join(PROJECT_ROOT, "compute", "RCPsiSquared.Propagate")
OUT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v4_global.txt")
CHECKPOINT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v4_checkpoint.json")

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


def evaluate_profile_csharp(gammas, n):
    gamma_str = ",".join(f"{g:.6f}" for g in gammas)
    cmd = ["dotnet", "run", "-c", "Release", "--",
           "profile", str(n), gamma_str]
    try:
        result = subprocess.run(
            cmd, cwd=CSHARP_DIR,
            capture_output=True, text=True, timeout=600
        )
    except subprocess.TimeoutExpired:
        return None
    for line in result.stdout.strip().split("\n"):
        if line.startswith("RESULT"):
            parts = dict(kv.split("=") for kv in line.split()[1:])
            return {k: float(v) for k, v in parts.items()}
    return None


if __name__ == "__main__":
    t_start = _time.time()
    N = 7
    gamma_base = 0.05
    NELDER_BEST = 0.143910  # from v3 deep run (converged)
    NELDER_PROFILE = [0.1296, 0.1219, 0.0455, 0.0500, 0.0010, 0.0010, 0.0010]

    log("=== RESONANT RETURN V4: GLOBAL OPTIMIZATION (Differential Evolution) ===")
    log(f"Started: {_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, gamma_base={gamma_base}")
    log(f"Nelder-Mead reference: SumMI={NELDER_BEST:.6f}")
    log(f"Nelder-Mead profile:   [{', '.join(f'{g:.4f}' for g in NELDER_PROFILE)}]")
    log()

    # ========================================================
    # Differential Evolution setup
    # ========================================================
    # Parametrize: N-1 free deviations from gamma_base.
    # delta_k in [-0.049, 0.20] (so gamma_k in [0.001, 0.25])
    # The Nth delta is computed as -sum(others) to maintain mean.
    #
    # Bounds: each delta in [-0.045, 0.15]
    # (conservative to keep all gammas physical)
    bounds = [(-0.045, 0.15)] * (N - 1)

    eval_count = [0]
    best_found = [0.0]
    best_gammas = [None]
    t0 = _time.time()

    def objective(x):
        eval_count[0] += 1
        delta_last = -np.sum(x)
        gammas = [gamma_base + x[k] for k in range(N - 1)] + [gamma_base + delta_last]

        # Hard constraint: all gammas must be physical
        for g in gammas:
            if g < 0.001 or g > 0.30:
                return 1.0  # penalty (DE minimizes, so positive = bad)

        result = evaluate_profile_csharp(gammas, N)
        if result is None:
            return 1.0

        mi = result.get("SumMI", 0.0)

        if mi > best_found[0]:
            best_found[0] = mi
            best_gammas[0] = list(gammas)

        if eval_count[0] % 20 == 0:
            elapsed = _time.time() - t0
            rate = elapsed / eval_count[0]
            log(f"  [{eval_count[0]:5d} evals  {elapsed:7.0f}s  {rate:.1f}s/eval]  "
                f"best SumMI = {best_found[0]:.6f}  "
                f"({'>' if best_found[0] > NELDER_BEST else '<='} Nelder)")

            ckpt = {
                "eval": eval_count[0],
                "best_mi": best_found[0],
                "best_gammas": best_gammas[0],
                "elapsed_s": round(elapsed, 1),
                "vs_nelder": round(best_found[0] / NELDER_BEST, 4) if NELDER_BEST > 0 else 0,
            }
            with open(CHECKPOINT_PATH, "w") as f:
                json.dump(ckpt, f, indent=2)

        return -mi

    log("  Differential Evolution parameters:")
    log(f"    Strategy: best1bin")
    log(f"    Population: 15 (= {15*(N-1)} total per generation)")
    log(f"    Max iterations: 30")
    log(f"    Bounds per delta: [{bounds[0][0]}, {bounds[0][1]}]")
    log(f"    Estimated evals: ~{15*(N-1)*30} (upper bound)")
    log()

    result = differential_evolution(
        objective,
        bounds,
        strategy='best1bin',
        maxiter=30,
        popsize=15,
        tol=1e-6,
        mutation=(0.5, 1.5),
        recombination=0.8,
        seed=42,
        disp=False,
        init='sobol',
    )

    total_time = _time.time() - t_start

    # Reconstruct final profile
    x_opt = result.x
    delta_last = -np.sum(x_opt)
    de_gammas = [gamma_base + x_opt[k] for k in range(N - 1)] + [gamma_base + delta_last]
    de_result = evaluate_profile_csharp(de_gammas, N)
    de_mi = de_result['SumMI'] if de_result else -result.fun

    # Also evaluate best_gammas found during search (may differ from final)
    best_result = evaluate_profile_csharp(best_gammas[0], N) if best_gammas[0] else None
    best_mi = best_result['SumMI'] if best_result else best_found[0]

    log()
    log("=" * 60)
    log("RESULTS")
    log("=" * 60)
    log()
    log(f"  Evals: {eval_count[0]}, Time: {total_time:.0f}s ({total_time/60:.1f} min)")
    log(f"  DE converged: {result.success} ({result.message})")
    log()

    # Report DE result
    log(f"  DE final profile:  [{', '.join(f'{g:.4f}' for g in de_gammas)}]")
    log(f"  DE final SumMI:    {de_mi:.6f}")
    log()

    # Report best found during search
    if best_gammas[0]:
        log(f"  Best found profile: [{', '.join(f'{g:.4f}' for g in best_gammas[0])}]")
        log(f"  Best found SumMI:   {best_mi:.6f}")
        if best_result:
            log(f"  Best SumMI@5:       {best_result.get('SumMI5', 0):.6f}")
            log(f"  Best PeakT:         {best_result.get('PeakT', 0):.2f}")
    log()

    # Compare with Nelder-Mead
    winner_mi = max(de_mi, best_mi)
    winner_gammas = de_gammas if de_mi >= best_mi else best_gammas[0]

    log(f"  COMPARISON WITH NELDER-MEAD:")
    log(f"    Nelder-Mead:  {NELDER_BEST:.6f}  [{', '.join(f'{g:.4f}' for g in NELDER_PROFILE)}]")
    log(f"    DE best:      {winner_mi:.6f}  [{', '.join(f'{g:.4f}' for g in winner_gammas)}]")
    ratio = winner_mi / NELDER_BEST if NELDER_BEST > 0 else 0
    log(f"    Ratio DE/NM:  {ratio:.4f}")
    log()

    if ratio > 1.02:
        log(f"  ** DE FOUND BETTER OPTIMUM: {winner_mi:.6f} > {NELDER_BEST:.6f} (+{(ratio-1)*100:.1f}%) **")
        log(f"  The Nelder-Mead result was a local optimum!")
    elif ratio > 0.98:
        log(f"  DE CONFIRMS NELDER-MEAD RESULT (within 2%)")
        log(f"  The sacrifice-zone at SumMI ~ {NELDER_BEST:.4f} is likely the GLOBAL optimum.")
    else:
        log(f"  DE found worse result - may need more iterations.")
    log()

    # Profile analysis
    opt_g = np.array(winner_gammas)
    mirror = opt_g[::-1]
    asym = np.linalg.norm(opt_g - mirror) / np.linalg.norm(opt_g)
    log(f"  PROFILE ANALYSIS:")
    log(f"    Min: {np.min(opt_g):.4f} at site {np.argmin(opt_g)}")
    log(f"    Max: {np.max(opt_g):.4f} at site {np.argmax(opt_g)}")
    log(f"    Asymmetry: {asym:.4f}")

    # Similarity to Nelder-Mead profile
    nm_g = np.array(NELDER_PROFILE)
    # Compare both orientations (original and mirror)
    dist_orig = np.linalg.norm(opt_g - nm_g) / np.linalg.norm(nm_g)
    dist_mirror = np.linalg.norm(opt_g - nm_g[::-1]) / np.linalg.norm(nm_g)
    log(f"    Distance to NM profile: {min(dist_orig, dist_mirror):.4f} "
        f"({'original' if dist_orig < dist_mirror else 'mirror'} orientation)")
    log()

    # V-shape comparison
    vshape_r = evaluate_profile_csharp([0.08, 0.07, 0.06, 0.05, 0.06, 0.07, 0.08], N)
    vshape_mi = vshape_r['SumMI'] if vshape_r else 0.002412
    log(f"  IMPROVEMENT:")
    log(f"    vs V-shape: {winner_mi/vshape_mi:.1f}x")
    log(f"    vs Mode 2:  {winner_mi/0.019501:.1f}x")
    log()

    log(f"Total runtime: {total_time:.0f}s ({total_time/60:.1f} min)")
    log(f"Results: {OUT_PATH}")
    _outf.close()
