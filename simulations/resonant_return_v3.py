#!/usr/bin/env python3
"""
Resonant Return V3: C# Backend Optimization
=============================================
Uses the C# RCPsiSquared.Propagate profile evaluator as backend for
scipy.optimize.minimize. Enables feasible optimization at N=7.

Test A: Individual modes at N=7 via C# (validation)
Test B: Optimizer at N=7 (key experiment)
Test C: Optimizer at N=5 (validation against v2 result)

Script:  simulations/resonant_return_v3.py
Output:  simulations/results/resonant_return_v3.txt
Depends: compute/RCPsiSquared.Propagate (profile mode)
"""

import numpy as np
from scipy.optimize import minimize
import subprocess
import json
import os
import sys
import time as _time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
CSHARP_DIR = os.path.join(PROJECT_ROOT, "compute", "RCPsiSquared.Propagate")
OUT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v3.txt")
CHECKPOINT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v3_checkpoint.json")

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# C# BACKEND
# ============================================================

def evaluate_profile_csharp(gammas, n, tmax=20.0, dt=0.05):
    """Call C# profile evaluator, return dict of metrics."""
    gamma_str = ",".join(f"{g:.6f}" for g in gammas)
    cmd = ["dotnet", "run", "-c", "Release", "--",
           "profile", str(n), gamma_str,
           "--tmax", f"{tmax}", "--dt", f"{dt}"]
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

    # Debug: print stderr if parsing failed
    sys.stderr.write(f"C# error: {result.stderr}\n")
    return None


def get_sum_mi(gammas, n, metric="SumMI"):
    """Get Sum-MI metric from C# backend."""
    result = evaluate_profile_csharp(gammas, n)
    if result is None:
        return 0.0
    return result.get(metric, 0.0)


# ============================================================
# OPTIMIZER
# ============================================================

def run_optimizer(n, gamma_base, x0_gammas, max_evals=500, metric="SumMI",
                  label="", checkpoint_path=None):
    """Run Nelder-Mead optimizer via C# backend.

    Parametrizes N-1 free deviations from gamma_base, with the Nth
    computed to maintain mean(gamma) = gamma_base exactly.
    """
    N = n
    # Adjust starting profile to have exact mean = gamma_base
    x0_arr = np.array(x0_gammas, dtype=float)
    delta_start = x0_arr - gamma_base
    delta_start -= np.mean(delta_start)
    x0 = delta_start[:N - 1]

    eval_count = [0]
    best_found = [0.0]
    best_gammas = [list(x0_gammas)]
    all_evals = []
    t0 = _time.time()

    def objective(x):
        eval_count[0] += 1
        delta_last = -np.sum(x)
        gammas = [gamma_base + x[k] for k in range(N - 1)] + [gamma_base + delta_last]

        # Penalty for gammas below threshold
        penalty = 0.0
        for g in gammas:
            if g < 0.001:
                penalty += 10000 * (0.001 - g) ** 2
        if penalty > 0:
            return penalty

        mi = get_sum_mi(gammas, N, metric)

        if mi > best_found[0]:
            best_found[0] = mi
            best_gammas[0] = list(gammas)

        record = {
            "eval": eval_count[0],
            "gammas": [round(g, 6) for g in gammas],
            "metric": mi,
            "best": best_found[0],
            "elapsed": round(_time.time() - t0, 1),
        }
        all_evals.append(record)

        if eval_count[0] % 10 == 0:
            elapsed = _time.time() - t0
            log(f"    [{eval_count[0]:4d} evals, {elapsed:6.1f}s] "
                f"best {metric} = {best_found[0]:.6f}")

            # Save checkpoint
            if checkpoint_path:
                ckpt = {
                    "label": label, "n": N, "metric": metric,
                    "eval_count": eval_count[0],
                    "best_mi": best_found[0],
                    "best_gammas": best_gammas[0],
                    "elapsed_s": round(elapsed, 1),
                    "history": all_evals[-50:],  # keep last 50
                }
                with open(checkpoint_path, "w") as f:
                    json.dump(ckpt, f, indent=2)

        return -mi

    log(f"  Running Nelder-Mead optimizer ({label}, N={N}, max {max_evals} evals)...")
    log(f"  Starting profile: [{', '.join(f'{g:.4f}' for g in x0_gammas)}]")
    log(f"  Starting {metric}: {get_sum_mi(x0_gammas, N, metric):.6f}")
    log()

    result = minimize(objective, x0, method='Nelder-Mead',
                      options={'maxfev': max_evals, 'xatol': 1e-5,
                               'fatol': 1e-8, 'adaptive': True})

    total_time = _time.time() - t0

    # Reconstruct final profile
    x_opt = result.x
    delta_last = -np.sum(x_opt)
    opt_gammas = [gamma_base + x_opt[k] for k in range(N - 1)] + [gamma_base + delta_last]

    # Re-evaluate to confirm
    final_result = evaluate_profile_csharp(opt_gammas, N)
    final_mi = final_result.get(metric, 0.0) if final_result else best_found[0]

    log()
    log(f"  Optimizer finished: {eval_count[0]} evals in {total_time:.1f}s "
        f"({total_time / eval_count[0]:.1f}s/eval)")
    log(f"  Converged: {result.success} ({result.message})")
    log(f"  Optimized profile: [{', '.join(f'{g:.4f}' for g in opt_gammas)}]")
    log(f"  Mean(gamma): {np.mean(opt_gammas):.6f}")
    log(f"  Optimized {metric}: {final_mi:.6f}")
    log(f"  Best found during search: {best_found[0]:.6f}")

    if final_result:
        log(f"  Full result: SumMI={final_result.get('SumMI', 0):.6f} "
            f"SumMI5={final_result.get('SumMI5', 0):.6f} "
            f"PeakMI={final_result.get('PeakMI', 0):.6f} "
            f"PeakT={final_result.get('PeakT', 0):.2f}")

    return {
        'opt_gammas': opt_gammas,
        'best_gammas': best_gammas[0],
        'final_mi': final_mi,
        'best_mi': best_found[0],
        'evals': eval_count[0],
        'time_s': total_time,
        'converged': result.success,
        'full_result': final_result,
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    t_start = _time.time()
    gamma_base = 0.05

    log("=== RESONANT RETURN V3: C# BACKEND OPTIMIZATION ===")
    log(f"Started: {_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"C# backend: {CSHARP_DIR}")
    log()

    # ========================================================
    # TEST C: Optimizer validation at N=5
    # ========================================================
    log("=" * 70)
    log("TEST C: OPTIMIZER VALIDATION (N=5, C# backend)")
    log("=" * 70)
    log()
    log("  Target: reproduce Python v2 result SumMI5 ~ 0.031071")
    log("  Using peak SumMI as objective (C# advantage over Python)")
    log()

    # N=5 mode 2 profile as starting point
    x0_n5 = [0.0703, 0.0474, 0.0181, 0.0474, 0.0703]

    res_c = run_optimizer(
        n=5, gamma_base=gamma_base, x0_gammas=x0_n5,
        max_evals=500, metric="SumMI", label="N=5 validation",
    )

    log()
    python_v2_mi5 = 0.031071
    c_mi5 = res_c['full_result'].get('SumMI5', 0) if res_c['full_result'] else 0
    log(f"  VALIDATION:")
    log(f"    Python v2 SumMI@5:  {python_v2_mi5:.6f}")
    log(f"    C# optimizer SumMI5: {c_mi5:.6f}")
    log(f"    C# optimizer SumMI (peak): {res_c['best_mi']:.6f}")
    match_pct = c_mi5 / python_v2_mi5 * 100 if python_v2_mi5 > 0 else 0
    log(f"    Match: {match_pct:.1f}%")
    if abs(match_pct - 100) < 10:
        log("    [PASS] Within 10% tolerance")
    else:
        log(f"    [NOTE] Outside 10% tolerance ({match_pct:.1f}%)")
    log()

    # ========================================================
    # TEST B: Optimizer at N=7 (the key experiment)
    # ========================================================
    log("=" * 70)
    log("TEST B: OPTIMIZER AT N=7 (C# backend)")
    log("=" * 70)
    log()
    log("  Key question: does N=7 also show asymmetric sacrifice-zone?")
    log("  Starting from: N=7 SVD mode 2 profile")
    log()

    # N=7 mode 2 profile
    x0_n7 = [0.0718, 0.0759, 0.0636, 0.0500, 0.0364, 0.0241, 0.0282]

    res_b = run_optimizer(
        n=7, gamma_base=gamma_base, x0_gammas=x0_n7,
        max_evals=500, metric="SumMI", label="N=7 optimizer",
        checkpoint_path=CHECKPOINT_PATH,
    )

    log()

    # Compare with baselines
    vshape_n7 = get_sum_mi([0.08, 0.07, 0.06, 0.05, 0.06, 0.07, 0.08], 7, "SumMI")
    mode2_n7 = get_sum_mi(x0_n7, 7, "SumMI")

    log(f"  COMPARISON (N=7, peak SumMI):")
    log(f"    V-shape:            {vshape_n7:.6f}")
    log(f"    Mode 2:             {mode2_n7:.6f} ({mode2_n7/vshape_n7:.1f}x vs V-shape)")
    log(f"    Optimizer:          {res_b['best_mi']:.6f} ({res_b['best_mi']/vshape_n7:.1f}x vs V-shape)")
    log(f"    Improvement over mode 2: {res_b['best_mi']/mode2_n7:.2f}x")
    log()

    # Analyze the optimal profile
    opt_g = np.array(res_b['best_gammas'])
    log(f"  PROFILE ANALYSIS:")
    log(f"    Profile: [{', '.join(f'{g:.4f}' for g in opt_g)}]")
    log(f"    Mean: {np.mean(opt_g):.4f}")
    log(f"    Min: {np.min(opt_g):.4f} at site {np.argmin(opt_g)}")
    log(f"    Max: {np.max(opt_g):.4f} at site {np.argmax(opt_g)}")
    log(f"    Std: {np.std(opt_g):.4f}")

    # Symmetry check: compare with mirror profile
    mirror = opt_g[::-1]
    asym = np.linalg.norm(opt_g - mirror) / np.linalg.norm(opt_g)
    log(f"    Asymmetry (||g - g_mirror|| / ||g||): {asym:.4f}")
    if asym > 0.1:
        log("    -> ASYMMETRIC (sacrifice-zone pattern)")
    else:
        log("    -> SYMMETRIC (palindromic pattern)")
    log()

    # ========================================================
    # SUMMARY
    # ========================================================
    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()

    log("  N=5 optimizer (C# backend):")
    log(f"    Best profile: [{', '.join(f'{g:.4f}' for g in res_c['best_gammas'])}]")
    log(f"    Peak SumMI: {res_c['best_mi']:.6f}")
    log(f"    SumMI@5: {res_c['full_result'].get('SumMI5', 0):.6f}" if res_c['full_result'] else "    SumMI@5: N/A")
    log(f"    Evals: {res_c['evals']}, Time: {res_c['time_s']:.0f}s")
    log()

    log("  N=7 optimizer (C# backend):")
    log(f"    Best profile: [{', '.join(f'{g:.4f}' for g in res_b['best_gammas'])}]")
    log(f"    Peak SumMI: {res_b['best_mi']:.6f}")
    log(f"    SumMI@5: {res_b['full_result'].get('SumMI5', 0):.6f}" if res_b['full_result'] else "    SumMI@5: N/A")
    log(f"    Evals: {res_b['evals']}, Time: {res_b['time_s']:.0f}s")
    log(f"    vs Mode 2: {res_b['best_mi']/mode2_n7:.2f}x")
    log(f"    vs V-shape: {res_b['best_mi']/vshape_n7:.1f}x")
    log()

    log(f"  KEY FINDING:")
    if asym > 0.1:
        log(f"    N=7 optimizer finds ASYMMETRIC profile (asymmetry={asym:.2f})")
        log(f"    The sacrifice-zone pattern is universal, not a small-N artifact.")
    else:
        log(f"    N=7 optimizer finds SYMMETRIC profile (asymmetry={asym:.2f})")
        log(f"    The N=5 asymmetry was a small-system effect.")
    log()

    total_time = _time.time() - t_start
    log(f"Total runtime: {total_time:.1f}s ({total_time / 60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
