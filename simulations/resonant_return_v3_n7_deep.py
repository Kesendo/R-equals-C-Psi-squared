#!/usr/bin/env python3
"""
N=7 Deep Optimizer Run (2000 evals) via C# backend.
Continuation of resonant_return_v3.py Step 3.
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
OUT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v3_n7_deep.txt")
CHECKPOINT_PATH = os.path.join(SCRIPT_DIR, "results", "resonant_return_v3_n7_deep_checkpoint.json")

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
    MAX_EVALS = 2000

    log("=== N=7 DEEP OPTIMIZER RUN (2000 evals, C# backend) ===")
    log(f"Started: {_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # Starting point: N=7 mode 2 profile
    x0_gammas = [0.0718, 0.0759, 0.0636, 0.0500, 0.0364, 0.0241, 0.0282]

    # Adjust to exact mean = gamma_base
    delta_start = np.array(x0_gammas) - gamma_base
    delta_start -= np.mean(delta_start)
    x0 = delta_start[:N - 1]

    eval_count = [0]
    best_found = [0.0]
    best_gammas = [list(x0_gammas)]
    t0 = _time.time()

    def objective(x):
        eval_count[0] += 1
        delta_last = -np.sum(x)
        gammas = [gamma_base + x[k] for k in range(N - 1)] + [gamma_base + delta_last]

        penalty = 0.0
        for g in gammas:
            if g < 0.001:
                penalty += 10000 * (0.001 - g) ** 2
        if penalty > 0:
            return penalty

        result = evaluate_profile_csharp(gammas, N)
        if result is None:
            return 1.0  # penalty for failed eval

        mi = result.get("SumMI", 0.0)

        if mi > best_found[0]:
            best_found[0] = mi
            best_gammas[0] = list(gammas)

        if eval_count[0] % 10 == 0:
            elapsed = _time.time() - t0
            rate = elapsed / eval_count[0]
            eta = rate * (MAX_EVALS - eval_count[0])
            log(f"  [{eval_count[0]:4d}/{MAX_EVALS}  {elapsed:7.0f}s  {rate:.1f}s/eval  "
                f"ETA {eta/60:.0f}min]  best SumMI = {best_found[0]:.6f}")

            # Checkpoint
            ckpt = {
                "eval": eval_count[0], "max_evals": MAX_EVALS,
                "best_mi": best_found[0],
                "best_gammas": best_gammas[0],
                "elapsed_s": round(elapsed, 1),
            }
            with open(CHECKPOINT_PATH, "w") as f:
                json.dump(ckpt, f, indent=2)

        return -mi

    start_result = evaluate_profile_csharp(x0_gammas, N)
    start_mi = start_result.get("SumMI", 0) if start_result else 0
    log(f"  Starting profile: [{', '.join(f'{g:.4f}' for g in x0_gammas)}]")
    log(f"  Starting SumMI: {start_mi:.6f}")
    log()

    result = minimize(objective, x0, method='Nelder-Mead',
                      options={'maxfev': MAX_EVALS, 'xatol': 1e-5,
                               'fatol': 1e-8, 'adaptive': True})

    total_time = _time.time() - t_start

    # Reconstruct final
    x_opt = result.x
    delta_last = -np.sum(x_opt)
    opt_gammas = [gamma_base + x_opt[k] for k in range(N - 1)] + [gamma_base + delta_last]
    final = evaluate_profile_csharp(opt_gammas, N)
    best_final = evaluate_profile_csharp(best_gammas[0], N)

    log()
    log("=" * 60)
    log("RESULTS")
    log("=" * 60)
    log()
    log(f"  Evals: {eval_count[0]}, Time: {total_time:.0f}s ({total_time/60:.1f} min)")
    log(f"  Converged: {result.success} ({result.message})")
    log()

    # Report best found (may differ from final simplex vertex)
    log(f"  Best profile: [{', '.join(f'{g:.4f}' for g in best_gammas[0])}]")
    log(f"  Mean(gamma): {np.mean(best_gammas[0]):.6f}")
    if best_final:
        log(f"  Best SumMI (peak): {best_final['SumMI']:.6f}")
        log(f"  Best SumMI@5:      {best_final.get('SumMI5', 0):.6f}")
        log(f"  PeakMI:            {best_final.get('PeakMI', 0):.6f}")
        log(f"  PeakT:             {best_final.get('PeakT', 0):.2f}")
    log()

    # Compare with baselines
    vshape_r = evaluate_profile_csharp([0.08, 0.07, 0.06, 0.05, 0.06, 0.07, 0.08], N)
    mode2_r = evaluate_profile_csharp(x0_gammas, N)
    vshape_mi = vshape_r['SumMI'] if vshape_r else 0
    mode2_mi = mode2_r['SumMI'] if mode2_r else 0
    best_mi = best_final['SumMI'] if best_final else best_found[0]

    log(f"  COMPARISON (peak SumMI):")
    log(f"    V-shape:      {vshape_mi:.6f}")
    log(f"    Mode 2:       {mode2_mi:.6f}  ({mode2_mi/vshape_mi:.1f}x vs V-shape)")
    log(f"    Optimized:    {best_mi:.6f}  ({best_mi/vshape_mi:.1f}x vs V-shape)")
    log(f"    Opt/Mode2:    {best_mi/mode2_mi:.2f}x")
    log()

    # Profile analysis
    opt_g = np.array(best_gammas[0])
    mirror = opt_g[::-1]
    asym = np.linalg.norm(opt_g - mirror) / np.linalg.norm(opt_g)
    log(f"  PROFILE ANALYSIS:")
    log(f"    Min: {np.min(opt_g):.4f} at site {np.argmin(opt_g)}")
    log(f"    Max: {np.max(opt_g):.4f} at site {np.argmax(opt_g)}")
    log(f"    Asymmetry: {asym:.4f}")
    if asym > 0.1:
        log("    -> ASYMMETRIC sacrifice-zone confirmed at N=7")
    else:
        log("    -> SYMMETRIC")

    # Check mirror profile (should give same SumMI by chain symmetry)
    mirror_r = evaluate_profile_csharp(list(mirror), N)
    if mirror_r:
        log(f"    Mirror SumMI:  {mirror_r['SumMI']:.6f} (should match original)")
    log()

    log(f"Total runtime: {total_time:.0f}s ({total_time/60:.1f} min)")
    log(f"Results: {OUT_PATH}")
    _outf.close()
