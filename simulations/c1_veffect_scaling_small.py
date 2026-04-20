#!/usr/bin/env python3
"""c1_veffect_scaling_small.py

Rerun of TASK_C1_VEFFECT_SCALING restricted to N=3..6, with incremental
JSON saves after each step so partial progress survives a crash.

Scope: Step 1+2 (endpoint c_1, psi_1+vac, N=3..6), Step 3 (full bond-profile
c_1 vector, N=3..6), Step 4 (endpoint c_1 for psi_2+vac, N=3..6). N=7 is
deferred (the earlier fullrun crashed trying to eigendecompose the d=16384
Liouvillian).

Reuses compute_c1_for_bond / propagate_via_eig from c1_veffect_scaling.py.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))

from c1_veffect_scaling import (
    v_effect,
    bonding_plus_vacuum,
    compute_c1_for_bond,
    DJ_PROBE,
    RESULTS_DIR,
)
from pi_pair_closure_investigation import (
    T_FINAL,
    N_STEPS,
    vacuum_ket,
    single_excitation_mode,
)

N_LIST = [3, 4, 5, 6]
OUT_JSON = RESULTS_DIR / "c1_vs_N_small.json"


def save(data):
    with open(OUT_JSON, "w") as f:
        json.dump(data, f, indent=2, default=str)


def main():
    print("=" * 70)
    print("c_1 vs N scaling (N=3..6), incremental saves")
    print("=" * 70)
    print(f"  Output: {OUT_JSON}")

    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)
    results = {"N_list": N_LIST, "dJ": DJ_PROBE}

    # ------------------------------------------------------------------
    # Step 1+2: endpoint c_1 for psi_1+vac
    # ------------------------------------------------------------------
    print(f"\n{'-'*70}\nStep 1+2: endpoint c_1, psi_1+vac, N={N_LIST}\n{'-'*70}")
    print(f"  {'N':>3} {'V(N)':>8} {'0.5V':>8} {'c_1':>11} "
          f"{'c_2':>10} {'ratio':>8} {'time':>8}")
    step12 = {}
    for N in N_LIST:
        state = bonding_plus_vacuum(N, k=1)
        t0 = time.time()
        r = compute_c1_for_bond(N, 0, state, times, DJ_PROBE)
        dt = time.time() - t0
        V = v_effect(N)
        print(f"  {N:>3d} {V:>8.4f} {V/2:>8.4f} {r['c_1']:>+11.5f} "
              f"{r['c_2_approx']:>+10.4f} {r['c_1']/(V/2):>8.4f} {dt:>7.1f}s")
        step12[N] = {
            "V_N": V,
            "c_1": r["c_1"],
            "c_2_approx": r["c_2_approx"],
            "ratio_c1_over_05V": r["c_1"] / (V / 2),
            "alpha_plus": r["entries"][1]["alpha"],
            "alpha_minus": r["entries"][0]["alpha"],
            "closure_plus": next(e["closure"] for e in r["entries"] if e["dJ"] > 0),
            "closure_minus": next(e["closure"] for e in r["entries"] if e["dJ"] < 0),
            "walltime_s": dt,
        }
    results["step12_psi1_endpoint"] = step12
    save(results)

    # Fits
    Ns = sorted(step12.keys())
    c1s = np.array([step12[N]["c_1"] for N in Ns])
    Vs = np.array([v_effect(N) for N in Ns])
    log_c1 = np.log(np.abs(c1s))
    log_N = np.log(Ns)
    log_V = np.log(Vs)
    p_N, logA_N = np.polyfit(log_N, log_c1, 1)
    p_V, logB_V = np.polyfit(log_V, log_c1, 1)
    print(f"\n  Power-law fits (log-log):")
    print(f"    c_1 ~ A * N^p:    p={p_N:.3f}, A={np.exp(logA_N):.4f}")
    print(f"    c_1 ~ B * V(N)^q: q={p_V:.3f}, B={np.exp(logB_V):.4f}")
    for N in Ns:
        predN = np.exp(logA_N) * N ** p_N
        predV = np.exp(logB_V) * v_effect(N) ** p_V
        print(f"      N={N}: c_1={step12[N]['c_1']:+.5f}, "
              f"pred(N^p)={predN:+.5f}, pred(V^q)={predV:+.5f}")
    results["step12_fits"] = {
        "c_1_vs_N_power": {"p": float(p_N), "A": float(np.exp(logA_N))},
        "c_1_vs_V_power": {"q": float(p_V), "B": float(np.exp(logB_V))},
    }
    save(results)

    # ------------------------------------------------------------------
    # Step 3: full c_1 vector over all bonds, psi_1+vac
    # ------------------------------------------------------------------
    print(f"\n{'-'*70}\nStep 3: full c_1 vector, psi_1+vac, N={N_LIST}\n{'-'*70}")
    step3 = {}
    for N in N_LIST:
        state = bonding_plus_vacuum(N, k=1)
        c1_vec = np.zeros(N - 1)
        t0 = time.time()
        print(f"\n  N={N}:")
        for b in range(N - 1):
            r = compute_c1_for_bond(N, b, state, times, DJ_PROBE)
            c1_vec[b] = r["c_1"]
            print(f"    bond ({b},{b+1}): c_1 = {r['c_1']:+.5f}")
        mirror = float(np.linalg.norm(c1_vec - c1_vec[::-1]))
        dt = time.time() - t0
        endpoint_c1 = float(c1_vec[0])
        interior_c1 = float(c1_vec[(N - 1) // 2])
        ratio = endpoint_c1 / interior_c1 if abs(interior_c1) > 1e-10 else float("nan")
        print(f"    mirror residual: {mirror:.2e}, "
              f"endpoint/interior: {ratio:+.3f}, walltime: {dt:.1f}s")
        step3[N] = {
            "c_1_vector": c1_vec.tolist(),
            "mirror_residual": mirror,
            "endpoint_over_interior": ratio,
            "walltime_s": dt,
        }
        results["step3_bond_profile_psi1"] = step3
        save(results)

    # ------------------------------------------------------------------
    # Step 4: endpoint c_1 for psi_2+vac
    # ------------------------------------------------------------------
    print(f"\n{'-'*70}\nStep 4: endpoint c_1, psi_2+vac, N={N_LIST}\n{'-'*70}")
    print(f"  {'N':>3} {'V(N)':>8} {'0.5V':>8} {'c_1':>11} {'ratio':>8}")
    step4 = {}
    for N in N_LIST:
        state = (vacuum_ket(N) + single_excitation_mode(N, k=2)) / np.sqrt(2.0)
        r = compute_c1_for_bond(N, 0, state, times, DJ_PROBE)
        V = v_effect(N)
        ratio = r["c_1"] / (V / 2) if V != 0 else float("nan")
        print(f"  {N:>3d} {V:>8.4f} {V/2:>8.4f} {r['c_1']:>+11.5f} {ratio:>8.4f}")
        step4[N] = {
            "V_N": V,
            "c_1": r["c_1"],
            "ratio_c1_over_05V": ratio,
        }
        results["step4_psi2_endpoint"] = step4
        save(results)

    print(f"\nFull JSON saved: {OUT_JSON}")


if __name__ == "__main__":
    main()
