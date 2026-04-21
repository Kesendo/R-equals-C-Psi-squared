#!/usr/bin/env python3
"""eq018_c1_purity_response.py

Alternative closure quantity that is STRICTLY BILINEAR in rho_0.

Definition:
    c_1_pr(N, b, rho_0; t_0) = Sum_i [dP_B(i, t_0) / dJ]_{dJ=0}

Since P_B(i, t) = Tr(rho_B(t, i)^2) is bilinear in rho_0 (the Liouvillian
evolution is linear in rho_0, and purity is bilinear in rho), the
first-order perturbation dP_B/dJ is ALSO bilinear in rho_0. No ratio is
introduced, so the kernel is strictly bilinear.

Tradeoff vs pointwise-c_1:
  - Loses the "alpha rescaling" interpretation of PTF closure.
  - Is NOT scale-invariant: c_1_pr[lambda rho] = lambda^2 c_1_pr[rho].
  - Has a fixed dimension (inverse time), which propagates to K dimensions.
  - On the other hand: kernel is universal (probe-independent).

This script repeats the core EQ-018 tests under c_1_pr at N=5, bond 0,
and compares empirical residuals with the pointwise and LSQ results.

Hypothesis: c_1_pr gives zero w-scan residual (strict bilinearity),
and pure-coherence-operator probe ratio is exactly 1.00.

Rules: UTF-8 stdout. No em-dashes.
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
    vacuum_ket, density_matrix,
    per_site_purity,
)
from c1_bilinearity_test import dicke_state
from eq018_kernel_extract import build_decomps, propagate

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_c1_purity_response"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DJ_EXTRACT = 0.01


def measure_c1_purity_response(rho_0, decomps, times, t_0, N, dj=DJ_EXTRACT):
    """c_1_pr(t_0) = Sum_i dP_B(t_0, i)/dJ via symmetric difference."""
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)
    t_0_idx = int(np.argmin(np.abs(times - t_0)))
    dP_B_dJ = (P_Bp[t_0_idx] - P_Bm[t_0_idx]) / (2 * dj)  # (N,) bilinear in rho_0
    return {"c_1": float(np.sum(dP_B_dJ)), "dP_B_dJ": dP_B_dJ.tolist()}


def extract_K_pr(N, bond, decomps, times, t_0):
    S = {n: dicke_state(N, n) for n in range(N + 1)}
    c_1_pure = {}
    for n in range(N + 1):
        c_1_pure[n] = measure_c1_purity_response(
            density_matrix(S[n]), decomps, times, t_0, N)["c_1"]

    K_DD_cross = {}
    K_CC = {}
    for n in range(N + 1):
        for m in range(n + 1, N + 1):
            rho_mix = 0.5 * (density_matrix(S[n]) + density_matrix(S[m]))
            c_mix = measure_c1_purity_response(rho_mix, decomps, times, t_0, N)["c_1"]
            pred_diag = 0.25 * (c_1_pure[n] + c_1_pure[m])
            K_DD_cross[(n, m)] = 2.0 * (c_mix - pred_diag)
            if m - n == 1:
                v = (S[n] + S[m]) / np.sqrt(2.0)
                c_coh = measure_c1_purity_response(
                    density_matrix(v), decomps, times, t_0, N)["c_1"]
                K_CC[(n, m)] = 2.0 * (c_coh - c_mix)
    return {"c_1_pure": c_1_pure, "K_DD_cross": K_DD_cross, "K_CC": K_CC}


def w_scan(N, bond, decomps, times, t_0, K_ref):
    S0 = dicke_state(N, 0); S1 = dicke_state(N, 1)
    w_vals = np.linspace(0, np.pi / 2, 11)
    rows = []
    for w in w_vals:
        c, s = float(np.cos(w)), float(np.sin(w))
        rho = density_matrix(c * S0 + s * S1)
        meas = measure_c1_purity_response(rho, decomps, times, t_0, N)["c_1"]
        pred = (c**4) * K_ref["DD_00"] + (s**4) * K_ref["DD_11"] \
               + 2 * (c**2) * (s**2) * K_ref["DD_01"] \
               + 2 * (c**2) * (s**2) * K_ref["CC_01"]
        rows.append({"w": float(w), "measured": meas, "predicted": pred,
                     "residual": meas - pred})
    return {"rows": rows, "max_residual": max(abs(r["residual"]) for r in rows)}


def pure_coh_probe(N, bond, decomps, times, t_0, K_CC_01):
    S0 = dicke_state(N, 0); S1 = dicke_state(N, 1)
    op = 0.5 * (np.outer(S0, S1.conj()) + np.outer(S1, S0.conj()))
    r = measure_c1_purity_response(op, decomps, times, t_0, N)
    ratio = r["c_1"] / (K_CC_01 / 2) if abs(K_CC_01) > 1e-10 else None
    return {"c_1_direct": r["c_1"], "K_CC_half": K_CC_01 / 2, "ratio": ratio}


def main():
    start = time.time()
    N, bond = 5, 0
    t_0 = 1.0 / GAMMA_0  # = 20
    print("=" * 78)
    print(f"EQ-018 c_1_purity_response (STRICT BILINEAR variant) at N={N}, bond={bond}, t_0={t_0}")
    print("=" * 78)

    decomps = build_decomps(N, bond)
    times = np.linspace(0, T_FINAL, N_STEPS + 1)

    K = extract_K_pr(N, bond, decomps, times, t_0)
    print(f"\n  K_DD[n,n]:")
    for n in range(N + 1):
        print(f"    K[{n},{n}] = {K['c_1_pure'][n]:+.5e}")
    print(f"\n  K_DD[n,m] cross:")
    for (n, m), v in sorted(K["K_DD_cross"].items()):
        print(f"    K[{n},{m}] = {v:+.5e}  (dN={m-n})")
    print(f"\n  K_CC[n,n+1]:")
    for (n, m), v in sorted(K["K_CC"].items()):
        print(f"    K_CC[{n},{m}] = {v:+.5e}")

    K_ref = {
        "DD_00": K["c_1_pure"][0], "DD_11": K["c_1_pure"][1],
        "DD_01": K["K_DD_cross"][(0, 1)], "CC_01": K["K_CC"][(0, 1)],
    }

    print(f"\n  w-scan bilinearity (predictions from 50/50 K above):")
    scan = w_scan(N, bond, decomps, times, t_0, K_ref)
    print(f"  {'w/pi':>8} {'measured':>14} {'predicted':>14} {'residual':>12}")
    for row in scan["rows"]:
        print(f"  {row['w']/np.pi:>8.4f} {row['measured']:>+14.5e}"
              f" {row['predicted']:>+14.5e} {row['residual']:>+12.2e}")
    print(f"  max |residual| = {scan['max_residual']:.2e}")

    pcp = pure_coh_probe(N, bond, decomps, times, t_0, K_ref["CC_01"])
    print(f"\n  Pure-coherence-operator probe:")
    print(f"    K_CC/2                = {K_ref['CC_01']/2:+.5e}")
    print(f"    c_1_pr[rho_coh_only] = {pcp['c_1_direct']:+.5e}")
    if pcp["ratio"] is not None:
        print(f"    ratio                = {pcp['ratio']:.6f}")
    else:
        print(f"    ratio                = undefined (K_CC too small)")

    out = {
        "N": N, "bond": bond, "t_0": t_0,
        "K_DD_diag": {str(n): K["c_1_pure"][n] for n in range(N + 1)},
        "K_DD_cross": {f"{n}_{m}": v for (n, m), v in K["K_DD_cross"].items()},
        "K_CC": {f"{n}_{m}": v for (n, m), v in K["K_CC"].items()},
        "wscan": scan,
        "pure_coh_probe": pcp,
    }
    path = RESULTS_DIR / f"c1_purity_response_t0_{t_0}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")
    print(f"Walltime: {time.time() - start:.1f} s")


def main_multi_t0():
    """Run purity_response at multiple t_0 to see how K values scale."""
    start = time.time()
    N, bond = 5, 0
    print("=" * 78)
    print(f"EQ-018 c_1_purity_response multi-t_0 at N={N}, bond={bond}")
    print("=" * 78)

    decomps = build_decomps(N, bond)
    times = np.linspace(0, T_FINAL, N_STEPS + 1)

    for t_0 in [2.0, 5.0, 10.0, 20.0]:
        print(f"\n--- t_0 = {t_0} ---")
        K = extract_K_pr(N, bond, decomps, times, t_0)
        K_ref = {
            "DD_00": K["c_1_pure"][0], "DD_11": K["c_1_pure"][1],
            "DD_01": K["K_DD_cross"][(0, 1)], "CC_01": K["K_CC"][(0, 1)],
        }
        print(f"  K_DD[1,1] = {K_ref['DD_11']:+.5e}")
        print(f"  K_DD[2,2] = {K['c_1_pure'][2]:+.5e}")
        print(f"  K_DD[0,1] = {K_ref['DD_01']:+.5e}")
        print(f"  K_CC[0,1] = {K_ref['CC_01']:+.5e}")

        scan = w_scan(N, bond, decomps, times, t_0, K_ref)
        print(f"  w-scan max residual = {scan['max_residual']:.2e}")

        pcp = pure_coh_probe(N, bond, decomps, times, t_0, K_ref["CC_01"])
        if pcp["ratio"] is not None:
            print(f"  pcp ratio = {pcp['ratio']:.6f}")
        else:
            print(f"  pcp ratio = undefined (K_CC tiny)")

    print(f"\nWalltime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "multi":
        main_multi_t0()
    else:
        main()
