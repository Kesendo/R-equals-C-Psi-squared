#!/usr/bin/env python3
"""eq018_c1_pointwise.py

Prototype c_1_pointwise: pointwise alpha match at a fixed reference time t_0,
in place of the LSQ alpha fit.

Definition
----------
At first order in dJ:

    P_B(i, t_0) = P_A(i, t_0) + dP_B/dJ(t_0, i) * dJ + O(dJ^2)
    delta_alpha_i = delta_P_B(t_0, i) / (t_0 * dP_A/dt(t_0, i))
    c_1_pw = Sum_i delta_alpha_i / dJ
          = Sum_i (dP_B/dJ)(t_0, i) / (t_0 * dP_A/dt(t_0, i))

Numerically via symmetric difference (dJ = +-0.01):

    dP_B/dJ_numerical = (P_B+(t_0) - P_B-(t_0)) / (2 dJ)

Expectation from Tom: pointwise gives cleaner K that respects bilinearity
across the w-scan, and directly verifies F72 (pure coherence probe should
give K_CC/2 exactly).

Tests
-----
1. K tables at N=5, bond 0 under pointwise, compared to LSQ.
2. w-scan: does pointwise match bilinear prediction at all w?
3. Pure-coherence-operator probe: is c_1_pw = K_CC[0,1]/2 exactly?
4. t_0 sensitivity sweep: {5, 10, 20, 40} at N=5 bond 0. How does K drift?
5. Stretch: N=4, 6 at bond 0 for N-scaling of K_CC[0,1].

Rules: UTF-8 stdout. No em-dashes. Hyphens only. sys.stdout.reconfigure on win32.
"""
from __future__ import annotations

import json
import sys
import time
from math import comb
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

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_c1_pointwise"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


DJ_EXTRACT = 0.01


def measure_c1_pointwise(rho_0, decomps, times, t_0, N,
                         dj=DJ_EXTRACT, return_per_site=False):
    """Pointwise c_1 at fixed reference time t_0.

    Returns dict with c_1 and diagnostics. Implements the symmetric-
    difference version.
    """
    (evA, VRA, VLA), (evBp, VRBp, VLBp), (evBm, VRBm, VLBm) = decomps
    rho_A = propagate(evA, VRA, VLA, rho_0, times)
    rho_Bp = propagate(evBp, VRBp, VLBp, rho_0, times)
    rho_Bm = propagate(evBm, VRBm, VLBm, rho_0, times)
    P_A = per_site_purity(rho_A, N)  # (T, N)
    P_Bp = per_site_purity(rho_Bp, N)
    P_Bm = per_site_purity(rho_Bm, N)

    t_0_idx = int(np.argmin(np.abs(times - t_0)))
    if t_0_idx == 0 or t_0_idx == len(times) - 1:
        raise ValueError(f"t_0 = {t_0} too close to edge of grid "
                         f"[{times[0]}, {times[-1]}]")
    dt = times[1] - times[0]

    # dP_A/dt at t_0 via central differences on the unperturbed trajectory
    dP_A_dt_t0 = (P_A[t_0_idx + 1] - P_A[t_0_idx - 1]) / (2 * dt)  # (N,)
    # Symmetric-difference (d/d dJ) P_B at t_0
    dP_B_dJ_t0 = (P_Bp[t_0_idx] - P_Bm[t_0_idx]) / (2 * dj)  # (N,)

    # per-site delta_alpha/dJ, avoiding divide-by-zero
    eps = 1e-14
    denom = t_0 * dP_A_dt_t0
    zero_mask = np.abs(denom) < eps
    # For nearly-stationary sites (vac only), dP_A/dt = 0 and the pointwise
    # rule is undefined. Return NaN there so we can see it in diagnostics.
    dalpha = np.where(zero_mask, np.nan, dP_B_dJ_t0 / denom)

    # Sum over sites ignoring NaNs (should be rare for generic rho_0)
    c_1 = float(np.nansum(dalpha))
    n_nan_sites = int(np.sum(zero_mask))

    out = {
        "c_1": c_1,
        "t_0": float(t_0),
        "t_0_idx": t_0_idx,
        "n_nan_sites": n_nan_sites,
    }
    if return_per_site:
        out["P_A_t0"] = P_A[t_0_idx].tolist()
        out["dP_A_dt_t0"] = dP_A_dt_t0.tolist()
        out["dP_B_dJ_t0"] = dP_B_dJ_t0.tolist()
        out["dalpha_per_site"] = [float(x) if np.isfinite(x) else None
                                   for x in dalpha]
    return out


def extract_K_pointwise(N, bond, decomps, times, t_0):
    """Extract K_DD[n, n], K_DD[n, m] and K_CC[n, n+1] under pointwise c_1."""
    S = {n: dicke_state(N, n) for n in range(N + 1)}

    c_1_pure = {}
    for n in range(N + 1):
        r = measure_c1_pointwise(density_matrix(S[n]), decomps, times, t_0, N)
        c_1_pure[n] = r["c_1"]

    K_DD_cross = {}
    K_CC = {}
    for n in range(N + 1):
        for m in range(n + 1, N + 1):
            rho_mix = 0.5 * (density_matrix(S[n]) + density_matrix(S[m]))
            c_mix = measure_c1_pointwise(rho_mix, decomps, times, t_0, N)["c_1"]
            pred_diag = 0.25 * (c_1_pure[n] + c_1_pure[m])
            K_DD_cross[(n, m)] = 2.0 * (c_mix - pred_diag)

            if m - n == 1:
                v = (S[n] + S[m]) / np.sqrt(2.0)
                rho_coh = density_matrix(v)
                c_coh = measure_c1_pointwise(rho_coh, decomps, times, t_0, N)["c_1"]
                K_CC[(n, m)] = 2.0 * (c_coh - c_mix)

    return {
        "c_1_pure": c_1_pure,
        "K_DD_cross": K_DD_cross,
        "K_CC": K_CC,
    }


def w_scan_bilinearity(N, bond, decomps, times, t_0, K_ref):
    """Re-run the w-scan check at pointwise c_1.

    K_ref = {"DD_00", "DD_11", "DD_01", "CC_01"} from pointwise at 50/50.
    Prediction: c_1(w) = cos^4 * K_DD_00 + sin^4 * K_DD_11
                         + 2 cos^2 sin^2 * K_DD_01 + 2 cos^2 sin^2 * K_CC_01
    """
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)
    w_values = np.linspace(0, np.pi / 2, 11)
    rows = []
    for w in w_values:
        c = float(np.cos(w)); s = float(np.sin(w))
        psi = c * S0 + s * S1
        rho = density_matrix(psi)
        meas = measure_c1_pointwise(rho, decomps, times, t_0, N)["c_1"]
        pred = (c**4) * K_ref["DD_00"] + (s**4) * K_ref["DD_11"] \
               + 2 * (c**2) * (s**2) * K_ref["DD_01"] \
               + 2 * (c**2) * (s**2) * K_ref["CC_01"]
        rows.append({"w": float(w), "c_1_measured": meas, "c_1_predicted": pred,
                     "residual": meas - pred})
    max_residual = max(abs(r["residual"]) for r in rows)
    return {"rows": rows, "max_residual": max_residual}


def pure_coherence_probe(N, bond, decomps, times, t_0, K_CC_01):
    """Pure coherence operator probe.

    rho_coh_only = (|vac><S_1| + |S_1><vac|)/2 (Hermitian, trace 0).

    If c_1_pw is strictly bilinear in rho_0, c_1_pw[rho_coh_only] = K_CC[0,1]/2.
    """
    S0 = dicke_state(N, 0)
    S1 = dicke_state(N, 1)
    rho_coh_only = 0.5 * (np.outer(S0, S1.conj()) + np.outer(S1, S0.conj()))
    r = measure_c1_pointwise(rho_coh_only, decomps, times, t_0, N,
                              return_per_site=True)
    ratio = r["c_1"] / (K_CC_01 / 2) if abs(K_CC_01) > 1e-10 else None
    return {
        "c_1_direct": r["c_1"],
        "K_CC_half": K_CC_01 / 2 if abs(K_CC_01) > 1e-10 else None,
        "ratio": ratio,
        "details": r,
    }


def main():
    start = time.time()
    print("=" * 78)
    print("EQ-018 c_1_pointwise prototype at N=5, bond=0")
    print("=" * 78)
    print(f"  gamma_0 = {GAMMA_0}, J = {J_UNIFORM}, dJ = +/- {DJ_EXTRACT}")
    print(f"  T_final = {T_FINAL}, N_steps = {N_STEPS}")

    N = 5
    bond = 0
    t_0_default = 1.0 / GAMMA_0  # = 20.0
    t_0_list = [5.0, 10.0, 20.0, 40.0]

    t0_build = time.time()
    decomps = build_decomps(N, bond)
    print(f"  L_A, L_B+/- built in {time.time() - t0_build:.1f} s")
    times = np.linspace(0.0, T_FINAL, N_STEPS + 1)

    # ----- Primary test at t_0 = 1/gamma_0 = 20 -----
    print(f"\n--- Pointwise at t_0 = {t_0_default} (1/gamma_0) ---")
    K = extract_K_pointwise(N, bond, decomps, times, t_0_default)
    print(f"  K_DD[n,n] (pure Dicke probe):")
    for n in range(N + 1):
        print(f"    K[{n},{n}] = {K['c_1_pure'][n]:+.5f}")
    print(f"\n  K_DD[n,m] cross (reliable only; full list):")
    for (n, m), v in sorted(K["K_DD_cross"].items()):
        print(f"    K[{n},{m}] = {v:+.5f}  (dN={m-n})")
    print(f"\n  K_CC[n,n+1]:")
    for (n, m), v in sorted(K["K_CC"].items()):
        print(f"    K_CC[{n},{m}] = {v:+.5f}")

    K_ref = {
        "DD_00": K["c_1_pure"][0],
        "DD_11": K["c_1_pure"][1],
        "DD_01": K["K_DD_cross"][(0, 1)],
        "CC_01": K["K_CC"][(0, 1)],
    }
    # ----- Bilinearity check (w-scan) -----
    print(f"\n  w-scan bilinearity (prediction from above K):")
    scan = w_scan_bilinearity(N, bond, decomps, times, t_0_default, K_ref)
    print(f"  {'w/pi':>8} {'c_1 meas':>12} {'c_1 pred':>12} {'residual':>12}")
    for row in scan["rows"]:
        print(f"  {row['w']/np.pi:>8.4f} {row['c_1_measured']:>+12.5f}"
              f" {row['c_1_predicted']:>+12.5f} {row['residual']:>+12.2e}")
    print(f"  max |residual| = {scan['max_residual']:.2e}")

    # ----- Pure-coherence probe -----
    print(f"\n  Pure-coherence-operator probe:")
    print(f"    K_CC[0,1]/2 from 50/50 extraction = {K_ref['CC_01']/2:+.5f}")
    pcp = pure_coherence_probe(N, bond, decomps, times, t_0_default,
                                K_ref["CC_01"])
    print(f"    c_1[rho_coh_only] direct         = {pcp['c_1_direct']:+.5f}")
    if pcp['ratio'] is not None:
        print(f"    ratio (direct / (K_CC/2))        = {pcp['ratio']:.6f}")
        print(f"    (1.0 = F72 direct; else extent of deviation)")

    # ----- t_0 sensitivity sweep -----
    print(f"\n--- t_0 sensitivity sweep at N=5, bond=0 ---")
    t0_sweep = {}
    for t_0 in t_0_list:
        Ksw = extract_K_pointwise(N, bond, decomps, times, t_0)
        pcp_sw = pure_coherence_probe(N, bond, decomps, times, t_0,
                                       Ksw["K_CC"][(0, 1)])
        Ksw_K_ref = {
            "DD_00": Ksw["c_1_pure"][0],
            "DD_11": Ksw["c_1_pure"][1],
            "DD_01": Ksw["K_DD_cross"][(0, 1)],
            "CC_01": Ksw["K_CC"][(0, 1)],
        }
        scan_sw = w_scan_bilinearity(N, bond, decomps, times, t_0, Ksw_K_ref)
        t0_sweep[t_0] = {
            "K_DD_diag_1": Ksw["c_1_pure"][1],
            "K_DD_diag_2": Ksw["c_1_pure"][2],
            "K_DD_cross_01": Ksw["K_DD_cross"][(0, 1)],
            "K_DD_cross_02": Ksw["K_DD_cross"][(0, 2)],
            "K_CC_01": Ksw["K_CC"][(0, 1)],
            "pcp_c_1_direct": pcp_sw["c_1_direct"],
            "pcp_ratio": pcp_sw["ratio"],
            "wscan_max_residual": scan_sw["max_residual"],
        }
        print(f"  t_0={t_0:>5.1f}: K_DD[1,1]={Ksw['c_1_pure'][1]:+.4f}"
              f"  K_DD[0,1]={Ksw['K_DD_cross'][(0, 1)]:+.4f}"
              f"  K_CC[0,1]={Ksw['K_CC'][(0, 1)]:+.4f}"
              f"  pcp_ratio={pcp_sw['ratio']:.4f}"
              f"  wscan_max_res={scan_sw['max_residual']:.2e}")

    # ----- Stretch: N=4, 6 at bond=0, just K_CC[0,1] under pointwise -----
    print(f"\n--- Stretch: N-scaling of K_CC[0,1] at bond=0 under pointwise (t_0=20) ---")
    N_scan = {}
    for N_ in [4, 6]:
        if N_ == 4:
            t_0_use = 20.0
        else:
            t_0_use = 20.0
        t1 = time.time()
        decs = build_decomps(N_, 0)
        # Reuse times (T_FINAL=80, N_STEPS=200) regardless of N
        Kn = extract_K_pointwise(N_, 0, decs, times, t_0_use)
        K_CC_01_n = Kn["K_CC"][(0, 1)]
        print(f"  N={N_}: K_CC[0,1]={K_CC_01_n:+.4f}  walltime={time.time()-t1:.1f}s")
        N_scan[N_] = {
            "K_CC_01": K_CC_01_n,
            "c_1_pure_S1": Kn["c_1_pure"][1],
            "K_DD_diag_values": Kn["c_1_pure"],
            "t_0": t_0_use,
        }

    # ----- Save everything -----
    out = {
        "N": N, "bond": bond, "gamma_0": GAMMA_0, "J": J_UNIFORM, "dJ": DJ_EXTRACT,
        "t_0_default": t_0_default,
        "primary": {
            "c_1_pure": K["c_1_pure"],
            "K_DD_cross": {f"{n}_{m}": v for (n, m), v in K["K_DD_cross"].items()},
            "K_CC": {f"{n}_{m}": v for (n, m), v in K["K_CC"].items()},
            "wscan": scan,
            "pure_coherence_probe": pcp,
        },
        "t_0_sweep": {f"{t_0}": d for t_0, d in t0_sweep.items()},
        "N_scan": N_scan,
    }
    path = RESULTS_DIR / "c1_pointwise.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {path}")

    # Summary
    print(f"\n{'=' * 78}")
    print(f"SUMMARY")
    print(f"{'=' * 78}")
    print(f"  Primary at t_0 = {t_0_default}:")
    print(f"    K_DD[1,1] = {K_ref['DD_11']:+.5f}")
    print(f"    K_DD[0,1] = {K_ref['DD_01']:+.5f}")
    print(f"    K_CC[0,1] = {K_ref['CC_01']:+.5f}")
    print(f"    w-scan max residual = {scan['max_residual']:.2e}")
    if pcp["ratio"] is not None:
        print(f"    pure-coh / (K_CC/2) ratio = {pcp['ratio']:.4f}")
    print(f"  t_0 sensitivity: K_CC[0,1] ranges from "
          f"{min(t0_sweep[t]['K_CC_01'] for t in t_0_list):+.3f} to "
          f"{max(t0_sweep[t]['K_CC_01'] for t in t_0_list):+.3f}")
    print(f"  N-scaling K_CC[0,1]: "
          f"N=4:{N_scan[4]['K_CC_01']:+.3f}, "
          f"N=5:{K['K_CC'][(0, 1)]:+.3f}, "
          f"N=6:{N_scan[6]['K_CC_01']:+.3f}")
    print(f"\n  Total walltime: {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
