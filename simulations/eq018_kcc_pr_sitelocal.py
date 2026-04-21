#!/usr/bin/env python3
"""eq018_kcc_pr_sitelocal.py

Site-local spectator probe comparison for K_CC[n, n+1]_pr.

The Dicke-style probe rho_cc_init = (|S_n><S_{n+1}| + h.c.)/2 contains
BOTH 1-site-differing and multi-site-differing computational-basis-pair
content (see eq018_kcc_pr_extension.py docstring). This script uses the
pure 1-site-differing probe

    rho_cc_site = (|0_i s><1_i s| + h.c.) / 2

where s is a fixed computational-basis spectator pattern of popcount n
on the N-1 sites other than i. This probe lives exactly in the
1-site-differing slice of the (n, n+1) sector block: no 3-site-differing
amplitude at t=0.

Comparison with the Dicke probe isolates:
  - If K_pr(Dicke) and K_pr(site-local) agree in sign and order of
    magnitude: the J-dependence is primarily "1-site-diff slice amplitude
    being redistributed inside the slice, measured asymmetrically" or
    "genuine 1-site-to-3-site leak".
  - If K_pr(site-local) is much smaller: the Dicke probe's K is driven by
    its 3-site-diff content feeding back into visible 1-site-diff under
    H, not by H-induced leak from 1-site-diff.
  - If K_pr(site-local) is much larger: Dicke's symmetric spread averages
    away asymmetric local responses.

Parameters mirror eq018_kcc_pr_extension.py. Default N = 5, all bonds,
two site-local configurations per (n, n+1):
  - site = 0 (boundary), spectator = n excitations packed leftmost on
    other sites
  - site = 2 (interior), same spectator packing
  - plus: site = 0, spectator = n excitations RIGHTmost (tests spectator
    placement sensitivity)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from eq018_kcc_pr_extension import (
    build_H_XY, build_liouvillian,
    rho_cc_site_local, rho_cc_dicke,
    spatial_sum_coh_purity, propagate,
    GAMMA_0, J_UNIFORM, DJ, T_MAX, N_TIMES, T_POINTWISE,
    S_zero_closed_form,
)

RESULTS_DIR = Path(__file__).parent / "results" / "eq018_kcc_pr_extension"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def spectator_pattern(n_excitations, n_qubits_spectator, packing="left"):
    """Generate a weight-n_excitations bit pattern on n_qubits_spectator sites.
    packing: 'left' (excitations on leftmost sites) or 'right' (rightmost).
    Returns a list of 0/1 values of length n_qubits_spectator.
    """
    bits = [0] * n_qubits_spectator
    if packing == "left":
        for k in range(n_excitations):
            bits[k] = 1
    elif packing == "right":
        for k in range(n_excitations):
            bits[n_qubits_spectator - 1 - k] = 1
    else:
        raise ValueError(f"unknown packing {packing}")
    return bits


def extract_S_K(rho_init, L_A, L_Bp, L_Bm, times, n_qubits, dJ):
    S_A = np.array([spatial_sum_coh_purity(rho, n_qubits)
                    for rho in propagate(L_A, rho_init, times)])
    S_Bp = np.array([spatial_sum_coh_purity(rho, n_qubits)
                     for rho in propagate(L_Bp, rho_init, times)])
    S_Bm = np.array([spatial_sum_coh_purity(rho, n_qubits)
                     for rho in propagate(L_Bm, rho_init, times)])
    K = (S_Bp - S_Bm) / (2.0 * dJ)
    return S_A, S_Bp, S_Bm, K


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=5)
    parser.add_argument("--bonds", type=str, default=None)
    parser.add_argument("--n-values", type=str, default=None)
    parser.add_argument("--gamma-0", type=float, default=GAMMA_0)
    parser.add_argument("--J", type=float, default=J_UNIFORM)
    parser.add_argument("--dJ", type=float, default=DJ)
    parser.add_argument("--t-max", type=float, default=T_MAX)
    parser.add_argument("--n-times", type=int, default=N_TIMES)
    parser.add_argument("--tag", type=str, default="sitelocal_n5")
    args = parser.parse_args()

    N = args.N
    bonds = list(range(N - 1)) if args.bonds is None else [int(x) for x in args.bonds.split(",")]
    if args.n_values is None:
        # focus on interior blocks (skip n=0 and n=N-1 which are F73 + particle-hole mirror)
        n_values = list(range(1, N - 1))
    else:
        n_values = [int(x) for x in args.n_values.split(",")]

    times = np.linspace(0.0, args.t_max, args.n_times)
    idx_pw = int(np.argmin(np.abs(times - T_POINTWISE)))
    t_pw = times[idx_pw]

    print("=" * 78)
    print(f"Site-local spectator probe comparison at N = {N}")
    print(f"  gamma_0 = {args.gamma_0}, J = {args.J}, dJ = {args.dJ}")
    print(f"  Pointwise t = {t_pw}")
    print("=" * 78)

    t_start = time.time()
    J_A = [args.J] * (N - 1)
    L_A = build_liouvillian(build_H_XY(J_A, N), args.gamma_0, N)
    print(f"L_A built in {time.time() - t_start:.1f} s")

    # For each bond, for each (n, n+1), run:
    #   - Dicke probe (reference)
    #   - site-local probe with spectator "left" at site 0
    #   - site-local probe with spectator "left" at site 2 (interior)
    #   - site-local probe with spectator "right" at site 0

    all_out = {}
    for bond in bonds:
        print(f"\n{'-' * 78}")
        print(f"Bond b = {bond}")
        print(f"{'-' * 78}")
        J_Bp = list(J_A); J_Bp[bond] = args.J + args.dJ
        J_Bm = list(J_A); J_Bm[bond] = args.J - args.dJ
        L_Bp = build_liouvillian(build_H_XY(J_Bp, N), args.gamma_0, N)
        L_Bm = build_liouvillian(build_H_XY(J_Bm, N), args.gamma_0, N)

        bond_out = {}
        for n_coh in n_values:
            print(f"\n  (n, n+1) = ({n_coh}, {n_coh+1}):")

            # Dicke probe (for reference)
            rho_dicke = rho_cc_dicke(n_coh, N)
            S_A_d, _, _, K_d = extract_S_K(rho_dicke, L_A, L_Bp, L_Bm, times, N, args.dJ)
            print(f"    Dicke probe:")
            print(f"      S(0) = {S_A_d[0]:.4f}   S(t={t_pw:.0f}) = {S_A_d[idx_pw]:.4e}   "
                  f"K(t={t_pw:.0f}) = {K_d[idx_pw]:+.4e}")

            # Site-local probes
            configs = [
                ("site0_specLeft", 0, spectator_pattern(n_coh, N - 1, "left")),
                ("site0_specRight", 0, spectator_pattern(n_coh, N - 1, "right")),
            ]
            if N >= 5:
                configs.append(
                    ("site2_specLeft", 2, spectator_pattern(n_coh, N - 1, "left"))
                )

            probe_results = {"dicke": {
                "S_A": S_A_d.tolist(), "K": K_d.tolist(),
                "S_A_pointwise": float(S_A_d[idx_pw]),
                "K_pointwise": float(K_d[idx_pw]),
                "S_0": float(S_A_d[0]),
                "S_0_closed_form_dicke": S_zero_closed_form(n_coh, N),
            }}

            for label, site_i, spec in configs:
                rho_sl = rho_cc_site_local(site_i, spec, N)
                S_A, _, _, K = extract_S_K(rho_sl, L_A, L_Bp, L_Bm, times, N, args.dJ)
                spec_str = "".join(str(b) for b in spec)
                print(f"    site-local {label} (site={site_i}, spec={spec_str}):")
                print(f"      S(0) = {S_A[0]:.4f}   S(t={t_pw:.0f}) = {S_A[idx_pw]:.4e}   "
                      f"K(t={t_pw:.0f}) = {K[idx_pw]:+.4e}")

                probe_results[label] = {
                    "site": site_i,
                    "spectator": spec,
                    "S_A": S_A.tolist(),
                    "K": K.tolist(),
                    "S_A_pointwise": float(S_A[idx_pw]),
                    "K_pointwise": float(K[idx_pw]),
                    "S_0": float(S_A[0]),
                }

            bond_out[f"n_{n_coh}"] = probe_results
        all_out[f"bond_{bond}"] = bond_out

    result = {
        "config": {
            "N": N, "gamma_0": args.gamma_0, "J": args.J, "dJ": args.dJ,
            "t_max": args.t_max, "n_times": args.n_times,
            "T_pointwise": float(t_pw),
            "bonds": bonds, "n_values": n_values,
        },
        "times": times.tolist(),
        "bonds": all_out,
    }

    out_path = RESULTS_DIR / f"kcc_pr_{args.tag}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved: {out_path}")
    print(f"Total walltime: {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
