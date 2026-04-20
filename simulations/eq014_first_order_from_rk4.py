#!/usr/bin/env python3
"""eq014_first_order_from_rk4.py

Extract the first-order coefficient Σ f_i = lim_{δJ→0} Σ ln(α_i(δJ)) / δJ
from exact RK4 evolution. This is the DIRECT test of whether the closure
law is a first-order theorem (Σ f_i = 0) or a higher-order approximation.

Method: run RK4 at δJ ∈ {0.001, 0.01, 0.1} and check scaling:
 - if Σ ln α / δJ → constant ≠ 0 as δJ → 0:  first-order coefficient is nonzero
 - if Σ ln α / δJ → 0:                      first-order coefficient vanishes
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps

sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, str(Path(__file__).parent))
from eq014_validate_groundtruth import rebuild_L_A, evolve_and_sample
from eq014_step4567_closure import (
    build_initial_states, build_V_L, fit_alpha,
)

N = 7
D = 128
J = 1.0
GAMMA = 0.05
TIMES = np.arange(401) * 0.2
T_FIT = 20.0


def main():
    print("=== EQ-014 first-order coefficient from exact RK4 ===")
    L_A = rebuild_L_A(N, J, GAMMA).tocsr()
    states = build_initial_states(N, D)

    deltas = [0.1, 0.01, 0.001]

    for state_name in list(states.keys())[:3] + ["plus_7"]:
        rho0 = states[state_name]
        for bond in [(0, 1)]:
            V_L = build_V_L(bond, N)
            print(f"\n### State: {state_name}, bond {bond} ###")
            print(f"{'δJ':<10} {'Σ_i ln(α_i)':<20} {'Σ f_i = Σ ln α / δJ':<20}")
            for dJ in deltas:
                L_B = (L_A + dJ * V_L).tocsr()
                # P_A only once but we need consistent grid for fit
                if dJ == deltas[0]:
                    P_A = evolve_and_sample(L_A, rho0, D, N, TIMES, dt_small=0.01)
                P_B = evolve_and_sample(L_B, rho0, D, N, TIMES, dt_small=0.01)
                alpha = fit_alpha(P_A, P_B, TIMES, t_max=T_FIT)
                valid = np.all(np.isfinite(alpha) & (alpha > 0))
                sum_ln = float(np.sum(np.log(alpha))) if valid else float('nan')
                f_coeff = sum_ln / dJ if np.isfinite(sum_ln) else float('nan')
                print(f"{dJ:<10} {sum_ln:+.6e}     {f_coeff:+.6e}")


if __name__ == "__main__":
    main()
