#!/usr/bin/env python3
"""EQ-016 verification: are the triple-Dicke maxima local maxima on the
FULL symmetric Dicke sphere (CP^N), or just on the 3-Dicke slice?

F69 doc claim: "no non-product local maxima on the symmetric Dicke sphere
at any tested N". My triple-Dicke maxima are stationary on the (D_i, D_j,
D_k) slice. On the full sphere, perturbations in OTHER Dicke directions
may ascend (saddle) or descend (local max on full sphere).

Test: take each top triple-Dicke maximum, perturb in EACH of the unused
Dicke directions (and small random perturbations), check if cpsi
increases. If any direction increases cpsi, the maximum is a saddle on
full sphere (consistent with F69 doc). If all directions decrease, it's
a true local max — F69 doc is wrong.
"""
from __future__ import annotations

import math
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from _eq016_n4_full_landscape import (
    dicke_state_vector, reduced_rho_AB, reduced_rho_A,
    pair_cpsi, single_qubit_purity, triple_dicke_optimum,
)


def perturb_unused_direction(N, i, j, k, c_i, c_j, c_k, perturbation_dirs, eps_list):
    """For each unused Dicke index l in {0..N} \\ {i, j, k}, perturb the
    optimal triple by epsilon in the |D_l> direction and check pair-CPsi."""
    Di = dicke_state_vector(N, i)
    Dj = dicke_state_vector(N, j)
    Dk = dicke_state_vector(N, k)
    psi0 = c_i * Di + c_j * Dj + c_k * Dk
    psi0 = psi0 / np.linalg.norm(psi0)
    cpsi0 = pair_cpsi(reduced_rho_AB(psi0, N))
    purity0 = single_qubit_purity(reduced_rho_AB(psi0, N))

    print(f"  Triple D_{i}+D_{j}+D_{k}: c=({c_i:+.4f},{c_j:+.4f},{c_k:+.4f}), "
          f"cpsi₀={cpsi0:.6f}, purity={purity0:.4f}")

    used = {i, j, k}
    unused = [l for l in range(N + 1) if l not in used]
    print(f"    Perturbing in unused Dicke directions: {unused}")

    n_ascent = 0
    n_descent = 0
    for l in unused:
        Dl = dicke_state_vector(N, l)
        for eps in eps_list:
            for sign in [+1, -1]:
                psi_pert = psi0 + sign * eps * Dl
                psi_pert = psi_pert / np.linalg.norm(psi_pert)
                cpsi_p = pair_cpsi(reduced_rho_AB(psi_pert, N))
                if cpsi_p > cpsi0 + 1e-10:
                    n_ascent += 1
                else:
                    n_descent += 1
                if abs(eps - eps_list[0]) < 1e-12 and sign > 0:
                    sign_str = "ASCEND" if cpsi_p > cpsi0 else "descend"
                    print(f"      D_{l}, +ε={eps:.0e}: Δcpsi={cpsi_p-cpsi0:+.4e}  ({sign_str})")
    print(f"    Total: {n_ascent} ascent / {n_descent} descent across all "
          f"unused directions × {len(eps_list)} eps × 2 signs")
    return n_ascent > 0


def random_perturbation_test(N, c_i, c_j, c_k, i, j, k, n_perturb=200, eps=1e-3, seed=0):
    """Random complex perturbation test on full symmetric Dicke sphere."""
    rng = np.random.default_rng(seed)
    Di = dicke_state_vector(N, i)
    Dj = dicke_state_vector(N, j)
    Dk = dicke_state_vector(N, k)
    psi0 = c_i * Di + c_j * Dj + c_k * Dk
    psi0 = psi0 / np.linalg.norm(psi0)
    cpsi0 = pair_cpsi(reduced_rho_AB(psi0, N))

    # Random Dicke perturbations: mix random combinations of all N+1 Dicke states
    n_ascent = 0
    n_descent = 0
    for _ in range(n_perturb):
        # Random complex coefficients on each Dicke
        deltas_re = rng.normal(size=N + 1)
        deltas_im = rng.normal(size=N + 1)
        deltas = deltas_re + 1j * deltas_im
        delta_psi = sum(deltas[m] * dicke_state_vector(N, m) for m in range(N + 1))
        # Project out parallel to psi0
        delta_psi = delta_psi - (psi0.conj() @ delta_psi) * psi0
        norm = np.linalg.norm(delta_psi)
        if norm < 1e-12:
            continue
        delta_psi = delta_psi / norm
        psi_pert = psi0 + eps * delta_psi
        psi_pert = psi_pert / np.linalg.norm(psi_pert)
        cpsi_p = pair_cpsi(reduced_rho_AB(psi_pert, N))
        if cpsi_p > cpsi0 + 1e-10:
            n_ascent += 1
        else:
            n_descent += 1
    return n_ascent, n_descent


def main():
    print("=" * 80)
    print("EQ-016: are triple-Dicke maxima saddles on full symmetric Dicke sphere?")
    print("=" * 80)

    eps_list = [1e-4, 1e-3, 1e-2, 1e-1]

    for N in [3, 4, 5, 6]:
        print(f"\n--- N = {N} (testing top triple-Dicke maxima on CP^{N}) ---")
        # Find top 3 triples by cpsi
        triples = []
        for i, j, k in combinations(range(N + 1), 3):
            r = triple_dicke_optimum(N, i, j, k)
            triples.append(r)
        triples.sort(key=lambda r: -r["cpsi_opt"])
        # Test top 3 unique-by-cpsi
        seen_cpsi = set()
        unique_top = []
        for r in triples:
            key = round(r["cpsi_opt"], 6)
            if key in seen_cpsi:
                continue
            seen_cpsi.add(key)
            unique_top.append(r)
            if len(unique_top) == 3:
                break

        for r in unique_top:
            i, j, k = r["i"], r["j"], r["k"]
            ci, cj, ck = r["c_i"], r["c_j"], r["c_k"]
            print()
            is_saddle = perturb_unused_direction(
                N, i, j, k, ci, cj, ck, None, eps_list,
            )
            n_a, n_d = random_perturbation_test(
                N, ci, cj, ck, i, j, k, n_perturb=300, eps=1e-3,
            )
            verdict = ("SADDLE on full sphere" if is_saddle or n_a > 0
                        else "LOCAL MAX on full sphere")
            print(f"    Random Dicke perturbations: {n_a} ascent / {n_d} descent")
            print(f"    Verdict: {verdict}")


if __name__ == "__main__":
    main()
