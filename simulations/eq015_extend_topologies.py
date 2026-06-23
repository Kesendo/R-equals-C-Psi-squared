#!/usr/bin/env python3
"""EQ-015: cavity-mode-exposure formula at N>=5 and non-chain topologies.

The F64 formula:

    rate(mode k) = 2 gamma_B * |v_k(B)|^2

where v_k is an eigenvector of the single-excitation coherence Liouvillian
L_coh = i H_1 - 2 gamma_B P_B and rate(k) = -Re(lambda_k).

**This is an exact algebraic identity, not a perturbative formula.**
For any L_coh eigenvector v_k with eigenvalue lambda_k:
  L_coh v_k = lambda_k v_k
  => i <v|H_1|v> - 2 gamma_B |v(B)|^2 = lambda_k
  Since H_1 is Hermitian, <v|H_1|v> is real, so:
  -Re(lambda_k) = 2 gamma_B |v(B)|^2.

For an XY (U(1)-conserving) Hamiltonian with single-site Z-dephasing on B,
the Liouvillian preserves the vac<->1exc coherence sector (Z preserves
particle number; vac is a fixed eigenstate). Within this sector, L acts as
the N x N operator L_coh = i H_1 - 2 gamma_B |B><B|. F64 is exact here.

Verified at N=3 (max rel err 1.8 %, chain) and N=4 (ratio 1.0000 +- 0.0003,
chain) — these earlier verifications used H_1 eigenvectors (perturbative
formulation) and the discrepancy reflects O((gamma_B/J)^2) corrections.
The L_coh formulation removes this discrepancy entirely.

EQ-015 extends to:
  - chains at N=5, 6, 7
  - ring, star, Y-junction, K_5 at N=5

Topologies with degenerate H_1 spectra (ring, star, Y, K_5) admit
"B-decoupled rotations": linear combinations within a degenerate subspace
that have |v(B)|^2 = 0 exactly. These are protected modes (rate = 0) and
F64 holds trivially with both sides equal to zero. This is the structurally
interesting case: F64 captures protection.

Output: print verification table; save JSON.
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results" / "eq015_extend_topologies"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def single_exc_H1(bonds, J, N):
    H1 = np.zeros((N, N), dtype=complex)
    for (i, j) in bonds:
        H1[i, j] += J
        H1[j, i] += J
    return H1


def L_coh_matrix(H1, B, gamma_B):
    """L_coh acting on coherence amplitudes: L_coh = i H_1 - 2 gamma_B |B><B|."""
    L = 1j * H1.copy()
    L[B, B] -= 2.0 * gamma_B
    return L


def verify_F64_per_mode(name, bonds, N, S, B, J, gamma_B, eps_aS=1e-9):
    """For each L_coh eigenvector with nonzero a_S, check rate == 2*gamma_B*|v(B)|^2."""
    H1 = single_exc_H1(bonds, J, N)
    L_coh = L_coh_matrix(H1, B, gamma_B)
    evals, V = np.linalg.eig(L_coh)
    # Normalize each eigenvector
    for k in range(V.shape[1]):
        V[:, k] /= np.linalg.norm(V[:, k])

    print(f"\n  {name},  N={N},  S={S}, B={B}, gamma_B={gamma_B}")
    print(f"    bonds = {bonds}")
    print(f"    {'k':>3} {'-Re(lambda)':>14} {'|v(S)|':>10} {'|v(B)|^2':>12} "
          f"{'F64 pred':>14} {'rel err':>12}")

    mode_results = []
    for k in range(N):
        v = V[:, k]
        rate = -evals[k].real
        a_S = abs(v[S])
        a_B2 = abs(v[B]) ** 2
        if a_S < eps_aS:
            continue
        pred = 2.0 * gamma_B * a_B2
        if rate > 1e-15:
            rel_err = abs(rate - pred) / rate
        elif pred < 1e-15:
            rel_err = 0.0
        else:
            rel_err = float('nan')
        mode_results.append({
            "k": k, "rate": rate, "a_S": a_S, "a_B2": a_B2,
            "pred": pred, "rel_err": rel_err,
        })
        print(f"    {k:>3} {rate:>14.6e} {a_S:>10.6f} {a_B2:>12.6f} "
              f"{pred:>14.6e} {rel_err:>12.2e}")

    rel_errs = [abs(r["rel_err"]) for r in mode_results
                  if not np.isnan(r["rel_err"])]
    max_err = max(rel_errs) if rel_errs else float('nan')
    print(f"    -> max rel err over all S-coherence modes: {max_err:.2e}")
    return {
        "name": name, "N": N, "bonds": bonds, "S": S, "B": B,
        "J": J, "gamma_B": gamma_B,
        "modes": mode_results, "max_rel_err": max_err,
    }


def main():
    print("=" * 80)
    print("EQ-015: per-mode F64 verification for chains N>=5 and non-chain topologies")
    print("=" * 80)
    print()
    print("F64 (mode-by-mode): rate(mode k) = 2*gamma_B*|v_k(B)|^2")
    print("  Tests every L_coh eigenmode with nonzero S-overlap.")
    print("  Protected modes (|v(B)|^2 = 0) verify the formula trivially")
    print("  with both sides zero.")

    J = 1.0
    gamma_B = 0.01  # good-cavity regime

    topologies = []
    for N in [5, 6, 7]:
        topologies.append((f"chain N={N}",
                            [(i, i + 1) for i in range(N - 1)], N, 0, N - 1))
    N = 5
    topologies.append(("ring N=5 (S=0,B=2)",
                        [(i, (i + 1) % N) for i in range(N)], N, 0, 2))
    topologies.append(("star N=5 (S=4,B=0 hub)",
                        [(0, i) for i in range(1, N)], N, 4, 0))
    topologies.append(("star N=5 (S=1,B=4 leaf)",
                        [(0, i) for i in range(1, N)], N, 1, 4))
    topologies.append(("Y N=5 (S=0,B=2)",
                        [(0, 1), (1, 2), (1, 3), (1, 4)], N, 0, 2))
    topologies.append(("K_5 N=5 (S=0,B=4)",
                        list(combinations(range(N), 2)), N, 0, 4))

    results = []
    for name, bonds, N, S, B in topologies:
        results.append(verify_F64_per_mode(name, bonds, N, S, B, J, gamma_B))

    # Summary
    print("\n" + "=" * 80)
    print(f"Summary at gamma_B = {gamma_B}")
    print("=" * 80)
    print(f"\n  {'topology':<35} {'#S-modes':>10} {'max rel err':>14}")
    print("  " + "-" * 60)
    for r in results:
        n_modes = len([m for m in r['modes']])
        print(f"  {r['name']:<35} {n_modes:>10d} {r['max_rel_err']:>14.2e}")

    # gamma_B independence: F64 is exact, so rel err stays at machine
    # precision regardless of gamma_B. Demonstrate.
    print()
    print("-" * 80)
    print("gamma_B independence: F64 is EXACT, max rel err independent of gamma_B")
    print("-" * 80)
    bonds5 = [(i, i + 1) for i in range(4)]
    convergence = []
    for gB in [0.5, 0.1, 0.01, 0.001]:
        r = verify_F64_per_mode(f"chain N=5 (gamma_B={gB})", bonds5, 5, 0, 4, J, gB)
        convergence.append(r)
    print()
    print(f"  Even at gamma_B/J = 0.5 (NOT good-cavity), F64 is exact for L_coh:")
    print(f"  {'gamma_B':>10} {'gamma_B/J':>12} {'max rel err':>14}")
    print("  " + "-" * 45)
    for r in convergence:
        gB = r['gamma_B']
        print(f"  {gB:>10.4f} {gB/J:>12.4f} {r['max_rel_err']:>14.2e}")

    out = {
        "J": J, "gamma_B_main": gamma_B,
        "stage1_results": results,
        "stage2_gamma_independence": convergence,
    }
    p = RESULTS_DIR / "verification.json"
    with open(p, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {p}")


if __name__ == "__main__":
    main()
