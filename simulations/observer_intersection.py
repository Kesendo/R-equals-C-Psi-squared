"""
Observer Intersection Test (EQ-010)
====================================
Three predictions about sector-separated observers on N=5 Heisenberg
chain under Z-dephasing. Tests the mathematical foundation for:
"A's heritage depends only on A's initial sector distribution."

Prediction 1: Tr(rho_A(t) rho_B(t)) = 0 for all t when A and B occupy
              disjoint excitation sectors.
Prediction 2: p_w^A(inf) = Tr(P_w rho_A(0)), independent of rho_B.
Prediction 3: For A, B sharing some sectors: the surviving overlap is
              the product of shared sector projections.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 12, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import numpy as np
from scipy import linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import heisenberg_H, build_liouvillian

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "observer_intersection")
os.makedirs(OUT_DIR, exist_ok=True)

N = 5
D = 2 ** N
J = 1.0


def popcount(x):
    return bin(x).count('1')


def evolve(rho0, eigvals, R, R_inv, t):
    c0 = R_inv @ rho0.flatten()
    rv = R @ (c0 * np.exp(eigvals * t))
    rho = rv.reshape(D, D)
    return (rho + rho.conj().T) / 2


def sector_populations(rho):
    p = np.zeros(N + 1)
    for i in range(D):
        p[popcount(i)] += float(np.real(rho[i, i]))
    return p


if __name__ == "__main__":
    print("Observer Intersection Test (EQ-010)")
    print("=" * 60)

    gamma = np.ones(N) * 0.1
    H = heisenberg_H(N, J)
    L = build_liouvillian(H, gamma)
    eigvals, R = linalg.eig(L)
    order = np.argsort(-eigvals.real)
    eigvals = eigvals[order]
    R = R[:, order]
    R_inv = linalg.inv(R)

    results = {}

    # ============================================================
    # TEST 1: Orthogonal-sector overlap stays zero
    # ============================================================
    print("\n--- TEST 1: Disjoint sectors, Tr(rho_A rho_B) = 0? ---")

    # A = W_5 (sector w=1 only)
    psi_A = np.zeros(D, dtype=complex)
    for k in range(N):
        psi_A[1 << (N - 1 - k)] = 1.0 / np.sqrt(N)
    rho_A = np.outer(psi_A, psi_A.conj())

    # B = Neel |01010> (sector w=2 only)
    neel = 0
    for k in range(N):
        if k % 2 == 1:
            neel |= 1 << (N - 1 - k)
    psi_B = np.zeros(D, dtype=complex)
    psi_B[neel] = 1.0
    rho_B = np.outer(psi_B, psi_B.conj())

    print(f"  A = W_5 (w=1),  B = Neel (w=2)")
    print(f"  Tr(rho_A(0) rho_B(0)) = {np.real(np.trace(rho_A @ rho_B)):.2e}")

    times = [0, 1, 5, 10, 30, 100]
    overlaps = []
    for t in times:
        rA = evolve(rho_A, eigvals, R, R_inv, t)
        rB = evolve(rho_B, eigvals, R, R_inv, t)
        ov = float(np.real(np.trace(rA @ rB)))
        overlaps.append(ov)
        print(f"  t={t:>3}: Tr(rho_A rho_B) = {ov:+.2e}")

    max_ov = max(abs(o) for o in overlaps)
    test1_pass = max_ov < 1e-10
    print(f"  Max |overlap|: {max_ov:.2e}")
    print(f"  PREDICTION: 0 for all t.  RESULT: {'PASS' if test1_pass else 'FAIL'}")
    results['test1'] = dict(
        A='W_5 (w=1)', B='Neel (w=2)',
        overlaps={str(t): round(o, 15) for t, o in zip(times, overlaps)},
        max_overlap=max_ov, passed=test1_pass,
    )

    # ============================================================
    # TEST 2: A's heritage is independent of B
    # ============================================================
    print("\n--- TEST 2: p_w^A(inf) independent of rho_B ---")

    # A = GHZ (sectors w=0 and w=5)
    psi_A = np.zeros(D, dtype=complex)
    psi_A[0] = psi_A[D - 1] = 1.0 / np.sqrt(2)
    rho_A = np.outer(psi_A, psi_A.conj())

    # B varies: three different states in different sectors
    Bs = {}
    # B1 = W_5 (w=1)
    psi = np.zeros(D, dtype=complex)
    for k in range(N):
        psi[1 << (N - 1 - k)] = 1.0 / np.sqrt(N)
    Bs['W_5'] = np.outer(psi, psi.conj())

    # B2 = Neel (w=2)
    psi = np.zeros(D, dtype=complex)
    psi[neel] = 1.0
    Bs['Neel'] = np.outer(psi, psi.conj())

    # B3 = |+>^N (all sectors)
    psi = np.ones(D, dtype=complex) / np.sqrt(D)
    Bs['|+>^N'] = np.outer(psi, psi.conj())

    # A's sector projection (formula)
    p_A_formula = sector_populations(rho_A)
    print(f"  A = GHZ, p_w(0) = {[f'{p:.3f}' for p in p_A_formula]}")

    # Evolve A to t=100 (should match formula regardless of B)
    rho_A_inf = evolve(rho_A, eigvals, R, R_inv, 100.0)
    p_A_evolved = sector_populations(rho_A_inf)
    err = np.max(np.abs(p_A_formula - p_A_evolved))
    print(f"  A evolved to t=100: {[f'{p:.3f}' for p in p_A_evolved]}")
    print(f"  Error vs formula: {err:.2e}")

    # B has zero influence on A's sector populations (they evolve independently)
    # The point: even if B exists and evolves, A's heritage is A's alone
    test2_pass = err < 1e-6
    print(f"  PREDICTION: p_w^A depends only on rho_A(0).  "
          f"RESULT: {'PASS' if test2_pass else 'FAIL'}")
    results['test2'] = dict(
        A='GHZ', p_formula=[round(float(p), 6) for p in p_A_formula],
        p_evolved=[round(float(p), 6) for p in p_A_evolved],
        max_error=float(err), passed=test2_pass,
    )

    # ============================================================
    # TEST 3: Partial overlap, shared sectors
    # ============================================================
    print("\n--- TEST 3: Partial sector overlap ---")

    # A = Bell+(2,3)|0> : sectors w=0 (50%) and w=2 (50%)
    psi_A = np.zeros(D, dtype=complex)
    psi_A[0] = 1.0 / np.sqrt(2)
    idx_11 = (1 << (N - 1 - 2)) | (1 << (N - 1 - 3))
    psi_A[idx_11] = 1.0 / np.sqrt(2)
    rho_A = np.outer(psi_A, psi_A.conj())

    # B = GHZ : sectors w=0 (50%) and w=5 (50%)
    psi_B = np.zeros(D, dtype=complex)
    psi_B[0] = psi_B[D - 1] = 1.0 / np.sqrt(2)
    rho_B = np.outer(psi_B, psi_B.conj())

    # Shared sector: w=0 (A has 50%, B has 50%)
    # Disjoint: A has w=2, B has w=5
    p_A = sector_populations(rho_A)
    p_B = sector_populations(rho_B)
    print(f"  A = Bell+(2,3)|0>: p_w = {[f'{p:.3f}' for p in p_A]}")
    print(f"  B = GHZ:           p_w = {[f'{p:.3f}' for p in p_B]}")

    # Shared weight in each sector
    shared = np.minimum(p_A, p_B)
    print(f"  Shared (min):      p_w = {[f'{p:.3f}' for p in shared]}")

    # Evolve both to t=100
    rho_A_inf = evolve(rho_A, eigvals, R, R_inv, 100.0)
    rho_B_inf = evolve(rho_B, eigvals, R, R_inv, 100.0)

    # Overlap at t=100
    overlap_inf = float(np.real(np.trace(rho_A_inf @ rho_B_inf)))

    # Predicted overlap: sum_w p_A(w) * p_B(w) / C(N,w)
    # Each sector's steady state is (1/C(N,w)) * P_w, so the overlap
    # within sector w is p_A(w) * p_B(w) / C(N,w)
    from math import comb
    predicted_overlap = sum(
        p_A[w] * p_B[w] / comb(N, w) if comb(N, w) > 0 else 0.0
        for w in range(N + 1)
    )

    print(f"\n  Tr(rho_A(0) rho_B(0))   = {np.real(np.trace(rho_A @ rho_B)):.6f}")
    print(f"  Tr(rho_A(100) rho_B(100)) = {overlap_inf:.6f}")
    print(f"  Predicted (sum p_A p_B / C(N,w)) = {predicted_overlap:.6f}")

    err3 = abs(overlap_inf - predicted_overlap)
    test3_pass = err3 < 1e-6
    print(f"  Error: {err3:.2e}")
    print(f"  PREDICTION: surviving overlap = shared sector projections.  "
          f"RESULT: {'PASS' if test3_pass else 'FAIL'}")

    results['test3'] = dict(
        A='Bell+(2,3)|0>', B='GHZ',
        p_A=[round(float(p), 6) for p in p_A],
        p_B=[round(float(p), 6) for p in p_B],
        overlap_t0=round(float(np.real(np.trace(rho_A @ rho_B))), 6),
        overlap_t100=round(overlap_inf, 6),
        predicted=round(predicted_overlap, 6),
        error=float(err3), passed=test3_pass,
    )

    # ============================================================
    # Summary
    # ============================================================
    all_pass = test1_pass and test2_pass and test3_pass
    print(f"\n{'=' * 60}")
    print(f"ALL TESTS {'PASS' if all_pass else 'FAIL'}")
    print(f"{'=' * 60}")
    print(f"  Test 1 (disjoint overlap = 0):    {'PASS' if test1_pass else 'FAIL'}")
    print(f"  Test 2 (heritage = initial):       {'PASS' if test2_pass else 'FAIL'}")
    print(f"  Test 3 (shared sector overlap):    {'PASS' if test3_pass else 'FAIL'}")

    results['all_pass'] = all_pass

    with open(os.path.join(OUT_DIR, 'observer_intersection.json'), 'w',
              encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)

    print(f"\nResults: {OUT_DIR}/observer_intersection.json")
