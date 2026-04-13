"""
U(1) Breaking: How robust is sector decoupling? (EQ-002)
=========================================================
Adds a transverse field ε Σ X_k to the Heisenberg Hamiltonian.
This breaks U(1) excitation-number conservation. Measures:

1. SE fraction of the slow mode as function of ε
2. Bell-pair dependence of c_slow (should appear at ε > 0)
3. Scaling: linear or quadratic in ε?

At ε = 0: sectors exact, c_slow independent of Bell pair (EQ-001, proven).
At ε > 0: sectors leak, coupling emerges.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 13, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import numpy as np
from scipy import linalg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H, build_liouvillian, partial_trace_to_pair,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "u1_breaking")
os.makedirs(OUT_DIR, exist_ok=True)

N = 5
D = 2 ** N
J = 1.0


def popcount(x):
    return bin(x).count('1')


def transverse_field(N):
    """Build Σ X_k on N qubits."""
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    H_tf = np.zeros((2**N, 2**N), dtype=complex)
    for k in range(N):
        ops = [I2] * N
        ops[k] = X
        term = ops[0]
        for op in ops[1:]:
            term = np.kron(term, op)
        H_tf += term
    return H_tf


def find_slow_se_mode_full(eigvals, R, R_inv, N):
    """Find the mode most localized in the SE sector (from full Liouvillian).

    Returns (index, eigenvalue, SE_fraction, left_eigvec, right_eigvec).
    """
    D = 2 ** N
    SE_IDX = [1 << (N - 1 - k) for k in range(N)]
    rates = -eigvals.real
    nonstat = np.where(rates > 1e-10)[0]
    order = np.argsort(rates[nonstat])

    for rank in range(min(50, len(order))):
        i = nonstat[order[rank]]
        V = R[:, i].reshape(D, D)
        block_sq = sum(abs(V[SE_IDX[m], SE_IDX[n]])**2
                       for m in range(N) for n in range(N))
        full_sq = np.sum(np.abs(V)**2)
        se_frac = np.sqrt(block_sq / full_sq) if full_sq > 1e-15 else 0.0
        if se_frac > 0.05:
            left_vec = R_inv[i, :]
            return i, eigvals[i], float(se_frac), left_vec, R[:, i]
    return None, None, 0.0, None, None


def track_mode(eigvals, R, R_inv, ref_right_vec):
    """Find the mode at current ε with highest overlap to reference mode.

    Uses right eigenvector overlap for tracking.
    Returns (index, eigenvalue, SE_fraction, left_eigvec).
    """
    D2 = len(eigvals)
    D = int(np.sqrt(D2))
    N_q = int(np.log2(D))
    SE_IDX = [1 << (N_q - 1 - k) for k in range(N_q)]

    # Overlap of each mode's right eigenvector with reference
    overlaps = np.abs(R.conj().T @ ref_right_vec)
    best = np.argmax(overlaps)

    V = R[:, best].reshape(D, D)
    block_sq = sum(abs(V[SE_IDX[m], SE_IDX[n]])**2
                   for m in range(N_q) for n in range(N_q))
    full_sq = np.sum(np.abs(V)**2)
    se_frac = np.sqrt(block_sq / full_sq) if full_sq > 1e-15 else 0.0

    left_vec = R_inv[best, :]
    return best, eigvals[best], float(se_frac), left_vec


def make_bell_exc_state(N, bell_a, bell_b, exc_k):
    """Bell+(a,b) x |1>_k x |0>_rest."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    i0 = 1 << (N - 1 - exc_k)
    i1 = (1 << (N - 1 - bell_a)) | (1 << (N - 1 - bell_b)) | (1 << (N - 1 - exc_k))
    psi[i0] = 1.0 / np.sqrt(2)
    psi[i1] = 1.0 / np.sqrt(2)
    return psi


if __name__ == "__main__":
    print("U(1) Breaking: Transverse Field Perturbation (EQ-002)")
    print("=" * 60)

    gamma = np.ones(N) * 0.1
    H_heis = heisenberg_H(N, J)
    H_tf = transverse_field(N)

    epsilons = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    results = {}
    ref_right_vec = None  # reference eigenvector for mode tracking

    print(f"\n--- SE fraction of slow mode vs ε (mode-tracked) ---\n")
    print(f"  {'ε':>8}  {'SE frac':>8}  {'Rate':>10}  {'Im':>10}  {'Overlap':>8}")
    print(f"  {'-' * 50}")

    for eps in epsilons:
        H = H_heis + eps * H_tf
        L = build_liouvillian(H, gamma)
        eigvals, R = linalg.eig(L)
        order = np.argsort(-eigvals.real)
        eigvals = eigvals[order]
        R = R[:, order]
        R_inv = linalg.inv(R)

        if ref_right_vec is None:
            # ε = 0: find the SE slow mode as reference
            idx, ev, se_frac, left_vec, ref_right_vec = \
                find_slow_se_mode_full(eigvals, R, R_inv, N)
            track_ov = 1.0
        else:
            # ε > 0: track the same mode by eigenvector overlap
            idx, ev, se_frac, left_vec = \
                track_mode(eigvals, R, R_inv, ref_right_vec)
            track_ov = float(np.max(np.abs(R.conj().T @ ref_right_vec)))

        print(f"  {eps:>8.4f}  {se_frac:>8.4f}  {-ev.real:>10.6f}  "
              f"{ev.imag:>+10.6f}  {track_ov:>8.4f}")

        # Measure Bell-pair dependence of c_slow
        exc_k = 4
        c_slow_values = {}
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4)]:
            if exc_k in (a, b):
                continue
            psi = make_bell_exc_state(N, a, b, exc_k)
            rho0 = np.outer(psi, psi.conj())
            c_slow = complex(left_vec @ rho0.flatten())
            c_slow_values[f"({a},{b})"] = round(abs(c_slow), 8)

        vals = list(c_slow_values.values())
        spread = max(vals) - min(vals) if len(vals) > 1 else 0.0

        results[str(eps)] = dict(
            se_fraction=round(se_frac, 6),
            rate=round(float(-ev.real), 6),
            tracking_overlap=round(track_ov, 6),
            c_slow_by_bell_pair=c_slow_values,
            spread=round(spread, 8),
        )

    # Print Bell-pair dependence table
    print(f"\n--- Bell-pair dependence of |c_slow| (exc_k=4) ---\n")
    print(f"  {'ε':>8}  {'Bell(0,1)':>10}  {'Bell(1,2)':>10}  {'Bell(2,3)':>10}  {'Spread':>10}")
    print(f"  {'-' * 55}")
    for eps in epsilons:
        r = results[str(eps)]
        cs = r['c_slow_by_bell_pair']
        vals = [cs.get("(0,1)", 0), cs.get("(1,2)", 0), cs.get("(2,3)", 0)]
        print(f"  {eps:>8.4f}  {vals[0]:>10.6f}  {vals[1]:>10.6f}  {vals[2]:>10.6f}  {r['spread']:>10.2e}")

    # Scaling analysis: how does spread grow with ε?
    print(f"\n--- Scaling: spread vs ε ---\n")
    eps_nonzero = [e for e in epsilons if e > 0 and results[str(e)]['spread'] > 1e-15]
    if len(eps_nonzero) >= 3:
        log_eps = np.log10([e for e in eps_nonzero])
        log_spread = np.log10([results[str(e)]['spread'] for e in eps_nonzero])
        # Linear fit in log-log
        coeffs = np.polyfit(log_eps, log_spread, 1)
        slope = coeffs[0]
        print(f"  Log-log slope: {slope:.2f}")
        if abs(slope - 1.0) < 0.2:
            print(f"  -> Linear scaling (slope ≈ 1)")
        elif abs(slope - 2.0) < 0.3:
            print(f"  -> Quadratic scaling (slope ≈ 2)")
        else:
            print(f"  -> Non-integer scaling (slope = {slope:.2f})")
        results['scaling_slope'] = round(slope, 4)

    # SE fraction scaling
    print(f"\n--- SE fraction departure from 1.0 ---\n")
    eps_nz = [e for e in epsilons if e > 0]
    departures = [1.0 - results[str(e)]['se_fraction'] for e in eps_nz]
    if min(departures) > 1e-15:
        log_dep = np.log10(departures)
        log_e = np.log10(eps_nz)
        coeffs2 = np.polyfit(log_e, log_dep, 1)
        slope2 = coeffs2[0]
        print(f"  SE departure scaling: slope = {slope2:.2f}")
        results['se_departure_slope'] = round(slope2, 4)

    # Save
    with open(os.path.join(OUT_DIR, 'u1_breaking_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n{'=' * 60}")
    print(f"Results saved to {OUT_DIR}/")
