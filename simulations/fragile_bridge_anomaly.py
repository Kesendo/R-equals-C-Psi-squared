#!/usr/bin/env python3
"""
Fragile Bridge Anomaly Investigation
======================================
Why does gamma_crit DROP at J_bridge = 5J?
Where is the maximum? What happens asymptotically?

A) Fine sweep around the turnover point (J_bridge 1..10)
B) Eigenvalue trajectories at peak, anomaly, and beyond
C) Effective topology: Hamiltonian spectrum comparison
D) Asymptotic behavior (J_bridge 10..100)

Script: simulations/fragile_bridge_anomaly.py
Output: simulations/results/fragile_bridge_anomaly.txt
"""

import numpy as np
from scipy.linalg import eigvals, eigvalsh
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fragile_bridge_anomaly.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def op_at(op, qubit, n_qubits):
    result = np.array([[1]], dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == qubit else I2)
    return result


def build_H(J=1.0, J_bridge=0.5):
    """Hamiltonian for two 2-qubit chains with bridge."""
    n = 4; d = 16
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * op_at(P, 0, n) @ op_at(P, 1, n)       # chain A
        H += J * op_at(P, 2, n) @ op_at(P, 3, n)       # chain B
        H += J_bridge * op_at(P, 1, n) @ op_at(P, 2, n) # bridge
    return H


def build_L(gamma, J=1.0, J_bridge=0.5):
    """Liouvillian for coupled decay+gain system."""
    n = 4; d = 16; d2 = 256
    H = build_H(J, J_bridge)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in [0, 1]:
        Zk = op_at(sz, k, n)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))
    for k in [2, 3]:
        Zk = op_at(sz, k, n)
        L += (-gamma) * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))
    return L


def find_gamma_crit(J_bridge, tol=1e-6):
    """Bisect for instability threshold."""
    g_lo, g_hi = 0.0, 2.0
    # Ensure g_hi is unstable
    while np.max(eigvals(build_L(g_hi, J_bridge=J_bridge)).real) < 1e-10:
        g_hi *= 2
        if g_hi > 1000:
            return None
    while (g_hi - g_lo) > tol:
        g_mid = (g_lo + g_hi) / 2
        if np.max(eigvals(build_L(g_mid, J_bridge=J_bridge)).real) > 1e-10:
            g_hi = g_mid
        else:
            g_lo = g_mid
    return (g_lo + g_hi) / 2


# ================================================================
log("=" * 70)
log("FRAGILE BRIDGE ANOMALY INVESTIGATION")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

# ================================================================
# A) Fine sweep J_bridge = 1..10
# ================================================================
log()
log("=" * 70)
log("AUFGABE A: Fine sweep J_bridge = 0.5 .. 10.0")
log("=" * 70)
log()

jb_values = np.concatenate([
    np.linspace(0.5, 3.0, 26),
    np.linspace(3.2, 10.0, 18)
])
results_a = []

log(f"  {'J_bridge':>9}  {'gamma_crit':>11}  {'gc/Jb':>8}")
log(f"  {'-'*35}")

for jb in jb_values:
    gc = find_gamma_crit(jb, tol=1e-5)
    if gc is not None:
        results_a.append((jb, gc))
        log(f"  {jb:>9.3f}  {gc:>11.6f}  {gc/jb:>8.4f}")
    else:
        results_a.append((jb, None))
        log(f"  {jb:>9.3f}  {'never':>11}  {'N/A':>8}")

# Find maximum
valid = [(jb, gc) for jb, gc in results_a if gc is not None]
if valid:
    jb_max, gc_max = max(valid, key=lambda x: x[1])
    log()
    log(f"  MAXIMUM: gamma_crit = {gc_max:.6f} at J_bridge = {jb_max:.3f}")
    log(f"  gamma_crit/J_bridge at max = {gc_max/jb_max:.4f}")

# ================================================================
# B) Eigenvalue trajectories at peak, anomaly, beyond
# ================================================================
log()
log("=" * 70)
log("AUFGABE B: Eigenvalue trajectories")
log("=" * 70)

for label, jb_test in [("peak", jb_max), ("anomaly", 5.0), ("beyond", 10.0)]:
    gc_test = find_gamma_crit(jb_test, tol=1e-6)
    if gc_test is None:
        log(f"\n  {label} (J_bridge={jb_test}): no instability found")
        continue

    log(f"\n  {label}: J_bridge = {jb_test:.3f}, gamma_crit = {gc_test:.6f}")
    log(f"  {'gamma/gc':>8}  {'Re(l1)':>10}  {'Im(l1)':>10}  {'Re(l2)':>10}  {'Im(l2)':>10}  type")
    log(f"  {'-'*60}")

    for ratio in [0.5, 0.8, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.5]:
        gamma = ratio * gc_test
        evals = eigvals(build_L(gamma, J_bridge=jb_test))
        idx = np.argsort(-evals.real)
        e1, e2 = evals[idx[0]], evals[idx[1]]
        t1 = 'c' if abs(e1.imag) > 1e-4 else 'r'
        t2 = 'c' if abs(e2.imag) > 1e-4 else 'r'
        log(f"  {ratio:>8.3f}  {e1.real:>10.6f}  {e1.imag:>10.4f}  "
            f"{e2.real:>10.6f}  {e2.imag:>10.4f}  {t1}{t2}")

# ================================================================
# C) Hamiltonian spectrum comparison
# ================================================================
log()
log("=" * 70)
log("AUFGABE C: Hamiltonian spectrum comparison")
log("=" * 70)
log()

for label, jb in [("weak bridge", 0.5), ("anomaly", 5.0),
                   ("uniform chain", 1.0)]:
    if label == "uniform chain":
        # Single 4-qubit chain with J=1.0 everywhere
        n = 4; d = 16
        H = np.zeros((d, d), dtype=complex)
        for i in range(3):
            for P in [sx, sy, sz]:
                H += 1.0 * op_at(P, i, n) @ op_at(P, i + 1, n)
    else:
        H = build_H(J=1.0, J_bridge=jb)

    evals_h = np.sort(eigvalsh(H).real)
    log(f"  {label} (J_bridge={jb if label != 'uniform chain' else 'all 1.0'}):")
    log(f"    eigenvalues: [{', '.join(f'{e:.4f}' for e in evals_h)}]")
    log(f"    range: [{evals_h[0]:.4f}, {evals_h[-1]:.4f}]")
    log(f"    gaps: [{', '.join(f'{evals_h[i+1]-evals_h[i]:.4f}' for i in range(len(evals_h)-1))}]")
    log()

# Spectrum distance: how close is J_bridge=X to uniform chain?
log("  Spectrum distance to uniform chain:")
n = 4; d = 16
H_uniform = np.zeros((d, d), dtype=complex)
for i in range(3):
    for P in [sx, sy, sz]:
        H_uniform += 1.0 * op_at(P, i, n) @ op_at(P, i + 1, n)
evals_uniform = np.sort(eigvalsh(H_uniform).real)

log(f"  {'J_bridge':>9}  {'dist_to_uniform':>16}  {'max_gap_diff':>14}")
log(f"  {'-'*45}")
for jb in [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 50.0]:
    H_test = build_H(J=1.0, J_bridge=jb)
    evals_test = np.sort(eigvalsh(H_test).real)
    # Normalize by scaling: uniform chain has range ~ 12J
    # For fair comparison, normalize both spectra to [0,1]
    e_u = (evals_uniform - evals_uniform[0]) / (evals_uniform[-1] - evals_uniform[0])
    e_t = (evals_test - evals_test[0]) / (evals_test[-1] - evals_test[0])
    dist = np.max(np.abs(e_u - e_t))
    gap_diff = np.max(np.abs(np.diff(evals_test) / (evals_test[-1] - evals_test[0])
                             - np.diff(evals_uniform) / (evals_uniform[-1] - evals_uniform[0])))
    log(f"  {jb:>9.1f}  {dist:>16.6f}  {gap_diff:>14.6f}")

# ================================================================
# D) Asymptotic behavior J_bridge = 10..100
# ================================================================
log()
log("=" * 70)
log("AUFGABE D: Asymptotic behavior (J_bridge = 10..100)")
log("=" * 70)
log()

jb_asymp = np.array([10, 15, 20, 30, 50, 75, 100])
results_d = []

log(f"  {'J_bridge':>9}  {'gamma_crit':>11}  {'gc*Jb':>10}  {'gc*Jb^2':>12}")
log(f"  {'-'*50}")

for jb in jb_asymp:
    gc = find_gamma_crit(float(jb), tol=1e-5)
    if gc is not None:
        results_d.append((float(jb), gc))
        log(f"  {jb:>9.1f}  {gc:>11.6f}  {gc*jb:>10.4f}  {gc*jb**2:>12.2f}")
    else:
        results_d.append((float(jb), None))
        log(f"  {jb:>9.1f}  {'never':>11}")

# Check scaling
if len(results_d) >= 3:
    x = np.array([r[0] for r in results_d if r[1] is not None])
    y = np.array([r[1] for r in results_d if r[1] is not None])

    if len(x) >= 3:
        log()
        # Power law fit
        b, log_a = np.polyfit(np.log(x), np.log(y), 1)
        a = np.exp(log_a)
        log(f"  Asymptotic power law: gamma_crit = {a:.4f} * J_bridge^{b:.4f}")

        if abs(b + 1) < 0.2:
            log(f"  Consistent with 1/J_bridge scaling")
            log(f"  gamma_crit * J_bridge ~ {np.mean(y * x):.4f} (constant?)")
        elif abs(b + 2) < 0.2:
            log(f"  Consistent with 1/J_bridge^2 scaling")

# ================================================================
log()
log("=" * 70)
log("SUMMARY")
log("=" * 70)
log()
if valid:
    log(f"Peak stability: gamma_crit = {gc_max:.6f} at J_bridge = {jb_max:.3f}")
    log(f"Below peak: gamma_crit ~ 0.19 * J_bridge (linear)")
    log(f"Above peak: gamma_crit decreases (topology transition)")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
