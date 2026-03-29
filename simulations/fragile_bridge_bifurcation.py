#!/usr/bin/env python3
"""
Fragile Bridge Bifurcation Analysis
=====================================
Two Heisenberg chains (decay + gain) coupled by a bridge.
When does the system become unstable?

A) gamma_crit(J_bridge) sweep for N=2 per chain
B) Functional form: power law vs linear fit
C) N-scaling: repeat for N=3 per chain
D) Eigenvalue trajectories at J_bridge=1.0

Script: simulations/fragile_bridge_bifurcation.py
Output: simulations/results/fragile_bridge_bifurcation.txt
"""

import numpy as np
from scipy.linalg import eigvals
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fragile_bridge_bifurcation.txt")
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
    """Place single-qubit operator on qubit k in n-qubit system."""
    result = np.array([[1]], dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == qubit else I2)
    return result


def build_coupled_liouvillian(n_per_chain, gamma, J=1.0, J_bridge=0.5):
    """
    Two N-qubit Heisenberg chains coupled by a bridge.

    Chain A: qubits 0..n-1, dephasing +gamma per site
    Chain B: qubits n..2n-1, dephasing -gamma per site
    Bridge: Heisenberg coupling J_bridge between qubit n-1 and qubit n
    """
    n_total = 2 * n_per_chain
    d = 2**n_total
    d2 = d * d

    # Hamiltonian
    H = np.zeros((d, d), dtype=complex)
    # Chain A: nearest-neighbor
    for i in range(n_per_chain - 1):
        for P in [sx, sy, sz]:
            H += J * op_at(P, i, n_total) @ op_at(P, i + 1, n_total)
    # Chain B: nearest-neighbor
    for i in range(n_per_chain, 2 * n_per_chain - 1):
        for P in [sx, sy, sz]:
            H += J * op_at(P, i, n_total) @ op_at(P, i + 1, n_total)
    # Bridge: last qubit of A to first qubit of B
    for P in [sx, sy, sz]:
        H += J_bridge * op_at(P, n_per_chain - 1, n_total) @ \
             op_at(P, n_per_chain, n_total)

    # Liouvillian
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    # Dephasing: +gamma on chain A, -gamma on chain B
    for k in range(n_per_chain):
        Zk = op_at(sz, k, n_total)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_at(sz, k, n_total)
        L += (-gamma) * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))

    return L


def is_unstable(L, threshold=1e-10):
    """Check if Liouvillian has any eigenvalue with Re > threshold."""
    evals = eigvals(L)
    return float(np.max(evals.real)) > threshold


def find_gamma_crit(n_per_chain, J_bridge, J=1.0, tol=1e-6):
    """Bisect to find gamma where max Re(lambda) first becomes > 0."""
    # First: check if system is ever unstable
    if not is_unstable(build_coupled_liouvillian(n_per_chain, 2.0, J, J_bridge)):
        if not is_unstable(build_coupled_liouvillian(n_per_chain, 20.0, J, J_bridge)):
            return None  # never unstable

    # Bisection between 0 and known-unstable gamma
    g_lo, g_hi = 0.0, 2.0
    # Ensure g_hi is unstable
    while not is_unstable(build_coupled_liouvillian(n_per_chain, g_hi, J, J_bridge)):
        g_hi *= 2

    while (g_hi - g_lo) > tol:
        g_mid = (g_lo + g_hi) / 2
        if is_unstable(build_coupled_liouvillian(n_per_chain, g_mid, J, J_bridge)):
            g_hi = g_mid
        else:
            g_lo = g_mid

    return (g_lo + g_hi) / 2


# ================================================================
log("=" * 70)
log("FRAGILE BRIDGE BIFURCATION ANALYSIS")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

# ================================================================
# A) gamma_crit sweep over J_bridge (N=2 per chain)
# ================================================================
log()
log("=" * 70)
log("AUFGABE A: gamma_crit(J_bridge) for N=2 per chain")
log("4 qubits total, Liouvillian 256x256")
log("=" * 70)
log()

j_bridges = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
results_a = []

log(f"  {'J_bridge':>9}  {'gamma_crit':>11}  {'max_Re_at_crit':>15}  {'time':>6}")
log(f"  {'-'*50}")

for jb in j_bridges:
    t0 = clock.time()
    gc = find_gamma_crit(2, jb, J=1.0, tol=1e-6)
    elapsed = clock.time() - t0

    if gc is not None:
        # Verify: max Re just above gamma_crit
        L = build_coupled_liouvillian(2, gc * 1.01, J=1.0, J_bridge=jb)
        mr = float(np.max(eigvals(L).real))
        results_a.append((jb, gc, mr))
        log(f"  {jb:>9.3f}  {gc:>11.6f}  {mr:>15.2e}  {elapsed:>5.1f}s")
    else:
        results_a.append((jb, None, None))
        log(f"  {jb:>9.3f}  {'never':>11}  {'N/A':>15}  {elapsed:>5.1f}s")

# ================================================================
# B) Functional form
# ================================================================
log()
log("=" * 70)
log("AUFGABE B: Functional form of gamma_crit(J_bridge)")
log("=" * 70)
log()

valid = [(jb, gc) for jb, gc, _ in results_a if gc is not None]
if len(valid) >= 3:
    x = np.array([v[0] for v in valid])
    y = np.array([v[1] for v in valid])

    # Power law fit: log(gamma_crit) = log(a) + b * log(J_bridge)
    log_x = np.log(x)
    log_y = np.log(y)
    b_pow, log_a = np.polyfit(log_x, log_y, 1)
    a_pow = np.exp(log_a)
    y_pred_pow = a_pow * x**b_pow
    ss_res_pow = np.sum((y - y_pred_pow)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2_pow = 1 - ss_res_pow / ss_tot

    log(f"  Power law: gamma_crit = {a_pow:.6f} * J_bridge^{b_pow:.4f}")
    log(f"  R^2 = {r2_pow:.8f}")
    log()

    # Linear fit: gamma_crit = a * J_bridge + c
    coeffs_lin = np.polyfit(x, y, 1)
    a_lin, c_lin = coeffs_lin
    y_pred_lin = a_lin * x + c_lin
    ss_res_lin = np.sum((y - y_pred_lin)**2)
    r2_lin = 1 - ss_res_lin / ss_tot

    log(f"  Linear: gamma_crit = {a_lin:.6f} * J_bridge + {c_lin:.6f}")
    log(f"  R^2 = {r2_lin:.8f}")
    log()

    # Linear through origin: gamma_crit = a * J_bridge
    a_origin = np.sum(x * y) / np.sum(x**2)
    y_pred_origin = a_origin * x
    ss_res_origin = np.sum((y - y_pred_origin)**2)
    r2_origin = 1 - ss_res_origin / ss_tot

    log(f"  Linear (origin): gamma_crit = {a_origin:.6f} * J_bridge")
    log(f"  R^2 = {r2_origin:.8f}")
    log()

    best = max([("Power law", r2_pow), ("Linear", r2_lin),
                ("Linear (origin)", r2_origin)], key=lambda x: x[1])
    log(f"  Best fit: {best[0]} (R^2 = {best[1]:.8f})")

    # Comparison table
    log()
    log(f"  {'J_bridge':>9}  {'gamma_crit':>11}  {'power_law':>11}  {'linear':>11}  {'origin':>11}")
    log(f"  {'-'*60}")
    for jb, gc in valid:
        yp = a_pow * jb**b_pow
        yl = a_lin * jb + c_lin
        yo = a_origin * jb
        log(f"  {jb:>9.3f}  {gc:>11.6f}  {yp:>11.6f}  {yl:>11.6f}  {yo:>11.6f}")

# ================================================================
# C) N-scaling
# ================================================================
log()
log("=" * 70)
log("AUFGABE C: N-scaling (N=3 per chain, 6 qubits, 4096x4096)")
log("=" * 70)
log()

j_bridges_c = [0.1, 0.2, 0.5, 1.0, 2.0]
results_c = []

log(f"  {'J_bridge':>9}  {'gamma_crit(N=2)':>16}  {'gamma_crit(N=3)':>16}  {'ratio':>8}  {'time':>6}")
log(f"  {'-'*65}")

for jb in j_bridges_c:
    t0 = clock.time()
    gc3 = find_gamma_crit(3, jb, J=1.0, tol=1e-5)
    elapsed = clock.time() - t0

    # Find matching N=2 result
    gc2 = None
    for jb2, gc, _ in results_a:
        if abs(jb2 - jb) < 1e-10:
            gc2 = gc
            break

    if gc2 is not None and gc3 is not None:
        ratio = gc3 / gc2
        results_c.append((jb, gc2, gc3, ratio))
        log(f"  {jb:>9.3f}  {gc2:>16.6f}  {gc3:>16.6f}  {ratio:>8.3f}  {elapsed:>5.1f}s")
    elif gc3 is not None:
        results_c.append((jb, gc2, gc3, None))
        log(f"  {jb:>9.3f}  {'N/A':>16}  {gc3:>16.6f}  {'N/A':>8}  {elapsed:>5.1f}s")
    else:
        results_c.append((jb, gc2, None, None))
        log(f"  {jb:>9.3f}  {gc2:>16.6f}  {'never':>16}  {'N/A':>8}  {elapsed:>5.1f}s")

# ================================================================
# D) Eigenvalue trajectories at J_bridge=1.0
# ================================================================
log()
log("=" * 70)
log("AUFGABE D: Eigenvalue trajectories at J_bridge=1.0 (N=2)")
log("=" * 70)
log()

jb_fixed = 1.0
# Find gamma_crit for this J_bridge
gc_ref = None
for jb, gc, _ in results_a:
    if abs(jb - jb_fixed) < 1e-10:
        gc_ref = gc
        break

if gc_ref is not None:
    log(f"  gamma_crit = {gc_ref:.6f} at J_bridge = {jb_fixed}")
    log()

    # Sweep gamma from 0 to 2*gamma_crit
    gammas_d = np.linspace(0, 2.0 * gc_ref, 41)
    n_top = 10  # track top 10 eigenvalues

    log(f"  {'gamma':>8}  {'gamma/gc':>8}  ", end="")
    log("  ".join([f"Re(l{i+1:d})" for i in range(min(n_top, 6))]))
    log(f"  {'-'*80}")

    coalescence_gamma = None
    prev_top_real = None

    for gamma in gammas_d:
        L = build_coupled_liouvillian(2, gamma, J=1.0, J_bridge=jb_fixed)
        evals = eigvals(L)

        # Sort by real part (descending)
        idx = np.argsort(-evals.real)
        top = evals[idx[:n_top]]

        # Check for coalescence: two real eigenvalues merging into complex pair
        top_real_only = top[np.abs(top.imag) < 1e-8]
        top_complex = top[np.abs(top.imag) >= 1e-8]

        if coalescence_gamma is None and len(top) > 0 and top[0].real > 1e-6:
            coalescence_gamma = gamma

        ratio = gamma / gc_ref if gc_ref > 0 else 0
        vals_str = "  ".join([f"{v.real:>+8.4f}" for v in top[:min(n_top, 6)]])
        im_str = "  ".join([f"{'c' if abs(v.imag)>1e-4 else 'r':>1}" for v in top[:min(n_top, 6)]])
        log(f"  {gamma:>8.4f}  {ratio:>8.3f}  {vals_str}  [{im_str}]")

    # Analyze the transition
    log()
    log("  r = real eigenvalue, c = complex (has imaginary part)")
    log()
    if coalescence_gamma is not None:
        log(f"  Instability onset at gamma = {coalescence_gamma:.4f}"
            f" (gamma/gamma_crit = {coalescence_gamma/gc_ref:.3f})")
    else:
        log("  No instability found in sweep range")

    # Check PT-symmetry breaking signature
    log()
    log("  PT-symmetry check:")
    log("  If two REAL eigenvalues coalesce and become a COMPLEX pair,")
    log("  that is classical PT-symmetry breaking.")
    log()

    # Detailed look near gamma_crit
    log("  Near gamma_crit (fine scan):")
    fine_gammas = np.linspace(0.9 * gc_ref, 1.1 * gc_ref, 11)
    log(f"  {'gamma':>8}  {'g/gc':>6}  {'Re(l1)':>10}  {'Im(l1)':>10}  {'Re(l2)':>10}  {'Im(l2)':>10}")
    log(f"  {'-'*60}")
    for gamma in fine_gammas:
        L = build_coupled_liouvillian(2, gamma, J=1.0, J_bridge=jb_fixed)
        evals = eigvals(L)
        idx = np.argsort(-evals.real)
        e1, e2 = evals[idx[0]], evals[idx[1]]
        ratio = gamma / gc_ref
        log(f"  {gamma:>8.5f}  {ratio:>6.3f}  {e1.real:>10.6f}  {e1.imag:>10.6f}"
            f"  {e2.real:>10.6f}  {e2.imag:>10.6f}")


# ================================================================
log()
log("=" * 70)
log("SUMMARY")
log("=" * 70)
log()

if len(valid) >= 3:
    log(f"gamma_crit(J_bridge): {best[0]} fit")
    if best[0] == "Power law":
        log(f"  gamma_crit = {a_pow:.6f} * J_bridge^{b_pow:.4f}")
    elif best[0] == "Linear":
        log(f"  gamma_crit = {a_lin:.6f} * J_bridge + {c_lin:.6f}")
    else:
        log(f"  gamma_crit = {a_origin:.6f} * J_bridge")
    log()

if results_c:
    ratios_c = [r[3] for r in results_c if r[3] is not None]
    if ratios_c:
        log(f"N-scaling: gamma_crit(N=3)/gamma_crit(N=2) = "
            f"{np.mean(ratios_c):.3f} +/- {np.std(ratios_c):.3f}")
        if abs(np.mean(ratios_c) - 1.0) < 0.1:
            log("  N-INDEPENDENT (within 10%)")
        else:
            log(f"  N-DEPENDENT (ratio != 1)")
    log()

log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
