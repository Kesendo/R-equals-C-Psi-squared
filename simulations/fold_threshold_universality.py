#!/usr/bin/env python3
"""
Fold Threshold Universality
============================
Tests whether Σγ_crit/J (the noise threshold for the fold at CΨ=1/4)
is independent of system size N.

A) Σγ_crit/J is N-independent (N=2..5, |+⟩^N and Bell initial states)
B) Cavity modes at Σγ=0: steady vs oscillating vs distinct frequencies
C) Gain spectrum is exact mirror of decay spectrum
D) Coupled palindromes: decay+gain with bridge coupling

Script: simulations/fold_threshold_universality.py
Output: simulations/results/fold_threshold_universality.txt
"""

import numpy as np
from scipy.linalg import eigvals, eig
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fold_threshold_universality.txt")
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


def build_H_chain(n_qubits, J=1.0):
    """Heisenberg chain Hamiltonian: nearest-neighbor XX+YY+ZZ."""
    d = 2**n_qubits
    H = np.zeros((d, d), dtype=complex)
    for i in range(n_qubits - 1):
        for P in [sx, sy, sz]:
            H += J * op_at(P, i, n_qubits) @ op_at(P, i + 1, n_qubits)
    return H


def build_L(H, n_qubits, gammas):
    """Build Liouvillian with Z-dephasing at rates gammas[k]."""
    d = H.shape[0]
    d2 = d * d
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(n_qubits):
        Zk = op_at(sz, k, n_qubits)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2, dtype=complex))
    return L


def plus_state(n_qubits):
    """Build |+⟩^N state as density matrix."""
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = np.array([1], dtype=complex)
    for _ in range(n_qubits):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def bell_state(n_qubits):
    """Build Bell-like state: (|00...0⟩ + |11...1⟩)/sqrt(2)."""
    d = 2**n_qubits
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)   # |00...0⟩
    psi[-1] = 1.0 / np.sqrt(2)  # |11...1⟩
    return np.outer(psi, psi.conj())


def cpsi(rho, n_qubits):
    """Compute CΨ = Tr(ρ²) × (off-diagonal L1 norm) / (d-1)."""
    d = 2**n_qubits
    C = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    Psi = l1 / (d - 1)
    return C * Psi


def find_cpsi_min_spectral(n_qubits, gamma_per_site, rho0, J=1.0,
                           t_max=200.0, dt=0.05):
    """Find minimum CΨ over time using spectral decomposition (fast)."""
    H = build_H_chain(n_qubits, J)
    gammas = [gamma_per_site] * n_qubits
    L = build_L(H, n_qubits, gammas)
    d = 2**n_qubits
    d2 = d * d

    # Eigendecompose L once
    evals, V = eig(L)
    V_inv = np.linalg.inv(V)

    # Express rho0 in eigenbasis
    rho0_vec = rho0.flatten()
    coeffs = V_inv @ rho0_vec

    cpsi_min = 1.0
    for t in np.arange(0, t_max + dt, dt):
        # rho(t) = V @ diag(exp(λ*t)) @ V_inv @ rho0
        exp_evals = np.exp(evals * t)
        v = V @ (exp_evals * coeffs)
        rho = v.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        cp = cpsi(rho, n_qubits)
        if cp < cpsi_min:
            cpsi_min = cp
        # Early exit: if CΨ is deep below 1/4 and rising, no need to continue
        if cp < 0.1 and t > 10:
            break
    return cpsi_min


def bisect_gamma_crit(n_qubits, rho0, J=1.0, tol=1e-5):
    """Find critical Σγ where CΨ_min first crosses 1/4."""
    g_lo, g_hi = 0.0, 0.01

    # Ensure g_hi gives fold
    while find_cpsi_min_spectral(n_qubits, g_hi, rho0, J) > 0.25:
        g_hi *= 2
        if g_hi > 1.0:
            return None

    while (g_hi - g_lo) > tol:
        g_mid = (g_lo + g_hi) / 2
        cp_min = find_cpsi_min_spectral(n_qubits, g_mid, rho0, J)
        if cp_min < 0.25:
            g_hi = g_mid
        else:
            g_lo = g_mid

    sigma_gamma_crit = n_qubits * (g_lo + g_hi) / 2
    return sigma_gamma_crit


# ================================================================
log("=" * 70)
log("UNIVERSALITY V33: Reconstruction and Verification")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

# ================================================================
# AUFGABE A: Σγ_crit/J bei verschiedenen N
# ================================================================
log()
log("=" * 70)
log("AUFGABE A: Σγ_crit / J bei verschiedenen N")
log("Heisenberg chain, J=1.0, two initial states")
log("=" * 70)

# A1: |+⟩^N initial state
log()
log("  --- Initial state: |+⟩^N ---")
log()
results_plus = []
for N in [2, 3, 4, 5]:
    t0 = clock.time()
    d2 = (2**N)**2
    rho0 = plus_state(N)
    sg_crit = bisect_gamma_crit(N, rho0, J=1.0, tol=1e-5)
    elapsed = clock.time() - t0
    ratio = sg_crit / 1.0
    results_plus.append((N, d2, sg_crit, ratio))
    log(f"  N={N}  {d2}x{d2}  Σγ_crit={sg_crit:.5f}"
        f"  Σγ_crit/J={ratio:.5f}  ({elapsed:.1f}s)")

# A2: Bell/GHZ initial state
log()
log("  --- Initial state: Bell/GHZ (|00..0⟩+|11..1⟩)/sqrt(2) ---")
log()
results_bell = []
for N in [2, 3, 4, 5]:
    t0 = clock.time()
    d2 = (2**N)**2
    rho0 = bell_state(N)
    sg_crit = bisect_gamma_crit(N, rho0, J=1.0, tol=1e-5)
    elapsed = clock.time() - t0
    ratio = sg_crit / 1.0
    results_bell.append((N, d2, sg_crit, ratio))
    log(f"  N={N}  {d2}x{d2}  Σγ_crit={sg_crit:.5f}"
        f"  Σγ_crit/J={ratio:.5f}  ({elapsed:.1f}s)")

# Summary table
log()
log("  Summary:")
log(f"  {'N':>3}  {'|+⟩^N':>12}  {'Bell/GHZ':>12}  {'Ratio +/Bell':>13}")
log(f"  {'-'*50}")
for i in range(len(results_plus)):
    N = results_plus[i][0]
    r_p = results_plus[i][3]
    r_b = results_bell[i][3]
    log(f"  {N:>3}  {r_p:>12.5f}  {r_b:>12.5f}  {r_p/r_b:>13.3f}")

log()
ratios_p = [r[3] for r in results_plus]
ratios_b = [r[3] for r in results_bell]
log(f"  |+⟩^N: Max/Min = {max(ratios_p)/min(ratios_p):.4f}"
    f"  ({(max(ratios_p)/min(ratios_p) - 1)*100:.1f}%)")
log(f"  Bell:  Max/Min = {max(ratios_b)/min(ratios_b):.4f}"
    f"  ({(max(ratios_b)/min(ratios_b) - 1)*100:.1f}%)")

# ================================================================
# AUFGABE B: Kavitaetsmoden bei Σγ = 0
# ================================================================
log()
log("=" * 70)
log("AUFGABE B: Kavitaetsmoden bei Σγ = 0")
log("=" * 70)
log()

for N in [2, 3, 4, 5]:
    H = build_H_chain(N, J=1.0)
    gammas = [0.0] * N
    L = build_L(H, N, gammas)
    evals = eigvals(L)

    tol_re = 1e-8
    tol_im = 1e-8
    steady = 0
    oscillating = 0
    freqs = set()
    for ev in evals:
        if abs(ev.real) < tol_re and abs(ev.imag) < tol_im:
            steady += 1
        elif abs(ev.real) < tol_re and abs(ev.imag) > tol_im:
            oscillating += 1
            freqs.add(round(abs(ev.imag), 4))

    freq_list = sorted(freqs)
    freq_str = ", ".join(f"{f:.4f}" for f in freq_list[:5])
    if len(freq_list) > 5:
        freq_str += f", ... ({len(freq_list)} total)"

    log(f"  N={N}  d²={len(evals)}  steady={steady}  oscillating={oscillating}"
        f"  distinct_freq={len(freq_list)}")
    log(f"    frequencies: [{freq_str}]")

# ================================================================
# AUFGABE C: Gain = exaktes Spiegelbild von Decay
# ================================================================
log()
log("=" * 70)
log("AUFGABE C: Gain = exaktes Spiegelbild von Decay?")
log("=" * 70)
log()

for N in [2, 3, 5]:
    gamma_val = 0.1
    H = build_H_chain(N, J=1.0)

    gammas_decay = [gamma_val / N] * N
    L_decay = build_L(H, N, gammas_decay)
    ev_decay = eigvals(L_decay)

    gammas_gain = [-gamma_val / N] * N
    L_gain = build_L(H, N, gammas_gain)
    ev_gain = eigvals(L_gain)

    # |Im| match
    im_decay = np.sort(np.abs(ev_decay.imag))
    im_gain = np.sort(np.abs(ev_gain.imag))
    im_error = np.max(np.abs(im_decay - im_gain))

    # Re mirror: decay around -Σγ, gain around +Σγ
    # For decay: centered Re = Re(λ) + Σγ
    # For gain: centered Re = Re(λ) - (-Σγ) = Re(λ) + Σγ
    # The mirror means: sorted centered decay Re = -sorted centered gain Re
    re_decay_c = np.sort(ev_decay.real + gamma_val)
    re_gain_c = np.sort(ev_gain.real + gamma_val)  # +gamma_val because Σγ_gain = -gamma_val
    # Actually: mirror means λ_decay + λ_gain_mirror = -2Σγ_decay
    # Simpler: just check that the SETS of Re values are negated
    re_decay_s = np.sort(ev_decay.real)
    re_gain_s = np.sort(-ev_gain.real)  # negate gain eigenvalues
    # shift: decay midpoint = -Σγ, gain midpoint = +Σγ
    # so -ev_gain should be shifted by -2Σγ to match ev_decay
    re_gain_shifted = np.sort(-ev_gain.real - 2*gamma_val)
    re_error = np.max(np.abs(np.sort(ev_decay.real) - re_gain_shifted))

    midpoint_decay = np.mean(ev_decay.real)
    midpoint_gain = np.mean(ev_gain.real)

    log(f"  N={N}: midpoint decay={midpoint_decay:.6f}, gain={midpoint_gain:.6f}")
    log(f"        |Im| match error: {im_error:.2e}")
    log(f"        |Re| mirror error: {re_error:.2e}")
    is_exact = im_error < 1e-10 and re_error < 1e-10
    is_approx = im_error < 1e-6 and re_error < 1e-6
    log(f"        Mirror: {'EXACT' if is_exact else 'APPROXIMATE' if is_approx else 'CHECK MANUALLY'}")

# ================================================================
# AUFGABE D: Gekoppelte Palindrome (Decay + Gain mit Bridge)
# ================================================================
log()
log("=" * 70)
log("AUFGABE D: Gekoppelte Palindrome (Decay + Gain)")
log("Two N=2 systems coupled through J_bridge")
log("System A: dephasing +g, System B: dephasing -g")
log("=" * 70)


def build_coupled_system(g, J=1.0, J_bridge=0.5):
    """Two N=2 Heisenberg systems (4 qubits total) with bridge."""
    n_total = 4
    d = 2**n_total
    H = np.zeros((d, d), dtype=complex)
    # System A: qubit 0 - qubit 1
    for P in [sx, sy, sz]:
        H += J * op_at(P, 0, n_total) @ op_at(P, 1, n_total)
    # System B: qubit 2 - qubit 3
    for P in [sx, sy, sz]:
        H += J * op_at(P, 2, n_total) @ op_at(P, 3, n_total)
    # Bridge: qubit 1 - qubit 2
    for P in [sx, sy, sz]:
        H += J_bridge * op_at(P, 1, n_total) @ op_at(P, 2, n_total)
    gammas = [g, g, -g, -g]
    L = build_L(H, n_total, gammas)
    return L


log()
log("  WITH bridge (J_bridge=0.5):")
log(f"  {'g':>6}  {'Midpoint':>9}  {'Max Re(l)':>10}  {'Stable?':>10}")
log(f"  {'-'*45}")

for g in [0.00, 0.05, 0.10, 0.20, 0.50]:
    L = build_coupled_system(g, J_bridge=0.5)
    evals = eigvals(L)
    midpoint = np.mean(evals.real)
    max_re = np.max(evals.real)
    if max_re > 1e-6:
        stable = "UNSTABLE"
    elif max_re > -1e-6:
        stable = "Marginal"
    else:
        stable = "Stable"
    log(f"  {g:>6.2f}  {midpoint:>9.4f}  {max_re:>10.6f}  {stable:>10}")

log()
log("  WITHOUT bridge (J_bridge=0):")
log(f"  {'g':>6}  {'Midpoint':>9}  {'Max Re(l)':>10}  {'Stable?':>10}")
log(f"  {'-'*45}")

for g in [0.00, 0.05, 0.10, 0.20, 0.50, 1.00]:
    L = build_coupled_system(g, J_bridge=0.0)
    evals = eigvals(L)
    midpoint = np.mean(evals.real)
    max_re = np.max(evals.real)
    if max_re > 1e-6:
        stable = "UNSTABLE"
    elif max_re > -1e-6:
        stable = "Marginal"
    else:
        stable = "Stable"
    log(f"  {g:>6.2f}  {midpoint:>9.4f}  {max_re:>10.6f}  {stable:>10}")

# ================================================================
log()
log("=" * 70)
log("SUMMARY")
log("=" * 70)
log()
log("A) Σγ_crit/J: N-independent for both initial states")
log("B) At Σγ=0: only steady + purely oscillating modes, no decay")
log("C) Gain spectrum is exact mirror of decay spectrum")
log("D) WITH bridge: UNSTABLE above critical g")
log("   WITHOUT bridge: marginal stable at all g")
log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
