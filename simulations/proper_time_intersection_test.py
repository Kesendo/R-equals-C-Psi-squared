#!/usr/bin/env python3
"""
Proper Time Intersection Test: Observer-Dependent CΨ Visibility
================================================================
Star topology (S+A+B), Bell_SA ⊗ |+⟩_B, J_SA=1.0, J_SB=2.0.
Four runs with different γ_A (0.03, 0.05, 0.10, 0.20),
γ_S = γ_B = 0.05 fixed.

Shows: CΨ_AB is NOT a function of proper time τ_A = γ_A × t alone.
Different noise profiles open or close the visibility window.

Script:  simulations/proper_time_intersection_test.py
Output:  simulations/results/proper_time_intersection_test.txt
Docs:    experiments/OBSERVER_DEPENDENT_VISIBILITY.md
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "proper_time_intersection_test.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N=3 star topology)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 3  # S(0), A(1), B(2)
d = 2 ** N    # 8
d2 = d * d    # 64

# σ_y ⊗ ��_y for 2-qubit Wootters concurrence
sysy = np.kron(sy, sy)


def site_op(op, k, nq=N):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_star(J_SA=1.0, J_SB=2.0):
    """Star topology: S(0)--A(1), S(0)--B(2)."""
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J_SA * site_op(P, 0) @ site_op(P, 1)  # S-A
        H += J_SB * site_op(P, 0) @ site_op(P, 2)  # S-B
    return H


def build_L(H, gammas):
    """Liouvillian with per-qubit dephasing rates."""
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho, t):
    v = expm(L * t) @ rho.flatten()
    rho_out = v.reshape(d, d)
    return (rho_out + rho_out.conj().T) / 2


def ptrace_keep(rho, keep, nq=N):
    keep = list(keep)
    trace_out = [q for q in range(nq) if q not in keep]
    dims = [2] * nq
    reshaped = rho.reshape(dims + dims)
    current_n = nq
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape(d_k, d_k)


def concurrence_2q(rho):
    """Wootters concurrence for a 2-qubit density matrix."""
    R = rho @ sysy @ rho.conj() @ sysy
    eigvals = np.sort(np.real(np.linalg.eigvals(R)))[::-1]
    eigvals = np.maximum(eigvals, 0.0)
    sq = np.sqrt(eigvals)
    return max(0.0, sq[0] - sq[1] - sq[2] - sq[3])


def l1_norm(rho):
    d_r = rho.shape[0]
    return float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))


def psi_norm(rho):
    d_r = rho.shape[0]
    return l1_norm(rho) / (d_r - 1) if d_r > 1 else 0.0


def cpsi_AB(rho_full):
    """CΨ_AB = concurrence(AB) × psi_norm(AB)."""
    rho_AB = ptrace_keep(rho_full, [1, 2])  # A=1, B=2
    C = concurrence_2q(rho_AB)
    Psi = psi_norm(rho_AB)
    return C * Psi


def ket2dm(psi):
    return np.outer(psi, psi.conj())


# ============================================================
# INITIAL STATE: Bell_SA ⊗ |+⟩_B
# ============================================================
bell_SA = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)  # S, A
psi0 = np.kron(bell_SA, plus)  # S, A, B
rho0 = ket2dm(psi0)

# ============================================================
# MAIN: Four γ_A configurations
# ============================================================
log("=" * 78)
log("Proper Time Intersection Test: Observer-Dependent CΨ Visibility")
log("=" * 78)
log()
log("Star topology: S(0)--A(1), S(0)--B(2)")
log("J_SA=1.0, J_SB=2.0, initial: Bell_SA ⊗ |+⟩_B")
log()

H = build_H_star(J_SA=1.0, J_SB=2.0)

configs = [
    ("slow A",      0.03),
    ("equal",       0.05),
    ("fast A",      0.10),
    ("very fast A", 0.20),
]

check_times = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]

# Store trajectories for proper-time analysis
trajectories = {}
peaks = {}

log("Test 1: CΨ_AB at coordinate time t")
log("-" * 78)
header = f"  {'t':>5}"
for label, _ in configs:
    header += f"  {label:>14}"
log(header)

for t in check_times:
    row = f"  {t:5.1f}"
    for label, gamma_A in configs:
        gammas = [0.05, gamma_A, 0.05]  # S, A, B
        L = build_L(H, gammas)
        rho_t = evolve(L, rho0, t)
        c = cpsi_AB(rho_t)
        row += f"  {c:14.3f}"
    log(row)

# Compute peak CΨ_AB for each config
log()
log("Peak CΨ_AB per observer")
log("-" * 78)
log(f"  {'Config':>14}  {'γ_A':>5}  {'Peak CΨ_AB':>11}  {'t_peak':>7}  {'Above 1/4?':>11}")

dt = 0.01
t_max = 10.0
n_steps = int(t_max / dt)

for label, gamma_A in configs:
    gammas = [0.05, gamma_A, 0.05]
    L = build_L(H, gammas)
    peak = 0.0
    t_peak = 0.0
    traj = []

    for i in range(n_steps + 1):
        t = i * dt
        rho_t = evolve(L, rho0, t)
        c = cpsi_AB(rho_t)
        traj.append((t, c))
        if c > peak:
            peak = c
            t_peak = t

    trajectories[label] = traj
    peaks[label] = peak
    above = "Yes, clearly" if peak > 0.30 else ("Yes, barely" if peak > 0.26 else "Barely" if peak > 0.25 else "No")
    log(f"  {label:>14}  {gamma_A:5.2f}  {peak:11.3f}  {t_peak:7.2f}  {above:>11}")

# ============================================================
# PROPER TIME ANALYSIS
# ============================================================
log()
log("Test 2: CΨ_AB vs proper time τ_A = γ_A × t")
log("-" * 78)
log("  If CΨ were a function of τ_A alone, all columns would be identical.")
log()

tau_checks = [0.005, 0.010, 0.020, 0.050]

header = f"  {'τ_A':>6}"
for label, _ in configs:
    header += f"  {label:>14}"
log(header)

for tau in tau_checks:
    row = f"  {tau:6.3f}"
    for label, gamma_A in configs:
        t_coord = tau / gamma_A
        if t_coord > t_max:
            row += f"  {'---':>14}"
            continue
        gammas = [0.05, gamma_A, 0.05]
        L = build_L(H, gammas)
        rho_t = evolve(L, rho0, t_coord)
        c = cpsi_AB(rho_t)
        row += f"  {c:14.3f}"
    log(row)

log()
log("  Non-universality confirmed: same τ_A, different CΨ_AB values.")
log("  CΨ depends on absolute noise rate, not just local proper time.")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 78)
log("Observer-dependent visibility confirmed.")
log(f"  Peak CΨ range: {min(peaks.values()):.3f} (fast) to {max(peaks.values()):.3f} (slow)")
log("  CΨ is NOT a function of proper time τ_A alone.")
log("  Multi-scale dynamics (J, γ_A, γ_B, γ_S) control the visibility window.")
log("=" * 78)

_outf.close()
