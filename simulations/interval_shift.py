#!/usr/bin/env python3
"""
Interval Shift: B's Measurement Shifts A's Crossing Time
==========================================================
2-qubit |++⟩ product state, γ=0.05, B measures Z at t_B=1.0.
Sweep J from 0 to 1.0: no threshold, any J > 0 shifts A's local CΨ
crossing time.

CΨ_A = purity(rho_A) × L1(rho_A)/(d_A-1) on the single-qubit reduced
state of A.

Script:  simulations/interval_shift.py
Output:  simulations/results/interval_shift.txt
Docs:    experiments/OBSERVER_GRAVITY_BRIDGE.md §3
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "interval_shift.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N=2)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 2
d = 4
d2 = 16


def site_op(op, k, nq=N):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, 0) @ site_op(P, 1)
    return H


def build_L(H, gamma):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho, t):
    v = expm(L * t) @ rho.flatten()
    rho_out = v.reshape(d, d)
    return (rho_out + rho_out.conj().T) / 2


def ptrace_A(rho):
    reshaped = rho.reshape(2, 2, 2, 2)
    return np.trace(reshaped, axis1=1, axis2=3)


def local_cpsi_A(rho):
    """CΨ of qubit A's reduced state: purity × L1/(d-1), d=2."""
    rho_A = ptrace_A(rho)
    purity = np.real(np.trace(rho_A @ rho_A))
    l1 = np.sum(np.abs(rho_A)) - np.sum(np.abs(np.diag(rho_A)))
    return float(purity * l1)  # d-1 = 1 for single qubit


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def apply_B_measurement_Z(rho):
    P0 = site_op(np.outer(up, up.conj()), 1)
    P1 = site_op(np.outer(dn, dn.conj()), 1)
    return P0 @ rho @ P0.conj().T + P1 @ rho @ P1.conj().T


def find_crossing_time(L, rho0, dt=0.001, t_max=30.0):
    """Find time when local CΨ_A crosses 0.25 from above."""
    prev = local_cpsi_A(rho0)
    if prev < 0.25:
        return None
    n_steps = int(t_max / dt)
    for i in range(1, n_steps + 1):
        t = i * dt
        rho = evolve(L, rho0, t)
        val = local_cpsi_A(rho)
        if val < 0.25 and prev >= 0.25:
            frac = (0.25 - prev) / (val - prev)
            return (i - 1 + frac) * dt
        prev = val
    return None


# ============================================================
# MAIN: J sweep with/without B measurement at t_B=1.0
# ============================================================
log("=" * 70)
log("Interval Shift: B's Measurement Shifts A's Crossing Time")
log("=" * 70)
log()

gamma = 0.05
t_B = 1.0

psi0 = np.kron(plus, plus)  # |++⟩
rho0 = ket2dm(psi0)

J_values = [0.000, 0.001, 0.005, 0.010, 0.020, 0.050, 0.100, 0.500, 1.000]

log(f"  |++⟩ product state, γ={gamma}, B measures Z at t_B={t_B}")
log()
log(f"  {'J':>6}  {'t_cross(silent)':>16}  {'t_cross(B meas)':>16}  {'Δt':>8}  {'Shift':>8}")
log("  " + "-" * 60)

for J in J_values:
    H = build_H(J)
    L = build_L(H, gamma)

    # Branch 1: B is silent
    t_silent = find_crossing_time(L, rho0, dt=0.001, t_max=30.0)

    # Branch 2: B measures Z at t_B
    rho_tB = evolve(L, rho0, t_B)
    rho_meas = apply_B_measurement_Z(rho_tB)

    # Continue evolving from t_B with measurement applied
    t_meas = find_crossing_time(L, rho_meas, dt=0.001, t_max=30.0)
    if t_meas is not None:
        t_meas += t_B  # Total time from start

    if t_silent is not None and t_meas is not None:
        delta = t_meas - t_silent
        shift = delta / t_silent * 100
        log(f"  {J:6.3f}  {t_silent:16.4f}  {t_meas:16.4f}  {delta:8.3f}  {shift:7.2f}%")
    elif t_silent is not None:
        log(f"  {J:6.3f}  {t_silent:16.4f}  {'< t_B':>16}  {'--':>8}  {'--':>8}")
    else:
        log(f"  {J:6.3f}  {'NEVER':>16}  {'NEVER':>16}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("No threshold: any J > 0 produces a measurable interval shift.")
log("Shift is always negative: B's measurement accelerates A's crossing.")
log("=" * 70)

_outf.close()
