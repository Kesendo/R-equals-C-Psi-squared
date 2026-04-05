#!/usr/bin/env python3
"""
Shift Mechanism: Why B's Measurement Accelerates A's Crossing
==============================================================
2-qubit |++⟩, J=0.5, γ=0.05. B measures Z at t_B=1.0.

Tracks:
1. Local coherence A, purity A, concurrence, nonlocal coherence
   before/after B's measurement
2. Aftermath: A's coherence decay comparison (silent vs measured)
3. Timing dependence: sweep t_B, measure damage
4. Regeneration: concurrence after B's measurement
5. Coupling vs isolation: J sweep for local crossing time

Key finding: coupling REDISTRIBUTES coherence (not protects). B's
measurement destroys the nonlocal reservoir, cutting the return flow.

Script:  simulations/shift_mechanism.py
Output:  simulations/results/shift_mechanism.txt
Docs:    experiments/OBSERVER_GRAVITY_BRIDGE.md §7
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "shift_mechanism.txt")
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

N = 2; d = 4; d2 = 16
sysy = np.kron(sy, sy)


def site_op(op, k, nq=N):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
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
    r = v.reshape(d, d)
    return (r + r.conj().T) / 2


def ptrace_A(rho):
    return np.trace(rho.reshape(2, 2, 2, 2), axis1=1, axis2=3)


def local_coherence_A(rho):
    rA = ptrace_A(rho)
    return float(np.sum(np.abs(rA)) - np.sum(np.abs(np.diag(rA))))


def local_purity_A(rho):
    rA = ptrace_A(rho)
    return float(np.real(np.trace(rA @ rA)))


def concurrence(rho):
    R = rho @ sysy @ rho.conj() @ sysy
    eigvals = np.sort(np.real(np.linalg.eigvals(R)))[::-1]
    eigvals = np.maximum(eigvals, 0.0)
    sq = np.sqrt(eigvals)
    return max(0.0, sq[0] - sq[1] - sq[2] - sq[3])


def nonlocal_coherence(rho):
    """Off-diagonal weight in the 2-qubit state beyond local contributions."""
    total_l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    rA = ptrace_A(rho)
    rB = np.trace(rho.reshape(2, 2, 2, 2), axis1=0, axis2=2)
    local_A = float(np.sum(np.abs(rA)) - np.sum(np.abs(np.diag(rA))))
    local_B = float(np.sum(np.abs(rB)) - np.sum(np.abs(np.diag(rB))))
    return max(0.0, total_l1 - local_A - local_B)


def local_cpsi_A(rho):
    rA = ptrace_A(rho)
    purity = float(np.real(np.trace(rA @ rA)))
    l1 = float(np.sum(np.abs(rA)) - np.sum(np.abs(np.diag(rA))))
    return purity * l1


def apply_B_measurement_Z(rho):
    P0 = site_op(np.outer(up, up.conj()), 1)
    P1 = site_op(np.outer(dn, dn.conj()), 1)
    return P0 @ rho @ P0.conj().T + P1 @ rho @ P1.conj().T


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def find_crossing_time(L, rho0, dt=0.001, t_max=30.0):
    prev = local_cpsi_A(rho0)
    if prev < 0.25: return None
    for i in range(1, int(t_max / dt) + 1):
        t = i * dt
        rho = evolve(L, rho0, t)
        val = local_cpsi_A(rho)
        if val < 0.25 and prev >= 0.25:
            frac = (0.25 - prev) / (val - prev)
            return (i - 1 + frac) * dt
        prev = val
    return None


# ============================================================
# SETUP
# ============================================================
log("=" * 70)
log("Shift Mechanism: Coherence Reservoir and B's Measurement")
log("=" * 70)
log()

J = 0.5
gamma = 0.05
t_B = 1.0
H = build_H(J)
L = build_L(H, gamma)

psi0 = np.kron(plus, plus)
rho0 = ket2dm(psi0)

# ============================================================
# TEST 1: At the moment of B's measurement
# ============================================================
log("Test 1: Observable changes at t_B = 1.0 (J=0.5, γ=0.05)")
log("-" * 70)

rho_tB = evolve(L, rho0, t_B)
rho_meas = apply_B_measurement_Z(rho_tB)

obs = [
    ("Local coherence A", local_coherence_A),
    ("Local purity A", local_purity_A),
    ("Concurrence", concurrence),
    ("Nonlocal coherence", nonlocal_coherence),
]

log(f"  {'Observable':>25}  {'Before B':>9}  {'After B':>9}  {'Δ':>10}")
log("  " + "-" * 56)

for name, fn in obs:
    before = fn(rho_tB)
    after = fn(rho_meas)
    delta = after - before
    log(f"  {name:>25}  {before:9.3f}  {after:9.3f}  {delta:10.3f}")

# ============================================================
# TEST 2: Aftermath — A's coherence decay comparison
# ============================================================
log()
log("Test 2: A's coherence after B's measurement vs silent")
log("-" * 70)
log(f"  {'t after B':>10}  {'Coh_A(silent)':>14}  {'Coh_A(meas)':>12}  {'Ratio':>7}")

for dt_after in [0.0, 0.5, 1.0, 2.0, 5.0]:
    rho_s = evolve(L, rho_tB, dt_after)
    rho_m = evolve(L, rho_meas, dt_after)
    coh_s = local_coherence_A(rho_s)
    coh_m = local_coherence_A(rho_m)
    ratio = coh_m / coh_s if coh_s > 1e-10 else 0
    log(f"  {dt_after:10.1f}  {coh_s:14.3f}  {coh_m:12.3f}  {ratio:7.2f}")

# ============================================================
# TEST 3: Timing dependence — sweep t_B
# ============================================================
log()
log("Test 3: Damage depends on WHEN B measures (J=0.5, γ=0.05)")
log("-" * 70)

# Baseline: A's crossing time without any B measurement
t_silent = find_crossing_time(L, rho0)
log(f"  Baseline t_cross(silent) = {t_silent:.4f}")
log()

t_B_values = [0.01, 0.10, 0.20, 0.50, 1.00, 2.00, 3.00, 5.00]
log(f"  {'t_B':>6}  {'NL coh':>8}  {'Remaining lifetime':>20}  {'Damage':>8}")

for tb in t_B_values:
    rho_tb = evolve(L, rho0, tb)
    nl_coh = nonlocal_coherence(rho_tb)
    rho_m = apply_B_measurement_Z(rho_tb)
    t_m = find_crossing_time(L, rho_m)
    if t_m is not None:
        total = tb + t_m
        remaining = total / t_silent * 100 if t_silent else 0
        damage = 100 - remaining
    else:
        remaining = 0
        damage = 100
    log(f"  {tb:6.2f}  {nl_coh:8.3f}  {remaining:19.1f}%  {damage:7.1f}%")

# ============================================================
# TEST 4: Regeneration after measurement
# ============================================================
log()
log("Test 4: Concurrence regeneration after B's measurement")
log("-" * 70)

log(f"  {'t after B':>10}  {'Concurrence':>12}")

for dt_after in [0.0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00]:
    rho_m = evolve(L, rho_meas, dt_after)
    c = concurrence(rho_m)
    arrow = "↑" if dt_after > 0 and c > concurrence(evolve(L, rho_meas, max(0, dt_after - 0.25))) else "↓" if c < concurrence(evolve(L, rho_meas, max(0, dt_after - 0.25))) else ""
    log(f"  {dt_after:10.2f}  {c:12.3f} {arrow}")

# ============================================================
# TEST 5: Coupling vs isolation
# ============================================================
log()
log("Test 5: Coupling accelerates local crossing vs isolated qubit")
log("-" * 70)

# Single |+⟩ qubit: build 1-qubit Liouvillian
H1 = np.zeros((2, 2), dtype=complex)
L1 = -1j * (np.kron(H1, I2) - np.kron(I2, H1.T))
L1 += gamma * (np.kron(sz, sz.conj()) - np.eye(4))
rho0_single = np.outer(plus, plus.conj())
# Find single-qubit crossing
prev = float(np.real(np.trace(rho0_single @ rho0_single)) *
             (np.sum(np.abs(rho0_single)) - np.sum(np.abs(np.diag(rho0_single)))))
t_single = None
for i in range(1, 30001):
    t = i * 0.001
    v = expm(L1 * t) @ rho0_single.flatten()
    rho_t = v.reshape(2, 2)
    rho_t = (rho_t + rho_t.conj().T) / 2
    pur = float(np.real(np.trace(rho_t @ rho_t)))
    l1 = float(np.sum(np.abs(rho_t)) - np.sum(np.abs(np.diag(rho_t))))
    val = pur * l1
    if val < 0.25 and prev >= 0.25:
        frac = (0.25 - prev) / (val - prev)
        t_single = (i - 1 + frac) * 0.001
        break
    prev = val

log(f"  {'System':>18}  {'t_cross_A':>10}  {'vs single':>10}")
log("  " + "-" * 42)
log(f"  {'Single |+⟩':>18}  {t_single:10.3f}  {'1.00x':>10}")

for J_val in [0.05, 0.10, 0.50, 1.00]:
    H_v = build_H(J_val)
    L_v = build_L(H_v, gamma)
    rho0_pp = ket2dm(np.kron(plus, up))  # |+,0⟩
    t_c = find_crossing_time(L_v, rho0_pp)
    if t_c and t_single:
        ratio = t_c / t_single
        log(f"  {'|+,0⟩ J=' + str(J_val):>18}  {t_c:10.3f}  {ratio:.2f}x")
    else:
        log(f"  {'|+,0⟩ J=' + str(J_val):>18}  {'NEVER':>10}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("Coupling REDISTRIBUTES coherence, not protects.")
log("B's measurement destroys nonlocal reservoir, cutting return flow.")
log("Damage is timing-dependent (oscillation phase matters).")
log("Regeneration peaks briefly, then dies.")
log("=" * 70)

_outf.close()
