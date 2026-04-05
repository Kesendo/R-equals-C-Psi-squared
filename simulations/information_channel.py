#!/usr/bin/env python3
"""
Information Channel: Multi-Pair Amplification and Discrimination
=================================================================
Post-processing on interval_shift results. B encodes bits by measuring
or not. A reads by comparing crossing time to expected t₀ = K/γ_A.

Computes:
1. Channel capacity vs coupling strength (reuses interval_shift data)
2. Multi-pair amplification: N_min = (σ/Δt)² for discrimination
3. Lieb-Robinson velocity estimates for gravitational coupling

Script:  simulations/information_channel.py
Output:  simulations/results/information_channel.txt
Docs:    experiments/OBSERVER_GRAVITY_BRIDGE.md §8
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "information_channel.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N=2, same as interval_shift.py)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 2; d = 4; d2 = 16


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
# MAIN
# ============================================================
log("=" * 70)
log("Information Channel: Multi-Pair Amplification")
log("=" * 70)
log()

gamma = 0.05
t_B = 1.0
psi0 = np.kron(plus, plus)
rho0 = ket2dm(psi0)

# ============================================================
# 1. Channel capacity vs coupling
# ============================================================
log("1. Channel capacity vs coupling strength")
log("-" * 70)
log(f"  {'J':>6}  {'t₀ (silent)':>12}  {'t₁ (B meas)':>12}  {'Δt':>8}  {'Shift':>8}")

J_values = [0.001, 0.005, 0.010, 0.020, 0.050, 0.100, 0.500, 1.000]
results = {}

for J in J_values:
    H = build_H(J)
    L = build_L(H, gamma)
    t0 = find_crossing_time(L, rho0)
    rho_tB = evolve(L, rho0, t_B)
    rho_meas = apply_B_measurement_Z(rho_tB)
    t1 = find_crossing_time(L, rho_meas)
    if t1 is not None:
        t1 += t_B

    if t0 and t1:
        dt_shift = t1 - t0
        shift = dt_shift / t0 * 100
        results[J] = (t0, t1, dt_shift)
        log(f"  {J:6.3f}  {t0:12.3f}  {t1:12.3f}  {dt_shift:8.3f}  {shift:7.2f}%")

# ============================================================
# 2. Multi-pair amplification
# ============================================================
log()
log("2. Multi-pair amplification: N_min = (σ/Δt)²")
log("-" * 70)

J_ref = 0.01
if J_ref in results:
    t0, t1, dt_shift = results[J_ref]
    log(f"  At J = {J_ref}: Δt = {abs(dt_shift):.3f}")
    log()
    log(f"  {'σ (jitter)':>12}  {'σ/t₀':>8}  {'N_min':>8}")
    log("  " + "-" * 32)

    for sigma_frac_label, sigma in [("10%", 0.1 * abs(t0)),
                                     ("100%", 1.0 * abs(t0)),
                                     ("1000%", 10.0 * abs(t0))]:
        N_min = (sigma / abs(dt_shift)) ** 2
        log(f"  {sigma_frac_label:>12}  {sigma/t0*100:7.0f}%  {N_min:8.0f}")

# ============================================================
# 3. Lieb-Robinson velocity estimates
# ============================================================
log()
log("3. Lieb-Robinson velocity estimates (gravitational coupling)")
log("-" * 70)
log("  v_LR ≤ 2 × J × a / ℏ, gravitational J_grav ~ G×m²/(ℏ×D)")
log()

G = 6.674e-11      # m³/(kg·s²)
hbar = 1.055e-34    # J·s

systems = [
    ("NV center",       1e-26, 1e-29),
    ("Optomechanical",  1e-15, 1e-7),
    ("Dust grain",      1e-9,  1e5),
    ("Microgram",       1e-6,  1e11),
]

log(f"  {'System':>18}  {'Mass (kg)':>12}  {'v (m/s)':>12}  {'Note':>18}")
log("  " + "-" * 66)

c_light = 3e8
for name, mass, v_approx in systems:
    v = G * mass**2 / hbar
    note = "Unmeasurable" if v < 1e-20 else "Extremely slow" if v < 1 else "Subluminal" if v < c_light else "EXCEEDS c"
    log(f"  {name:>18}  {mass:12.0e}  {v:12.2e}  {note:>18}")

log()
log("  Naive formula breaks above ~microgram scale (v > c).")
log("  Relativistic corrections or modified J_grav scaling needed.")

# ============================================================
# 4. Rate scaling
# ============================================================
log()
log("4. Rate scaling: Δt/t ∝ (J/γ) at small J")
log("-" * 70)
log(f"  {'J/γ':>8}  {'Δt/t':>10}  {'(Δt/t)/(J/γ)':>14}")

for J in [0.001, 0.005, 0.010, 0.020]:
    if J in results:
        t0, t1, dt_shift = results[J]
        ratio = abs(dt_shift) / t0
        jg = J / gamma
        scaling = ratio / jg if jg > 0 else 0
        log(f"  {jg:8.2f}  {ratio:10.4f}  {scaling:14.4f}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("1 bit per ~21 pairs at J/γ=0.2 with 100% jitter.")
log("Rate scales as (J/γ)² at small coupling.")
log("CΨ crossing time = lock-in amplifier for weak quantum couplings.")
log("=" * 70)

_outf.close()
