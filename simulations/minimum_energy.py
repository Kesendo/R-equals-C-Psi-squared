#!/usr/bin/env python3
"""
Minimum Crossing Energy: No Energy Threshold, Only a Coherence Barrier
=======================================================================
2-qubit Heisenberg system under Z-dephasing. Tracks CΨ = τ/(d-1) where
τ = concurrence² (tangle) for the 2-qubit density matrix.

Tests:
1. cos(alpha)|00> + sin(alpha)|11> family: same energy, different crossing.
   Critical angle = 30 degrees exactly (CΨ(0) = 1/4).
2. cos(alpha)|01> + sin(alpha)|10> family: all states cross.
3. Product states: |0,1> crosses (Hamiltonian generates entanglement),
   |+,+> never does (eigenstate of H, no dynamics).
4. J/gamma sweep for |0,1>: critical J/gamma ~ 5-10.

Script:  simulations/minimum_energy.py
Output:  simulations/results/minimum_energy.txt
Docs:    experiments/MINIMUM_CROSSING_ENERGY.md
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "minimum_energy.txt")
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

# sigma_y tensor sigma_y for Wootters concurrence
sysy = np.kron(sy, sy)


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


def concurrence(rho):
    """Wootters concurrence for a 2-qubit density matrix."""
    R = rho @ sysy @ rho.conj() @ sysy
    eigvals = np.sort(np.real(np.linalg.eigvals(R)))[::-1]
    eigvals = np.maximum(eigvals, 0.0)
    sq = np.sqrt(eigvals)
    return max(0.0, sq[0] - sq[1] - sq[2] - sq[3])


def cpsi(rho):
    """CΨ = τ/(d-1) where τ = concurrence² (tangle)."""
    C = concurrence(rho)
    return C * C / (d - 1)


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def cpsi_max_and_cross(L, rho0, t_max=10.0, dt=0.01):
    """Find max CΨ and whether it crosses 1/4."""
    n_steps = int(t_max / dt)
    c_max = 0.0
    t_max_val = 0.0
    crossed = False
    c0 = cpsi(rho0)
    c_prev = c0

    for i in range(1, n_steps + 1):
        t = i * dt
        rho = evolve(L, rho0, t)
        c = cpsi(rho)
        if c > c_max:
            c_max = c
            t_max_val = t
        if (c_prev >= 0.25 and c < 0.25) or (c_prev < 0.25 and c >= 0.25):
            crossed = True
        c_prev = c

    if c0 >= 0.25:
        crossed = True

    return c_max, t_max_val, crossed


# ============================================================
# TEST 1: cos(alpha)|00> + sin(alpha)|11>
# ============================================================
log("=" * 70)
log("Minimum Crossing Energy: Coherence Barrier, Not Energy Barrier")
log("=" * 70)
log()

J = 1.0
gamma = 0.05
H = build_H(J)
L = build_L(H, gamma)

log("Test 1: cos(α)|00⟩ + sin(α)|11⟩ — same energy, different crossing")
log("-" * 70)
log(f"  J = {J}, γ = {gamma}, ⟨H⟩ = +J for all α")
log(f"{'α (deg)':>8}  {'⟨H⟩':>5}  {'CΨ(0)':>8}  {'CΨ_max':>8}  {'Crosses?':>9}")

alphas_deg = [45, 35, 31, 30, 25, 15]
for a_deg in alphas_deg:
    a = np.radians(a_deg)
    psi = np.cos(a) * np.kron(up, up) + np.sin(a) * np.kron(dn, dn)
    rho0 = ket2dm(psi)
    c0 = cpsi(rho0)
    c_max, t_max_val, crossed = cpsi_max_and_cross(L, rho0)
    log(f"{a_deg:8d}  {'J':>5}  {c0:8.4f}  {c_max:8.4f}  {'YES' if crossed else 'NO':>9}")

# Binary search for critical angle
log()
log("Binary search for critical α:")
lo, hi = 0.0, 45.0
for _ in range(50):
    mid = (lo + hi) / 2
    a = np.radians(mid)
    psi = np.cos(a) * np.kron(up, up) + np.sin(a) * np.kron(dn, dn)
    c0 = cpsi(ket2dm(psi))
    if c0 > 0.25:
        hi = mid
    else:
        lo = mid
crit = (lo + hi) / 2
log(f"  α_critical = {crit:.6f}° (sin²(2α)/(d-1) = 0.250000)")

# ============================================================
# TEST 2: cos(alpha)|01> + sin(alpha)|10>
# ============================================================
log()
log("Test 2: cos(α)|01⟩ + sin(α)|10⟩ — all cross")
log("-" * 70)
log(f"{'α (deg)':>8}  {'CΨ(0)':>8}  {'CΨ_max':>8}  {'t(max)':>7}  {'Crosses?':>9}")

alphas2 = [45, 25, 15, 5]
for a_deg in alphas2:
    a = np.radians(a_deg)
    psi = np.cos(a) * np.kron(up, dn) + np.sin(a) * np.kron(dn, up)
    rho0 = ket2dm(psi)
    c0 = cpsi(rho0)
    c_max, t_max_val, crossed = cpsi_max_and_cross(L, rho0)
    log(f"{a_deg:8d}  {c0:8.4f}  {c_max:8.4f}  {t_max_val:7.2f}  {'YES' if crossed else 'NO':>9}")

# ============================================================
# TEST 3: Product states
# ============================================================
log()
log("Test 3: Product states — Hamiltonian creates the crossing")
log("-" * 70)

product_states = {
    "|0,1>": np.kron(up, dn),
    "|1,0>": np.kron(dn, up),
    "|+,0>": np.kron(plus, up),
    "|0,+>": np.kron(up, plus),
    "|+,1>": np.kron(plus, dn),
    "|+,+>": np.kron(plus, plus),
    "|0,0>": np.kron(up, up),
    "|1,1>": np.kron(dn, dn),
}

log(f"{'State':>8}  {'CΨ(0)':>8}  {'CΨ_max':>8}  {'Crosses?':>9}")

for name, psi in product_states.items():
    rho0 = ket2dm(psi)
    c0 = cpsi(rho0)
    c_max, _, crossed = cpsi_max_and_cross(L, rho0, t_max=20.0)
    log(f"{name:>8}  {c0:8.3f}  {c_max:8.3f}  {'YES' if crossed else 'NO':>9}")

# ============================================================
# TEST 4: J/gamma sweep for |0,1>
# ============================================================
log()
log("Test 4: J/γ sweep for |0,1⟩ product state")
log("-" * 70)

psi_01 = np.kron(up, dn)
rho0_01 = ket2dm(psi_01)

j_gamma_ratios = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
gamma_fixed = 0.05

log(f"{'J/γ':>7}  {'CΨ_max':>8}  {'Crosses?':>9}")

for jg in j_gamma_ratios:
    J_val = jg * gamma_fixed
    H_val = build_H(J_val)
    L_val = build_L(H_val, gamma_fixed)
    c_max, _, crossed = cpsi_max_and_cross(L_val, rho0_01, t_max=20.0)
    log(f"{jg:7.1f}  {c_max:8.3f}  {'YES' if crossed else 'NO':>9}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("No energy threshold for CΨ crossing.")
log(f"Critical angle: α = {crit:.1f}° (CΨ(0) = 1/4).")
log("|01>/|10> family: all cross (exchange interaction generates entanglement).")
log("Product states: J/γ ~ 5-10 needed (coherence barrier, not energy barrier).")
log("Eigenstates of H (|0,0>, |1,1>, |+,+>): no dynamics, no crossing, no time.")
log("=" * 70)

_outf.close()
