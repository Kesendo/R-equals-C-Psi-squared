#!/usr/bin/env python3
"""
No-Signalling Boundary: CΨ Sees a Regime Change That No Local Observer Can Detect
==================================================================================
Bell+ pair. B measures Z (averaged over outcomes). A's reduced state rho_A
is unchanged (no-signalling). But CΨ = C × Ψ drops from 0.500 to 0.250,
exactly onto the 1/4 boundary.

Convention: C = Tr(rho_AB^2) (global purity), Ψ = max eigenvalue of rho_A.

Script:  simulations/test2_no_signalling.py
Output:  simulations/results/test2_no_signalling.txt
Docs:    experiments/NO_SIGNALLING_BOUNDARY.md
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "test2_no_signalling.txt")
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
    """Partial trace over B, keeping A."""
    reshaped = rho.reshape(2, 2, 2, 2)
    return np.trace(reshaped, axis1=1, axis2=3)


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def purity(rho):
    return np.real(np.trace(rho @ rho))


def von_neumann_entropy(rho):
    eigvals = np.real(np.linalg.eigvalsh(rho))
    eigvals = eigvals[eigvals > 1e-15]
    return -np.sum(eigvals * np.log2(eigvals))


def apply_B_measurement_Z(rho):
    """B measures in Z basis (averaged over outcomes)."""
    P0 = site_op(np.outer(up, up.conj()), 1)
    P1 = site_op(np.outer(dn, dn.conj()), 1)
    return P0 @ rho @ P0.conj().T + P1 @ rho @ P1.conj().T


def cpsi(rho_AB):
    """CΨ = Tr(rho_AB^2) × max eigenvalue of rho_A."""
    C = purity(rho_AB)
    rho_A = ptrace_A(rho_AB)
    Psi = np.max(np.real(np.linalg.eigvalsh(rho_A)))
    return C * Psi


def regime(val):
    if val > 0.2501:
        return "quantum"
    elif val > 0.2499:
        return "boundary"
    else:
        return "classical"


# ============================================================
# TEST 1: Static (no evolution) — B measures Z on Bell+
# ============================================================
log("=" * 70)
log("No-Signalling Boundary: CΨ Regime Change")
log("=" * 70)
log()

bell_plus = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
rho_before = ket2dm(bell_plus)
rho_after = apply_B_measurement_Z(rho_before)

rho_A_before = ptrace_A(rho_before)
rho_A_after = ptrace_A(rho_after)
delta = np.linalg.norm(rho_A_before - rho_A_after)

log("Test 1: Bell+ pair, B measures Z (averaged over outcomes)")
log("-" * 70)
log(f"  ||Δρ_A|| = {delta:.10f}")
log()

# Full observable table
log(f"  {'Quantity':>20}  {'Before':>10}  {'After':>10}  {'Changed?':>9}  {'A sees?':>8}")
log("  " + "-" * 62)

observables = [
    ("rho_A", "I/2", "I/2", delta < 1e-12, False),
]

# Pauli expectations on A
for name, P in [("⟨σx⟩_A", sx), ("⟨σy⟩_A", sy), ("⟨σz⟩_A", sz)]:
    before = np.real(np.trace(rho_A_before @ P))
    after = np.real(np.trace(rho_A_after @ P))
    changed = abs(before - after) > 1e-12
    log(f"  {name:>20}  {before:10.3f}  {after:10.3f}  {'YES' if changed else 'NO':>9}  {'-':>8}")

# Purity of A
pur_A_before = purity(rho_A_before)
pur_A_after = purity(rho_A_after)
log(f"  {'Purity(rho_A)':>20}  {pur_A_before:10.3f}  {pur_A_after:10.3f}  "
    f"{'YES' if abs(pur_A_before - pur_A_after) > 1e-12 else 'NO':>9}  {'-':>8}")

# Entropy of A
S_before = von_neumann_entropy(rho_A_before)
S_after = von_neumann_entropy(rho_A_after)
log(f"  {'S(rho_A)':>20}  {S_before:10.3f}  {S_after:10.3f}  "
    f"{'YES' if abs(S_before - S_after) > 1e-12 else 'NO':>9}  {'-':>8}")

# Ψ = max eigenvalue of rho_A
Psi_before = np.max(np.real(np.linalg.eigvalsh(rho_A_before)))
Psi_after = np.max(np.real(np.linalg.eigvalsh(rho_A_after)))
log(f"  {'Ψ (max eig rho_A)':>20}  {Psi_before:10.3f}  {Psi_after:10.3f}  "
    f"{'YES' if abs(Psi_before - Psi_after) > 1e-12 else 'NO':>9}  {'-':>8}")

# C = Tr(rho_AB^2)
C_before = purity(rho_before)
C_after = purity(rho_after)
log(f"  {'C = Tr(rho_AB²)':>20}  {C_before:10.3f}  {C_after:10.3f}  "
    f"{'YES' if abs(C_before - C_after) > 1e-12 else 'NO':>9}  {'NO':>8}")

# CΨ
cpsi_before = C_before * Psi_before
cpsi_after = C_after * Psi_after
log(f"  {'CΨ':>20}  {cpsi_before:10.3f}  {cpsi_after:10.3f}  "
    f"{'YES' if abs(cpsi_before - cpsi_after) > 1e-12 else 'NO':>9}  {'NO':>8}")

# Regime
log(f"  {'Regime':>20}  {regime(cpsi_before):>10}  {regime(cpsi_after):>10}  "
    f"{'YES':>9}  {'NO':>8}")

log()
if delta < 1e-12:
    log("PASS: No-signalling holds exactly. rho_A unchanged.")
    log(f"  CΨ: {cpsi_before:.3f} → {cpsi_after:.3f} ({regime(cpsi_before)} → {regime(cpsi_after)})")
else:
    log(f"FAIL: ||Δρ_A|| = {delta:.2e}")

# ============================================================
# TEST 2: Time evolution with dephasing
# ============================================================
log()
log("Test 2: Bell+ under dephasing, B measures Z at t=2")
log("-" * 70)

gamma = 0.05
J = 1.0
H = build_H(J)
L = build_L(H, gamma)

rho0 = ket2dm(bell_plus)

# Evolve to several time points, apply B measurement, compare rho_A
log(f"  {'t':>6}  {'CΨ(no meas)':>12}  {'CΨ(B→Z)':>10}  {'||Δρ_A||':>12}  {'Regime(no)':>11}  {'Regime(Z)':>10}")

check_times = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
for t in check_times:
    rho_t = evolve(L, rho0, t)
    rho_Bz = apply_B_measurement_Z(rho_t)

    rho_A_no = ptrace_A(rho_t)
    rho_A_Bz = ptrace_A(rho_Bz)
    diff = np.linalg.norm(rho_A_no - rho_A_Bz)

    cpsi_no = cpsi(rho_t)
    cpsi_Bz = cpsi(rho_Bz)

    log(f"  {t:6.1f}  {cpsi_no:12.4f}  {cpsi_Bz:10.4f}  {diff:12.2e}  "
        f"{regime(cpsi_no):>11}  {regime(cpsi_Bz):>10}")

# ============================================================
# TEST 3: CΨ trajectory comparison (with/without measurement)
# ============================================================
log()
log("Test 3: CΨ trajectory — measurement at t=2, then continue evolving")
log("-" * 70)

rho_t2 = evolve(L, rho0, 2.0)
rho_t2_Bz = apply_B_measurement_Z(rho_t2)

log(f"  {'t':>6}  {'CΨ(no meas)':>12}  {'CΨ(meas@t=2)':>14}  {'ΔCΨ':>8}")

for dt_after in [0.0, 0.5, 1.0, 2.0, 5.0]:
    rho_no = evolve(L, rho_t2, dt_after)
    rho_Bz = evolve(L, rho_t2_Bz, dt_after)

    c_no = cpsi(rho_no)
    c_Bz = cpsi(rho_Bz)

    log(f"  {2.0 + dt_after:6.1f}  {c_no:12.4f}  {c_Bz:14.4f}  {c_Bz - c_no:8.4f}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("CΨ sees the regime change (0.500 → 0.250). A cannot see it.")
log("C drops (1.0 → 0.5). Ψ stays (0.5 → 0.5). rho_A unchanged.")
log("Dynamic bridge eliminated: B cannot send new information to A.")
log("=" * 70)

_outf.close()
