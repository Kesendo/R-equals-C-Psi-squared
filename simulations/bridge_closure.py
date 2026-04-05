#!/usr/bin/env python3
"""
Bridge Closure: Pre-Shared Entanglement Without a Channel = Shared Randomness
==============================================================================
Verifies the no-signalling consequence: after separation (J=0), A's reduced
state rho_A is identical regardless of B's actions.

Tests:
1. Bell+ pair evolved under independent dephasing, then branch into three
   scenarios (B does nothing, B measures Z, B measures X). rho_A must be
   identical across all branches at every time step.
2. Product states |++>, |+0>, |+-> all have identical rho_A = |+><+|.

Expected: max ||Delta rho_A|| = 0.0 at machine precision for all branches.

Script:  simulations/bridge_closure.py
Output:  simulations/results/bridge_closure.txt
Docs:    experiments/BRIDGE_CLOSURE.md
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "bridge_closure.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ============================================================
# INFRASTRUCTURE (N=2, J=0)
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)
minus = (up - dn) / np.sqrt(2)

N = 2
d = 4   # 2^2
d2 = 16  # d^2


def site_op(op, k, nq=N):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_L_no_coupling(gamma):
    """Liouvillian with J=0 (no Hamiltonian), Z-dephasing on both qubits."""
    H = np.zeros((d, d), dtype=complex)
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
    """Partial trace over B (qubit 1), keeping A (qubit 0)."""
    reshaped = rho.reshape(2, 2, 2, 2)
    return np.trace(reshaped, axis1=1, axis2=3)


def ptrace_B(rho):
    """Partial trace over A (qubit 0), keeping B (qubit 1)."""
    reshaped = rho.reshape(2, 2, 2, 2)
    return np.trace(reshaped, axis1=0, axis2=2)


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def apply_B_measurement_Z(rho):
    """B measures in Z basis (averaged over outcomes): rho -> sum_k P_k rho P_k."""
    P0 = site_op(np.outer(up, up.conj()), 1)
    P1 = site_op(np.outer(dn, dn.conj()), 1)
    return P0 @ rho @ P0.conj().T + P1 @ rho @ P1.conj().T


def apply_B_measurement_X(rho):
    """B measures in X basis (averaged over outcomes)."""
    Pp = site_op(np.outer(plus, plus.conj()), 1)
    Pm = site_op(np.outer(minus, minus.conj()), 1)
    return Pp @ rho @ Pp.conj().T + Pm @ rho @ Pm.conj().T


# ============================================================
# TEST 1: Bell+ under three B-action scenarios
# ============================================================
log("=" * 70)
log("Bridge Closure: Pre-Shared Entanglement = Shared Randomness")
log("=" * 70)
log()

gamma = 0.05
L = build_L_no_coupling(gamma)

bell_plus = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
rho0 = ket2dm(bell_plus)

log("Test 1: Bell+ pair, J=0, gamma=0.05, three B-action scenarios")
log("-" * 70)
log(f"{'t_after':>8}  {'rho_A(nothing)':>15}  {'rho_A(B->Z)':>12}  {'rho_A(B->X)':>12}  {'max ||Delta||':>14}")

check_times = [0.0, 0.5, 1.0, 2.0, 5.0]
max_diff_all = 0.0

for t_check in check_times:
    # Evolve to t_check first
    rho_t = evolve(L, rho0, t_check)

    # Branch 1: B does nothing, continue evolving
    rho_nothing = rho_t

    # Branch 2: B measures Z at t_check
    rho_Bz = apply_B_measurement_Z(rho_t)

    # Branch 3: B measures X at t_check
    rho_Bx = apply_B_measurement_X(rho_t)

    # Compare rho_A across all three
    rA_nothing = ptrace_A(rho_nothing)
    rA_Bz = ptrace_A(rho_Bz)
    rA_Bx = ptrace_A(rho_Bx)

    d1 = np.linalg.norm(rA_nothing - rA_Bz)
    d2 = np.linalg.norm(rA_nothing - rA_Bx)
    d3 = np.linalg.norm(rA_Bz - rA_Bx)
    max_d = max(d1, d2, d3)
    max_diff_all = max(max_diff_all, max_d)

    purity_A = np.real(np.trace(rA_nothing @ rA_nothing))
    log(f"{t_check:8.1f}  {purity_A:15.6f}  {purity_A:12.6f}  {purity_A:12.6f}  {max_d:14.2e}")

log()
if max_diff_all < 1e-12:
    log("PASS: rho_A identical across all B-actions at machine precision.")
else:
    log(f"FAIL: max difference = {max_diff_all:.2e}")

# ============================================================
# TEST 2: Product states with identical rho_A
# ============================================================
log()
log("Test 2: Product states |++>, |+0>, |+-> all yield rho_A = |+><+|")
log("-" * 70)

states = {
    "|++>": np.kron(plus, plus),
    "|+0>": np.kron(plus, up),
    "|+->": np.kron(plus, minus),
}

log(f"{'State':>8}  {'rho_A[0,0]':>10}  {'rho_A[0,1]':>10}  {'rho_A[1,0]':>10}  {'rho_A[1,1]':>10}  {'== |+><+|?':>12}")

rho_plus = np.outer(plus, plus.conj())
all_match = True

for name, psi in states.items():
    rho = ket2dm(psi)
    rA = ptrace_A(rho)
    diff = np.linalg.norm(rA - rho_plus)
    match = diff < 1e-12
    all_match = all_match and match
    log(f"{name:>8}  {rA[0,0].real:10.6f}  {rA[0,1].real:10.6f}  {rA[1,0].real:10.6f}  {rA[1,1].real:10.6f}  {'YES' if match else 'NO':>12}")

log()
if all_match:
    log("PASS: A cannot distinguish product states differing only in B's qubit.")
else:
    log("FAIL: Some product states are distinguishable on A alone.")

# ============================================================
# TEST 3: Evolve and branch (extended, like reproduction code)
# ============================================================
log()
log("Test 3: Evolve Bell+ to t=2, branch, continue evolving")
log("-" * 70)

rho_t2 = evolve(L, rho0, 2.0)

# Branch from t=2: B does nothing vs B measures Z
check_steps = [0, 50, 100, 200, 400]
dt = 0.02
max_diff = 0.0

log(f"{'step':>6}  {'t_after_branch':>15}  {'||rho_A(nothing) - rho_A(B->Z)||':>35}")

for step in check_steps:
    t_after = step * dt
    rho_nothing = evolve(L, rho_t2, t_after)
    rho_Bz = apply_B_measurement_Z(rho_t2)
    rho_Bz_evol = evolve(L, rho_Bz, t_after)

    diff = np.linalg.norm(ptrace_A(rho_nothing) - ptrace_A(rho_Bz_evol))
    max_diff = max(max_diff, diff)
    log(f"{step:6d}  {t_after:15.2f}  {diff:35.2e}")

log()
if max_diff < 1e-12:
    log("PASS: All zero. Bridge is dead.")
else:
    log(f"FAIL: max difference = {max_diff:.2e}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("Bridge permanently closed for J=0.")
log("CΨ fingerprints require access to joint state rho_AB.")
log("A's qubit carries at most 1 bit; the schedule carries log2(N) bits.")
log("=" * 70)

_outf.close()
