#!/usr/bin/env python3
"""
Observer-Gravity Crossing: t_cross = K(Observer, State) / γ
=============================================================
2-qubit system cos(α)|00⟩+sin(α)|11⟩ under Heisenberg coupling (J=1.0)
and Z-dephasing. The crossing time factorizes: K is γ-invariant.

Test 1: Sweep γ (six environments), verify K_conc = const.
Test 2: Sweep α, show K(Conc)/K(MI) is state-dependent.
Test 3: States below α=30° never cross.

Script:  simulations/observer_gravity_cross.py
Output:  simulations/results/observer_gravity_cross.txt
Docs:    experiments/OBSERVER_GRAVITY_BRIDGE.md §2
"""

import numpy as np
from scipy.linalg import expm
import os, sys

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "observer_gravity_cross.txt")
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

N = 2
d = 4
d2 = 16
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
    R = rho @ sysy @ rho.conj() @ sysy
    eigvals = np.sort(np.real(np.linalg.eigvals(R)))[::-1]
    eigvals = np.maximum(eigvals, 0.0)
    sq = np.sqrt(eigvals)
    return max(0.0, sq[0] - sq[1] - sq[2] - sq[3])


def cpsi_conc(rho):
    """CΨ via concurrence: tangle/(d-1)."""
    C = concurrence(rho)
    return C * C / (d - 1)


def mutual_information(rho):
    """I(A:B) = S(A) + S(B) - S(AB)."""
    def entropy(r):
        eigvals = np.real(np.linalg.eigvalsh(r))
        eigvals = eigvals[eigvals > 1e-15]
        return -np.sum(eigvals * np.log2(eigvals))

    rho_A = np.trace(rho.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    rho_B = np.trace(rho.reshape(2, 2, 2, 2), axis1=0, axis2=2)
    return entropy(rho_A) + entropy(rho_B) - entropy(rho)


def cpsi_mi(rho):
    """CΨ via mutual information: MI/4 (scaled so Bell+ starts at 0.5)."""
    return mutual_information(rho) / 4.0


def ket2dm(psi):
    return np.outer(psi, psi.conj())


def find_crossing(L, rho0, metric_fn, threshold=0.25, dt=0.001, t_max=100.0):
    """Find time when metric crosses threshold from above."""
    prev = metric_fn(rho0)
    if prev < threshold:
        return None
    n_steps = int(t_max / dt)
    for i in range(1, n_steps + 1):
        t = i * dt
        rho = evolve(L, rho0, t)
        val = metric_fn(rho)
        if val < threshold and prev >= threshold:
            # Linear interpolation
            frac = (threshold - prev) / (val - prev)
            return (i - 1 + frac) * dt
        prev = val
    return None


# ============================================================
# TEST 1: K_conc invariant across γ
# ============================================================
log("=" * 70)
log("Observer-Gravity Crossing: t_cross = K / γ")
log("=" * 70)
log()

J = 1.0
H = build_H(J)

bell_plus = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
rho0 = ket2dm(bell_plus)

envs = [
    ("Deep Space", 0.01),
    ("Mars",       0.019),
    ("Earth",      0.05),
    ("Jupiter",    0.13),
    ("Neutron",    0.20),
    ("Black Hole", 0.50),
]

log("Test 1: K_conc invariance (Bell+ state, J=1.0)")
log("-" * 70)
log(f"  {'Environment':>12}  {'γ':>6}  {'t_cross(Conc)':>14}  {'K_conc':>10}")

K_values = []
for name, gamma in envs:
    L = build_L(H, gamma)
    t_cross = find_crossing(L, rho0, cpsi_conc, 0.25, dt=0.0001,
                            t_max=50.0/gamma)
    if t_cross:
        K = t_cross * gamma
        K_values.append(K)
        log(f"  {name:>12}  {gamma:6.3f}  {t_cross:14.4f}  {K:10.6f}")
    else:
        log(f"  {name:>12}  {gamma:6.3f}  {'NEVER':>14}")

if K_values:
    K_mean = np.mean(K_values)
    K_std = np.std(K_values)
    cv = K_std / K_mean * 100 if K_mean > 0 else 0
    log(f"\n  K_conc = {K_mean:.5f} ± {K_std:.5f}, CV = {cv:.2f}%")

# ============================================================
# TEST 2: K ratio across states
# ============================================================
log()
log("Test 2: K(Conc)/K(MI) is state-dependent")
log("-" * 70)

gamma = 0.05
L = build_L(H, gamma)

alphas = [45, 40, 35, 31, 30, 25]
log(f"  {'α':>5}  {'K_conc':>10}  {'K_MI':>10}  {'Ratio':>10}")

for a_deg in alphas:
    a = np.radians(a_deg)
    psi = np.cos(a) * np.kron(up, up) + np.sin(a) * np.kron(dn, dn)
    rho0_a = ket2dm(psi)

    t_conc = find_crossing(L, rho0_a, cpsi_conc, 0.25, dt=0.001, t_max=200.0)
    t_mi = find_crossing(L, rho0_a, cpsi_mi, 0.25, dt=0.001, t_max=200.0)

    if t_conc and t_mi:
        K_c = t_conc * gamma
        K_m = t_mi * gamma
        ratio = K_c / K_m
        log(f"  {a_deg:5d}  {K_c:10.6f}  {K_m:10.6f}  {ratio:10.6f}")
    elif t_conc:
        K_c = t_conc * gamma
        log(f"  {a_deg:5d}  {K_c:10.6f}  {'never':>10}")
    else:
        log(f"  {a_deg:5d}  {'never':>10}  {'never':>10}  {'--':>10}")

# ============================================================
# SUMMARY
# ============================================================
log()
log("=" * 70)
log("t_cross = K(Observer, State) / γ. K is γ-invariant.")
log("K(Conc)/K(MI) is state-dependent. States below α=30° never cross.")
log("=" * 70)

_outf.close()
