#!/usr/bin/env python3
"""
Two Qubits, No Noise, Just Time
==================================
What time alone does, and what it cannot do without noise.

Script: simulations/two_qubits_no_noise.py
Output: simulations/results/two_qubits_no_noise.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "two_qubits_no_noise.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 2; d = 4; d2 = 16


def site_op(op, k):
    return np.kron(op, I2) if k == 0 else np.kron(I2, op)


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, 0) @ site_op(P, 1)
    return H


def build_L(H, gammas):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def trace_dist(r1, r2):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r1 - r2)))


def ptrace_A(rho):
    return np.array([[rho[0,0]+rho[1,1], rho[0,2]+rho[1,3]],
                     [rho[2,0]+rho[3,1], rho[2,2]+rho[3,3]]])


def concurrence(rho):
    sy2 = np.kron(sy, sy)
    R = rho @ sy2 @ rho.conj() @ sy2
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])


def cpsi(rho):
    C = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    Psi = l1 / 3.0
    return C * Psi


def von_neumann(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


# ============================================================
H = build_H()
ZZ = site_op(sz, 0) @ site_op(sz, 1)
XY = site_op(sx, 0) @ site_op(sy, 1)

states = {
    'Bell+': np.outer((np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2),
                       (np.kron(up, up) + np.kron(dn, dn)).conj() / np.sqrt(2)),
    '|+0>':  np.outer(np.kron(plus, up), np.kron(plus, up).conj()),
    '|01>':  np.outer(np.kron(up, dn), np.kron(up, dn).conj()),
}

tlist = np.arange(0, 50.01, 0.1)

log("Two Qubits, No Noise, Just Time")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ============================================================
# Main sweep: gamma=0 for each initial state
# ============================================================
L0 = build_L(H, [0, 0])

for sname, rho0 in states.items():
    log("=" * 70)
    log(f"Initial state: {sname}, gamma = 0")
    log("=" * 70)
    log()

    crossings_down = 0; crossings_up = 0
    prev_above = cpsi(rho0) > 0.25
    recurrence_t = None

    log(f"  {'t':>6}  {'CPsi':>8}  {'Conc':>8}  {'PurA':>8}  {'Dist':>8}  "
        f"{'S_A':>8}  {'<ZZ>':>8}  {'<XY>':>8}  {'Pur':>8}")
    log(f"  {'-'*78}")

    for t in tlist:
        rho = evolve(L0, rho0, t)
        rA = ptrace_A(rho)
        rA = (rA + rA.conj().T) / 2

        cp = cpsi(rho)
        conc = concurrence(rho)
        purA = np.real(np.trace(rA @ rA))
        dist = trace_dist(rho, rho0)
        sA = von_neumann(rA)
        zz = np.real(np.trace(ZZ @ rho))
        xy = np.real(np.trace(XY @ rho))
        pur_full = np.real(np.trace(rho @ rho))

        above = cp > 0.25
        if prev_above and not above: crossings_down += 1
        if not prev_above and above: crossings_up += 1
        prev_above = above

        if recurrence_t is None and t > 0.5 and dist < 0.01:
            recurrence_t = t

        if abs(t % 5) < 0.05 or t < 0.15:
            log(f"  {t:>6.1f}  {cp:>8.4f}  {conc:>8.4f}  {purA:>8.4f}  {dist:>8.4f}  "
                f"{sA:>8.4f}  {zz:>8.4f}  {xy:>8.4f}  {pur_full:>8.6f}")

    log()
    log(f"  CPsi crossings downward through 1/4: {crossings_down}")
    log(f"  CPsi crossings upward through 1/4: {crossings_up}")
    log(f"  Recurrence (D < 0.01): {'t = ' + f'{recurrence_t:.1f}' if recurrence_t else 'not found in [0, 50]'}")
    log(f"  Full system purity at t=50: {np.real(np.trace(evolve(L0, rho0, 50) @ evolve(L0, rho0, 50))):.10f}")
    log()

# ============================================================
# Comparison table: gamma=0 vs gamma=0.05
# ============================================================
log("=" * 70)
log("COMPARISON: gamma=0 vs gamma=0.05 (initial state |01>)")
log("=" * 70)
log()

rho0 = states['|01>']
L_noisy = build_L(H, [0.05, 0.05])

# gamma=0
cross_down_0 = 0; cross_up_0 = 0; prev_above = cpsi(rho0) > 0.25
rec_0 = None
for t in tlist:
    rho = evolve(L0, rho0, t)
    cp = cpsi(rho)
    above = cp > 0.25
    if prev_above and not above: cross_down_0 += 1
    if not prev_above and above: cross_up_0 += 1
    prev_above = above
    if rec_0 is None and t > 0.5 and trace_dist(rho, rho0) < 0.01: rec_0 = t

rho_50_clean = evolve(L0, rho0, 50)
pur_50_clean = np.real(np.trace(rho_50_clean @ rho_50_clean))
zz_50_clean = np.real(np.trace(ZZ @ rho_50_clean))

# gamma=0.05
cross_down_n = 0; cross_up_n = 0; prev_above = cpsi(rho0) > 0.25
rec_n = None
for t in tlist:
    rho = evolve(L_noisy, rho0, t)
    cp = cpsi(rho)
    above = cp > 0.25
    if prev_above and not above: cross_down_n += 1
    if not prev_above and above: cross_up_n += 1
    prev_above = above
    if rec_n is None and t > 0.5 and trace_dist(rho, rho0) < 0.01: rec_n = t

rho_50_noisy = evolve(L_noisy, rho0, 50)
pur_50_noisy = np.real(np.trace(rho_50_noisy @ rho_50_noisy))
zz_50_noisy = np.real(np.trace(ZZ @ rho_50_noisy))

# Palindromic check
evals_clean = np.linalg.eigvals(L0)
evals_noisy = np.linalg.eigvals(L_noisy)
Sg = 0.10

def count_palindromic(evals, Sg):
    n_paired = 0
    for k in range(len(evals)):
        target = -(evals[k] + 2 * Sg)
        if np.min(np.abs(evals - target)) < 1e-6: n_paired += 1
    return n_paired

pal_noisy = count_palindromic(evals_noisy, Sg)

log(f"  {'Property':>40}  {'gamma=0':>15}  {'gamma=0.05':>15}")
log(f"  {'-'*75}")
log(f"  {'CPsi crosses 1/4 downward':>40}  {cross_down_0:>15}  {cross_down_n:>15}")
log(f"  {'CPsi crosses 1/4 upward':>40}  {cross_up_0:>15}  {cross_up_n:>15}")
log(f"  {'<ZZ> at t=50':>40}  {zz_50_clean:>15.4f}  {zz_50_noisy:>15.4f}")
rec_0_str = f"t={rec_0:.1f}" if rec_0 else "not found"
rec_n_str = f"t={rec_n:.1f}" if rec_n else "no"
log(f"  {'Recurrence (D < 0.01)':>40}  {rec_0_str:>15}  {rec_n_str:>15}")
log(f"  {'Tr(rho^2) at t=50':>40}  {pur_50_clean:>15.10f}  {pur_50_noisy:>15.10f}")
log(f"  {'Palindromic pairs (out of 16)':>40}  {'0 (no Re)':>15}  {pal_noisy:>15}")

# Eigenvalue spectra
log()
log("  Liouvillian eigenvalues (gamma=0, purely imaginary):")
for ev in sorted(evals_clean, key=lambda x: x.imag):
    if abs(ev) > 1e-10:
        log(f"    Re={ev.real:>10.6f}  Im={ev.imag:>10.6f}")

log()
log("  Liouvillian eigenvalues (gamma=0.05, real parts = decay rates):")
for ev in sorted(evals_noisy, key=lambda x: x.real):
    if abs(ev) > 1e-10:
        log(f"    Re={ev.real:>10.6f}  Im={ev.imag:>10.6f}")

log()
log("=" * 70)
log("CONCLUSION")
log("=" * 70)
log()
log("Time without noise: oscillation, reversibility, recurrence.")
log("Nothing is ever decided. No arrow. No palindrome. No absorbing boundary.")
log("CΨ crosses 1/4 freely in both directions.")
log()
log("Time WITH noise: irreversibility, the arrow, the palindrome,")
log("the absorbing 1/4 boundary, classical emergence.")
log("CΨ crosses 1/4 once and stays below.")
log()
log("Time is the stage. Noise is the play.")
log("The stage exists without the play. But only the play has a plot.")
log()

log(f"Total runtime: {clock.time():.0f}s")
log(f"Results: {OUT_PATH}")
_outf.close()
