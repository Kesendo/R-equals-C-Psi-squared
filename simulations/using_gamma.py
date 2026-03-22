#!/usr/bin/env python3
"""
Using the Mediator: Practical γ Control
=========================================
Exp 1: DD simulation (reduced effective γ on specific qubits)
Exp 2: γ gradient (shaped noise profile)
Exp 3: γ modulation (AC on the gate)
Exp 4: Time-resolved decoder (detecting γ changes)
Exp 5: γ feedback loop (self-tuning transistor)

Script: simulations/using_gamma.py
Output: simulations/results/using_gamma.txt
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
import os, sys, time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "using_gamma.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

N = 5
d = 2 ** N
d2 = d * d


def site_op(op, k, nq=N):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H(J=1.0, bonds=None):
    if bonds is None:
        bonds = [(i, i+1) for i in range(N-1)]
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, j)
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


def ptrace_keep(rho, nq, keep):
    dim = 2 ** nq; nk = len(keep); dk = 2 ** nk
    traced = [k for k in range(nq) if k not in keep]
    rho_r = np.zeros((dk, dk), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            bi = [(i >> (nq-1-k)) & 1 for k in range(nq)]
            bj = [(j >> (nq-1-k)) & 1 for k in range(nq)]
            if all(bi[k] == bj[k] for k in traced):
                ki = sum(bi[keep[m]] << (nk-1-m) for m in range(nk))
                kj = sum(bj[keep[m]] << (nk-1-m) for m in range(nk))
                rho_r[ki, kj] += rho[i, j]
    return rho_r


def von_neumann(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


def mutual_info(rho, nq, kA, kB):
    kAB = sorted(set(kA) | set(kB))
    rA = ptrace_keep(rho, nq, kA)
    rB = ptrace_keep(rho, nq, kB)
    rAB = ptrace_keep(rho, nq, kAB)
    rA = (rA + rA.conj().T) / 2
    rB = (rB + rB.conj().T) / 2
    rAB = (rAB + rAB.conj().T) / 2
    return von_neumann(rA) + von_neumann(rB) - von_neumann(rAB)


def bell_A_0M_pp_B():
    psi_A = (np.kron(np.array([1,0]), np.array([1,0])) +
             np.kron(np.array([0,1]), np.array([0,1]))) / np.sqrt(2)
    plus = np.array([1, 1]) / np.sqrt(2)
    psi = np.kron(np.kron(psi_A, np.array([1, 0])), np.kron(plus, plus))
    return np.outer(psi, psi.conj())


# ============================================================
# EXPERIMENT 2: γ Gradient
# ============================================================
def run_exp_2():
    log("=" * 70)
    log("EXPERIMENT 2: gamma GRADIENT")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_A_0M_pp_B()
    t_measure = 5.0

    profiles = {
        'uniform':    [0.05, 0.05, 0.05, 0.05, 0.05],
        'V-shape':    [0.01, 0.03, 0.05, 0.03, 0.01],
        'inv-V':      [0.05, 0.03, 0.01, 0.03, 0.05],
        'ramp-up':    [0.01, 0.02, 0.03, 0.04, 0.05],
        'ramp-down':  [0.05, 0.04, 0.03, 0.02, 0.01],
        'gate-open':  [0.05, 0.05, 0.01, 0.05, 0.05],
        'gate-closed':[0.05, 0.05, 0.20, 0.05, 0.05],
        'quiet-recv': [0.05, 0.05, 0.05, 0.01, 0.01],
    }

    log(f"  {'Profile':>14}  {'MI(A:B)':>10}  {'MI(0,1:3,4)':>12}  {'Sg':>8}")
    log(f"  {'-'*50}")

    for name, gammas in profiles.items():
        L = build_L(H, gammas)
        rho_t = evolve(L, rho0, t_measure)
        mi = mutual_info(rho_t, N, [0, 1], [3, 4])
        Sg = sum(gammas)
        log(f"  {name:>14}  {mi:>10.6f}  {mi:>12.6f}  {Sg:>8.3f}")

    log()


# ============================================================
# EXPERIMENT 1: DD Simulation
# ============================================================
def run_exp_1():
    log("=" * 70)
    log("EXPERIMENT 1: DYNAMICAL DECOUPLING SIMULATION")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_A_0M_pp_B()
    t_measure = 5.0
    gamma_base = 0.05

    strategies = {
        'no DD (baseline)':   [1.0, 1.0, 1.0, 1.0, 1.0],
        'DD on M only':       [1.0, 1.0, 0.1, 1.0, 1.0],
        'DD on receiver':     [1.0, 1.0, 1.0, 0.1, 0.1],
        'DD on M+receiver':   [1.0, 1.0, 0.1, 0.1, 0.1],
        'DD on all but sender':[1.0, 1.0, 0.1, 0.1, 0.1],
        'DD everywhere':      [0.1, 0.1, 0.1, 0.1, 0.1],
    }

    log(f"  {'Strategy':>22}  {'MI(A:B)':>10}  {'gamma_eff profile':>30}")
    log(f"  {'-'*65}")

    for name, factors in strategies.items():
        gammas = [gamma_base * f for f in factors]
        L = build_L(H, gammas)
        rho_t = evolve(L, rho0, t_measure)
        mi = mutual_info(rho_t, N, [0, 1], [3, 4])
        profile = '[' + ', '.join(f'{g:.3f}' for g in gammas) + ']'
        log(f"  {name:>22}  {mi:>10.6f}  {profile:>30}")

    log()


# ============================================================
# EXPERIMENT 3: γ Modulation (AC on gate)
# ============================================================
def run_exp_3():
    log("=" * 70)
    log("EXPERIMENT 3: gamma MODULATION (AC ON GATE)")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_A_0M_pp_B()
    gamma_base = 0.05
    gamma_amp = 0.04  # modulation amplitude
    dt = 0.05
    t_max = 10.0
    n_steps = int(t_max / dt)

    log(f"  {'freq':>8}  {'MI(A:B) at t=5':>15}  {'MI(A:B) at t=10':>16}")
    log(f"  {'-'*42}")

    for freq in [0, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0]:
        rho = rho0.copy()
        mi_5 = 0; mi_10 = 0

        for step in range(n_steps):
            t = step * dt
            gM = gamma_base + gamma_amp * np.sin(2 * np.pi * freq * t) if freq > 0 else gamma_base
            gM = max(0.005, gM)  # floor
            gammas = [gamma_base, gamma_base, gM, gamma_base, gamma_base]
            L = build_L(H, gammas)
            rho = evolve(L, rho, dt)

            if abs(t + dt - 5.0) < dt / 2:
                mi_5 = mutual_info(rho, N, [0, 1], [3, 4])
            if abs(t + dt - 10.0) < dt / 2:
                mi_10 = mutual_info(rho, N, [0, 1], [3, 4])

        log(f"  {freq:>8.1f}  {mi_5:>15.6f}  {mi_10:>16.6f}")

    log()


# ============================================================
# EXPERIMENT 4: Time-Resolved Decoder
# ============================================================
def run_exp_4():
    log("=" * 70)
    log("EXPERIMENT 4: TIME-RESOLVED DECODER")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_A_0M_pp_B()
    dt = 0.1
    t_switch = 5.0
    t_max = 10.0

    gammas_before = [0.05] * N
    gammas_after = [0.05, 0.05, 0.10, 0.05, 0.05]  # qubit 2 doubles

    log(f"  At t={t_switch}: gamma on qubit 2 switches from 0.05 to 0.10")
    log()
    log(f"  {'t':>6}  {'MI(A:B)':>10}  {'Purity':>10}  {'Phase':>8}")
    log(f"  {'-'*38}")

    rho = rho0.copy()
    for step in range(int(t_max / dt)):
        t = step * dt
        gammas = gammas_before if t < t_switch else gammas_after
        L = build_L(H, gammas)
        rho = evolve(L, rho, dt)

        if step % 5 == 0:
            mi = mutual_info(rho, N, [0, 1], [3, 4])
            pur = np.real(np.trace(rho @ rho))
            phase = "before" if t < t_switch else "AFTER"
            log(f"  {t+dt:>6.1f}  {mi:>10.6f}  {pur:>10.6f}  {phase:>8}")

    log()


# ============================================================
# EXPERIMENT 5: γ Feedback Loop
# ============================================================
def run_exp_5():
    log("=" * 70)
    log("EXPERIMENT 5: gamma FEEDBACK LOOP")
    log("=" * 70)
    log()

    H = build_H()
    rho0 = bell_A_0M_pp_B()
    gamma_base = 0.05
    dt = 0.05
    t_max = 10.0
    n_steps = int(t_max / dt)

    # O_int = Z_1 Z_2 (interaction observable on the boundary)
    O_int = site_op(sz, 1) @ site_op(sz, 2)

    log(f"  {'kappa':>8}  {'MI(t=5)':>10}  {'MI(t=10)':>10}  {'gamma_M range':>18}")
    log(f"  {'-'*50}")

    for kappa in [0, 0.1, 0.3, 0.5, 0.7, 1.0]:
        rho = rho0.copy()
        mi_5 = 0; mi_10 = 0
        gM_min = 1e10; gM_max = 0

        for step in range(n_steps):
            t = step * dt
            O_val = np.real(np.trace(O_int @ rho))
            gM = gamma_base * (1 + kappa * abs(O_val))
            gM = max(0.001, min(0.5, gM))  # clamp
            gM_min = min(gM_min, gM)
            gM_max = max(gM_max, gM)

            gammas = [gamma_base, gamma_base, gM, gamma_base, gamma_base]
            L = build_L(H, gammas)
            rho = evolve(L, rho, dt)

            if abs(t + dt - 5.0) < dt / 2:
                mi_5 = mutual_info(rho, N, [0, 1], [3, 4])
            if abs(t + dt - 10.0) < dt / 2:
                mi_10 = mutual_info(rho, N, [0, 1], [3, 4])

        gM_range = f"[{gM_min:.4f}, {gM_max:.4f}]"
        log(f"  {kappa:>8.1f}  {mi_5:>10.6f}  {mi_10:>10.6f}  {gM_range:>18}")

    log()


# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Using the Mediator: Practical gamma Control")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, d={d}")
    log()

    t1 = time.time()
    run_exp_2()
    log(f"[Exp 2 completed in {time.time()-t1:.1f}s]")
    log()

    t2 = time.time()
    run_exp_1()
    log(f"[Exp 1 completed in {time.time()-t2:.1f}s]")
    log()

    t3 = time.time()
    run_exp_3()
    log(f"[Exp 3 completed in {time.time()-t3:.1f}s]")
    log()

    t4 = time.time()
    run_exp_4()
    log(f"[Exp 4 completed in {time.time()-t4:.1f}s]")
    log()

    t5 = time.time()
    run_exp_5()
    log(f"[Exp 5 completed in {time.time()-t5:.1f}s]")
    log()

    log(f"Total runtime: {time.time()-t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
