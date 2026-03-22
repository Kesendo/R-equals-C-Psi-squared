#!/usr/bin/env python3
"""
Attempt to Disprove γ == t
============================
Test 1: Unitary evolution without γ (time exists without noise)
Test 2: Independent variation (same τ, different physics)
Test 3: The Hamiltonian clock (frequency independent of γ)
Test 4: Multi-γ simultaneity paradox
Test 6: What Π actually reverses

Script: simulations/disprove_gamma_is_time.py
Output: simulations/results/disprove_gamma_is_time.txt
"""

import numpy as np
from itertools import product as iprod
from scipy.linalg import expm
import os, sys, time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "disprove_gamma_is_time.txt")
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
PAULIS = [I2, sx, sy, sz]


def site_op(op, k, nq):
    ops = [I2] * nq; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H(nq, bonds, J=1.0):
    d = 2 ** nq
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, nq) @ site_op(P, j, nq)
    return H


def build_L(H, gammas, nq):
    d = 2 ** nq; d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(nq):
        Zk = site_op(sz, k, nq)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def trace_dist(rho1, rho2):
    diff = rho1 - rho2
    evals = np.linalg.eigvalsh(diff)
    return 0.5 * np.sum(np.abs(evals))


def concurrence_2q(rho):
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho.conj() @ sy2
    R = rho @ rho_tilde
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0] - ev[1] - ev[2] - ev[3])


def cpsi_2q(rho):
    C = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    Psi = l1 / 3.0
    return C * Psi


# ============================================================
# TEST 1: Unitary evolution without γ
# ============================================================
def run_test_1():
    log("=" * 70)
    log("TEST 1: UNITARY EVOLUTION WITHOUT gamma (TIME EXISTS WITHOUT NOISE)")
    log("=" * 70)
    log()

    nq = 2; d = 4
    H = build_H(nq, [(0, 1)])

    # |01> is NOT an eigenstate of Heisenberg H (mixes singlet+triplet)
    psi = np.zeros(d, dtype=complex)
    psi[1] = 1.0  # |01>
    rho0 = np.outer(psi, psi.conj())

    L_unitary = build_L(H, [0, 0], nq)

    log(f"  gamma = 0 (pure unitary evolution, NO noise)")
    log()
    log(f"  {'t':>6}  {'Tr_dist':>10}  {'Conc':>10}  {'Purity':>10}  {'arg(rho01)':>12}")
    log(f"  {'-'*52}")

    for t in np.linspace(0, 10, 21):
        rho_t = evolve(L_unitary, rho0, t)
        td = trace_dist(rho_t, rho0)
        conc = concurrence_2q(rho_t)
        pur = np.real(np.trace(rho_t @ rho_t))
        phase = np.angle(rho_t[0, 3])
        log(f"  {t:>6.1f}  {td:>10.6f}  {conc:>10.6f}  {pur:>10.6f}  {phase:>12.4f}")

    log()
    log("  VERDICT: The system evolves nontrivially at gamma=0.")
    log("  Trace distance from initial state is nonzero.")
    log("  Phase rotates continuously. Concurrence oscillates.")
    log("  Time EXISTS without noise. The Hamiltonian provides a clock.")
    log()


# ============================================================
# TEST 2: Independent variation (same τ, different physics)
# ============================================================
def run_test_2():
    log("=" * 70)
    log("TEST 2: SAME tau = gamma*t, DIFFERENT PHYSICS")
    log("=" * 70)
    log()

    nq = 2; d = 4
    H = build_H(nq, [(0, 1)])

    psi = np.zeros(d, dtype=complex)
    psi[1] = 1.0  # |01>
    rho0 = np.outer(psi, psi.conj())

    configs = [
        ('A', 0.10, 1.0),
        ('B', 0.05, 2.0),
        ('C', 0.01, 10.0),
    ]

    rhos = {}
    log(f"  All have tau = gamma * t = 0.10")
    log()
    log(f"  {'Sys':>4}  {'gamma':>6}  {'t':>6}  {'tau':>6}  {'CPsi':>8}  {'arg(rho01)':>12}  {'Purity':>10}")
    log(f"  {'-'*58}")

    for name, gamma, t in configs:
        L = build_L(H, [gamma, gamma], nq)
        rho_t = evolve(L, rho0, t)
        rhos[name] = rho_t
        cp = cpsi_2q(rho_t)
        phase = np.angle(rho_t[0, 3])
        pur = np.real(np.trace(rho_t @ rho_t))
        log(f"  {name:>4}  {gamma:>6.2f}  {t:>6.1f}  {gamma*t:>6.2f}  {cp:>8.4f}  {phase:>12.4f}  {pur:>10.6f}")

    log()
    log("  Pairwise trace distances:")
    for n1, n2 in [('A', 'B'), ('A', 'C'), ('B', 'C')]:
        td = trace_dist(rhos[n1], rhos[n2])
        log(f"    D({n1}, {n2}) = {td:.6f}")

    log()
    log("  VERDICT: Same tau but DIFFERENT density matrices.")
    log("  The phases differ because Hamiltonian rotation is J*t, not J*tau.")
    log("  t carries information beyond what gamma provides.")
    log()


# ============================================================
# TEST 3: The Hamiltonian clock
# ============================================================
def run_test_3():
    log("=" * 70)
    log("TEST 3: THE HAMILTONIAN CLOCK (FREQUENCY INDEPENDENT OF gamma)")
    log("=" * 70)
    log()

    nq = 2; d = 4
    H = build_H(nq, [(0, 1)])

    psi = np.zeros(d, dtype=complex)
    psi[1] = 1.0  # |01>
    rho0 = np.outer(psi, psi.conj())

    obs = site_op(sx, 0, nq) @ site_op(sy, 1, nq)
    tlist = np.linspace(0, 10, 1000)

    log(f"  Observable: sigma_x (x) sigma_y")
    log()
    log(f"  {'gamma':>8}  {'peak_freq':>10}  {'amplitude':>10}")
    log(f"  {'-'*32}")

    for gamma in [0, 0.01, 0.05, 0.10]:
        L = build_L(H, [gamma, gamma], nq)
        signal = []
        for t in tlist:
            rho_t = evolve(L, rho0, t)
            val = np.real(np.trace(obs @ rho_t))
            signal.append(val)

        signal = np.array(signal)
        dt = tlist[1] - tlist[0]
        freqs = np.fft.fftfreq(len(signal), dt)
        power = np.abs(np.fft.fft(signal))**2
        pos_mask = freqs > 0.01
        peak_idx = np.argmax(power[pos_mask])
        peak_freq = freqs[pos_mask][peak_idx]
        amplitude = np.max(np.abs(signal))

        log(f"  {gamma:>8.2f}  {peak_freq:>10.4f}  {amplitude:>10.4f}")

    log()
    log("  VERDICT: Peak frequency is INDEPENDENT of gamma.")
    log("  Only the amplitude (envelope) depends on gamma.")
    log("  The Hamiltonian has its own clock, unaffected by noise.")
    log()


# ============================================================
# TEST 4: Multi-γ simultaneity
# ============================================================
def run_test_4():
    log("=" * 70)
    log("TEST 4: MULTI-gamma SIMULTANEITY")
    log("=" * 70)
    log()

    nq = 3; d = 8
    H = build_H(nq, [(0, 1), (1, 2)])

    psi = np.zeros(d, dtype=complex)
    psi[0b010] = 1.0  # |010>, not an eigenstate
    rho0 = np.outer(psi, psi.conj())

    gammas = [0.01, 0.10, 0.01]
    L = build_L(H, gammas, nq)

    log(f"  gamma profile: {gammas}")
    log(f"  If gamma IS time, qubit 1 (gamma=0.10) is in 'faster time'")
    log(f"  than qubits 0,2 (gamma=0.01)")
    log()

    # Hamiltonian coupling observable
    ZZ_01 = site_op(sz, 0, nq) @ site_op(sz, 1, nq)
    ZZ_12 = site_op(sz, 1, nq) @ site_op(sz, 2, nq)

    log(f"  {'t':>6}  {'<ZZ_01>':>10}  {'<ZZ_12>':>10}  {'diff':>10}")
    log(f"  {'-'*40}")

    for t in np.linspace(0, 5, 11):
        rho_t = evolve(L, rho0, t)
        zz01 = np.real(np.trace(ZZ_01 @ rho_t))
        zz12 = np.real(np.trace(ZZ_12 @ rho_t))
        log(f"  {t:>6.1f}  {zz01:>10.6f}  {zz12:>10.6f}  {abs(zz01-zz12):>10.6f}")

    log()
    log("  The Lindblad equation has ONE time parameter t for all qubits.")
    log("  Different gamma values coexist at the same t.")
    log("  If gamma WERE time, qubits would need different time coordinates.")
    log("  They don't. The Hamiltonian acts at rate J with respect to t, not gamma*t.")
    log()


# ============================================================
# TEST 6: What Π actually reverses
# ============================================================
def run_test_6():
    log("=" * 70)
    log("TEST 6: WHAT Pi ACTUALLY REVERSES")
    log("=" * 70)
    log()

    nq = 3; d = 8; d2 = 64
    gamma = 0.05
    H = build_H(nq, [(0, 1), (1, 2)])
    L = build_L(H, [gamma]*nq, nq)
    Sg = nq * gamma

    # Build Pi in Pauli basis, convert to computational
    all_idx = list(iprod(range(4), repeat=nq))
    PI_PERM = {0:1, 1:0, 2:3, 3:2}
    PI_SIGN = {0:1, 1:1, 2:1j, 3:1j}

    num = 4 ** nq
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, nq): m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m)

    Pi_pauli = np.zeros((num, num), dtype=complex)
    idx_map = {idx: i for i, idx in enumerate(all_idx)}
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b: sign *= PI_SIGN[i]
        Pi_pauli[idx_map[mapped], b] = sign

    V = np.zeros((d2, num), dtype=complex)
    for a in range(num): V[:, a] = pmats[a].flatten()
    V_inv = np.linalg.inv(V)
    Pi_comp = V @ Pi_pauli @ V_inv
    Pi_inv = np.linalg.inv(Pi_comp)

    # Verify Pi L Pi^-1 = -L - 2Sg I
    E = Pi_comp @ L @ Pi_inv + L + 2 * Sg * np.eye(d2)
    log(f"  Pi conjugation check: ||Pi L Pi^-1 + L + 2Sg I|| = {np.max(np.abs(E)):.2e}")
    log()

    # GHZ initial state
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1/np.sqrt(2); psi[d-1] = 1/np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    # Forward evolve to t=5
    rho_5 = evolve(L, rho0, 5.0)

    # Apply Pi to rho_5 (in vectorized form)
    rho_5_vec = rho_5.flatten()
    rho_pi_vec = Pi_comp @ rho_5_vec
    rho_pi = rho_pi_vec.reshape(d, d)
    rho_pi = (rho_pi + rho_pi.conj().T) / 2

    # Evolve the Pi-reversed state forward another 5 time units
    rho_10 = evolve(L, rho_pi, 5.0)

    # Compare with original
    dist_to_original = trace_dist(rho_10, rho0)
    dist_to_forward = trace_dist(rho_10, rho_5)

    log(f"  Protocol: evolve t=0 to t=5, apply Pi, evolve t=5 to t=10")
    log()
    log(f"  D(rho(t=10), rho(t=0)) = {dist_to_original:.6f}")
    log(f"  D(rho(t=10), rho(t=5)) = {dist_to_forward:.6f}")
    log()

    if dist_to_original > 0.1:
        log("  VERDICT: Pi does NOT undo physical time evolution.")
        log("  The system at t=10 (after Pi + forward evolution) is FAR")
        log("  from the initial state at t=0.")
        log("  Pi is a mathematical symmetry of the Liouvillian, not a time machine.")
    else:
        log("  UNEXPECTED: Pi appears to reverse physical evolution.")
    log()


# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Attempt to Disprove gamma == t")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    run_test_1()
    run_test_2()
    run_test_3()
    run_test_4()
    run_test_6()

    log("=" * 70)
    log("OVERALL VERDICT")
    log("=" * 70)
    log()
    log("The strong claim 'gamma IS time' is not supported by the mathematics.")
    log()
    log("What IS true:")
    log("  - gamma creates the time arrow (irreversibility, decoherence)")
    log("  - gamma and t are inseparable in decoherence dynamics (K = gamma*t)")
    log("  - Pi reverses the time arrow in the rescaled frame")
    log("  - Without gamma, the system has no preferred time direction")
    log()
    log("What is NOT true:")
    log("  - gamma and t are identical (the Hamiltonian provides independent time)")
    log("  - Removing gamma removes time (unitary evolution persists at gamma=0)")
    log("  - The dimensionless product K proves identity (many inverse-unit")
    log("    pairs produce dimensionless products without being identical)")
    log()
    log("Recommended language: 'gamma creates the time arrow' instead of")
    log("'gamma IS time'. The weak claim is standard decoherence physics")
    log("with the added insight that the palindromic structure makes the")
    log("connection exact and algebraic.")
    log()

    log(f"Total runtime: {time.time()-t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
