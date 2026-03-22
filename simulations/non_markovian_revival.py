#!/usr/bin/env python3
"""
Non-Markovian CΨ Revival Test
===============================
Can non-Markovian dynamics push CΨ permanently back above 1/4?

Test 1: Markovian baseline (constant γ, J>0) — do oscillations cross 1/4?
Test 2: Structured bath — 2 system qubits coupled to 1 bath qubit
Test 3: Pulsed γ — dephasing on/off periodically
Test 4: Oscillating γ — γ(t) = γ_0 · (1 + A·sin(ωt))
Test 5: Worst-case sweep — find parameters that maximize revival

Script:  simulations/non_markovian_revival.py
Output:  simulations/results/non_markovian_revival.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

# ====================================================================
# Output
# ====================================================================

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "non_markovian_revival.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ====================================================================
# Pauli matrices
# ====================================================================

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)


# ====================================================================
# Infrastructure
# ====================================================================

def site_op(op, k, nq):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def heisenberg_bond(i, j, nq, J=1.0):
    d = 2 ** nq
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, i, nq) @ site_op(P, j, nq)
    return H


def build_liouvillian(H, collapse_ops):
    """Build full Liouvillian superoperator."""
    d = H.shape[0]
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for C_op in collapse_ops:
        Cd = C_op.conj().T
        CdC = Cd @ C_op
        L += np.kron(C_op, C_op.conj()) - 0.5 * (
            np.kron(CdC, Id) + np.kron(Id, CdC.T))
    return L


def evolve_expm(L, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def lindblad_rhs(rho, H, L_ops, LdL_ops):
    drho = -1j * (H @ rho - rho @ H)
    for L, LdL in zip(L_ops, LdL_ops):
        drho += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return drho


def rk4_step(rho, H, L_ops, LdL_ops, dt):
    k1 = lindblad_rhs(rho, H, L_ops, LdL_ops)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, L_ops, LdL_ops)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, L_ops, LdL_ops)
    k4 = lindblad_rhs(rho + dt * k3, H, L_ops, LdL_ops)
    rho_new = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    rho_new = 0.5 * (rho_new + rho_new.conj().T)
    tr = np.trace(rho_new).real
    if abs(tr) < 1e-14:
        raise FloatingPointError("Trace collapsed.")
    rho_new /= tr
    return rho_new


def partial_trace_keep(rho, keep, n_qubits):
    keep = list(keep)
    trace_out = [q for q in range(n_qubits) if q not in keep]
    dims = [2] * n_qubits
    reshaped = rho.reshape(dims + dims)
    current_n = n_qubits
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_keep = 2 ** len(keep)
    return reshaped.reshape((d_keep, d_keep))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1) if d > 1 else 0.0


def cpsi(rho):
    """CΨ = purity × psi_norm (the product that crosses 1/4)."""
    return purity(rho) * psi_norm(rho)


def analyze_crossings(t_arr, cpsi_arr, threshold=0.25):
    """Analyze CΨ trajectory relative to 1/4 boundary."""
    cross_down = 0
    cross_up = 0
    first_below_t = None
    max_revival = 0.0
    max_revival_t = 0.0
    sustained_above = 0.0  # longest continuous time above 1/4 after first crossing below

    above = cpsi_arr[0] > threshold
    below_since = None
    above_since = None

    for i in range(1, len(t_arr)):
        now_above = cpsi_arr[i] > threshold
        if above and not now_above:
            cross_down += 1
            if first_below_t is None:
                first_below_t = t_arr[i]
            below_since = t_arr[i]
            if above_since is not None and first_below_t is not None:
                dur = t_arr[i] - above_since
                if dur > sustained_above and above_since > first_below_t:
                    sustained_above = dur
            above_since = None
        if not above and now_above:
            cross_up += 1
            above_since = t_arr[i]
            if first_below_t is not None:
                revival = cpsi_arr[i]
                if revival > max_revival:
                    max_revival = revival
                    max_revival_t = t_arr[i]
        above = now_above

    # Final sustained period
    if above_since is not None and first_below_t is not None:
        dur = t_arr[-1] - above_since
        if dur > sustained_above and above_since > first_below_t:
            sustained_above = dur

    return {
        'cross_down': cross_down,
        'cross_up': cross_up,
        'first_below_t': first_below_t,
        'max_revival': max_revival,
        'max_revival_t': max_revival_t,
        'sustained_above': sustained_above,
        'final_cpsi': cpsi_arr[-1],
    }


def log_analysis(info, label=""):
    if label:
        log(f"  {label}:")
    fb = info['first_below_t']
    log(f"    Crossings: {info['cross_down']} down, {info['cross_up']} up")
    log(f"    First below 1/4: t={fb:.2f}" if fb else "    Never crossed below 1/4")
    if info['cross_up'] > 0 and fb:
        log(f"    Max revival after first drop: CΨ={info['max_revival']:.4f} at t={info['max_revival_t']:.2f}")
        log(f"    Longest sustained above 1/4 (after first drop): {info['sustained_above']:.2f}")
    log(f"    Final CΨ: {info['final_cpsi']:.6f}")


# ====================================================================
# Test 1: Markovian baseline
# ====================================================================

def test_1_markovian():
    log("=" * 70)
    log("TEST 1: MARKOVIAN BASELINE (constant γ, J>0)")
    log("=" * 70)
    log()

    nq = 2; d = 4; d2 = 16
    J = 1.0

    # Bell+ initial state
    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    # |01> initial state
    psi01 = np.kron(up, dn)
    rho0_01 = np.outer(psi01, psi01.conj())

    H = heisenberg_bond(0, 1, nq, J)
    t_max = 60.0
    dt = 0.02
    n_steps = int(round(t_max / dt))
    t_arr = np.array([i * dt for i in range(n_steps + 1)])

    for gamma, state_name, rho_init in [
        (0.05, "Bell+", rho0),
        (0.05, "|01>", rho0_01),
        (0.01, "Bell+", rho0),
        (0.10, "Bell+", rho0),
    ]:
        c_ops = [np.sqrt(gamma) * site_op(sz, k, nq) for k in range(nq)]
        L = build_liouvillian(H, c_ops)

        cpsi_arr = np.zeros(n_steps + 1)
        for i, t in enumerate(t_arr):
            rho = evolve_expm(L, rho_init, t)
            cpsi_arr[i] = cpsi(rho)

        info = analyze_crossings(t_arr, cpsi_arr)
        log(f"  γ={gamma}, J={J}, state={state_name}, t_max={t_max}")
        log_analysis(info)
        log()


# ====================================================================
# Test 2: Structured bath (2 system + 1 bath qubit)
# ====================================================================

def test_2_structured_bath():
    log("=" * 70)
    log("TEST 2: STRUCTURED BATH (S1-S2 entangled, B coupled to S2)")
    log("=" * 70)
    log()

    nq = 3; d = 8
    t_max = 100.0; dt = 0.02
    n_steps = int(round(t_max / dt))
    t_arr = np.array([i * dt for i in range(n_steps + 1)])

    # S1-S2 in Bell+, B in |0>
    bell_12 = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    psi_init = np.kron(bell_12, up)  # S1 S2 B
    rho0 = np.outer(psi_init, psi_init.conj())

    J_sys = 1.0  # S1-S2 coupling
    H_sys = heisenberg_bond(0, 1, nq, J_sys)

    log(f"  System: S1-S2 Bell+, Bath: B=|0>")
    log(f"  J_sys={J_sys}, t_max={t_max}")
    log()

    # Sweep J_SB (system-bath coupling) and γ_B (bath dephasing)
    best_revival = 0.0
    best_params = {}

    configs = [
        # J_SB, γ_B, γ_sys, label
        (0.5,  0.1, 0.0,  "J_SB=0.5, γ_B=0.1, γ_S=0"),
        (1.0,  0.1, 0.0,  "J_SB=1.0, γ_B=0.1, γ_S=0"),
        (2.0,  0.1, 0.0,  "J_SB=2.0, γ_B=0.1, γ_S=0"),
        (1.0,  0.05, 0.0, "J_SB=1.0, γ_B=0.05, γ_S=0"),
        (1.0,  0.5, 0.0,  "J_SB=1.0, γ_B=0.5, γ_S=0"),
        (1.0,  0.01, 0.0, "J_SB=1.0, γ_B=0.01, γ_S=0"),
        (0.5,  0.1, 0.01, "J_SB=0.5, γ_B=0.1, γ_S=0.01"),
        (1.0,  0.1, 0.01, "J_SB=1.0, γ_B=0.1, γ_S=0.01"),
    ]

    for J_SB, gamma_B, gamma_S, label in configs:
        H = H_sys + heisenberg_bond(1, 2, nq, J_SB)

        c_ops = []
        if gamma_B > 0:
            c_ops.append(np.sqrt(gamma_B) * site_op(sz, 2, nq))
        if gamma_S > 0:
            for k in range(2):
                c_ops.append(np.sqrt(gamma_S) * site_op(sz, k, nq))

        L = build_liouvillian(H, c_ops)

        cpsi_arr = np.zeros(n_steps + 1)
        for i, t in enumerate(t_arr):
            rho_full = evolve_expm(L, rho0, t)
            rho_sys = partial_trace_keep(rho_full, [0, 1], nq)
            cpsi_arr[i] = cpsi(rho_sys)

        info = analyze_crossings(t_arr, cpsi_arr)
        log(f"  {label}")
        log_analysis(info)

        if info['max_revival'] > best_revival:
            best_revival = info['max_revival']
            best_params = {'label': label, 'info': info}
        log()

    log(f"  BEST REVIVAL: CΨ={best_revival:.4f} ({label})")
    if best_revival > 0.25:
        log(f"  *** CΨ CROSSED BACK ABOVE 1/4! ***")
    else:
        log(f"  CΨ never returned above 1/4")
    log()


# ====================================================================
# Test 3: Pulsed γ (on/off)
# ====================================================================

def test_3_pulsed():
    log("=" * 70)
    log("TEST 3: PULSED γ (dephasing on/off periodically)")
    log("=" * 70)
    log()

    nq = 2; d = 4
    J = 1.0
    t_max = 100.0; dt = 0.02
    n_steps = int(round(t_max / dt))

    # Bell+ initial state
    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    H = heisenberg_bond(0, 1, nq, J)

    # Pre-build Liouvillians for γ=0 and γ>0
    gamma_on = 0.05
    c_ops_on = [np.sqrt(gamma_on) * site_op(sz, k, nq) for k in range(nq)]
    LdL_on = [c.conj().T @ c for c in c_ops_on]
    c_ops_off = []
    LdL_off = []

    configs = [
        # T_on, T_off, label
        (2.0,  2.0,  "T_on=2, T_off=2"),
        (1.0,  3.0,  "T_on=1, T_off=3"),
        (3.0,  1.0,  "T_on=3, T_off=1"),
        (0.5,  4.5,  "T_on=0.5, T_off=4.5"),
        (5.0,  5.0,  "T_on=5, T_off=5"),
        (10.0, 10.0, "T_on=10, T_off=10"),
    ]

    for T_on, T_off, label in configs:
        period = T_on + T_off
        rho = rho0.copy()
        t_arr = []
        cpsi_arr = []

        for step in range(n_steps + 1):
            t = step * dt
            t_arr.append(t)
            cpsi_arr.append(cpsi(rho))

            if step < n_steps:
                phase = t % period
                if phase < T_on:
                    rho = rk4_step(rho, H, c_ops_on, LdL_on, dt)
                else:
                    rho = rk4_step(rho, H, c_ops_off, LdL_off, dt)

        t_arr = np.array(t_arr)
        cpsi_arr = np.array(cpsi_arr)

        info = analyze_crossings(t_arr, cpsi_arr)
        log(f"  γ={gamma_on} (on/off), J={J}, {label}")
        log_analysis(info)
        log()


# ====================================================================
# Test 4: Oscillating γ(t)
# ====================================================================

def test_4_oscillating():
    log("=" * 70)
    log("TEST 4: OSCILLATING γ(t) = γ_0 * max(0, 1 + A·sin(ωt))")
    log("=" * 70)
    log()

    nq = 2; d = 4
    J = 1.0
    gamma_0 = 0.05
    t_max = 100.0; dt = 0.02
    n_steps = int(round(t_max / dt))

    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    H = heisenberg_bond(0, 1, nq, J)

    configs = [
        # A (amplitude), omega, label
        (0.9,  1.0,  "A=0.9, ω=1.0"),
        (0.9,  0.5,  "A=0.9, ω=0.5"),
        (0.9,  2.0,  "A=0.9, ω=2.0"),
        (1.0,  1.0,  "A=1.0, ω=1.0 (γ touches 0)"),
        (1.0,  0.3,  "A=1.0, ω=0.3 (slow, γ touches 0)"),
        (1.5,  1.0,  "A=1.5, ω=1.0 (γ clips at 0)"),
    ]

    for A, omega, label in configs:
        rho = rho0.copy()
        t_arr = []
        cpsi_arr = []

        for step in range(n_steps + 1):
            t = step * dt
            t_arr.append(t)
            cpsi_arr.append(cpsi(rho))

            if step < n_steps:
                gamma_t = gamma_0 * max(0.0, 1.0 + A * np.sin(omega * t))
                if gamma_t > 0:
                    c_ops = [np.sqrt(gamma_t) * site_op(sz, k, nq) for k in range(nq)]
                    LdL = [c.conj().T @ c for c in c_ops]
                else:
                    c_ops = []
                    LdL = []
                rho = rk4_step(rho, H, c_ops, LdL, dt)

        t_arr = np.array(t_arr)
        cpsi_arr = np.array(cpsi_arr)

        info = analyze_crossings(t_arr, cpsi_arr)
        log(f"  γ_0={gamma_0}, {label}")
        log_analysis(info)
        log()


# ====================================================================
# Test 5: Worst-case sweep (structured bath)
# ====================================================================

def test_5_worst_case():
    log("=" * 70)
    log("TEST 5: WORST-CASE SWEEP (maximize revival above 1/4)")
    log("=" * 70)
    log()

    nq = 3; d = 8
    t_max = 80.0; dt = 0.02
    n_steps = int(round(t_max / dt))

    # S1-S2 in Bell+, B in |+>  (bath starts with coherence)
    bell_12 = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    plus_state = (up + dn) / np.sqrt(2)
    psi_init_plus = np.kron(bell_12, plus_state)
    rho0_plus = np.outer(psi_init_plus, psi_init_plus.conj())

    # Also try B in |0>
    psi_init_0 = np.kron(bell_12, up)
    rho0_0 = np.outer(psi_init_0, psi_init_0.conj())

    J_sys = 1.0
    H_sys = heisenberg_bond(0, 1, nq, J_sys)

    best_overall = 0.0
    best_label = ""

    # Wide sweep
    J_SBs = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]
    gamma_Bs = [0.01, 0.05, 0.1, 0.5]

    log(f"  Sweeping J_SB × γ_B × bath_state ({len(J_SBs)}×{len(gamma_Bs)}×2 = {len(J_SBs)*len(gamma_Bs)*2} configs)")
    log()

    header = f"  {'J_SB':>6}  {'γ_B':>6}  {'Bath':>5}  {'↓':>3}  {'↑':>3}  {'MaxRev':>8}  {'Sust':>6}  {'Final':>8}"
    log(header)
    log("  " + "-" * 62)

    for J_SB in J_SBs:
        for gamma_B in gamma_Bs:
            for bath_label, rho0_bath in [("  |0>", rho0_0), ("  |+>", rho0_plus)]:
                H = H_sys + heisenberg_bond(1, 2, nq, J_SB)
                c_ops = [np.sqrt(gamma_B) * site_op(sz, 2, nq)]
                L = build_liouvillian(H, c_ops)

                cpsi_arr = np.zeros(n_steps + 1)
                for i in range(n_steps + 1):
                    t = i * dt
                    rho_full = evolve_expm(L, rho0_bath, t)
                    rho_sys = partial_trace_keep(rho_full, [0, 1], nq)
                    cpsi_arr[i] = cpsi(rho_sys)

                t_arr = np.array([i * dt for i in range(n_steps + 1)])
                info = analyze_crossings(t_arr, cpsi_arr)

                rev = info['max_revival']
                tag = " ***" if rev > 0.25 else ""
                log(f"  {J_SB:6.1f}  {gamma_B:6.2f}  {bath_label}  "
                    f"{info['cross_down']:3d}  {info['cross_up']:3d}  "
                    f"{rev:8.4f}  {info['sustained_above']:6.1f}  "
                    f"{info['final_cpsi']:8.4f}{tag}")

                if rev > best_overall:
                    best_overall = rev
                    best_label = f"J_SB={J_SB}, γ_B={gamma_B}, bath={bath_label.strip()}"

    log()
    log(f"  BEST REVIVAL: CΨ = {best_overall:.4f}")
    log(f"  Parameters: {best_label}")
    if best_overall > 0.25:
        log(f"  *** NON-MARKOVIAN REVIVAL ABOVE 1/4 FOUND! ***")
        log(f"  The 1/4 boundary is NOT absorbing under non-Markovian dynamics.")
    else:
        log(f"  1/4 boundary holds: no revival above 1/4 found.")
    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Non-Markovian CΨ Revival Test")
    log("=" * 70)
    log(f"Question: Can non-Markovian dynamics push CΨ back above 1/4?")
    log()

    test_1_markovian()
    test_2_structured_bath()
    test_3_pulsed()
    test_4_oscillating()
    test_5_worst_case()

    log("=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log("If ANY test shows sustained CΨ > 1/4 after first crossing below:")
    log("  → 1/4 is NOT an absorbing boundary under non-Markovian dynamics")
    log("  → Layer 5 remains OPEN")
    log()
    log("If NO test shows revival above 1/4:")
    log("  → Strong evidence that 1/4 is absorbing even non-Markovianly")
    log("  → Layer 5 can be considered closed (pending formal proof)")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
