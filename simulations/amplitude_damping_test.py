#!/usr/bin/env python3
"""
Amplitude Damping: Direct CΨ Trajectory
==========================================
The non-unital channel. L = √γ |0⟩⟨1| (spontaneous decay to |0⟩).
Fixed point: |00⟩ (pure state with CΨ = 0).

Test 1: Pure amplitude damping (no H) - Bell+ trajectory
Test 2: Amplitude damping + Heisenberg coupling
Test 3: K-invariance - does K = γ·t_cross hold?
Test 4: Combined amplitude damping + Z-dephasing
Test 5: Asymmetric rates and various initial states
Test 6: Monotonicity check - is CΨ non-increasing above 1/4?

Script:  simulations/amplitude_damping_test.py
Output:  simulations/results/amplitude_damping_test.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

# ====================================================================
# Output
# ====================================================================

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "amplitude_damping_test.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ====================================================================
# Pauli matrices and states
# ====================================================================

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
# Lowering operator |0><1|
sm = np.array([[0, 1], [0, 0]], dtype=complex)  # σ_-
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

N = 2; d = 4; d2 = 16


# ====================================================================
# Infrastructure
# ====================================================================

def site_op(op, k):
    return np.kron(op, I2) if k == 0 else np.kron(I2, op)


def build_H(J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J * site_op(P, 0) @ site_op(P, 1)
    return H


def build_L(H, collapse_ops):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for C_op in collapse_ops:
        Cd = C_op.conj().T
        CdC = Cd @ C_op
        L += np.kron(C_op, C_op.conj()) - 0.5 * (
            np.kron(CdC, Id) + np.kron(Id, CdC.T))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def concurrence(rho):
    try:
        sy2 = np.kron(sy, sy)
        rho_tilde = sy2 @ rho.conj() @ sy2
        R = rho @ rho_tilde
        evals = np.linalg.eigvals(R)
        evals = np.real(np.sqrt(np.maximum(evals, 0.0)))
        evals = np.sort(evals)[::-1]
        return float(max(0.0, evals[0] - evals[1] - evals[2] - evals[3]))
    except np.linalg.LinAlgError:
        return 0.0


def make_state(name):
    if name == "Bell+":
        psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    elif name == "|01>":
        psi = np.kron(up, dn)
    elif name == "|11>":
        psi = np.kron(dn, dn)
    elif name == "|+,+>":
        psi = np.kron(plus, plus)
    elif name == "Bell-":
        psi = (np.kron(up, up) - np.kron(dn, dn)) / np.sqrt(2)
    elif name == "Psi+":
        psi = (np.kron(up, dn) + np.kron(dn, up)) / np.sqrt(2)
    elif name == "Psi-":
        psi = (np.kron(up, dn) - np.kron(dn, up)) / np.sqrt(2)
    return np.outer(psi, psi.conj())


def ad_ops(gamma_ad, gamma_z=0.0):
    """Build collapse operators: amplitude damping + optional Z-dephasing per qubit."""
    ops = []
    for k in range(2):
        if gamma_ad > 0:
            ops.append(np.sqrt(gamma_ad) * site_op(sm, k))
        if gamma_z > 0:
            ops.append(np.sqrt(gamma_z) * site_op(sz, k))
    return ops


def ad_ops_asym(gamma_ad_0, gamma_ad_1, gamma_z=0.0):
    """Asymmetric amplitude damping: different rates per qubit."""
    ops = []
    if gamma_ad_0 > 0:
        ops.append(np.sqrt(gamma_ad_0) * site_op(sm, 0))
    if gamma_ad_1 > 0:
        ops.append(np.sqrt(gamma_ad_1) * site_op(sm, 1))
    for k in range(2):
        if gamma_z > 0:
            ops.append(np.sqrt(gamma_z) * site_op(sz, k))
    return ops


def find_crossing_detailed(L, rho0, t_max=100.0, dt=0.02, threshold=0.25):
    """Track CΨ trajectory and find crossing details."""
    n_steps = int(round(t_max / dt))
    t_arr = np.zeros(n_steps + 1)
    cpsi_arr = np.zeros(n_steps + 1)
    pur_arr = np.zeros(n_steps + 1)
    conc_arr = np.zeros(n_steps + 1)

    for i in range(n_steps + 1):
        t = i * dt
        t_arr[i] = t
        rho = evolve(L, rho0, t)
        cpsi_arr[i] = cpsi(rho)
        pur_arr[i] = purity(rho)
        conc_arr[i] = concurrence(rho)

    # Find crossing
    t_cross = None
    n_inc_above = 0
    max_inc = 0.0
    for i in range(1, len(t_arr)):
        if cpsi_arr[i-1] > threshold and cpsi_arr[i] <= threshold and t_cross is None:
            frac = (cpsi_arr[i-1] - threshold) / (cpsi_arr[i-1] - cpsi_arr[i]) if cpsi_arr[i-1] != cpsi_arr[i] else 0.5
            t_cross = (i - 1 + frac) * dt
        if cpsi_arr[i] > cpsi_arr[i-1] and cpsi_arr[i-1] > threshold:
            n_inc_above += 1
            inc = cpsi_arr[i] - cpsi_arr[i-1]
            if inc > max_inc:
                max_inc = inc

    return {
        't_cross': t_cross,
        'n_inc_above': n_inc_above,
        'max_inc': max_inc,
        'cpsi_final': cpsi_arr[-1],
        'pur_final': pur_arr[-1],
        'conc_final': conc_arr[-1],
        'cpsi_max': np.max(cpsi_arr),
        'cpsi_min_after_cross': np.min(cpsi_arr[cpsi_arr < threshold]) if np.any(cpsi_arr < threshold) else None,
    }


# ====================================================================
# Test 1: Pure amplitude damping (no Hamiltonian)
# ====================================================================

def test_1_pure_ad():
    log("=" * 70)
    log("TEST 1: PURE AMPLITUDE DAMPING (no Hamiltonian)")
    log("=" * 70)
    log()
    log("  L = √γ |0⟩⟨1| on each qubit, H = 0")
    log("  Fixed point: |00⟩ (CΨ = 0)")
    log()

    H = np.zeros((d, d), dtype=complex)
    rho0 = make_state("Bell+")

    log(f"  {'γ_AD':>6}  {'t_cross':>8}  {'CΨ(0)':>8}  {'CΨ_final':>8}  {'Pur_f':>6}  {'Conc_f':>6}  {'Mono':>5}")
    log("  " + "-" * 58)

    for gamma in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
        c_ops = ad_ops(gamma)
        L = build_L(H, c_ops)
        info = find_crossing_detailed(L, rho0)

        tc = f"{info['t_cross']:.3f}" if info['t_cross'] else "NEVER"
        mono = "YES" if info['n_inc_above'] == 0 else f"NO({info['n_inc_above']})"
        log(f"  {gamma:6.2f}  {tc:>8}  {0.3333:8.4f}  {info['cpsi_final']:8.6f}  "
            f"{info['pur_final']:6.4f}  {info['conc_final']:6.4f}  {mono:>5}")

    log()


# ====================================================================
# Test 2: Amplitude damping + Heisenberg coupling
# ====================================================================

def test_2_ad_with_H():
    log("=" * 70)
    log("TEST 2: AMPLITUDE DAMPING + HEISENBERG COUPLING")
    log("=" * 70)
    log()

    rho0 = make_state("Bell+")

    log(f"  {'J':>4}  {'γ_AD':>6}  {'t_cross':>8}  {'CΨ_final':>8}  {'Pur_f':>6}  {'#↑':>4}  {'Mono':>5}")
    log("  " + "-" * 50)

    for J in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]:
        for gamma in [0.01, 0.05, 0.1, 0.5]:
            H = build_H(J) if J > 0 else np.zeros((d, d), dtype=complex)
            c_ops = ad_ops(gamma)
            L = build_L(H, c_ops)
            info = find_crossing_detailed(L, rho0)

            tc = f"{info['t_cross']:.3f}" if info['t_cross'] else "NEVER"
            mono = "YES" if info['n_inc_above'] == 0 else "NO"
            log(f"  {J:4.1f}  {gamma:6.2f}  {tc:>8}  {info['cpsi_final']:8.6f}  "
                f"{info['pur_final']:6.4f}  {info['n_inc_above']:4d}  {mono:>5}")

    log()


# ====================================================================
# Test 3: K-invariance for amplitude damping
# ====================================================================

def test_3_K_invariance():
    log("=" * 70)
    log("TEST 3: K-INVARIANCE - K = γ_AD × t_cross")
    log("=" * 70)
    log()

    rho0 = make_state("Bell+")

    for J, label in [(0.0, "J=0 (pure AD)"), (1.0, "J=1.0 (AD + Heisenberg)")]:
        H = build_H(J) if J > 0 else np.zeros((d, d), dtype=complex)
        log(f"  {label}")
        log(f"    {'γ_AD':>6}  {'t_cross':>8}  {'K':>8}")
        log("    " + "-" * 28)

        K_vals = []
        for gamma in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
            c_ops = ad_ops(gamma)
            L = build_L(H, c_ops)
            info = find_crossing_detailed(L, rho0, t_max=200.0)

            if info['t_cross']:
                K = gamma * info['t_cross']
                K_vals.append(K)
                log(f"    {gamma:6.3f}  {info['t_cross']:8.3f}  {K:8.4f}")
            else:
                log(f"    {gamma:6.3f}  {'NEVER':>8}  {'---':>8}")

        if len(K_vals) >= 2:
            K_m = np.mean(K_vals)
            K_s = np.std(K_vals)
            cv = K_s / K_m * 100 if K_m > 0 else 0
            log(f"    K = {K_m:.4f} ± {K_s:.4f} (CV={cv:.1f}%)")
            log(f"    K-invariant: {'YES' if cv < 5 else 'NO'}")
        log()


# ====================================================================
# Test 4: Combined amplitude damping + Z-dephasing
# ====================================================================

def test_4_combined():
    log("=" * 70)
    log("TEST 4: COMBINED - Amplitude damping + Z-dephasing")
    log("=" * 70)
    log()

    H = build_H(J=1.0)
    rho0 = make_state("Bell+")

    log(f"  {'γ_AD':>6}  {'γ_Z':>6}  {'γ_tot':>6}  {'t_cross':>8}  {'CΨ_f':>8}  {'Pur_f':>6}  {'Mono':>5}")
    log("  " + "-" * 56)

    for gamma_ad in [0.0, 0.01, 0.05, 0.1]:
        for gamma_z in [0.0, 0.01, 0.05, 0.1]:
            if gamma_ad == 0 and gamma_z == 0:
                continue
            c_ops = ad_ops(gamma_ad, gamma_z)
            L = build_L(H, c_ops)
            info = find_crossing_detailed(L, rho0)

            tc = f"{info['t_cross']:.3f}" if info['t_cross'] else "NEVER"
            mono = "YES" if info['n_inc_above'] == 0 else "NO"
            g_tot = gamma_ad + gamma_z
            log(f"  {gamma_ad:6.2f}  {gamma_z:6.2f}  {g_tot:6.2f}  {tc:>8}  "
                f"{info['cpsi_final']:8.6f}  {info['pur_final']:6.4f}  {mono:>5}")

    log()


# ====================================================================
# Test 5: Various initial states
# ====================================================================

def test_5_states():
    log("=" * 70)
    log("TEST 5: VARIOUS INITIAL STATES")
    log("=" * 70)
    log()

    H = build_H(J=1.0)
    gamma = 0.05

    states = ["Bell+", "Bell-", "Psi+", "Psi-", "|01>", "|11>", "|+,+>"]

    log(f"  γ_AD={gamma}, J=1.0")
    log(f"  {'State':>8}  {'CΨ(0)':>8}  {'t_cross':>8}  {'CΨ_f':>8}  {'Pur_f':>6}  {'#↑':>4}  {'Mono':>5}")
    log("  " + "-" * 56)

    for sname in states:
        rho0 = make_state(sname)
        c_ops = ad_ops(gamma)
        L = build_L(H, c_ops)
        info = find_crossing_detailed(L, rho0)

        cpsi_0 = cpsi(rho0)
        tc = f"{info['t_cross']:.3f}" if info['t_cross'] else "NEVER"
        mono = "YES" if info['n_inc_above'] == 0 else f"NO({info['n_inc_above']})"
        log(f"  {sname:>8}  {cpsi_0:8.4f}  {tc:>8}  {info['cpsi_final']:8.6f}  "
            f"{info['pur_final']:6.4f}  {info['n_inc_above']:4d}  {mono:>5}")

    log()
    log("  Note: states starting below 1/4 (|01>, |11>, |+,+>) show 'NEVER'")
    log("  because they never crossed DOWN through 1/4")
    log()


# ====================================================================
# Test 6: Trajectory detail - CΨ over time for key cases
# ====================================================================

def test_6_trajectory():
    log("=" * 70)
    log("TEST 6: TRAJECTORY DETAIL - CΨ(t) snapshots")
    log("=" * 70)
    log()

    rho0 = make_state("Bell+")

    cases = [
        ("AD only, γ=0.05",       0.0, 0.05, 0.0),
        ("AD + H, γ=0.05, J=1",   1.0, 0.05, 0.0),
        ("AD + Z, γ_AD=0.05, γ_Z=0.05", 1.0, 0.05, 0.05),
        ("Z only, γ=0.05",        1.0, 0.0,  0.05),
    ]

    for label, J, gamma_ad, gamma_z in cases:
        H = build_H(J) if J > 0 else np.zeros((d, d), dtype=complex)
        c_ops = ad_ops(gamma_ad, gamma_z)
        L = build_L(H, c_ops)

        log(f"  {label}")
        log(f"  {'t':>6}  {'CΨ':>8}  {'Purity':>8}  {'Conc':>8}  {'Ψ':>8}  {'Above':>5}")
        log("  " + "-" * 48)

        for t in [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
            rho = evolve(L, rho0, t)
            c = cpsi(rho)
            p = purity(rho)
            co = concurrence(rho)
            pn = psi_norm(rho)
            above = ">" if c > 0.25 else "<"
            log(f"  {t:6.1f}  {c:8.4f}  {p:8.4f}  {co:8.4f}  {pn:8.4f}  {above:>5}")

        log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Amplitude Damping: Direct CΨ Trajectory")
    log("=" * 70)
    log("L = √γ |0⟩⟨1| (spontaneous decay, non-unital)")
    log("Fixed point: |00⟩ (CΨ = 0)")
    log()

    test_1_pure_ad()
    test_2_ad_with_H()
    test_3_K_invariance()
    test_4_combined()
    test_5_states()
    test_6_trajectory()

    log("=" * 70)
    log("VERDICT")
    log("=" * 70)
    log()
    log("Amplitude damping is the NON-UNITAL channel (fixed point = |00⟩, not I/d).")
    log("If CΨ crosses 1/4 for all tested configs: Conjecture 5.1 extends to non-unital.")
    log("If monotonic for Bell+: Conjecture 5.2 extends to non-unital Markovian channels.")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
