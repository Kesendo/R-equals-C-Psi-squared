#!/usr/bin/env python3
"""
Generalized Pauli Channels: 1/4 Boundary Test
================================================
ℰ(ρ) = (1-p)ρ + p_x σ_x ρ σ_x + p_y σ_y ρ σ_y + p_z σ_z ρ σ_z

Test 1: Continuous Lindblad with (γ_x, γ_y, γ_z) sweep — does CΨ always cross 1/4?
Test 2: Monotonicity check — is CΨ non-increasing above 1/4?
Test 3: Palindromic eigenvalue pairing — does it survive non-Z noise?
Test 4: K-invariance — does K = γ_eff × t_cross hold across noise types?
Test 5: Discrete channel application — repeated ℰ^n

Script:  simulations/generalized_pauli_channels.py
Output:  simulations/results/generalized_pauli_channels.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time
from itertools import product as iprod

# ====================================================================
# Output
# ====================================================================

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "generalized_pauli_channels.txt")
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
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)

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


def build_L(H, noise_params):
    """Build Liouvillian for generalized Pauli noise on both qubits.
    noise_params = [(gamma_x, gamma_y, gamma_z)] per qubit, or single tuple for uniform.
    """
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))

    if isinstance(noise_params, tuple):
        noise_params = [noise_params, noise_params]

    paulis = [sx, sy, sz]
    for k, (gx, gy, gz) in enumerate(noise_params):
        for g, P in zip([gx, gy, gz], paulis):
            if g > 0:
                Pk = site_op(P, k)
                L += g * (np.kron(Pk, Pk.conj()) - np.eye(d2))
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


def make_bell_plus():
    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    return np.outer(psi, psi.conj())


def make_01():
    psi = np.kron(up, dn)
    return np.outer(psi, psi.conj())


def find_crossing(L, rho0, threshold=0.25, t_max=200.0, dt=0.02):
    """Find the time CΨ first crosses below threshold.
    Returns (t_cross, max_increase, n_increases, cpsi_trajectory_sample).
    """
    t_cross = None
    prev_cpsi = cpsi(rho0)
    max_increase = 0.0
    n_increases = 0
    n_steps = int(round(t_max / dt))

    for step in range(1, n_steps + 1):
        t = step * dt
        rho = evolve(L, rho0, t)
        c = cpsi(rho)

        if prev_cpsi > threshold and c <= threshold and t_cross is None:
            # Linear interpolation
            frac = (prev_cpsi - threshold) / (prev_cpsi - c) if prev_cpsi != c else 0.5
            t_cross = (step - 1 + frac) * dt

        if c > prev_cpsi and prev_cpsi > threshold:
            inc = c - prev_cpsi
            if inc > max_increase:
                max_increase = inc
            n_increases += 1

        prev_cpsi = c

    return t_cross, max_increase, n_increases, prev_cpsi


def eigenvalue_pairing(L):
    """Check palindromic eigenvalue pairing of the Liouvillian."""
    evals = np.linalg.eigvals(L)
    rates = np.sort(np.real(evals))[::-1]  # Descending
    # Remove the zero eigenvalue (steady state)
    rates_nonzero = rates[np.abs(rates) > 1e-10]

    if len(rates_nonzero) == 0:
        return 0, 0, 0.0

    # Check pairing: r_k + r_{N-1-k} = const?
    n = len(rates_nonzero)
    sums = []
    for i in range(n // 2):
        sums.append(rates_nonzero[i] + rates_nonzero[n - 1 - i])

    if len(sums) == 0:
        return 0, 0, 0.0

    mean_sum = np.mean(sums)
    max_err = max(abs(s - mean_sum) for s in sums)
    n_paired = sum(1 for s in sums if abs(s - mean_sum) < 0.01 * abs(mean_sum + 1e-15))

    return n_paired, len(sums), max_err


# ====================================================================
# Test 1: Sweep (γ_x, γ_y, γ_z) — does CΨ always cross 1/4?
# ====================================================================

def test_1_sweep():
    log("=" * 70)
    log("TEST 1: GENERALIZED PAULI SWEEP — Does CΨ always cross 1/4?")
    log("=" * 70)
    log()

    H = build_H(J=1.0)
    rho0 = make_bell_plus()

    # Sweep: each gamma from 0 to 0.2 in steps
    gamma_vals = [0.0, 0.01, 0.05, 0.1, 0.2]
    total = 0
    crossed = 0
    never_crossed = []

    header = f"  {'γ_x':>6}  {'γ_y':>6}  {'γ_z':>6}  {'γ_eff':>6}  {'t_cross':>8}  {'K=γ·t':>8}  {'Mono':>5}  {'Final':>8}"
    log(header)
    log("  " + "-" * 68)

    for gx, gy, gz in iprod(gamma_vals, repeat=3):
        if gx == 0 and gy == 0 and gz == 0:
            continue  # Skip no-noise case
        total += 1
        gamma_eff = gx + gy + gz

        L = build_L(H, (gx, gy, gz))
        t_cross, max_inc, n_inc, final_c = find_crossing(L, rho0)

        mono = "YES" if n_inc == 0 else f"NO({n_inc})"

        if t_cross is not None:
            crossed += 1
            K = gamma_eff * t_cross
            log(f"  {gx:6.2f}  {gy:6.2f}  {gz:6.2f}  {gamma_eff:6.2f}  "
                f"{t_cross:8.3f}  {K:8.4f}  {mono:>5}  {final_c:8.6f}")
        else:
            never_crossed.append((gx, gy, gz, final_c))
            log(f"  {gx:6.2f}  {gy:6.2f}  {gz:6.2f}  {gamma_eff:6.2f}  "
                f"{'NEVER':>8}  {'---':>8}  {mono:>5}  {final_c:8.6f}")

    log()
    log(f"  Total configs: {total}")
    log(f"  Crossed 1/4: {crossed}/{total}")
    if never_crossed:
        log(f"  NEVER crossed: {len(never_crossed)}")
        for gx, gy, gz, fc in never_crossed:
            log(f"    γ=({gx},{gy},{gz}), final CΨ={fc:.6f}")
    else:
        log(f"  ALL configurations crossed 1/4. Conjecture 5.1 holds.")
    log()


# ====================================================================
# Test 2: Monotonicity — is CΨ non-increasing above 1/4?
# ====================================================================

def test_2_monotonicity():
    log("=" * 70)
    log("TEST 2: MONOTONICITY — Is CΨ non-increasing above 1/4?")
    log("=" * 70)
    log()

    H = build_H(J=1.0)
    rho0_bell = make_bell_plus()
    rho0_01 = make_01()

    # Representative noise types
    configs = [
        ("Pure Z",      (0.0, 0.0, 0.05)),
        ("Pure X",      (0.05, 0.0, 0.0)),
        ("Pure Y",      (0.0, 0.05, 0.0)),
        ("Depolarizing", (0.017, 0.017, 0.017)),
        ("X+Z",         (0.025, 0.0, 0.025)),
        ("Y+Z",         (0.0, 0.025, 0.025)),
        ("X+Y",         (0.025, 0.025, 0.0)),
        ("Asymmetric",  (0.1, 0.01, 0.03)),
        ("Strong X",    (0.2, 0.0, 0.0)),
        ("Strong Z",    (0.0, 0.0, 0.2)),
    ]

    for state_name, rho0 in [("Bell+", rho0_bell), ("|01>", rho0_01)]:
        log(f"  State: {state_name}")
        log(f"  {'Type':>14}  {'γ_eff':>6}  {'t_cross':>8}  {'#↑>1/4':>6}  {'MaxΔ':>8}  {'Mono':>5}")
        log("  " + "-" * 54)

        for name, (gx, gy, gz) in configs:
            gamma_eff = gx + gy + gz
            L = build_L(H, (gx, gy, gz))
            t_cross, max_inc, n_inc, final_c = find_crossing(L, rho0)

            mono = "YES" if n_inc == 0 else "NO"
            tc_str = f"{t_cross:.3f}" if t_cross else "NEVER"
            log(f"  {name:>14}  {gamma_eff:6.3f}  {tc_str:>8}  {n_inc:6d}  "
                f"{max_inc:8.6f}  {mono:>5}")

        log()

    log("  Note: '#↑>1/4' counts timesteps where CΨ increased while above 1/4")
    log("  If #↑>1/4 = 0 for all: Conjecture 5.2 (Markovian monotonicity) holds")
    log()


# ====================================================================
# Test 3: Palindromic eigenvalue pairing
# ====================================================================

def test_3_palindrome():
    log("=" * 70)
    log("TEST 3: PALINDROMIC PAIRING — Does it survive non-Z noise?")
    log("=" * 70)
    log()

    H = build_H(J=1.0)

    configs = [
        ("Pure Z (γ=0.05)",   (0.0, 0.0, 0.05)),
        ("Pure X (γ=0.05)",   (0.05, 0.0, 0.0)),
        ("Pure Y (γ=0.05)",   (0.0, 0.05, 0.0)),
        ("Depolarizing",      (0.017, 0.017, 0.017)),
        ("X+Z symmetric",     (0.025, 0.0, 0.025)),
        ("Y+Z symmetric",     (0.0, 0.025, 0.025)),
        ("X+Y symmetric",     (0.025, 0.025, 0.0)),
        ("X+Y+Z uniform",     (0.05, 0.05, 0.05)),
        ("Asymmetric",        (0.1, 0.01, 0.03)),
        ("Pure Z (γ=0.2)",    (0.0, 0.0, 0.2)),
        ("Pure X (γ=0.2)",    (0.2, 0.0, 0.0)),
    ]

    log(f"  {'Type':>20}  {'Paired':>8}  {'Total':>6}  {'MaxErr':>10}  {'Status':>10}")
    log("  " + "-" * 60)

    for name, (gx, gy, gz) in configs:
        L = build_L(H, (gx, gy, gz))
        n_paired, n_total, max_err = eigenvalue_pairing(L)

        if n_total > 0:
            pct = n_paired / n_total * 100
            status = "PALINDROME" if pct == 100 else f"{pct:.0f}%"
        else:
            status = "N/A"

        log(f"  {name:>20}  {n_paired:>4}/{n_total:<3}  {n_total:6d}  "
            f"{max_err:10.2e}  {status:>10}")

    log()
    log("  Palindromic pairing requires: r_k + r_{N-1-k} = const for all k")
    log("  If only Z-noise produces 100% pairing: palindrome is Z-specific")
    log()


# ====================================================================
# Test 4: K-invariance across noise types
# ====================================================================

def test_4_K_invariance():
    log("=" * 70)
    log("TEST 4: K-INVARIANCE — Does K = γ_eff × t_cross hold?")
    log("=" * 70)
    log()

    H = build_H(J=1.0)
    rho0 = make_bell_plus()

    # For each noise type, sweep overall strength and check K = const
    noise_types = [
        ("Pure Z",       lambda g: (0.0, 0.0, g)),
        ("Pure X",       lambda g: (g, 0.0, 0.0)),
        ("Pure Y",       lambda g: (0.0, g, 0.0)),
        ("Depolarizing", lambda g: (g/3, g/3, g/3)),
        ("X+Z equal",    lambda g: (g/2, 0.0, g/2)),
    ]

    gamma_strengths = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

    for type_name, noise_fn in noise_types:
        log(f"  Noise type: {type_name}")
        log(f"    {'γ_eff':>6}  {'t_cross':>8}  {'K':>8}")
        log("    " + "-" * 28)

        K_values = []
        for g in gamma_strengths:
            gx, gy, gz = noise_fn(g)
            gamma_eff = gx + gy + gz
            L = build_L(H, (gx, gy, gz))
            t_cross, _, _, _ = find_crossing(L, rho0)

            if t_cross is not None:
                K = gamma_eff * t_cross
                K_values.append(K)
                log(f"    {gamma_eff:6.3f}  {t_cross:8.3f}  {K:8.4f}")
            else:
                log(f"    {gamma_eff:6.3f}  {'NEVER':>8}  {'---':>8}")

        if len(K_values) >= 2:
            K_mean = np.mean(K_values)
            K_std = np.std(K_values)
            K_cv = K_std / K_mean * 100 if K_mean > 0 else 0
            log(f"    K = {K_mean:.4f} ± {K_std:.4f} (CV={K_cv:.1f}%)")
            log(f"    K-invariant: {'YES' if K_cv < 5 else 'NO'}")
        log()


# ====================================================================
# Test 5: Discrete channel application ℰ^n
# ====================================================================

def test_5_discrete():
    log("=" * 70)
    log("TEST 5: DISCRETE CHANNEL — Repeated ℰ^n application")
    log("=" * 70)
    log()

    rho0 = make_bell_plus()

    # Build discrete Pauli channel
    def apply_pauli_channel(rho, px, py, pz):
        """Apply ℰ(ρ) = (1-p)ρ + px X ρ X + py Y ρ Y + pz Z ρ Z on each qubit."""
        p0 = 1.0 - px - py - pz
        paulis = [(p0, I2), (px, sx), (py, sy), (pz, sz)]
        rho_out = np.zeros_like(rho)
        for p1, P1 in paulis:
            for p2, P2 in paulis:
                K = np.kron(P1, P2)
                rho_out += p1 * p2 * (K @ rho @ K.conj().T)
        return rho_out

    configs = [
        ("Pure Z (p=0.05)",      0.0,   0.0,   0.05),
        ("Pure X (p=0.05)",      0.05,  0.0,   0.0),
        ("Pure Y (p=0.05)",      0.0,   0.05,  0.0),
        ("Depolarizing (p=0.05)", 0.017, 0.017, 0.017),
        ("Strong Z (p=0.2)",     0.0,   0.0,   0.2),
        ("Strong X (p=0.2)",     0.2,   0.0,   0.0),
        ("Asymmetric",           0.1,   0.02,  0.05),
        ("Near-depol (p=0.3)",   0.1,   0.1,   0.1),
    ]

    n_max = 500

    log(f"  {'Channel':>24}  {'n_cross':>7}  {'CΨ(0)':>8}  {'CΨ(n_c)':>8}  {'CΨ({n_max})':>8}  {'Mono':>5}")
    log("  " + "-" * 68)

    for name, px, py, pz in configs:
        rho = rho0.copy()
        crossed = False
        n_cross = None
        prev_c = cpsi(rho)
        n_inc = 0

        for n in range(1, n_max + 1):
            rho = apply_pauli_channel(rho, px, py, pz)
            c = cpsi(rho)

            if c > prev_c and prev_c > 0.25:
                n_inc += 1

            if not crossed and c <= 0.25:
                crossed = True
                n_cross = n

            prev_c = c

        c_init = cpsi(rho0)
        mono = "YES" if n_inc == 0 else f"NO({n_inc})"
        nc_str = str(n_cross) if n_cross else "NEVER"
        log(f"  {name:>24}  {nc_str:>7}  {c_init:8.4f}  "
            f"{0.25:8.4f}  {prev_c:8.6f}  {mono:>5}")

    log()
    log("  If all channels cross 1/4: Conjecture 5.1 holds for Pauli channels")
    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Generalized Pauli Channels: 1/4 Boundary Test")
    log("=" * 70)
    log("ℰ(ρ) = (1-p)ρ + p_x σ_x ρ σ_x + p_y σ_y ρ σ_y + p_z σ_z ρ σ_z")
    log()

    test_1_sweep()
    test_2_monotonicity()
    test_3_palindrome()
    test_4_K_invariance()
    test_5_discrete()

    log("=" * 70)
    log("OVERALL VERDICT")
    log("=" * 70)
    log()
    log("Layer 5 (Channel Independence) status for generalized Pauli channels:")
    log("  - Does CΨ always cross 1/4? → Test 1")
    log("  - Is the trajectory monotonic above 1/4? → Test 2")
    log("  - Does palindromic pairing survive non-Z noise? → Test 3")
    log("  - Does K-invariance hold across noise types? → Test 4")
    log("  - Does the discrete channel also work? → Test 5")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
