#!/usr/bin/env python3
"""
Gravity Chain Test: Schwarzschild Self-Consistency
====================================================
Phase A: Reproduce SELF_CONSISTENCY_SCHWARZSCHILD tables (uncoupled 2-qubit)
Phase B: 8-qubit coupled chain with spatial gamma profiles
Phase C: Comparison and verdict

Tests: does Schwarzschild gamma profile produce unbounded R-concentration
while alternatives plateau?

Script:  simulations/gravity_chain_test.py
Output:  simulations/results/gravity_chain_test.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

# ====================================================================
# Output setup
# ====================================================================

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "gravity_chain_test.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ====================================================================
# Pauli matrices and basic states
# ====================================================================

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)
plus = (up + dn) / np.sqrt(2)

# ====================================================================
# Parameters
# ====================================================================

GAMMA_0 = 0.05
J_COUPLING = 1.0
R_S = 1.0               # Schwarzschild radius (normalized)
DT = 0.01               # RK4 timestep
T_GAMMA_0_VALUES = [0.1, 0.2, 0.3, 0.5, 0.7]

# Phase A positions (match SELF_CONSISTENCY table)
R_PHASE_A = R_S * np.array([1.001, 1.01, 1.05, 1.10, 1.50, 3.0, 10.0, 50.0])

# Phase B: 8 qubits, log-spaced
N_CHAIN = 8
R_PHASE_B = R_S * np.logspace(np.log10(1.1), np.log10(50.0), N_CHAIN)


# ====================================================================
# Section 1: Infrastructure
# ====================================================================

def site_op(op, k, nq):
    """Place operator on qubit k in nq-qubit system."""
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def chain_H(nq, J=1.0):
    """Nearest-neighbor Heisenberg chain Hamiltonian."""
    d = 2 ** nq
    H = np.zeros((d, d), dtype=complex)
    for i in range(nq - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, nq) @ site_op(P, i + 1, nq)
    return H


def build_L_ops(gammas, nq):
    """Z-dephasing Lindblad operators, one per qubit."""
    ops = []
    LdL = []
    for k, g in enumerate(gammas):
        if g > 0:
            Lk = np.sqrt(g) * site_op(sz, k, nq)
            ops.append(Lk)
            LdL.append(Lk.conj().T @ Lk)
    return ops, LdL


def lindblad_rhs(rho, H, L_ops, LdL_ops):
    """Lindblad master equation RHS with pre-computed L†L."""
    drho = -1j * (H @ rho - rho @ H)
    for L, LdL in zip(L_ops, LdL_ops):
        drho += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return drho


def rk4_step(rho, H, L_ops, LdL_ops, dt):
    """RK4 integrator with numerical hygiene."""
    k1 = lindblad_rhs(rho, H, L_ops, LdL_ops)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, L_ops, LdL_ops)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, L_ops, LdL_ops)
    k4 = lindblad_rhs(rho + dt * k3, H, L_ops, LdL_ops)
    rho_new = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    # Hermitianize and renormalize
    rho_new = 0.5 * (rho_new + rho_new.conj().T)
    tr = np.trace(rho_new).real
    if abs(tr) < 1e-14:
        raise FloatingPointError("Trace collapsed.")
    rho_new /= tr
    return rho_new


def partial_trace_keep(rho, keep, n_qubits):
    """Partial trace keeping specified qubits."""
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


def l1_coherence(rho):
    return float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))


def psi_norm(rho):
    d = rho.shape[0]
    if d <= 1:
        return 0.0
    return l1_coherence(rho) / (d - 1)


def concurrence_2q(rho4):
    """Wootters concurrence for 4x4 density matrix."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho4.conj() @ sy2
    R = rho4 @ rho_tilde
    evals = np.linalg.eigvals(R)
    evals = np.real(np.sqrt(np.maximum(evals, 0.0)))
    evals = np.sort(evals)[::-1]
    return float(max(0.0, evals[0] - evals[1] - evals[2] - evals[3]))


def compute_R(rho):
    """R = purity * psi_norm^2 (SELF_CONSISTENCY convention)."""
    return purity(rho) * psi_norm(rho) ** 2


# ====================================================================
# Section 2: Metric profile functions
# ====================================================================

def f_schwarzschild(r, r_s):
    return np.sqrt(1.0 - r_s / r)


def f_inverse(r, r_s):
    return r / (r + r_s)


def f_inv_square(r, r_s):
    return r**2 / (r**2 + r_s**2)


def f_uniform(r, r_s):
    return 1.0


METRICS = [
    ("Schwarzschild", f_schwarzschild),
    ("Inverse",       f_inverse),
    ("Inv-Square",    f_inv_square),
    ("Uniform",       f_uniform),
]


def gamma_profile(positions, r_s, gamma_0, f_metric):
    """Compute per-position gamma: gamma_i = gamma_0 * f(r_i)."""
    return np.array([gamma_0 * f_metric(r, r_s) for r in positions])


# ====================================================================
# Section 3: Phase A — Analytical Baseline (2-qubit Bell+)
# ====================================================================

def run_phase_A():
    log("=" * 70)
    log("PHASE A: ANALYTICAL BASELINE (Uncoupled 2-qubit Bell+)")
    log("=" * 70)
    log()
    log(f"Setup: Bell+ initial state, J={J_COUPLING}, gamma={GAMMA_0}")
    log(f"Method: Liouvillian expm (d²=16)")
    log()

    t0 = _time.time()

    # Build 2-qubit system
    d = 4; d2 = 16
    H = np.zeros((d, d), dtype=complex)
    for P in [sx, sy, sz]:
        H += J_COUPLING * np.kron(P, I2) @ np.kron(I2, P)

    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(2):
        Zk = np.kron(sz, I2) if k == 0 else np.kron(I2, sz)
        L += GAMMA_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))

    # Bell+ initial state
    psi = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(psi, psi.conj())

    # Evolve and record R(t) at fine resolution
    T_max = max(T_GAMMA_0_VALUES) / GAMMA_0  # coordinate time
    dt_fine = 0.01
    n_steps = int(T_max / dt_fine) + 1
    t_arr = np.linspace(0, T_max, n_steps)
    R_arr = np.zeros(n_steps)

    for i, t in enumerate(t_arr):
        v = expm(L * t) @ rho0.flatten()
        rho_t = v.reshape(d, d)
        rho_t = (rho_t + rho_t.conj().T) / 2
        R_arr[i] = compute_R(rho_t)

    log(f"Universal curve computed: {n_steps} points, t=0..{T_max:.1f}")
    log(f"R(0) = {R_arr[0]:.4f}, R({T_max:.0f}) = {R_arr[-1]:.6f}")
    log()

    # ---- R(r) profile at T*gamma_0 = 0.5 ----
    log("--- R(r) Profile at T*gamma_0 = 0.5 ---")
    log()
    header = f"  {'r/r_s':>8}"
    for name, _ in METRICS:
        header += f"  {name:>14}"
    log(header)
    log("  " + "-" * (8 + 16 * len(METRICS)))

    Tg = 0.5
    T_coord = Tg / GAMMA_0
    for r in R_PHASE_A:
        row = f"  {r/R_S:8.3f}"
        for name, f_met in METRICS:
            tau = T_coord * f_met(r, R_S)
            R_val = np.interp(tau, t_arr, R_arr)
            row += f"  {R_val:14.4f}"
        log(row)
    log()

    # ---- Concentration ratios ----
    log("--- Concentration Ratios (R_near / R_far) ---")
    log(f"  r_near = {R_PHASE_A[0]/R_S:.3f}*r_s, r_far = {R_PHASE_A[-1]/R_S:.1f}*r_s")
    log()

    header = f"  {'T*g0':>6}"
    for name, _ in METRICS:
        header += f"  {name:>14}"
    log(header)
    log("  " + "-" * (6 + 16 * len(METRICS)))

    ratios_A = {}  # Store for Phase C comparison
    for Tg in T_GAMMA_0_VALUES:
        T_coord = Tg / GAMMA_0
        row = f"  {Tg:6.1f}"
        for name, f_met in METRICS:
            tau_near = T_coord * f_met(R_PHASE_A[0], R_S)
            tau_far = T_coord * f_met(R_PHASE_A[-1], R_S)
            R_near = np.interp(tau_near, t_arr, R_arr)
            R_far = np.interp(tau_far, t_arr, R_arr)
            ratio = R_near / R_far if R_far > 1e-15 else float('inf')
            row += f"  {ratio:14.1f}"
            ratios_A[(name, Tg)] = ratio
        log(row)
    log()

    # Compare to SELF_CONSISTENCY table
    log("--- Comparison to SELF_CONSISTENCY_SCHWARZSCHILD ---")
    expected = {
        ("Schwarzschild", 0.1): 2.4, ("Schwarzschild", 0.2): 4.6,
        ("Schwarzschild", 0.3): 7.5, ("Schwarzschild", 0.5): 15.3,
        ("Schwarzschild", 0.7): 25.7,
        ("Inverse", 0.1): 1.6, ("Inverse", 0.2): 2.1,
        ("Inverse", 0.3): 2.5, ("Inverse", 0.5): 3.2,
        ("Inverse", 0.7): 3.8,
        ("Inv-Square", 0.1): 1.6, ("Inv-Square", 0.2): 2.1,
        ("Inv-Square", 0.3): 2.6, ("Inv-Square", 0.5): 3.3,
        ("Inv-Square", 0.7): 3.9,
    }
    log(f"  {'Metric':>14}  {'T*g0':>5}  {'Got':>8}  {'Expected':>8}  {'Match':>8}")
    log("  " + "-" * 52)
    all_match = True
    for (name, Tg), exp_val in sorted(expected.items()):
        got = ratios_A.get((name, Tg), 0)
        pct = abs(got - exp_val) / exp_val * 100
        ok = "YES" if pct < 20 else "NO"
        if pct >= 20:
            all_match = False
        log(f"  {name:>14}  {Tg:5.1f}  {got:8.1f}  {exp_val:8.1f}  {ok:>5} ({pct:.0f}%)")
    log()
    log(f"Phase A validation: {'PASS' if all_match else 'FAIL'}")
    log(f"Runtime: {_time.time() - t0:.1f}s")
    log()

    return ratios_A, t_arr, R_arr


# ====================================================================
# Section 4: Phase B — Coupled 8-qubit chain
# ====================================================================

def run_phase_B():
    log("=" * 70)
    log("PHASE B: COUPLED CHAIN (8-qubit Lindblad simulation)")
    log("=" * 70)
    log()
    log(f"N={N_CHAIN}, J={J_COUPLING}, gamma_0={GAMMA_0}, dt={DT}")
    log(f"Initial state: |+>^{N_CHAIN}")
    log(f"Positions (r/r_s): [{', '.join(f'{r/R_S:.2f}' for r in R_PHASE_B)}]")
    log()

    t0_total = _time.time()
    d = 2 ** N_CHAIN

    # Build chain Hamiltonian (once for all metrics)
    H = chain_H(N_CHAIN, J_COUPLING)

    # Initial state: |+>^N
    psi_plus_n = plus
    for _ in range(N_CHAIN - 1):
        psi_plus_n = np.kron(psi_plus_n, plus)
    rho0 = np.outer(psi_plus_n, psi_plus_n.conj())

    # Time snapshots
    T_coords = [Tg / GAMMA_0 for Tg in T_GAMMA_0_VALUES]
    T_max = max(T_coords)

    results_B = {}

    for met_name, f_met in METRICS:
        t0_met = _time.time()
        gammas = gamma_profile(R_PHASE_B, R_S, GAMMA_0, f_met)
        log(f"--- Metric: {met_name} ---")
        log(f"  gammas: [{', '.join(f'{g:.4f}' for g in gammas)}]")

        L_ops, LdL_ops = build_L_ops(gammas, N_CHAIN)
        rho = rho0.copy()

        # Evolve with RK4, snapshot at each T_coord
        snapshot_idx = 0
        t_current = 0.0
        n_steps = int(round(T_max / DT))

        # Header for single-qubit R
        log()
        log(f"  {'T*g0':>6}  " + "  ".join(f"  R_q{i:d}" for i in range(N_CHAIN))
            + "  Ratio")
        log("  " + "-" * (6 + 8 * N_CHAIN + 8))

        for step in range(n_steps + 1):
            t_current = step * DT

            # Check if we hit a snapshot
            if snapshot_idx < len(T_coords) and t_current >= T_coords[snapshot_idx] - DT/2:
                Tg = T_GAMMA_0_VALUES[snapshot_idx]

                # Single-qubit R at each site
                R_sites = []
                for q in range(N_CHAIN):
                    rho_q = partial_trace_keep(rho, [q], N_CHAIN)
                    R_sites.append(compute_R(rho_q))

                ratio = R_sites[0] / R_sites[-1] if R_sites[-1] > 1e-15 else float('inf')
                results_B[(met_name, Tg)] = {
                    'R_sites': R_sites[:],
                    'ratio': ratio,
                }

                row = f"  {Tg:6.1f}  " + "  ".join(f"{R:7.4f}" for R in R_sites)
                row += f"  {ratio:6.1f}"
                log(row)

                snapshot_idx += 1

            # RK4 step (skip on last iteration)
            if step < n_steps:
                rho = rk4_step(rho, H, L_ops, LdL_ops, DT)

        elapsed = _time.time() - t0_met
        log(f"  Runtime: {elapsed:.1f}s")
        log()

    log(f"Phase B total runtime: {_time.time() - t0_total:.1f}s")
    log()

    return results_B


# ====================================================================
# Section 5: Phase C — Comparison + Verdict
# ====================================================================

def run_phase_C(ratios_A, results_B):
    log("=" * 70)
    log("PHASE C: COMPARISON (Coupled vs Uncoupled)")
    log("=" * 70)
    log()
    log("Note: Phase A uses 2-qubit Bell+ at 8 discrete positions (r/r_s = 1.001..50)")
    log("      Phase B uses 8-qubit chain with |+>^8 at log-spaced positions (r/r_s = 1.1..50)")
    log("      Comparison is qualitative (same trend?) not quantitative (same values)")
    log()

    header = f"  {'Metric':>14}  {'T*g0':>5}  {'Uncoupled':>10}  {'Coupled':>10}  {'Trend':>8}"
    log(header)
    log("  " + "-" * 55)

    for Tg in T_GAMMA_0_VALUES:
        for met_name, _ in METRICS:
            r_a = ratios_A.get((met_name, Tg), 0)
            b_data = results_B.get((met_name, Tg), {})
            r_b = b_data.get('ratio', 0)
            if r_a > 0 and r_b > 0:
                trend = "SAME" if abs(r_b - r_a) / max(r_a, 1) < 0.5 else (
                    "HIGHER" if r_b > r_a else "LOWER")
            else:
                trend = "N/A"
            log(f"  {met_name:>14}  {Tg:5.1f}  {r_a:10.1f}  {r_b:10.1f}  {trend:>8}")

    log()
    log("=" * 70)
    log("VERDICT")
    log("=" * 70)
    log()

    # Test 1: Schwarzschild grows monotonically (Phase B)
    schw_ratios = []
    for Tg in T_GAMMA_0_VALUES:
        d = results_B.get(("Schwarzschild", Tg), {})
        schw_ratios.append(d.get('ratio', 0))

    monotonic = all(schw_ratios[i] < schw_ratios[i+1]
                    for i in range(len(schw_ratios)-1))
    log(f"1. Schwarzschild concentration grows monotonically: "
        f"{'PASS' if monotonic else 'FAIL'}")
    log(f"   Ratios: [{', '.join(f'{r:.1f}' for r in schw_ratios)}]")

    # Test 2: Alternatives grow slower
    for alt_name in ["Inverse", "Inv-Square", "Uniform"]:
        alt_ratios = []
        for Tg in T_GAMMA_0_VALUES:
            d = results_B.get((alt_name, Tg), {})
            alt_ratios.append(d.get('ratio', 0))

        # Check growth rate: Schwarzschild should pull ahead
        if len(schw_ratios) >= 2 and len(alt_ratios) >= 2:
            schw_growth = schw_ratios[-1] / max(schw_ratios[0], 0.01)
            alt_growth = alt_ratios[-1] / max(alt_ratios[0], 0.01)
            passes = schw_growth > alt_growth * 1.5
        else:
            passes = False

        log(f"2. {alt_name} grows slower than Schwarzschild: "
            f"{'PASS' if passes else 'FAIL'}")
        log(f"   Ratios: [{', '.join(f'{r:.1f}' for r in alt_ratios)}]")
        log(f"   Growth factor: Schw={schw_growth:.1f}x, {alt_name}={alt_growth:.1f}x")

    # Test 3: Coupling effect
    log()
    log("3. Coupling effect on metric discrimination:")
    for met_name in ["Schwarzschild", "Inverse", "Inv-Square"]:
        Tg = 0.5
        r_a = ratios_A.get((met_name, Tg), 0)
        r_b = results_B.get((met_name, Tg), {}).get('ratio', 0)
        if r_a > 0:
            pct = (r_b - r_a) / r_a * 100
            log(f"   {met_name}: uncoupled={r_a:.1f}, coupled={r_b:.1f} ({pct:+.0f}%)")

    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Gravity Chain Test: Schwarzschild Self-Consistency")
    log("=" * 70)
    log(f"N={N_CHAIN}, J={J_COUPLING}, gamma_0={GAMMA_0}, r_s={R_S}")
    log(f"dt={DT}, T*gamma_0 = {T_GAMMA_0_VALUES}")
    log()

    ratios_A, t_arr, R_arr = run_phase_A()
    results_B = run_phase_B()
    run_phase_C(ratios_A, results_B)

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log()
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
