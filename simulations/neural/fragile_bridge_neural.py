#!/usr/bin/env python3
"""
Fragile Bridge -- Neural Validation (Wilson-Cowan)
===================================================
Does a Wilson-Cowan E-I model show the same three-regime structure
as the quantum gain-loss system?

A) Hopf threshold P_crit as function of E-I coupling scale s
B) Fine sweep around critical points
C) Comparison with quantum fragile bridge
D) Nonlinear validation (full simulation)
E) Epilepsy test (strong coupling regime)

Script: simulations/neural/fragile_bridge_neural.py
Output: simulations/results/fragile_bridge_neural.txt
"""

import numpy as np
import os, sys, time as clock

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "results", "fragile_bridge_neural.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# === Wilson-Cowan primitives ===

def sigmoid(x, a=1.3, theta=4.0):
    arg = np.clip(-a * (x - theta), -500, 500)
    return 1.0 / (1.0 + np.exp(arg))


def dsigmoid(x, a=1.3, theta=4.0):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def find_fixed_point(s, P, tau_E=8.0, tau_I=18.0,
                     w_EE=16.0, w_EI_base=12.0, w_IE_base=15.0, w_II=3.0,
                     a_E=1.3, a_I=2.0, theta_E=4.0, theta_I=3.7,
                     n_iter=10000, damping=0.01):
    """Find fixed point (E*, I*) by iteration."""
    w_EI = s * w_EI_base
    w_IE = s * w_IE_base
    E, I = 0.2, 0.4
    for _ in range(n_iter):
        inp_E = w_EE * E - w_EI * I + P
        inp_I = w_IE * E - w_II * I
        E_new = sigmoid(inp_E, a_E, theta_E)
        I_new = sigmoid(inp_I, a_I, theta_I)
        E = E + damping * (E_new - E)
        I = I + damping * (I_new - I)
    return E, I


def jacobian_eigenvalues(s, P, tau_E=8.0, tau_I=18.0,
                         w_EE=16.0, w_EI_base=12.0, w_IE_base=15.0, w_II=3.0,
                         a_E=1.3, a_I=2.0, theta_E=4.0, theta_I=3.7):
    """Compute Jacobian eigenvalues at fixed point."""
    w_EI = s * w_EI_base
    w_IE = s * w_IE_base
    E_star, I_star = find_fixed_point(s, P)

    inp_E = w_EE * E_star - w_EI * I_star + P
    inp_I = w_IE * E_star - w_II * I_star

    fE = dsigmoid(inp_E, a_E, theta_E)
    fI = dsigmoid(inp_I, a_I, theta_I)

    # Jacobian: dE/dt = (-E + sigmoid(inp_E)) / tau_E
    #           dI/dt = (-I + sigmoid(inp_I)) / tau_I
    J = np.array([
        [(-1 + w_EE * fE) / tau_E,   (-w_EI * fE) / tau_E],
        [(w_IE * fI) / tau_I,         (-1 - w_II * fI) / tau_I]
    ])

    evals = np.linalg.eigvals(J)
    return evals, E_star, I_star


def find_P_crit(s, P_lo=0.0, P_hi=10.0, tol=1e-5):
    """Bisect to find P where fixed point becomes unstable (Hopf)."""
    # Check if already unstable at P_lo
    evals, _, _ = jacobian_eigenvalues(s, P_lo)
    if np.max(evals.real) > 0:
        return 0.0  # always unstable

    # Check if stable even at P_hi
    evals, _, _ = jacobian_eigenvalues(s, P_hi)
    if np.max(evals.real) <= 0:
        return None  # always stable (no Hopf)

    while (P_hi - P_lo) > tol:
        P_mid = (P_lo + P_hi) / 2
        evals, _, _ = jacobian_eigenvalues(s, P_mid)
        if np.max(evals.real) > 0:
            P_hi = P_mid
        else:
            P_lo = P_mid

    return (P_lo + P_hi) / 2


def simulate_wc(s, P, dt=0.1, t_max=500.0,
                tau_E=8.0, tau_I=18.0,
                w_EE=16.0, w_EI_base=12.0, w_IE_base=15.0, w_II=3.0,
                a_E=1.3, a_I=2.0, theta_E=4.0, theta_I=3.7):
    """Full nonlinear Wilson-Cowan simulation, single node."""
    w_EI = s * w_EI_base
    w_IE = s * w_IE_base
    n_steps = int(t_max / dt)
    E, I = 0.2, 0.4
    E_trace = np.zeros(n_steps)
    I_trace = np.zeros(n_steps)

    for t in range(n_steps):
        inp_E = w_EE * E - w_EI * I + P
        inp_I = w_IE * E - w_II * I
        dE = (-E + sigmoid(inp_E, a_E, theta_E)) / tau_E
        dI = (-I + sigmoid(inp_I, a_I, theta_I)) / tau_I
        E += dE * dt
        I += dI * dt
        E_trace[t] = E
        I_trace[t] = I

    return E_trace, I_trace


# ================================================================
log("=" * 70)
log("FRAGILE BRIDGE -- NEURAL VALIDATION (Wilson-Cowan)")
log(f"Started: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 70)

# ================================================================
# A) P_crit as function of E-I coupling scale s
# ================================================================
log()
log("=" * 70)
log("AUFGABE A: P_crit(s) -- Hopf threshold vs E-I coupling")
log("w_EI = s*12.0, w_IE = s*15.0, w_EE=16.0, w_II=3.0 fixed")
log("=" * 70)
log()

s_values = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]
results_a = []

log(f"  {'s':>5}  {'P_crit':>8}  {'Re(l)':>10}  {'Im(l)':>10}  {'freq_Hz':>8}  {'type':>6}  {'E*':>6}  {'I*':>6}")
log(f"  {'-'*75}")

for s in s_values:
    pc = find_P_crit(s)

    if pc is not None and pc > 0:
        # Eigenvalues at P just above P_crit
        evals, E_star, I_star = jacobian_eigenvalues(s, pc * 1.01)
        idx = np.argmax(evals.real)
        ev = evals[idx]
        freq = abs(ev.imag) / (2 * np.pi) * 1000  # convert from 1/ms to Hz
        bif_type = "Hopf" if abs(ev.imag) > 1e-4 else "Saddle"
        results_a.append((s, pc, ev.real, ev.imag, freq, bif_type))
        log(f"  {s:>5.1f}  {pc:>8.4f}  {ev.real:>10.6f}  {ev.imag:>10.4f}  {freq:>8.1f}  {bif_type:>6}  {E_star:>6.3f}  {I_star:>6.3f}")
    elif pc == 0.0:
        evals, E_star, I_star = jacobian_eigenvalues(s, 0.0)
        idx = np.argmax(evals.real)
        ev = evals[idx]
        freq = abs(ev.imag) / (2 * np.pi) * 1000
        bif_type = "Hopf" if abs(ev.imag) > 1e-4 else "Saddle"
        results_a.append((s, 0.0, ev.real, ev.imag, freq, bif_type))
        log(f"  {s:>5.1f}  {'always':>8}  {ev.real:>10.6f}  {ev.imag:>10.4f}  {freq:>8.1f}  {bif_type:>6}  {E_star:>6.3f}  {I_star:>6.3f}")
    else:
        evals, E_star, I_star = jacobian_eigenvalues(s, 5.0)
        results_a.append((s, None, None, None, None, "stable"))
        log(f"  {s:>5.1f}  {'never':>8}  {'stable':>10}  {'':>10}  {'':>8}  {'':>6}  {E_star:>6.3f}  {I_star:>6.3f}")

# Find maximum P_crit
valid_pc = [(s, pc) for s, pc, *_ in results_a if pc is not None and pc > 0]
if valid_pc:
    s_max, pc_max = max(valid_pc, key=lambda x: x[1])
    log()
    log(f"  MAXIMUM P_crit = {pc_max:.4f} at s = {s_max:.1f}")

# ================================================================
# B) Fine sweep around critical points
# ================================================================
log()
log("=" * 70)
log("AUFGABE B: Fine sweep")
log("=" * 70)
log()

# Determine sweep range from A results
s_fine = np.linspace(0.1, 10.0, 100)
results_b = []

log(f"  {'s':>6}  {'P_crit':>8}  {'Im(l)':>10}  {'freq_Hz':>8}")
log(f"  {'-'*40}")

for s in s_fine:
    pc = find_P_crit(s, tol=1e-4)
    if pc is not None and pc > 0:
        evals, _, _ = jacobian_eigenvalues(s, pc * 1.01)
        idx = np.argmax(evals.real)
        ev = evals[idx]
        freq = abs(ev.imag) / (2 * np.pi) * 1000
        results_b.append((s, pc, ev.imag, freq))
        if abs(s - round(s)) < 0.06 or pc > 0.95 * pc_max:
            log(f"  {s:>6.2f}  {pc:>8.4f}  {ev.imag:>10.4f}  {freq:>8.1f}")
    else:
        results_b.append((s, 0.0 if pc == 0 else None, None, None))

# Find peak in fine sweep
valid_fine = [(s, pc) for s, pc, *_ in results_b if pc is not None and pc > 0]
if valid_fine:
    s_peak, pc_peak = max(valid_fine, key=lambda x: x[1])
    log()
    log(f"  FINE PEAK: P_crit = {pc_peak:.4f} at s = {s_peak:.2f}")

# ================================================================
# C) Comparison with quantum
# ================================================================
log()
log("=" * 70)
log("AUFGABE C: Normalized comparison (quantum vs neural)")
log("=" * 70)
log()

# Quantum data from fragile_bridge_bifurcation.py
quantum_data = [
    (0.01, 0.001628), (0.02, 0.003279), (0.05, 0.008364),
    (0.1, 0.017292), (0.2, 0.036861), (0.5, 0.093820),
    (1.0, 0.187310), (2.0, 0.383837), (5.0, 0.131967)
]
q_jb = np.array([d[0] for d in quantum_data])
q_gc = np.array([d[1] for d in quantum_data])
q_peak_idx = np.argmax(q_gc)
q_jb_norm = q_jb / q_jb[q_peak_idx]
q_gc_norm = q_gc / q_gc[q_peak_idx]

# Neural data
n_s = np.array([s for s, pc, *_ in results_b if pc is not None and pc > 0])
n_pc = np.array([pc for s, pc, *_ in results_b if pc is not None and pc > 0])
if len(n_s) > 0 and len(valid_fine) > 0:
    n_s_norm = n_s / s_peak
    n_pc_norm = n_pc / pc_peak

    log("  Normalized comparison (x = coupling/peak_coupling, y = threshold/peak):")
    log()
    log(f"  {'x_quantum':>10}  {'y_quantum':>10}  |  {'x_neural':>10}  {'y_neural':>10}")
    log(f"  {'-'*50}")

    # Sample neural at comparable normalized x values
    for qx, qy in zip(q_jb_norm, q_gc_norm):
        # Find closest neural point
        if len(n_s_norm) > 0:
            idx = np.argmin(np.abs(n_s_norm - qx))
            if abs(n_s_norm[idx] - qx) < 0.5:
                log(f"  {qx:>10.3f}  {qy:>10.3f}  |  {n_s_norm[idx]:>10.3f}  {n_pc_norm[idx]:>10.3f}")
            else:
                log(f"  {qx:>10.3f}  {qy:>10.3f}  |  {'---':>10}  {'---':>10}")

# ================================================================
# D) Nonlinear validation
# ================================================================
log()
log("=" * 70)
log("AUFGABE D: Nonlinear validation (full simulation)")
log("=" * 70)
log()

test_cases = []
if valid_fine:
    # Weak coupling
    s_weak = 0.5
    pc_weak = find_P_crit(s_weak)
    if pc_weak and pc_weak > 0:
        test_cases.append(("weak", s_weak, pc_weak))

    # Optimal coupling (peak)
    pc_opt = find_P_crit(s_peak)
    if pc_opt and pc_opt > 0:
        test_cases.append(("optimal", s_peak, pc_opt))

    # Strong coupling
    s_strong = min(5.0, s_peak * 3)
    pc_strong = find_P_crit(s_strong)
    if pc_strong and pc_strong > 0:
        test_cases.append(("strong", s_strong, pc_strong))

for label, s, pc in test_cases:
    log(f"  {label}: s={s:.2f}, P_crit={pc:.4f}")

    for P_ratio, desc in [(0.9, "below"), (1.1, "above")]:
        P_test = pc * P_ratio
        E_trace, I_trace = simulate_wc(s, P_test, t_max=500.0)

        # Analyze last 200ms
        E_last = E_trace[-2000:]
        amplitude = np.max(E_last) - np.min(E_last)
        mean_E = np.mean(E_last)

        # Frequency from zero crossings
        E_centered = E_last - mean_E
        crossings = np.where(np.diff(np.sign(E_centered)))[0]
        if len(crossings) >= 2:
            period_ms = 2 * np.mean(np.diff(crossings)) * 0.1  # dt=0.1ms
            freq_hz = 1000.0 / period_ms if period_ms > 0 else 0
        else:
            freq_hz = 0

        oscillating = "YES" if amplitude > 0.01 else "no"
        log(f"    P={P_test:.4f} ({desc}): amp={amplitude:.4f}  "
            f"osc={oscillating}  freq={freq_hz:.1f} Hz")

    # Jacobian frequency for comparison
    evals, _, _ = jacobian_eigenvalues(s, pc * 1.01)
    jac_freq = abs(evals[np.argmax(evals.real)].imag) / (2 * np.pi) * 1000
    log(f"    Jacobian predicted frequency: {jac_freq:.1f} Hz")
    log()

# ================================================================
# E) Epilepsy test
# ================================================================
log()
log("=" * 70)
log("AUFGABE E: Strong coupling (epilepsy regime)")
log("=" * 70)
log()

for s_epi in [5.0, 10.0]:
    # Check if always unstable
    evals_0, E0, I0 = jacobian_eigenvalues(s_epi, 0.0)
    evals_p, Ep, Ip = jacobian_eigenvalues(s_epi, 1.5)

    log(f"  s = {s_epi:.1f}:")
    log(f"    At P=0.0: Re(l)={np.max(evals_0.real):.4f}  "
        f"Im(l)={evals_0[np.argmax(evals_0.real)].imag:.4f}  "
        f"{'UNSTABLE' if np.max(evals_0.real) > 0 else 'stable'}")
    log(f"    At P=1.5: Re(l)={np.max(evals_p.real):.4f}  "
        f"Im(l)={evals_p[np.argmax(evals_p.real)].imag:.4f}  "
        f"{'UNSTABLE' if np.max(evals_p.real) > 0 else 'stable'}")

    # Simulate at P=1.5
    E_trace, I_trace = simulate_wc(s_epi, 1.5, t_max=500.0)
    E_last = E_trace[-2000:]
    amp = np.max(E_last) - np.min(E_last)
    mean_E = np.mean(E_last)
    E_centered = E_last - mean_E
    crossings = np.where(np.diff(np.sign(E_centered)))[0]
    if len(crossings) >= 2:
        period = 2 * np.mean(np.diff(crossings)) * 0.1
        freq = 1000.0 / period if period > 0 else 0
    else:
        freq = 0

    log(f"    Simulation (P=1.5): amp={amp:.4f}  freq={freq:.1f} Hz  "
        f"E range=[{np.min(E_last):.3f}, {np.max(E_last):.3f}]")

    if 3 <= freq <= 30:
        log(f"    -> In epileptic range (3-30 Hz)")
    elif freq > 30:
        log(f"    -> Above typical epileptic range")
    elif freq > 0:
        log(f"    -> Below typical epileptic range")
    log()

# ================================================================
log()
log("=" * 70)
log("SUMMARY")
log("=" * 70)
log()
if valid_fine:
    log(f"Peak P_crit = {pc_peak:.4f} at s = {s_peak:.2f}")
    log(f"Quantum peak at J_bridge/J ~ 1.9")
    log(f"Neural peak at s ~ {s_peak:.2f}")
    log()

    # Check hypotheses
    has_three_regimes = pc_peak > 0 and any(pc == 0 for s, pc, *_ in results_a)
    log(f"H1 (three regimes): {'CONFIRMED' if has_three_regimes else 'CHECK MANUALLY'}")

    bif_types = set(bt for _, _, _, _, _, bt in results_a if bt != "stable")
    log(f"H3 (Hopf everywhere): {'CONFIRMED' if bif_types == {'Hopf'} else 'MIXED: ' + str(bif_types)}")

log()
log(f"Completed: {clock.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"Results: {OUT_PATH}")
_outf.close()
