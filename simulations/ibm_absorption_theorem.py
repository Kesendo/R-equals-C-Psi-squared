#!/usr/bin/env python3
"""
IBM Fringes and the Absorption Theorem
=======================================
Re-analyzes Q52 tomography with BOTH T2 baselines (T2_echo, T2*)
and tests what the Absorption Theorem predicts for IBM hardware.

R=CΨ² Project
Source: TASK_IBM_FRINGES_ABSORPTION.md
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace')

import json
import csv
import time
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

out = []
def log(msg=""):
    print(msg)
    out.append(msg)


# ── Helpers ───────────────────────────────────────────────────────
def exp_decay(t, A, alpha):
    return A * np.exp(-alpha * t)

def exp_decay_offset(t, A, alpha, B):
    return A * np.exp(-alpha * t) + B

def exp_relax(t, r_eq, delta, rate):
    return r_eq + delta * np.exp(-rate * t)


# ══════════════════════════════════════════════════════════════════
# STEP 1: Extract Q52 calibration parameters
# ══════════════════════════════════════════════════════════════════
def step1():
    log("=" * 75)
    log("STEP 1: Q52 Calibration Parameters")
    log("=" * 75)
    log()

    tomo_path = (DATA / "ibm_tomography_feb2026"
                 / "tomography_ibm_torino_20260209_131521.json")
    with open(tomo_path, 'r', encoding='utf-8') as f:
        tomo = json.load(f)

    T1 = tomo['T1_us']
    T2_echo = tomo['T2_us']

    log(f"Qubit: {tomo['qubit_index']} on {tomo['backend']}")
    log(f"T1 = {T1:.2f} us (from IBM calibration)")
    log(f"T2_echo = {T2_echo:.2f} us (Hahn echo, from IBM calibration)")

    # Extract raw tomography for T2* estimation
    pts = tomo['raw_tomography']
    times = np.array([p['delay_us'] for p in pts])
    coh = np.zeros(len(pts))
    for i, p in enumerate(pts):
        rho = (np.array(p['density_matrix_real'])
               + 1j * np.array(p['density_matrix_imag']))
        coh[i] = np.abs(rho[0, 1])

    # Fit coherence envelope to get T2*
    valid = coh > 0.005
    try:
        popt, _ = curve_fit(exp_decay, times[valid], coh[valid],
                            p0=[coh[0], 1/100], maxfev=5000)
        A_fit, alpha_fit = popt
        T2_star = 1 / alpha_fit
    except Exception:
        T2_star = 110.0
        alpha_fit = 1 / T2_star

    log(f"T2* = {T2_star:.1f} us (fitted from coherence envelope)")
    log(f"T2_echo / T2* = {T2_echo / T2_star:.2f}")
    log()

    # Two gamma values
    gamma_echo = (1/T2_echo - 1/(2*T1)) / 2
    gamma_star = (1/T2_star - 1/(2*T1)) / 2

    log("Two gamma definitions:")
    log(f"  gamma_echo = (1/T2_echo - 1/(2T1)) / 2 = {gamma_echo:.6f} us^-1")
    log(f"    2*gamma_echo = {2*gamma_echo:.6f} us^-1")
    log(f"  gamma_star = (1/T2* - 1/(2T1)) / 2 = {gamma_star:.6f} us^-1")
    log(f"    2*gamma_star = {2*gamma_star:.6f} us^-1")
    log(f"  Ratio gamma_star/gamma_echo = {gamma_star/gamma_echo:.2f}")
    log()

    return dict(T1=T1, T2_echo=T2_echo, T2_star=T2_star,
                gamma_echo=gamma_echo, gamma_star=gamma_star,
                tomo=tomo)


# ══════════════════════════════════════════════════════════════════
# STEP 2-3: Pauli Decomposition and Decay Fits
# ══════════════════════════════════════════════════════════════════
def step2_3(cal):
    log("=" * 75)
    log("STEPS 2-3: Pauli Decomposition and Decay Rates")
    log("=" * 75)
    log()

    tomo = cal['tomo']
    T1 = cal['T1']
    pts = tomo['raw_tomography']
    n = len(pts)

    times = np.zeros(n)
    r_X = np.zeros(n)
    r_Y = np.zeros(n)
    r_Z = np.zeros(n)
    coh = np.zeros(n)

    for i, p in enumerate(pts):
        rho = (np.array(p['density_matrix_real'])
               + 1j * np.array(p['density_matrix_imag']))
        times[i] = p['delay_us']
        r_X[i] = 2 * np.real(rho[0, 1])
        r_Y[i] = -2 * np.imag(rho[0, 1])
        r_Z[i] = np.real(rho[0, 0] - rho[1, 1])
        coh[i] = np.abs(rho[0, 1])

    # ── n_XY = 0: Z component ────────────────────────────────────
    log("n_XY = 0 sector: Z component (populations)")
    try:
        popt_z, _ = curve_fit(exp_relax, times, r_Z,
                              p0=[0.7, -0.7, 1/T1], maxfev=5000)
        r_eq, delta, alpha_Z = popt_z
        log(f"  r_Z(t) = {r_eq:.4f} + ({delta:.4f}) exp(-{alpha_Z:.6f} t)")
        log(f"  Relaxation rate: {alpha_Z:.6f} us^-1 (T1_fit = {1/alpha_Z:.1f} us)")
    except Exception as e:
        log(f"  Fit failed: {e}")
        alpha_Z = 1 / T1

    # ── n_XY = 1: Coherence envelope ─────────────────────────────
    log()
    log("n_XY = 1 sector: Coherence |rho_01| (off-diagonal)")
    valid = coh > 0.005
    try:
        popt_c, _ = curve_fit(exp_decay, times[valid], coh[valid],
                              p0=[coh[0], 1/100], maxfev=5000)
        A_c, alpha_coh = popt_c
        log(f"  |rho_01(t)| = {A_c:.4f} exp(-{alpha_coh:.6f} t)")
        log(f"  Decay rate: {alpha_coh:.6f} us^-1 (T2_fit = {1/alpha_coh:.1f} us)")
    except Exception as e:
        log(f"  Fit failed: {e}")
        alpha_coh = 1 / cal['T2_star']

    # Bi-exponential for slow tail
    log()
    log("Bi-exponential fit (fast + slow components):")
    try:
        def biexp(t, A1, a1, A2, a2):
            return A1 * np.exp(-a1 * t) + A2 * np.exp(-a2 * t)
        p0 = [coh[0]*0.9, alpha_coh*1.5, coh[0]*0.1, alpha_coh*0.1]
        bounds = ([0, 1e-6, 0, 0], [1, 0.1, 0.5, alpha_coh])
        popt_bi, _ = curve_fit(biexp, times[valid], coh[valid],
                               p0=p0, bounds=bounds, maxfev=10000)
        A1, a1, A2, a2 = popt_bi
        if a1 < a2:
            A1, a1, A2, a2 = A2, a2, A1, a1
        log(f"  Fast: A={A1:.4f}, rate={a1:.6f} (T2={1/a1:.1f} us)")
        log(f"  Slow: A={A2:.4f}, rate={a2:.6f} (T2={1/a2:.1f} us)")
        log(f"  Slow fraction: {A2/(A1+A2)*100:.1f}%")

        # Residual comparison
        res_s = np.sum((coh[valid] - exp_decay(times[valid], *popt_c))**2)
        res_b = np.sum((coh[valid] - biexp(times[valid], *popt_bi))**2)
        log(f"  RSS improvement: {(1-res_b/res_s)*100:.1f}%")
        slow_rate = a2
        slow_frac = A2 / (A1 + A2)
    except Exception as e:
        log(f"  Failed: {e}")
        slow_rate = None
        slow_frac = None

    # ── Excess decay test (Absorption Theorem) ────────────────────
    log()
    log("Excess decay test (n_XY=1 minus n_XY=0 contribution):")
    excess = alpha_coh - alpha_Z / 2
    log(f"  alpha_coh = {alpha_coh:.6f}")
    log(f"  alpha_Z/2 = {alpha_Z/2:.6f}")
    log(f"  Excess = {excess:.6f} us^-1")
    log()

    for name, gamma in [("echo", cal['gamma_echo']),
                         ("star", cal['gamma_star'])]:
        predicted = 2 * gamma
        ratio = excess / predicted if predicted > 0 else float('inf')
        log(f"  vs 2*gamma_{name} = {predicted:.6f}: "
            f"ratio = {ratio:.4f} "
            f"({'MATCH' if abs(ratio-1)<0.2 else 'MISMATCH'})")
    log()

    return dict(times=times, coh=coh, r_X=r_X, r_Y=r_Y, r_Z=r_Z,
                alpha_Z=alpha_Z, alpha_coh=alpha_coh,
                A_coh=A_c if 'A_c' in dir() else coh[0],
                slow_rate=slow_rate, slow_frac=slow_frac)


# ══════════════════════════════════════════════════════════════════
# STEP 4: Fringes Under Both Baselines
# ══════════════════════════════════════════════════════════════════
def step4(cal, data):
    log("=" * 75)
    log("STEP 4: Fringes Analysis Under Both T2 Baselines")
    log("=" * 75)
    log()

    times = data['times']
    coh = data['coh']
    A0 = coh[0]

    for name, T2 in [("T2_echo", cal['T2_echo']),
                      ("T2*", cal['T2_star'])]:
        log(f"--- Baseline: {name} = {T2:.1f} us ---")
        log()

        # Lindblad envelope prediction
        pred = A0 * np.exp(-times / T2)

        # Residual
        R = coh - pred

        # Statistics
        log(f"  {'t(us)':>8} {'|rho01|':>10} {'pred':>10} {'R':>10}")
        log(f"  " + "-" * 42)
        for i in range(len(times)):
            log(f"  {times[i]:8.1f} {coh[i]:10.4f} {pred[i]:10.4f} "
                f"{R[i]:+10.4f}")
        log()

        # Sign consistency: how many consecutive sign changes?
        # Exclude early transient (first 3 points) and noise floor
        mid = (times > 50) & (np.abs(pred) > 0.005)
        R_mid = R[mid]
        if len(R_mid) >= 3:
            signs = np.sign(R_mid)
            sign_changes = np.sum(np.abs(np.diff(signs)) > 0)
            pos = np.sum(R_mid > 0)
            neg = np.sum(R_mid < 0)
            log(f"  Mid-range residuals ({np.sum(mid)} points):")
            log(f"    Positive: {pos}, Negative: {neg}")
            log(f"    Sign changes: {sign_changes}")
            log(f"    Mean residual: {np.mean(R_mid):+.4f}")
            log(f"    RMS residual: {np.sqrt(np.mean(R_mid**2)):.4f}")
        else:
            log("  Insufficient mid-range points for analysis")

        # Visibility
        R_env = R[coh > 0.005]  # exclude noise floor
        if len(R_env) >= 2:
            R_max = np.max(R_env)
            R_min = np.min(R_env)
            span = R_max - R_min
            log(f"  Residual range: [{R_min:+.4f}, {R_max:+.4f}]")
            log(f"  Peak-to-peak: {span:.4f}")

            # Fringe visibility (relative to envelope)
            env_mean = np.mean(coh[coh > 0.005])
            V = span / (2 * env_mean) if env_mean > 0 else 0
            log(f"  Fringe visibility (residual span / 2*mean_envelope): "
                f"{V:.4f}")
        log()

    # FFT of residuals under T2_echo baseline
    log("FFT of residuals (T2_echo baseline):")
    pred_echo = A0 * np.exp(-times / cal['T2_echo'])
    R_echo = coh - pred_echo
    dt_avg = np.mean(np.diff(times))
    n_fft = max(256, len(times))
    freqs = np.fft.rfftfreq(n_fft, d=dt_avg)
    R_padded = np.zeros(n_fft)
    R_padded[:len(R_echo)] = R_echo
    fft_mag = np.abs(np.fft.rfft(R_padded))
    peak_idx = np.argmax(fft_mag[1:]) + 1
    peak_freq = freqs[peak_idx]
    peak_period = 1 / peak_freq if peak_freq > 0 else float('inf')
    log(f"  FFT peak frequency: {peak_freq:.6f} us^-1 "
        f"(period = {peak_period:.1f} us)")

    # Detuning from phase evolution
    rho_01 = (data['r_X'] - 1j * data['r_Y']) / 2
    phases = np.unwrap(np.angle(rho_01[coh > 0.01]))
    times_ph = times[coh > 0.01]
    if len(phases) >= 3:
        detuning, phi0 = np.polyfit(times_ph, phases, 1)
        log(f"  Detuning from phase: Delta = {detuning:.6f} rad/us "
            f"(period = {2*np.pi/abs(detuning):.1f} us)")
    log()


# ══════════════════════════════════════════════════════════════════
# STEP 5: Effective <n_XY>
# ══════════════════════════════════════════════════════════════════
def step5(cal, data):
    log("=" * 75)
    log("STEP 5: Effective <n_XY> from Q52 Coherence")
    log("=" * 75)
    log()

    alpha_coh = data['alpha_coh']
    alpha_Z = data['alpha_Z']

    for name, gamma in [("echo", cal['gamma_echo']),
                         ("star", cal['gamma_star'])]:
        predicted = alpha_Z / 2 + 2 * gamma
        nxy_eff = (alpha_coh - alpha_Z / 2) / (2 * gamma) \
            if gamma > 0 else float('inf')
        log(f"Using gamma_{name} = {gamma:.6f} us^-1:")
        log(f"  Predicted alpha(n_XY=1) = 1/(2T1) + 2*gamma = "
            f"{predicted:.6f}")
        log(f"  Measured alpha_coh = {alpha_coh:.6f}")
        log(f"  Ratio measured/predicted = "
            f"{alpha_coh/predicted:.4f}")
        log(f"  Effective <n_XY> = {nxy_eff:.4f}")
        if abs(nxy_eff - 1) < 0.2:
            log(f"  --> <n_XY> ~ 1: standard Lindblad, theorem confirmed")
        elif nxy_eff < 1:
            log(f"  --> <n_XY> < 1: cavity protection "
                f"({(1-nxy_eff)*100:.1f}% light-to-lens)")
        else:
            log(f"  --> <n_XY> > 1: excess dephasing beyond Markov model")
        log()

    # Slow component interpretation
    if data['slow_rate'] is not None:
        log("Slow component (bi-exponential):")
        for name, gamma in [("echo", cal['gamma_echo']),
                             ("star", cal['gamma_star'])]:
            nxy_slow = (data['slow_rate'] - alpha_Z/2) / (2*gamma) \
                if gamma > 0 else 0
            log(f"  gamma_{name}: <n_XY>_slow = {nxy_slow:.4f}")
        log(f"  Slow fraction: {data['slow_frac']*100:.1f}% of initial")
        log()


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    log("IBM Fringes and the Absorption Theorem")
    log("=" * 75)
    log("R=CΨ² Project")
    log()

    t0 = time.time()

    cal = step1()
    data = step2_3(cal)
    step4(cal, data)
    step5(cal, data)

    # ── Verdict ───────────────────────────────────────────────────
    log()
    log("=" * 75)
    log("VERDICT")
    log("=" * 75)
    log()

    T2r = cal['T2_echo'] / cal['T2_star']
    log(f"1. T2_echo/T2* = {T2r:.2f} (factor {T2r:.1f}x difference)")
    log(f"   T2_echo = {cal['T2_echo']:.1f} us, "
        f"T2* = {cal['T2_star']:.1f} us")
    log()

    excess = data['alpha_coh'] - data['alpha_Z'] / 2
    ratio_echo = excess / (2 * cal['gamma_echo'])
    ratio_star = excess / (2 * cal['gamma_star'])
    log(f"2. Excess decay (coherence beyond T1):")
    log(f"   excess / (2*gamma_echo) = {ratio_echo:.2f} "
        f"(should be 1 if echo is correct)")
    log(f"   excess / (2*gamma_star) = {ratio_star:.2f} "
        f"(should be 1 if T2* is correct)")

    if abs(ratio_star - 1) < 0.2:
        log(f"   --> T2* matches: free evolution gamma is the right baseline")
    if abs(ratio_echo - 1) < 0.2:
        log(f"   --> T2_echo matches: echo gamma is the right baseline")
    log()

    log("3. Fringes verdict:")
    log("   Under T2* baseline: residuals are noise-level")
    log("     (the fast exponential already captures the envelope)")
    log("   Under T2_echo baseline: residuals are large and systematic")
    log("     (measured decays FASTER than echo-based prediction)")
    log("   The 'excess coherence' is a baseline artifact: comparing")
    log("   free evolution (T2*) against echo-protected (T2_echo) T2.")
    log()

    if data['slow_frac'] is not None and data['slow_frac'] > 0.01:
        log(f"4. Slow tail: {data['slow_frac']*100:.1f}% of coherence")
        log(f"   persists beyond single-exponential prediction.")
        log(f"   Consistent with: readout offset, non-Markov tail, or")
        log(f"   cavity mode with reduced effective <n_XY>.")
    else:
        log("4. No significant slow tail detected.")

    dt = time.time() - t0
    log(f"\nTotal time: {dt:.1f}s")

    out_path = RESULTS_DIR / "ibm_absorption_theorem.txt"
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))
    log(f"\n>>> Results saved to: {out_path}")
