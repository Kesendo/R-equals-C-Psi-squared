#!/usr/bin/env python3
"""
SHADOW RETRODICTION: IBM Q80 late-time phase from detuning

The shadow framework predicts: at late times, only shadow-protected
modes survive. For a single qubit with detuning delta_omega under
Z-dephasing, the off-diagonal element evolves as:

    rho_01(t) = rho_01(0) * exp(-t/T2) * exp(-i * delta_omega * t)

The Lindblad model (no detuning) predicts phase = 0. The shadow
(detuning) produces a phase that grows linearly with time. At late
times (t >> T2), the exponential has killed the amplitude, but the
RATIO Im/Re = tan(delta_omega * t) reveals the shadow.

This script:
1. Extracts delta_omega from Q80's phase data (fit)
2. Retrodicts each data point's phase from the fitted detuning
3. Compares retrodiction with measured phase
"""

import numpy as np
import json
import os

REPO = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
DATA_DIR = os.path.join(REPO, "data", "ibm_shadow_march2026")

# Load Q80 hardware data
with open(os.path.join(DATA_DIR, "shadow_hardware_q80_20260309_181852.json")) as f:
    q80 = json.load(f)

print("=" * 70)
print("SHADOW RETRODICTION: IBM Torino Q80")
print("=" * 70)
print(f"Qubit: {q80['qubit']}, Backend: {q80['backend']}")
print(f"T1 = {q80['T1_us']:.1f} us, T2 = {q80['T2_us']:.1f} us, "
      f"T2* = {q80['T2_star_us']:.1f} us")
print()

# Extract late-time points (t/T2 >= 0.5)
points = q80["points"]
late = [p for p in points if p["t_over_T2"] >= 0.5]

print(f"Late-time points: {len(late)} (t/T2 >= 0.5)")
print()

# ============================================================
# METHOD 1: Fit delta_omega from phase vs time
# ============================================================
# Phase(t) = delta_omega * t (in radians)
# Measured phase in degrees -> convert to radians

times_us = np.array([p["delay_us"] for p in late])
phases_rad = np.array([p["rho01_phase_deg"] for p in late]) * np.pi / 180

# Linear fit: phase = delta_omega * t (no intercept, shadow starts at t=0)
# Use least squares: delta_omega = sum(t * phase) / sum(t^2)
delta_omega_fit = np.sum(times_us * phases_rad) / np.sum(times_us ** 2)
delta_omega_kHz = delta_omega_fit * 1000 / (2 * np.pi)  # convert rad/us to kHz

print("METHOD 1: Linear fit phase(t) = delta_omega * t")
print(f"  delta_omega = {delta_omega_fit:.6f} rad/us")
print(f"  delta_omega = {delta_omega_kHz:.2f} kHz")
print(f"  (= {abs(delta_omega_fit * 180/np.pi):.2f} deg/us)")
print()

# ============================================================
# METHOD 2: Fit with intercept (phase = delta_omega * t + phi_0)
# ============================================================
A = np.vstack([times_us, np.ones(len(times_us))]).T
result = np.linalg.lstsq(A, phases_rad, rcond=None)
delta_omega_fit2, phi_0 = result[0]
delta_omega_kHz2 = delta_omega_fit2 * 1000 / (2 * np.pi)

print("METHOD 2: Linear fit phase(t) = delta_omega * t + phi_0")
print(f"  delta_omega = {delta_omega_fit2:.6f} rad/us = {delta_omega_kHz2:.2f} kHz")
print(f"  phi_0 = {phi_0 * 180/np.pi:.2f} deg")
print()

# ============================================================
# RETRODICTION: compare predicted vs measured
# ============================================================
print("=" * 70)
print("RETRODICTION (Method 2)")
print("=" * 70)
print(f"  {'t/T2':>6}  {'t [us]':>8}  {'Meas [deg]':>10}  "
      f"{'Pred [deg]':>10}  {'Error [deg]':>10}")
print(f"  {'-' * 50}")

errors = []
for p in late:
    t = p["delay_us"]
    meas = p["rho01_phase_deg"]
    pred = (delta_omega_fit2 * t + phi_0) * 180 / np.pi
    err = meas - pred
    errors.append(err)
    print(f"  {p['t_over_T2']:>6.1f}  {t:>8.1f}  {meas:>10.1f}  "
          f"{pred:>10.1f}  {err:>10.1f}")

rms = np.sqrt(np.mean(np.array(errors) ** 2))
print(f"\n  RMS error: {rms:.1f} deg")
print(f"  Mean |error|: {np.mean(np.abs(errors)):.1f} deg")

# ============================================================
# FULL MODEL: rho_01(t) with detuning
# ============================================================
print()
print("=" * 70)
print("FULL SHADOW MODEL: rho_01(t) = rho_01(0) * exp(-t/T2) * exp(-i*dw*t)")
print("=" * 70)

T2 = q80["T2_us"]
rho01_0 = complex(points[0]["rho01_re"], points[0]["rho01_im"])

print(f"  rho_01(0) = {rho01_0.real:.4f} + {rho01_0.imag:.4f}i")
print(f"  T2 = {T2:.1f} us")
print(f"  delta_omega = {delta_omega_fit2:.6f} rad/us ({delta_omega_kHz2:.2f} kHz)")
print()
print(f"  {'t/T2':>6}  {'Re meas':>8}  {'Re pred':>8}  "
      f"{'Im meas':>8}  {'Im pred':>8}  {'|err|':>8}")
print(f"  {'-' * 55}")

full_errors = []
for p in points[1:]:  # skip t=0 (reference)
    t = p["delay_us"]
    # Shadow model: exp(-t/T2) * exp(-i*dw*t)
    # Convention: rho_01 phase = -dw*t, so positive measured phase means dw < 0
    # The fit gives positive dw from positive phases, so we negate for the model
    decay = np.exp(-t / T2)
    dw_physical = -delta_omega_fit2  # sign correction
    phase = np.exp(-1j * dw_physical * t)
    pred = rho01_0 * decay * phase
    meas = complex(p["rho01_re"], p["rho01_im"])
    err = abs(meas - pred)
    full_errors.append(err)

    # Lindblad (no detuning) for comparison
    lindblad = rho01_0.real * decay  # pure real decay

    print(f"  {p['t_over_T2']:>6.1f}  {meas.real:>8.4f}  {pred.real:>8.4f}  "
          f"{meas.imag:>8.4f}  {pred.imag:>8.4f}  {err:>8.4f}")

print(f"\n  Mean |error| (shadow model): {np.mean(full_errors):.4f}")

# Compare with Lindblad (no detuning)
lindblad_errors = []
for p in points[1:]:
    t = p["delay_us"]
    decay = np.exp(-t / T2)
    pred_lindblad = rho01_0.real * decay
    meas = complex(p["rho01_re"], p["rho01_im"])
    err = abs(meas - complex(pred_lindblad, 0))
    lindblad_errors.append(err)

print(f"  Mean |error| (Lindblad, no detuning): {np.mean(lindblad_errors):.4f}")
print(f"  Improvement factor: {np.mean(lindblad_errors)/np.mean(full_errors):.1f}x")

# ============================================================
# METHOD 3: Free fit of T_eff and dw to complex rho_01(t)
# ============================================================
print()
print("=" * 70)
print("FREE FIT: rho_01(t) = A * exp(-t/T_eff) * exp(i*phi(t))")
print("  Fitting T_eff and dw_phys simultaneously")
print("=" * 70)

from scipy.optimize import minimize

all_points = points[1:]  # skip t=0 reference
t_arr = np.array([p["delay_us"] for p in all_points])
re_arr = np.array([p["rho01_re"] for p in all_points])
im_arr = np.array([p["rho01_im"] for p in all_points])

def model_error(params):
    T_eff, dw = params
    if T_eff <= 0:
        return 1e10
    pred = rho01_0 * np.exp(-t_arr / T_eff) * np.exp(-1j * dw * t_arr)
    return np.sum((re_arr - pred.real)**2 + (im_arr - pred.imag)**2)

# Try multiple starting points
best = None
for T_start in [T2, q80["T2_star_us"], 15, 20]:
    for dw_start in [-0.01, -0.005, 0.005, 0.01]:
        res = minimize(model_error, [T_start, dw_start], method="Nelder-Mead")
        if best is None or res.fun < best.fun:
            best = res

T_eff_fit, dw_fit = best.x
dw_kHz_fit = dw_fit * 1000 / (2 * np.pi)

print(f"  T_eff = {T_eff_fit:.2f} us (T2 = {T2:.1f}, T2* = {q80['T2_star_us']:.1f})")
print(f"  dw = {dw_fit:.6f} rad/us = {dw_kHz_fit:.2f} kHz")
print()

# Full retrodiction with free fit
print(f"  {'t/T2':>6}  {'Re meas':>8}  {'Re pred':>8}  "
      f"{'Im meas':>8}  {'Im pred':>8}  {'|err|':>8}")
print(f"  {'-' * 55}")

free_errors = []
for p in all_points:
    t = p["delay_us"]
    pred = rho01_0 * np.exp(-t / T_eff_fit) * np.exp(-1j * dw_fit * t)
    meas = complex(p["rho01_re"], p["rho01_im"])
    err = abs(meas - pred)
    free_errors.append(err)
    print(f"  {p['t_over_T2']:>6.1f}  {meas.real:>8.4f}  {pred.real:>8.4f}  "
          f"{meas.imag:>8.4f}  {pred.imag:>8.4f}  {err:>8.4f}")

print(f"\n  Mean |error| (free fit): {np.mean(free_errors):.4f}")
print(f"  Mean |error| (Lindblad): {np.mean(lindblad_errors):.4f}")
print(f"  Improvement: {np.mean(lindblad_errors)/np.mean(free_errors):.1f}x")

# ============================================================
# VERDICT
# ============================================================
print()
print("=" * 70)
print("VERDICT")
print("=" * 70)
improvement = np.mean(lindblad_errors) / np.mean(full_errors)
print(f"  The shadow model (Lindblad + detuning = {delta_omega_kHz2:.1f} kHz)")
print(f"  reduces the mean prediction error by {improvement:.1f}x compared to")
print(f"  standard Lindblad (no detuning).")
print()
if improvement > 2:
    print("  The shadow is real. The detuning is the Hamiltonian's fingerprint")
    print("  in the shadow of the dephasing light. It becomes visible when the")
    print("  illuminated modes (exponential decay) have died and only the")
    print("  shadow-protected phase rotation remains.")
print()
print("=" * 70)
