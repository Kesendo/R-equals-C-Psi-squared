#!/usr/bin/env python3
"""
EXPLORATORY Q80 IN-SAMPLE DETUNING-COMPATIBLE FIT

This script is Q80-only. It compares two descriptive models with the same
Q80 record used to fit them.
It is not a held-out prediction and does not identify a Q52 mechanism.

The fixed-Hahn-T2 descriptive curve is

    rho_fit(t) = A0 * exp(-t/T2) * exp(i * (phi_0 + m*t)),

with A0 = |rho_01(0)|. Under the Hamiltonian convention
rho_01(t) proportional to exp(-i*delta_omega*t), delta_omega = -m.

This script:
1. fits a linear Q80 phase model and an intercept-only comparator on the same
   eight late rows, then scores both fixed-Hahn-T2 curves on the same nine
   nonzero-time rows;
2. fits a free (T_eff, delta_omega) curve and a nested no-detuning envelope
   comparator on the same nine rows, then scores both there.

The sampled delays lie on one time grid. Frequencies are therefore reported
as near-zero representatives modulo that grid's alias spacing; no physical
prior in this script selects an absolute detuning.

The JSON field named T2_star_us is a scheduling proxy T2_echo/2.5, not a
measured Ramsey T2*. It is retained only as an optimizer start value.
"""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = REPO / "data" / "ibm_shadow_march2026"
ALIAS_SCOPE = (
    "All quoted frequencies are near-zero representatives of alias classes; "
    "the sampled complex trajectory does not identify an absolute detuning."
)
Q80_SOURCE_IDENTITY_ERROR = (
    "standalone Q80 record does not match the combined March Q80 row"
)

# Load Q80 hardware data
with (DATA_DIR / "shadow_hardware_q80_20260309_181852.json").open(
    encoding="utf-8"
) as f:
    q80 = json.load(f)
with (DATA_DIR / "shadow_hardware_combined_20260309_181852.json").open(
    encoding="utf-8"
) as f:
    combined = json.load(f)

combined_q80_rows = [
    row for row in combined["qubit_results"] if int(row["qubit"]) == 80
]
if len(combined_q80_rows) != 1:
    raise ValueError("combined March record does not contain exactly one Q80 row")
combined_q80 = combined_q80_rows[0]
if (
    combined.get("experiment") != "shadow_hunt_hardware"
    or combined.get("mode") != "hardware"
    or combined.get("backend") != "ibm_torino"
    or combined.get("timestamp") != "20260309_181852"
    or int(q80["qubit"]) != 80
    or q80.get("backend") != "ibm_torino"
    or q80.get("timestamp") != "20260309_181852"
    or combined["backend"] != q80["backend"]
    or combined["timestamp"] != q80["timestamp"]
    or combined_q80["points"] != q80["points"]
    or combined_q80["verdict"] != q80["verdict"]
):
    raise ValueError(Q80_SOURCE_IDENTITY_ERROR)
for field in ("T1_us", "T2_us", "T2_star_us"):
    if float(combined_q80["verdict"][field]) != float(q80[field]):
        raise ValueError(Q80_SOURCE_IDENTITY_ERROR)

t2_echo_us = float(q80["T2_us"])
schedule_proxy_factor = float(combined["T2_star_factor"])
schedule_proxy_us = t2_echo_us / schedule_proxy_factor
if not np.isclose(
    float(q80["T2_star_us"]), schedule_proxy_us, rtol=0.0, atol=1e-12
):
    raise ValueError("Q80 stored scheduling proxy is not T2_echo/T2_star_factor")

print("=" * 70)
print("Q80 IN-SAMPLE DETUNING-COMPATIBLE FIT: IBM Torino Q80")
print("=" * 70)
print(f"Qubit: {q80['qubit']}, Backend: {q80['backend']}")
print(f"T1 = {q80['T1_us']:.1f} us, T2_echo = {t2_echo_us:.1f} us")
print(
    f"Scheduling proxy T2_echo/{schedule_proxy_factor:.1f} = "
    f"{schedule_proxy_us:.1f} us (not a measured Ramsey T2*)"
)
print()

# Extract late-time points (t/T2 >= 0.5)
points = q80["points"]
for point in points:
    delay_us = float(point["delay_us"])
    if not np.isclose(
        float(point["t_over_T2"]), delay_us / t2_echo_us, rtol=0.0, atol=1e-15
    ):
        raise ValueError("Q80 t_over_T2 disagrees with delay_us/T2_echo")
    if not np.isclose(
        float(point["t_over_T2star"]),
        delay_us / schedule_proxy_us,
        rtol=0.0,
        atol=1e-15,
    ):
        raise ValueError(
            "Q80 t_over_T2star disagrees with delay_us/scheduling_proxy"
        )
late = [p for p in points if float(p["delay_us"]) / t2_echo_us >= 0.5]
positive_delays = np.array(
    [float(p["delay_us"]) for p in points if float(p["delay_us"]) > 0.0]
)
sample_grid_us = float(np.min(positive_delays))
sample_indices = np.rint(positive_delays / sample_grid_us)
if not np.allclose(
    positive_delays,
    sample_indices * sample_grid_us,
    rtol=0.0,
    atol=1e-12,
):
    raise ValueError("Q80 delays do not lie on the inferred sampling grid")
alias_spacing_khz = 1000.0 / sample_grid_us


def near_zero_alias_representative(frequency_khz):
    return (
        (frequency_khz + alias_spacing_khz / 2.0) % alias_spacing_khz
        - alias_spacing_khz / 2.0
    )

print(f"Late-time points: {len(late)} (t/T2 >= 0.5)")
print(
    f"Sampling grid: {sample_grid_us:.6f} us; frequency alias spacing: "
    f"{alias_spacing_khz:.6f} kHz"
)
print(ALIAS_SCOPE)
print()

times_us = np.array([p["delay_us"] for p in late])
measured_late = np.array(
    [complex(p["rho01_re"], p["rho01_im"]) for p in late], dtype=complex
)
phases_rad = np.unwrap(np.angle(measured_late))
stored_phases_rad = np.deg2rad(
    np.array([p["rho01_phase_deg"] for p in late], dtype=float)
)
phase_coordinate_error = np.angle(np.exp(1j * (stored_phases_rad - phases_rad)))
if not np.allclose(phase_coordinate_error, 0.0, rtol=0.0, atol=1e-15):
    raise ValueError("Q80 stored phases disagree with the raw complex samples")

# ============================================================
# FIXED-T2 PHASE-LINE FIT
# ============================================================
A = np.vstack([times_us, np.ones(len(times_us))]).T
result = np.linalg.lstsq(A, phases_rad, rcond=None)
phase_slope_fit, phi_0 = result[0]
phase_slope_khz = near_zero_alias_representative(
    phase_slope_fit * 1000 / (2 * np.pi)
)
detuning_fit = -phase_slope_fit
detuning_khz = near_zero_alias_representative(
    detuning_fit * 1000 / (2 * np.pi)
)
intercept_only_phase = float(np.mean(phases_rad))

print("FIXED-T2 PHASE LINE: phase(t) = phi_0 + m*t")
print(
    f"  phase slope representative m = {phase_slope_fit:.6f} rad/us "
    f"= {phase_slope_khz:.2f} kHz"
)
print(f"  detuning representative delta_omega = -m = {detuning_khz:.2f} kHz")
print(f"  phi_0 = {phi_0 * 180/np.pi:.2f} deg")
print(f"  intercept-only phase = {intercept_only_phase * 180/np.pi:.2f} deg")
print()

# ============================================================
# IN-SAMPLE PHASE FIT: compare fitted vs measured phase
# ============================================================
print("=" * 70)
print("IN-SAMPLE PHASE-LINE FIT")
print("=" * 70)
print(f"  {'t/T2':>6}  {'t [us]':>8}  {'Meas [deg]':>10}  "
      f"{'Pred [deg]':>10}  {'Error [deg]':>10}")
print(f"  {'-' * 50}")

errors = []
for p in late:
    t = p["delay_us"]
    meas = p["rho01_phase_deg"]
    pred = (phase_slope_fit * t + phi_0) * 180 / np.pi
    err = meas - pred
    errors.append(err)
    print(f"  {p['t_over_T2']:>6.1f}  {t:>8.1f}  {meas:>10.1f}  "
          f"{pred:>10.1f}  {err:>10.1f}")

rms = np.sqrt(np.mean(np.array(errors) ** 2))
print(f"\n  RMS error: {rms:.1f} deg")
print(f"  Mean |error|: {np.mean(np.abs(errors)):.1f} deg")

# ============================================================
# FIXED-T2 COMPLEX SCORE
# ============================================================
print()
print("=" * 70)
print("FIXED-T2 COMPLEX SCORE: |rho_01(0)| * exp(-t/T2) with fitted phase")
print("=" * 70)

T2 = t2_echo_us
rho01_0 = complex(points[0]["rho01_re"], points[0]["rho01_im"])

print(f"  rho_01(0) = {rho01_0.real:.4f} + {rho01_0.imag:.4f}i")
print(f"  T2 = {T2:.1f} us")
print(
    f"  phase slope representative m = {phase_slope_fit:.6f} rad/us "
    f"({phase_slope_khz:.2f} kHz)"
)
print(f"  detuning representative delta_omega = -m = {detuning_khz:.2f} kHz")
print()
print(f"  {'t/T2':>6}  {'Re meas':>8}  {'Re pred':>8}  "
      f"{'Im meas':>8}  {'Im pred':>8}  {'|err|':>8}")
print(f"  {'-' * 55}")

fixed_t2_errors = []
intercept_only_errors = []
for p in points[1:]:  # skip t=0 (reference)
    t = p["delay_us"]
    decay = np.exp(-t / T2)
    fitted_phase = phase_slope_fit * t + phi_0
    pred = abs(rho01_0) * decay * np.exp(1j * fitted_phase)
    meas = complex(p["rho01_re"], p["rho01_im"])
    err = abs(meas - pred)
    fixed_t2_errors.append(err)
    intercept_pred = (
        abs(rho01_0) * decay * np.exp(1j * intercept_only_phase)
    )
    intercept_only_errors.append(abs(meas - intercept_pred))

    print(f"  {p['t_over_T2']:>6.1f}  {meas.real:>8.4f}  {pred.real:>8.4f}  "
          f"{meas.imag:>8.4f}  {pred.imag:>8.4f}  {err:>8.4f}")

fixed_error = float(np.mean(fixed_t2_errors))
intercept_error = float(np.mean(intercept_only_errors))
fixed_ratio = intercept_error / fixed_error
print(f"\n  Mean |error| (fixed-T2 phase-line): {fixed_error:.4f}")
print(f"  Mean |error| (fixed-T2 intercept-only): {intercept_error:.4f}")
print(f"  In-sample comparator/error ratio: {fixed_ratio:.1f}x")

# ============================================================
# METHOD 3: Free fit of T_eff and dw to complex rho_01(t)
# ============================================================
print()
print("=" * 70)
print(
    "FREE FIT: rho_01(t) = rho_01(0) * exp(-t/T_eff) "
    "* exp(-i*delta_omega*t)"
)
print("  Fitting T_eff and delta_omega simultaneously")
print("=" * 70)

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
for T_start in [T2, schedule_proxy_us, 15, 20]:
    for dw_start in [-0.01, -0.005, 0.005, 0.01]:
        res = minimize(model_error, [T_start, dw_start], method="Nelder-Mead")
        if best is None or res.fun < best.fun:
            best = res

T_eff_fit, dw_fit = best.x
dw_kHz_fit = near_zero_alias_representative(
    dw_fit * 1000 / (2 * np.pi)
)
dw_rep_fit = dw_kHz_fit * 2.0 * np.pi / 1000.0

print(
    f"  T_eff = {T_eff_fit:.2f} us (T2_echo = {T2:.1f}, "
    f"scheduling proxy = {schedule_proxy_us:.1f})"
)
print(
    f"  near-zero dw representative = {dw_rep_fit:.6f} rad/us "
    f"= {dw_kHz_fit:.2f} kHz"
)
print()


def envelope_only_error(T_eff):
    if T_eff <= 0:
        return 1e10
    pred = rho01_0 * np.exp(-t_arr / T_eff)
    return np.sum((re_arr - pred.real)**2 + (im_arr - pred.imag)**2)


envelope = minimize_scalar(
    envelope_only_error, bounds=(1.0, 200.0), method="bounded"
)
envelope_T_eff = float(envelope.x)

# Same-record complex score with the free fit
print(f"  {'t/T2':>6}  {'Re meas':>8}  {'Re pred':>8}  "
      f"{'Im meas':>8}  {'Im pred':>8}  {'|err|':>8}")
print(f"  {'-' * 55}")

free_errors = []
envelope_errors = []
for p in all_points:
    t = p["delay_us"]
    pred = rho01_0 * np.exp(-t / T_eff_fit) * np.exp(-1j * dw_fit * t)
    meas = complex(p["rho01_re"], p["rho01_im"])
    err = abs(meas - pred)
    free_errors.append(err)
    envelope_pred = rho01_0 * np.exp(-t / envelope_T_eff)
    envelope_errors.append(abs(meas - envelope_pred))
    print(f"  {p['t_over_T2']:>6.1f}  {meas.real:>8.4f}  {pred.real:>8.4f}  "
          f"{meas.imag:>8.4f}  {pred.imag:>8.4f}  {err:>8.4f}")

free_error = float(np.mean(free_errors))
envelope_error = float(np.mean(envelope_errors))
free_ratio = envelope_error / free_error
print(f"\n  Mean |error| (free complex fit): {free_error:.4f}")
print(
    "  Mean |error| (fitted no-detuning envelope): "
    f"{envelope_error:.4f} (T_eff = {envelope_T_eff:.2f} us)"
)
print(f"  In-sample comparator/error ratio: {free_ratio:.1f}x")

# ============================================================
# SCOPE VERDICT
# ============================================================
print()
print("=" * 70)
print("EXPLORATORY IN-SAMPLE VERDICT")
print("=" * 70)
print(
    f"  Fixed-T2 phase line ({detuning_khz:.1f} kHz near-zero detuning "
    "representative): "
    f"{fixed_ratio:.1f}x smaller mean error than the intercept-only comparator."
)
print(
    f"  Free (T_eff, delta_omega) fit: {free_ratio:.1f}x smaller mean error "
    "than its fitted no-detuning envelope comparator."
)
print(
    "  Both comparisons reuse the Q80 record; neither is a Q52 mechanism "
    "fit or a held-out prediction."
)
print(
    f"RESULT phase_slope_rep_khz={phase_slope_khz:.6f} "
    f"detuning_rep_khz={detuning_khz:.6f} "
    f"phase_slope_rep_rad_per_us={phase_slope_fit:.6f} "
    f"detuning_rep_rad_per_us={detuning_fit:.6f} "
    f"phase_intercept_rad={phi_0:.6f} "
    f"intercept_only_phase_rad={intercept_only_phase:.6f} "
    f"phase_rms_deg={rms:.6f} "
    f"phase_mae_deg={np.mean(np.abs(errors)):.6f} "
    f"sample_grid_us={sample_grid_us:.6f} "
    f"alias_spacing_khz={alias_spacing_khz:.6f} "
    f"t2_echo_us={T2:.6f} "
    f"schedule_proxy_factor={schedule_proxy_factor:.6f} "
    f"schedule_proxy_us={schedule_proxy_us:.6f} "
    f"fixed_error={fixed_error:.6f} "
    f"intercept_error={intercept_error:.6f} "
    f"fixed_ratio={fixed_ratio:.6f} "
    f"free_T_eff_us={T_eff_fit:.6f} "
    f"free_detuning_rep_khz={dw_kHz_fit:.6f} "
    f"free_detuning_rep_rad_per_us={dw_rep_fit:.6f} "
    f"free_error={free_error:.6f} "
    f"envelope_T_eff_us={envelope_T_eff:.6f} "
    f"envelope_error={envelope_error:.6f} "
    f"free_ratio={free_ratio:.6f}"
)
print()
print("=" * 70)
