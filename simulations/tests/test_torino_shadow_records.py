"""Truth gates for the Torino shadow records.

The records: the Q52 tomography of 2026-02-09 (data/ibm_tomography_feb2026/), the March 9 shadow run on
Q80 and Q102 with its simulator null (data/ibm_shadow_march2026/), and the Ramsey fits of March 12 and 18
(data/ibm_run3_march2026/). Recomputed here from the raw JSON: the numbers that
experiments/FIXED_POINT_SHADOW.md and experiments/RESIDUAL_ANALYSIS.md rest on; the two-fits ratio, the
bi-exponential slow tail and the early rotation of experiments/IBM_ABSORPTION_THEOREM.md; the crossing
entry's t*, 9.1%, CPsi(0) and late purity and the absorption entry's 1.03 and 2.8% in the two Torino registry
entries; the figure readings visualizations/README.md gives for the two crossing figures; and the lines of
simulations/results/cockpit_validation.txt.

How a value is compared, and why:
  * An exact identity is compared exactly: the stored axis delay_over_T2 = t/T2_echo, Psi = 2|rho01| and
    CPsi = C*Psi, the uniform delay grid, the scheduling proxy T2_star_us = T2_us/2.5, the raw count
    ratios behind the stored Pauli readings, the two-fits absorption ratio with one T1 on both sides.
  * A fit or optimizer output is compared at the rounding the documents print. The solvers converge at
    least two decades below the last printed digit, so what can fail is the printed digit itself.
  * A Monte Carlo value is compared with its own standard error: the cited value must lie within four
    standard errors of the estimate, plus half a unit of the cited value's last digit.

No producer is imported. The cockpit producer runs in a subprocess with its output redirected to a
scratch file; that output must equal the committed results file, line endings normalized (Python writes
the platform's, git checks out the configured ones).

Run with explicit paths only (python -m pytest simulations/tests/test_torino_shadow_records.py).
"""
from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy import stats
from scipy.optimize import brentq, curve_fit, least_squares, minimize, minimize_scalar

ROOT = Path(os.environ.get("RCPSI_REPO_ROOT") or Path(__file__).resolve().parents[2])
Q52_JSON = "data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json"
FIXTURE_JSON = "data/ibm_tomography_feb2026/simulator_test_20260209_125106.json"
MARCH_JSON = "data/ibm_shadow_march2026/shadow_hardware_combined_20260309_181852.json"
MARCH_SIM_JSON = "data/ibm_shadow_march2026/shadow_simulate_20260309_181709.json"
RAMSEY_MAR12 = "data/ibm_run3_march2026/ramsey_march12_20260312.json"
RAMSEY_MAR18 = "data/ibm_run3_march2026/ramsey_sameday_20260318.json"
RUN3_JSON = "data/ibm_run3_march2026/palindrome_ibm_torino_20260318_191348.json"
SHOT_SIGMA = 1.0 / (2.0 * np.sqrt(8192))  # shot noise of one component of rho01 near zero
Q52_LATE_MEAN = 0.0185  # the 15-sample late mean at t/T2_echo >= 1.25, as the pages quote it


def load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def assert_mc(estimate: float, standard_error: float, cited: float, decimals: int, label: str) -> None:
    """A Monte Carlo estimate against a cited value: four standard errors plus half a printed unit."""
    bound = 4.0 * standard_error + 0.5 * 10.0 ** (-decimals)
    assert abs(estimate - cited) <= bound, f"{label}: {estimate} vs cited {cited} (bound {bound:.2e})"


def sem(values) -> float:
    values = np.asarray(values)
    return float(values.std(ddof=1) / np.sqrt(len(values)))


# ------------------------------------------------------------------------------------------ Q52 record
class Q52:
    def __init__(self, record: dict):
        self.record = record
        self.T1 = float(record["T1_us"])
        self.T2e = float(record["T2_us"])
        raw = record["raw_tomography"]
        self.t = np.array([float(r["delay_us"]) for r in raw])
        rho = np.array([np.array(r["density_matrix_real"]) + 1j * np.array(r["density_matrix_imag"]) for r in raw])
        self.r01 = rho[:, 0, 1]
        self.a = np.abs(self.r01)
        self.X = 2.0 * self.r01.real
        self.Y = -2.0 * self.r01.imag
        self.Z = (rho[:, 0, 0] - rho[:, 1, 1]).real
        self.x = self.t / self.T2e
        self.late = self.x >= 1.0


def late_direction_counts(record: Q52) -> tuple[int, int, int]:
    """(late samples, how many with Re > 0, how many with Im < 0)."""
    late = record.r01[record.late]
    return int(late.size), int(np.sum(late.real > 0)), int(np.sum(late.imag < 0))


def late_phase_drift(record: Q52):
    """Straight-line drift of the unwrapped late phases, in degrees per microsecond."""
    return stats.linregress(record.t[record.late], np.degrees(np.unwrap(np.angle(record.r01[record.late]))))


def envelope(t):
    """The early coherence envelope of the record (the fit above the 0.005 noise floor)."""
    return 0.4785 * np.exp(-np.asarray(t) / 110.73)


def mean_direction(values) -> tuple[float, float]:
    """Direction of the mean of complex samples and its standard error, in degrees (first-order propagation with
    the Re/Im covariance neglected; on the late Q52 samples that is the cautious side, 5.5 deg against 4.7 with it)."""
    values = np.asarray(values)
    z = values.mean()
    se = np.sqrt((z.imag * sem(values.real)) ** 2 + (z.real * sem(values.imag)) ** 2) / abs(z) ** 2
    return float(np.degrees(np.angle(z))), float(np.degrees(se))


def fixed_point_r_minus(record) -> complex:
    """R- = (1 - 2 C Psi - i sqrt(4 C Psi - 1)) / (2 C) at the last sample before the crossing, as proposed."""
    analysis = record.record["analysis"]
    k = max(i for i, a in enumerate(analysis) if a["cpsi_measured"] >= 0.25)
    c, psi = analysis[k]["C_measured"], analysis[k]["psi_measured"]
    return (1 - 2 * c * psi - 1j * np.sqrt(4 * c * psi - 1)) / (2 * c)


def quadrant(z: complex) -> int:
    return 1 if z.real > 0 and z.imag >= 0 else 2 if z.real <= 0 < z.imag else 3 if z.real < 0 and z.imag <= 0 else 4


@pytest.fixture(scope="module")
def q52():
    return Q52(load(Q52_JSON))


def test_q52_axis_grid_and_stored_derived_quantities_are_exact(q52):
    analysis = q52.record["analysis"]
    assert len(q52.t) == len(analysis) == 25 and q52.record["shots"] == 8192
    # the stored axis is t over the calibration echo time, exactly
    assert all(a["delay_over_T2"] == q52.t[k] / q52.T2e for k, a in enumerate(analysis))
    # the 25 delays are one uniform grid of T2_echo/8, not logarithmic
    assert np.array_equal(q52.t, np.arange(25) * (q52.T2e / 8))
    assert f"{q52.T2e / 8:.2f}" == "37.28" and f"{1e3 / (q52.T2e / 8):.6f}" == "26.823361"
    # Psi = 2|rho01| and CPsi = C*Psi exactly, so r(|rho01|, 1/4 - CPsi) reuses one quantity
    psi = np.array([a["psi_measured"] for a in analysis])
    cpsi = np.array([a["cpsi_measured"] for a in analysis])
    c = np.array([a["C_measured"] for a in analysis])
    assert np.array_equal(psi, 2.0 * q52.a) and np.array_equal(cpsi, c * psi)
    r = np.corrcoef(q52.a[q52.late], 0.25 - cpsi[q52.late])[0, 1]
    assert f"{r:.4f}" == "-0.9955"
    # T2_echo is 2.7 free-decay times on this qubit; the late window starts where 7% of the coherence is left
    assert (f"{q52.T2e / 110.7:.1f}", f"{np.exp(-q52.T2e / 110.7):.3f}") == ("2.7", "0.068")
    # the stored Pauli readings are raw count ratios, 2n/8192 - 1 with n an integer: no readout mitigation
    counts = (np.concatenate([q52.X, q52.Y, q52.Z]) + 1.0) * 4096.0
    assert np.max(np.abs(counts - np.rint(counts))) < 1e-9


def test_q52_late_direction_is_static_with_equal_x_and_y_offsets(q52):
    late = q52.late
    assert late_direction_counts(q52) == (17, 17, 17)
    mean = q52.r01[late].mean()
    assert (f"{mean.real:+.5f}", f"{mean.imag:+.5f}") == ("+0.01137", "-0.01279")
    assert (f"{np.degrees(np.angle(mean)):.1f}", f"{np.degrees(np.angle(q52.r01[late])).mean():.0f}") == ("-48.4", "-44")
    assert (f"{q52.r01[late].real.std():.5f}", f"{q52.r01[late].imag.std():.5f}") == ("0.00499", "0.00805")
    x_late, y_late = q52.X[late], q52.Y[late]
    assert (f"{x_late.mean():.4f}", f"{sem(x_late):.4f}") == ("0.0227", "0.0025")
    assert (f"{y_late.mean():.4f}", f"{sem(y_late):.4f}") == ("0.0256", "0.0040")
    assert (f"{x_late.mean() / sem(x_late):.1f}", f"{y_late.mean() / sem(y_late):.1f}") == ("9.1", "6.4")
    # <X> = <Y> at every cut, and the mean direction with its standard error
    expected = {1.0: (17, "-0.0028", "0.0041", "0.7", "-48.4", "5.5"),
                1.25: (15, "-0.0045", "0.0044", "1.0", "-50.2", "5.8"),
                1.5: (13, "-0.0063", "0.0049", "1.3", "-51.9", "5.6")}
    for cut, (n, diff, diff_se, sigmas, angle, angle_se) in expected.items():
        m = q52.x >= cut
        d = (q52.X - q52.Y)[m]
        got_angle, got_error = mean_direction(q52.r01[m])
        got = (int(m.sum()), f"{d.mean():+.4f}", f"{sem(d):.4f}", f"{abs(d.mean()) / sem(d):.1f}",
               f"{got_angle:.1f}", f"{got_error:.1f}")
        assert got == (n, diff, diff_se, sigmas, angle, angle_se), (cut, got)
    # the late phase does not rotate
    drift = late_phase_drift(q52)
    assert f"{drift.slope:.3f}" == "-0.040" and f"{drift.pvalue:.2f}" == "0.18"
    # <Z> climbs across the late samples as the qubit relaxes
    assert (f"{q52.Z[0]:.3f}", f"{q52.Z[late].min():.2f}", f"{q52.Z[late].max():.2f}") == ("0.005", "0.50", "0.71")
    # 9 of the 15 phases at t/T2_echo >= 1.25 fall in the -60 deg bin, 4 in the -30 deg bin
    bins = np.rint(np.degrees(np.angle(q52.r01[q52.x >= 1.25])) / 30.0) * 30
    assert (int(np.sum(bins == -60)), int(np.sum(bins == -30))) == (9, 4)


@pytest.fixture(scope="module")
def late_fits(q52):
    """The three late-only offset fits the pages quote: a constant from 300 us, a constant plus the decaying
    rotation from 300 us and from 250 us."""
    return [offset_fit(q52.t, q52.r01, 300.0, False), offset_fit(q52.t, q52.r01, 300.0, True),
            offset_fit(q52.t, q52.r01, 250.0, True)]


def late_direction_distance(record, r_minus: complex) -> float:
    """Standard errors between the late mean direction and the fixed point's."""
    angle, error = mean_direction(record.r01[record.late])
    return (np.degrees(np.angle(r_minus)) - angle) / error


def test_q52_neither_tracked_nor_froze_at_the_fixed_point(q52, late_fits):
    """The proposal, tested on its own record: the phase at the crossing and after it, and the late direction
    against the fixed point's. (The pre-registered question is cross-qubit; one record cannot test it.)"""
    t_star = float(q52.record["crossing_us"])
    k = int(np.searchsorted(q52.t, t_star))
    assert (f"{q52.t[k - 1]:.1f}", f"{t_star:.1f}", f"{q52.t[k]:.1f}") == ("111.8", "114.7", "149.1")
    phases = np.degrees(np.angle(q52.r01))
    assert [f"{p:+.1f}" for p in phases[k - 1:k + 4]] == ["+79.4", "+89.0", "-48.1", "-146.9", "+100.0"]
    assert (f"{q52.t[k + 3]:.1f}", f"{q52.t[k + 4]:.1f}") == ("261.0", "298.2")
    # the phase is still turning 150 us after the crossing (the last turning sample, 261 us), and the fixed
    # direction starts at the next sample, t = T2_echo = 298 us
    assert round(q52.t[k + 3] - t_star, -1) == 150.0 and q52.x[k + 4] == 1.0
    # sampled to the March depth (the envelope above 9% of its start) the record turns through all four quadrants
    deep = envelope(q52.t) / envelope(0.0) >= 0.09
    assert int(deep.sum()) == 8 and f"{q52.t[deep][-1]:.1f}" == "261.0"
    assert {quadrant(z) for z in q52.r01[deep]} == {1, 2, 3, 4}
    # the late direction is not the fixed point's: 4.4 standard errors from R-, 0.6 from -45 deg
    r_minus = fixed_point_r_minus(q52)
    angle, error = mean_direction(q52.r01[q52.late])
    assert (f"{np.degrees(np.angle(r_minus)):.1f}", f"{angle:.1f}", f"{error:.1f}") == ("-24.1", "-48.4", "5.5")
    assert (f"{late_direction_distance(q52, r_minus):.1f}", f"{abs(angle + 45.0) / error:.1f}") == ("4.4", "0.6")
    # with the Re/Im covariance (sample r = -0.27) the error is 4.7 deg and the distance 5.1: 4.4 is the cautious figure
    late = q52.r01[q52.late]
    mean = late.mean()
    cov = np.cov(np.vstack([late.real, late.imag])) / late.size
    grad = np.array([-mean.imag, mean.real]) / abs(mean) ** 2
    full = np.degrees(np.sqrt(grad @ cov @ grad))
    assert (f"{np.corrcoef(late.real, late.imag)[0, 1]:.2f}", f"{full:.1f}",
            f"{(np.degrees(np.angle(r_minus)) - angle) / full:.1f}") == ("-0.27", "4.7", "5.1")
    # each late-only fit sits 4.0 to 4.3 standard errors from R- and within 0.6 to 1.5 of -45 deg
    from_r = [(np.degrees(np.angle(r_minus)) - np.degrees(np.angle(c))) / np.degrees(se) for c, _, se in late_fits]
    from_45 = [abs(np.degrees(np.angle(c)) + 45.0) / np.degrees(se) for c, _, se in late_fits]
    assert (f"{min(from_r):.1f}", f"{max(from_r):.1f}") == ("4.0", "4.3")
    assert (f"{min(from_45):.1f}", f"{max(from_45):.1f}") == ("0.6", "1.5")
    # a static component of about 0.017, which an exponential envelope reaches near 3.3 T2*; it showed from 2.7 T2*
    offset = abs(q52.r01[q52.late].mean())
    assert f"{offset:.3f}" == "0.017"
    assert (f"{np.log(0.4785 / offset):.1f}", f"{q52.t[k + 4] / 110.7:.1f}") == ("3.3", "2.7")


def test_q52_early_component_rotates_and_the_470_us_line_mixes_two_components(q52):
    phases = np.degrees(np.angle(q52.r01[:5]))
    assert [f"{p:+.1f}" for p in phases] == ["+1.1", "-62.4", "-146.7", "+79.4", "+89.0"]
    slope = np.polyfit(q52.t[:5], np.unwrap(np.angle(q52.r01[:5])), 1)[0]
    assert (f"{slope:.6f}", f"{slope / (2 * np.pi) * 1e3:.3f}") == ("-0.035696", "-5.681")
    assert (f"{np.degrees(slope * q52.T2e / 8):.1f}", f"{2 * np.pi / abs(slope):.0f}") == ("-76.2", "176")
    early = q52.t < 300
    slope300 = np.polyfit(q52.t[early], np.unwrap(np.angle(q52.r01[early])), 1)[0]
    assert f"{slope300 / (2 * np.pi) * 1e3:.1f}" == "-6.9"
    steps = np.degrees(np.angle(q52.r01[1:9] / q52.r01[:8]))
    assert [f"{s:+.1f}" for s in steps] == ["-63.5", "-84.3", "-133.9", "+9.6", "-137.0", "-98.9", "-113.1", "-142.8"]
    # the producer's "Detuning from phase": one line through the 24 samples above 0.01, 16 of them static
    above = q52.a > 0.01
    line = np.polyfit(q52.t[above], np.unwrap(np.angle(q52.r01[above])), 1)[0]
    assert above.sum() == 24 and (above & q52.late).sum() == 16
    assert (f"{line:.4f}", f"{2 * np.pi / abs(line):.0f}") == ("-0.0134", "470")
    # the early magnitude falls as 0.48 exp(-t/110 us)
    popt = curve_fit(lambda tt, A, T: A * np.exp(-tt / T), q52.t[q52.a > 0.005], q52.a[q52.a > 0.005],
                     p0=[q52.a[0], 100])[0]
    assert (f"{popt[0]:.2f}", f"{popt[1]:.1f}") == ("0.48", "110.7")


def test_q52_tail_slope_depends_on_the_cut(q52):
    expected = {  # cut: (samples, slope per T2_echo, two-sided p, 95% interval)
        ">=1.0": (17, "+0.00699", "0.014", ("+0.00164", "+0.01234")),
        ">=1.25": (15, "+0.00878", "0.015", ("+0.00203", "+0.01554")),
        ">=1.5": (13, "+0.00819", "0.053", ("-0.00013", "+0.01651")),
        ">1.5": (12, "+0.01041", "0.033", ("+0.00103", "+0.01978")),
    }
    masks = {">=1.0": q52.x >= 1.0, ">=1.25": q52.x >= 1.25, ">=1.5": q52.x >= 1.5, ">1.5": q52.x > 1.5}
    for cut, (n, slope, p, (lo, hi)) in expected.items():
        mask = masks[cut]
        fit = stats.linregress(q52.x[mask], q52.a[mask])
        tcrit = stats.t.ppf(0.975, mask.sum() - 2)
        got = (int(mask.sum()), f"{fit.slope:+.5f}", f"{fit.pvalue:.3f}",
               (f"{fit.slope - tcrit * fit.stderr:+.5f}", f"{fit.slope + tcrit * fit.stderr:+.5f}"))
        assert got == (n, slope, p, (lo, hi)), (cut, got)
    fit = stats.linregress(q52.x[masks[">=1.5"]], q52.a[masks[">=1.5"]])
    assert f"{fit.slope / q52.T2e:.2e}" == "2.75e-05" and f"{fit.rvalue:.4f}" == "0.5469"
    tail = q52.a[masks[">=1.5"]]
    assert not np.all(np.diff(tail) > 0)  # the amplitudes are non-monotone
    # the late local maxima sit on the sampling grid of 0.125 T2_echo
    xl, al = q52.x[q52.late], q52.a[q52.late]
    peaks = [i for i in range(1, len(al) - 1) if al[i] > al[i - 1] and al[i] > al[i + 1] and xl[i] >= 1.5]
    assert [f"{xl[i]:.3f}" for i in peaks] == ["1.500", "1.750", "2.000", "2.625"]
    # two roundings (the delay product and the quotient by T2_echo) separate x*8 from an integer
    assert np.max(np.abs(xl * 8 - np.rint(xl * 8))) < 1e-12


def run_null(t, T2s, runs, rng, offset=0j, amplitude=0.5):
    """The recorded null: amplitude exp(-t/T2*) at one random phase per run, binomial X and Y shots at 8192."""
    phi = rng.uniform(0.0, 2.0 * np.pi, size=(runs, 1))
    rho01 = amplitude * np.exp(-t / T2s)[None, :] * np.exp(1j * phi) + offset
    px = np.clip((1.0 + 2.0 * rho01.real) / 2.0, 0.0, 1.0)
    py = np.clip((1.0 - 2.0 * rho01.imag) / 2.0, 0.0, 1.0)
    xm = 2.0 * rng.binomial(8192, px) / 8192 - 1.0
    ym = 2.0 * rng.binomial(8192, py) / 8192 - 1.0
    return np.abs((xm - 1j * ym) / 2.0)


def test_q52_excess_over_the_rebuilt_null_is_real_and_robust(q52):
    rng = np.random.default_rng(20260209)
    selected = q52.t > 370  # the 15 samples at t/T2_echo >= 1.25
    assert selected.sum() == 15 and np.array_equal(selected, q52.x >= 1.25)
    observed = q52.a[selected].mean()
    assert f"{observed:.5f}" == "0.01852" and f"{observed:.4f}" == f"{Q52_LATE_MEAN:.4f}"
    samples = run_null(q52.t, 110.2, 10000, rng)
    late_means = samples[:, selected].mean(axis=1)
    m, s = late_means.mean(), late_means.std(ddof=1)
    assert_mc(m, s / np.sqrt(10000), 0.00859, 5, "rebuilt null mean")
    assert_mc(m, s / np.sqrt(10000), 0.00861, 5, "recorded null mean")
    assert_mc(s, s / np.sqrt(2 * 9999), 0.00105, 5, "null spread")
    assert int((late_means >= observed).sum()) == 0
    z = (observed - m) / s
    z_error = 1.0 / np.sqrt(10000) + z / np.sqrt(2 * 9999)
    assert abs(z - 9.5) <= 4.0 * z_error + 0.05, z
    # 9 of 25 samples above the per-sample 99th percentile, two early and seven late; one of them
    # (t/T2_echo = 2.125) sits inside the percentile's own rank uncertainty (four binomial standard
    # deviations of the rank), every other sample clears it
    above = q52.a > np.percentile(samples, 99, axis=0)
    assert [float(v) for v in np.round(q52.x[above], 3)] == [0.25, 0.375, 2.0, 2.125, 2.375, 2.5, 2.625, 2.875, 3.0]
    ordered = np.sort(samples, axis=0)
    rank_sd = np.sqrt(10000 * 0.99 * 0.01)
    low = ordered[int(9900 - 4 * rank_sd)]
    high = ordered[int(9900 + 4 * rank_sd)]
    edge = (q52.a > low) & (q52.a <= high)
    assert [float(v) for v in np.round(q52.x[edge], 3)] == [2.125]
    at_2625 = np.flatnonzero(np.isclose(q52.x, 2.625))[0]
    assert low[at_2625] <= 0.0167 <= high[at_2625] and f"{q52.a[at_2625] / 0.0167:.1f}" == "2.2"
    # robust to the decay constant; only 150 us, far off the fitted 110, lets a handful through
    for T2s in (100.0, 110.0, 120.0, 130.0):
        assert int((run_null(q52.t, T2s, 4000, rng)[:, selected].mean(axis=1) >= observed).sum()) == 0, T2s
    at_150 = int((run_null(q52.t, 150.0, 4000, rng)[:, selected].mean(axis=1) >= observed).sum())
    assert 0 < at_150 <= 20
    # with the fitted amplitude 0.48 instead of the ideal 1/2 the null's late mean drops to 0.00849
    fitted = run_null(q52.t, 110.2, 10000, rng, amplitude=0.4785)[:, selected].mean(axis=1)
    assert_mc(fitted.mean(), fitted.std(ddof=1) / np.sqrt(10000), 0.00849, 5, "null at the fitted amplitude")
    assert int((fitted >= observed).sum()) == 0
    # an offset of the late samples' own size (delta = 0.024) added to the null reproduces the late mean
    with_offset = run_null(q52.t, 110.2, 4000, rng, offset=(0.024 - 0.024j) / 2)[:, selected].mean(axis=1)
    mo, so = with_offset.mean(), with_offset.std(ddof=1)
    assert_mc(mo, so / np.sqrt(4000), 0.0186, 4, "null plus the late offset")
    assert_mc(so, so / np.sqrt(2 * 3999), 0.0028, 4, "null plus the late offset, spread")
    assert f"{(q52.X[q52.late].mean() + q52.Y[q52.late].mean()) / 2:.3f}" == "0.024"


def basis_change_readings(eta):
    """Exact (sympy): the Pauli coefficients of what a Z measurement reads after the X and after the Y basis change,
    with the basis-change pulse over-rotated by eta. The X change is rz(pi/2) sx rz(pi/2), the IBM form of a
    Hadamard; the Y change puts S-dagger in front of it. Also returns the Hadamard itself."""
    pauli_x = sp.Matrix([[0, 1], [1, 0]])
    pauli_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    pauli_z = sp.diag(1, -1)
    rz = sp.diag(sp.exp(-sp.I * sp.pi / 4), sp.exp(sp.I * sp.pi / 4))
    theta = sp.pi / 2 + eta
    sx = sp.cos(theta / 2) * sp.eye(2) - sp.I * sp.sin(theta / 2) * pauli_x
    hadamard = rz * sx * rz

    def reads(u):
        observable = u.H * pauli_z * u
        return [sp.simplify(sp.expand_complex((observable * p).trace() / 2)) for p in (pauli_x, pauli_y, pauli_z)]

    return reads(hadamard), reads(hadamard * sp.diag(1, -sp.I)), hadamard, pauli_z


def test_q52_the_record_does_not_rank_the_two_measurement_forms(q52):
    # the leak: one pulse error leaks <Z> into the X and the Y reading with the same coefficient, -sin(eta), because
    # the S-dagger that makes the Y change out of the X one commutes with Z; exact, and exactly the ideal at eta = 0
    eta = sp.symbols("eta", real=True)
    x_reads, y_reads, hadamard, pauli_z = basis_change_readings(eta)
    assert [sp.simplify(c - e) for c, e in zip(x_reads, (sp.cos(eta), 0, -sp.sin(eta)))] == [0, 0, 0]
    assert [sp.simplify(c - e) for c, e in zip(y_reads, (0, sp.cos(eta), -sp.sin(eta)))] == [0, 0, 0]
    # a |+> made with the same pulse starts at <Z> = -sin(eta), the leak coefficient itself
    prepared = hadamard * sp.Matrix([1, 0])
    assert sp.simplify(sp.expand_complex((prepared.H * pauli_z * prepared)[0]) + sp.sin(eta)) == 0
    # a pulse 0.04 rad short reads a coherence-free state with <Z> = 0.71 as <X> = <Y> = +0.028
    assert f"{float(-sp.sin(-0.04) * 0.71):+.3f}" == "+0.028"
    # the two forms on the 34 late readings, one parameter each (shot-noise chi2, 33 degrees of freedom each)
    sigma = 1.0 / np.sqrt(8192)  # shot noise of one Pauli expectation near zero
    readings = np.concatenate([q52.X[q52.late], q52.Y[q52.late]])
    z_late = np.concatenate([q52.Z[q52.late], q52.Z[q52.late]])
    delta = readings.mean()
    eps = float(readings @ z_late / (z_late @ z_late))
    chi_delta = np.sum((readings - delta) ** 2) / sigma ** 2
    chi_eps = np.sum((readings - eps * z_late) ** 2) / sigma ** 2
    assert (readings.size, f"{delta:.3f}", f"{chi_delta:.1f}", f"{eps:.3f}", f"{chi_eps:.1f}") == \
        (34, "0.024", "50.5", "0.039", "45.4")
    # the chi2 pair leans toward the leak (likelihood ratio ~13, p 0.03 against 0.07) without deciding
    assert (f"{np.exp((chi_delta - chi_eps) / 2):.0f}", f"{stats.chi2.sf(chi_delta, 33):.2f}",
            f"{stats.chi2.sf(chi_eps, 33):.2f}") == ("13", "0.03", "0.07")
    # the late slopes on <Z> separate neither form: <X> sits on the leak's 0.039 and 0.8 sigma from 0, <Y> is
    # steeper than both (1.8 sigma above the leak, 2.4 sigma above 0)
    sx_fit = stats.linregress(q52.Z[q52.late], q52.X[q52.late])
    sy_fit = stats.linregress(q52.Z[q52.late], q52.Y[q52.late])
    assert (f"{sx_fit.slope:.2f}", f"{sx_fit.stderr:.2f}", f"{sy_fit.slope:.2f}", f"{sy_fit.stderr:.2f}") == \
        ("0.04", "0.04", "0.15", "0.06")
    assert abs(sx_fit.slope - eps) < sx_fit.stderr and abs(sy_fit.slope - eps) < 2 * sy_fit.stderr
    assert sx_fit.slope < 2 * sx_fit.stderr and sy_fit.slope < 3 * sy_fit.stderr
    # the first sample splits: <Y> stands against the asymmetry (a 2.6 deg phase error inside that reading), not
    # against the leak; <Z> stands against the leak only if the preparation shares the pulse
    assert (f"{q52.Y[0]:+.3f}", f"{sigma:.3f}", f"{(q52.Y[0] - 0.024) / sigma:.1f}") == ("-0.018", "0.011", "-3.8")
    assert f"{(q52.Y[0] - eps * q52.Z[0]) / sigma:.1f}" == "-1.6"
    phase_error = np.degrees(np.arctan2(q52.Y[0] - 0.024, q52.X[0] - 0.024))
    assert f"{phase_error:.1f}" == "-2.6"
    assert (f"{q52.Z[0]:+.3f}", f"{(q52.Z[0] - 0.024) / sigma:.1f}", f"{(eps - q52.Z[0]) / sigma:.1f}") == \
        ("+0.005", "-1.7", "3.0")
    # the <Z> relaxation intercept leans to the asymmetry, with the relaxation time fitted (241 us) ...
    f = lambda tt, req, d, rate: req + d * np.exp(-rate * tt)
    popt, pcov = curve_fit(f, q52.t, q52.Z, p0=[0.7, -0.7, 1 / q52.T1])
    grad = np.array([1.0, 1.0, 0.0])
    intercept, intercept_error = popt[0] + popt[1], np.sqrt(grad @ pcov @ grad)
    assert (f"{intercept:+.3f}", f"{intercept_error:.3f}", f"{1 / popt[2]:.0f}") == ("+0.021", "0.008", "241")
    assert (f"{(eps - intercept) / intercept_error:.1f}", f"{intercept / intercept_error:.1f}") == ("2.2", "2.5")
    # ... and with the calibration T1 held fixed it reads +0.006 +- 0.007, a fit worse by 14 in shot-noise chi2
    g = lambda tt, req, d: req + d * np.exp(-tt / q52.T1)
    p_cal, c_cal = curve_fit(g, q52.t, q52.Z, p0=[0.7, -0.7])
    z_sigma = np.sqrt(1.0 - q52.Z ** 2) / np.sqrt(8192)
    worse = np.sum(((q52.Z - g(q52.t, *p_cal)) / z_sigma) ** 2) - np.sum(((q52.Z - f(q52.t, *popt)) / z_sigma) ** 2)
    assert (f"{p_cal[0] + p_cal[1]:+.3f}", f"{np.sqrt(np.ones(2) @ c_cal @ np.ones(2)):.3f}", f"{worse:.0f}") == \
        ("+0.006", "0.007", "14")


def offset_fit(t, r01, t_min, with_rotation):
    """Late-only offset fit: a constant, or a constant plus the decaying rotation on the fixed early envelope."""
    m = t >= t_min
    tt, rr = t[m], r01[m]
    if not with_rotation:
        c = rr.mean()
        chi2 = np.sum(np.abs(rr - c) ** 2) / SHOT_SIGMA ** 2
        se = np.sqrt((c.imag * sem(rr.real)) ** 2 + (c.real * sem(rr.imag)) ** 2) / abs(c) ** 2
        return c, chi2 / (2 * m.sum() - 2), se

    def residual(p):
        ph, w, cr, ci = p
        d = rr - (envelope(tt) * np.exp(1j * (ph + w * tt)) + cr + 1j * ci)
        return np.concatenate([d.real, d.imag]) / SHOT_SIGMA

    best = min((least_squares(residual, [ph0, w0, 0.01, -0.01], bounds=([-np.pi, -0.09, -0.1, -0.1], [np.pi, 0.09, 0.1, 0.1]))
                for w0 in np.linspace(-0.08, 0.08, 65) for ph0 in np.linspace(-np.pi, np.pi, 9)), key=lambda s: s.cost)
    dof = 2 * m.sum() - 4
    cov = np.linalg.inv(best.jac.T @ best.jac) * (2 * best.cost / dof)
    cr, ci = best.x[2], best.x[3]
    grad = np.array([0.0, 0.0, -ci, cr]) / (cr ** 2 + ci ** 2)
    return cr + 1j * ci, 2 * best.cost / dof, float(np.sqrt(grad @ cov @ grad))


def test_q52_the_offset_is_there_and_its_direction_is_only_suggestive(q52, late_fits):
    popt, pcov = curve_fit(lambda tt, A, T, c: A * np.exp(-tt / T) + c, q52.t, q52.a, p0=[0.47, 110.0, 0.01])
    assert (f"{popt[2]:.4f}", f"{np.sqrt(pcov[2, 2]):.4f}") == ("0.0123", "0.0051")

    def residual(p, with_offset):
        model = p[0] * np.exp(-q52.t / p[1]) * np.exp(1j * (p[2] + p[3] * q52.t))
        if with_offset:
            model = model + p[4] + 1j * p[5]
        d = q52.r01 - model
        return np.concatenate([d.real, d.imag]) / SHOT_SIGMA

    best = {}
    for with_offset in (False, True):
        extra = [0.01, -0.01] if with_offset else []
        lo = [0, 10, -np.pi, -0.09] + ([-0.1, -0.1] if with_offset else [])
        hi = [1, 1000, np.pi, 0.09] + ([0.1, 0.1] if with_offset else [])
        fits = [least_squares(residual, [0.48, 110.0, 0.0, w0] + extra, args=(with_offset,), bounds=(lo, hi))
                for w0 in np.linspace(-0.08, 0.08, 33)]
        best[with_offset] = min(fits, key=lambda f: f.cost)
    delta = 2.0 * best[False].cost - 2.0 * best[True].cost
    misfit = 2.0 * best[True].cost / 44
    assert (f"{delta:.0f}", f"{misfit:.0f}") == ("117", "34")
    # scaled by the early model's misfit the offset's improvement is not significant
    assert f"{stats.f.sf((delta / 2) / misfit, 2, 44):.2f}" == "0.19"
    # fits of the late samples alone, where a simple model fits, put the offset between -49 and -56 deg
    assert [f"{np.degrees(np.angle(c)):.0f}" for c, _, _ in late_fits] == ["-49", "-55", "-56"]
    assert [f"{chi:.1f}" for _, chi, _ in late_fits] == ["1.7", "1.5", "2.0"]
    assert all(abs(np.degrees(np.angle(c)) + 45.0) <= 1.5 * np.degrees(se) for c, _, se in late_fits)


def test_q52_proposed_fixed_point_and_the_absorption_two_fits_ratio(q52):
    analysis = q52.record["analysis"]
    k = max(i for i, a in enumerate(analysis) if a["cpsi_measured"] >= 0.25)
    c, psi = analysis[k]["C_measured"], analysis[k]["psi_measured"]
    r_minus = (1 - 2 * c * psi - 1j * np.sqrt(4 * c * psi - 1)) / (2 * c)
    assert (f"{q52.t[k]:.1f}", f"{c:.3f}", f"{psi:.3f}") == ("111.8", "0.624", "0.419")
    assert (f"{r_minus.real:.4f}", f"{r_minus.imag:.4f}", f"{np.degrees(np.angle(r_minus)):.1f}") == ("0.3824", "-0.1712", "-24.1")
    # the absorption ratio: coherence fit and population fit of one record, against the calibration T1
    floor = q52.a > 0.005
    alpha_coh = curve_fit(lambda tt, A, a: A * np.exp(-a * tt), q52.t[floor], q52.a[floor],
                          p0=[q52.a[0], 1 / 100], maxfev=5000)[0][1]
    alpha_z = curve_fit(lambda tt, req, d, a: req + d * np.exp(-a * tt), q52.t, q52.Z,
                        p0=[0.7, -0.7, 1 / q52.T1], maxfev=5000)[0][2]
    assert (f"{1 / alpha_coh:.1f}", f"{1 / alpha_z:.1f}") == ("110.7", "241.4")
    excess = alpha_coh - alpha_z / 2
    two_gamma_star = alpha_coh - 1 / (2 * q52.T1)
    assert (f"{excess:.6f}", f"{two_gamma_star:.6f}", f"{excess / two_gamma_star:.2f}") == ("0.006960", "0.006771", "1.03")
    # the whole 3% is the T1 difference: with the population-fit T1 on both sides the ratio is 1 (exact algebra)
    a, az, t1 = sp.symbols("alpha alpha_Z T1", positive=True)
    assert sp.simplify(((a - az / 2) / (a - 1 / (2 * t1))).subs(t1, 1 / az)) == 1
    two_gamma_echo = 1 / q52.T2e - 1 / (2 * q52.T1)
    assert f"{excess / two_gamma_echo:.2f}" == "6.37"


def generalized_crossing(T2s: float, T1: float) -> float:
    r = T2s / T1
    return -np.log(brentq(lambda b: (1 - b ** r + b ** (2 * r) / 2 + b * b / 2) * b - 0.25, 1e-9, 1 - 1e-12))


def test_q52_crossing_compares_one_t2star_on_both_sides(q52):
    floor = q52.a > 0.005
    fit_floor = curve_fit(lambda tt, A, T: A * np.exp(-tt / T), q52.t[floor], q52.a[floor], p0=[q52.a[0], 100])[0][1]
    fit_all = curve_fit(lambda tt, A, T: A * np.exp(-tt / T), q52.t, q52.a, p0=[q52.a[0], 100])[0][1]
    assert (f"{fit_floor:.1f}", f"{fit_all:.1f}") == ("110.7", "110.1")
    t_star = float(q52.record["crossing_us"])
    assert f"{t_star:.1f}" == "114.7"
    rows = {}
    # the registry's T2*(FID) = 110.7 us as it prints it, and the all-samples fit of IBM_QUANTUM_TOMOGRAPHY
    for label, T2s in (("registry", 110.7), ("all samples", fit_all)):
        predicted = generalized_crossing(T2s, q52.T1)
        rows[label] = (f"{T2s / q52.T1:.3f}", f"{predicted:.3f}", f"{t_star / T2s:.3f}", f"{100 * (t_star / T2s / predicted - 1):.1f}")
    assert rows["registry"] == ("0.500", "0.950", "1.036", "9.1")
    assert rows["all samples"] == ("0.498", "0.949", "1.041", "9.7")
    # r = 0.456 would need T2* = 100.9 us, a different T2* than the measured ratio uses
    assert f"{0.456 * q52.T1:.1f}" == "100.9"
    assert f"{generalized_crossing(0.456 * q52.T1, q52.T1):.3f}" == "0.936"


def test_q52_slow_tail_registry_crossing_numbers_and_figure_readings(q52):
    # IBM_ABSORPTION_THEOREM section 4: the bi-exponential of the envelope above the 0.005 floor, fitted the way
    # simulations/ibm_absorption_theorem.py fits it (its start and bounds copied, the producer not imported)
    floor = q52.a > 0.005
    single = curve_fit(lambda tt, A, a: A * np.exp(-a * tt), q52.t[floor], q52.a[floor],
                       p0=[q52.a[0], 1 / 100], maxfev=5000)[0][1]
    fit = curve_fit(lambda tt, A1, a1, A2, a2: A1 * np.exp(-a1 * tt) + A2 * np.exp(-a2 * tt), q52.t[floor], q52.a[floor],
                    p0=[q52.a[0] * 0.9, single * 1.5, q52.a[0] * 0.1, single * 0.1],
                    bounds=([0, 1e-6, 0, 0], [1, 0.1, 0.5, single]), maxfev=10000)[0]
    (fast, fast_rate), (slow, slow_rate) = sorted([(fit[0], fit[1]), (fit[2], fit[3])], key=lambda c: -c[1])
    total = fast + slow
    assert (f"{fast:.3f}", f"{fast_rate:.5f}", f"{1 / fast_rate:.0f}", f"{100 * fast / total:.1f}") == \
        ("0.470", "0.00984", "102", "97.2")
    assert (f"{slow:.3f}", f"{slow_rate:.5f}", f"{100 * slow / total:.1f}") == ("0.013", "0.00000", "2.8")
    # the crossing entry: CPsi(0) = 0.885, and the late purity 0.740 is the mean of the last five samples
    analysis = q52.record["analysis"]
    purity = np.array([a["C_measured"] for a in analysis])
    assert (f"{analysis[0]['cpsi_measured']:.3f}", f"{purity[-5:].mean():.3f}") == ("0.885", "0.740")
    # visualizations/README.md, the crossing figure: the purity falls to 0.58 at 149 us and climbs back to about
    # 0.75 by the end of the record; the figure's axis is t/T2_echo, where the crossing reads 0.384 and the
    # generalized line, drawn on the echo time, reads 1.250 at r = 1.35
    k = int(np.argmin(purity))
    assert (f"{purity[k]:.2f}", f"{q52.t[k]:.0f}", f"{purity[-1]:.2f}") == ("0.58", "149", "0.75")
    assert f"{float(q52.record['crossing_us']) / q52.T2e:.3f}" == "0.384"
    stored = q52.record["analytical_prediction_generalized"]
    assert (f"{q52.T2e / q52.T1:.2f}", f"{generalized_crossing(q52.T2e, q52.T1):.3f}") == ("1.35", "1.250")
    assert (f"{stored['r']:.2f}", f"{stored['t_star_over_T2']:.3f}") == ("1.35", "1.250")
    # the generalized-crossing figure: its five markers carry values of the curve itself (at r = 1/10, 1/6, 2/3,
    # 3/4 and 1), not measurements
    marked = [generalized_crossing(r, 1.0) for r in (0.1, 1 / 6, 2 / 3, 0.75, 1.0)]
    assert [f"{v:.3f}" for v in marked] == ["0.863", "0.870", "1.010", "1.043", "1.141"]


def test_simulator_fixture_is_a_cross_fixture_comparison(q52):
    fixture = load(FIXTURE_JSON)
    assert (fixture["T1_us"], fixture["T2_us"]) == (200.0, 150.0)
    assert "raw_tomography" not in fixture  # magnitudes only: no phases to compare
    xs = np.array([a["delay_over_T2"] for a in fixture["analysis"]])
    ts = np.array([a["delay_us"] for a in fixture["analysis"]])
    assert np.array_equal(xs, ts / 150.0)
    late = np.array([a["populations"]["rho_01_abs"] for a in fixture["analysis"]])[xs >= 1.25]
    assert len(late) == 15 and f"{late.mean():.4f}" == "0.0726"


# --------------------------------------------------------------------------------------- March records
@pytest.fixture(scope="module")
def march():
    record = load(MARCH_JSON)
    return {q["qubit"]: q for q in record["qubit_results"]}, record


def phase_line(qubit_result):
    t = np.array([p["delay_us"] for p in qubit_result["points"]])
    r = np.array([complex(p["rho01_re"], p["rho01_im"]) for p in qubit_result["points"]])
    unwrapped = np.unwrap(np.angle(r))
    slope, intercept = np.polyfit(t, unwrapped, 1)
    rms = np.degrees(np.sqrt(np.mean((unwrapped - (slope * t + intercept)) ** 2)))
    return t, r, slope / (2 * np.pi) * 1e3, rms


def march_quadrants(qubit_result):
    return [quadrant(complex(p["rho01_re"], p["rho01_im"])) for p in qubit_result["points"][2:]]


def stored_t2star_is_proxy(verdict) -> bool:
    return verdict["T2_star_us"] == verdict["T2_us"] / 2.5


def detuned_decay_fits(t, r):
    """A detuned decay, and a detuned decay plus a constant offset, fitted to one March record."""
    def residual(p, with_offset):
        model = p[0] * np.exp(-t / p[1]) * np.exp(1j * (p[2] + p[3] * t))
        if with_offset:
            model = model + p[4] + 1j * p[5]
        d = r - model
        return np.concatenate([d.real, d.imag]) / SHOT_SIGMA

    best = {}
    for with_offset in (False, True):
        extra = [0.0, 0.0] if with_offset else []
        lo = [0, 1, -np.pi, -1.0] + ([-0.2, -0.2] if with_offset else [])
        hi = [1, 500, np.pi, 1.0] + ([0.2, 0.2] if with_offset else [])
        best[with_offset] = min((least_squares(residual, [abs(r[0]), 25.0, 0.0, w0] + extra, args=(with_offset,), bounds=(lo, hi))
                                 for w0 in np.linspace(-0.3, 0.3, 61)), key=lambda s: s.cost)
    return best[False], best[True]


def test_march_scheduling_proxy_grids_and_directions(march):
    qubits, record = march
    assert set(qubits) == {80, 102} and record["backend"] == "ibm_torino"
    assert all(stored_t2star_is_proxy(q["verdict"]) for q in qubits.values())  # a scheduling proxy, exactly
    assert f"{qubits[80]['verdict']['T2_star_us']:.2f}" == "10.83"
    assert (f"{qubits[102]['verdict']['T2_us'] / qubits[102]['verdict']['T1_us']:.3f}",
            f"{qubits[80]['verdict']['T2_us'] / qubits[80]['verdict']['T1_us']:.3f}") == ("0.159", "0.170")
    assert record["delay_multiples"][-1] == 5.0  # five times the proxy T2_echo/2.5 is 2 T2_echo
    assert march_quadrants(qubits[80]) == [1] * 8
    assert march_quadrants(qubits[102]) == [3, 2, 1, 1, 4, 3, 3, 2]
    assert f"{qubits[80]['verdict']['mean_phase_deg']:.0f}" == "29"
    for q in qubits.values():
        assert all(p["residual_re"] < 0 for p in q["points"][2:])
        t2 = q["verdict"]["T2_us"]
        assert all(p["rho01_re_theory"] == 0.5 * np.exp(-p["delay_us"] / t2) for p in q["points"])
    floor = 1.0 / np.sqrt(record["shots"])
    assert [sum(p["residual_abs"] > floor for p in qubits[k]["points"][2:]) for k in (80, 102)] == [8, 7]
    for qid, grid, alias in ((102, "6.60", "151.5"), (80, "5.41", "184.7")):
        t = np.array([p["delay_us"] for p in qubits[qid]["points"]])
        step = t[t > 0].min()
        assert np.allclose(t / step, np.rint(t / step), rtol=0, atol=1e-9)
        assert (f"{step:.2f}", f"{1e3 / step:.1f}") == (grid, alias)


def test_march_phases_are_steady_rotations_and_ramsey_agrees(march):
    qubits, _ = march
    _, _, f102, rms102 = phase_line(qubits[102])
    _, r80, f80, rms80 = phase_line(qubits[80])
    assert (f"{f102:.1f}", f"{rms102:.1f}") == ("-25.7", "6.5")
    assert (f"{f80:.1f}", f"{rms80:.1f}") == ("1.9", "7.8")
    phases80 = np.degrees(np.angle(r80))
    assert np.all(np.diff(phases80[:7]) > 0) and not np.all(np.diff(phases80) > 0)
    ramsey = {}
    for label, path in (("mar12", RAMSEY_MAR12), ("mar18", RAMSEY_MAR18)):
        for q in load(path)["results"]:
            ramsey[(label, q["qubit"])] = (q["analysis"]["Frequency"]["value"], q["analysis"]["Frequency"]["stderr"])
    f12, s12 = ramsey[("mar12", 102)]
    f18, s18 = ramsey[("mar18", 102)]
    assert (f"{f12 / 1e3:.1f}", f"{s12 / 1e3:.1f}", f"{f18 / 1e3:.1f}", f"{s18 / 1e3:.1f}") == ("19.4", "1.3", "26.3", "0.5")
    assert f12 / s12 > 10 and f18 / s18 > 10  # Q102: tens of kHz, far outside zero on both days
    assert f12 / 1e3 <= abs(f102) <= f18 / 1e3  # the shadow run's own rotation sits between them in magnitude
    g12, t12 = ramsey[("mar12", 80)]
    g18, t18 = ramsey[("mar18", 80)]
    assert (f"{g12 / 1e3:.1f}", f"{t12 / 1e3:.1f}") == ("3.8", "8.1") and abs(g12) < t12  # Q80: consistent with zero
    assert abs(g18) < 1.0 and t18 != t18  # the March 18 fit ended at 0 Hz with no error estimate (NaN)
    assert abs(f80 * 1e3) < t12  # Q80's +1.9 kHz drift sits inside what its Ramsey fit resolves
    t2star = {q["qubit"]: q["T2star_us"] for q in load(RAMSEY_MAR12)["results"]}
    assert f"{t2star[80]:.2f}" == "11.01"
    assert f"{load(RAMSEY_MAR18)['results'][0]['T2star_us']:.2f}" == "17.36"


def test_march_window_closed_before_a_q52_sized_offset_could_show(march, q52):
    qubits, _ = march
    rows = {}
    for qid in (80, 102):
        t = np.array([p["delay_us"] for p in qubits[qid]["points"]])
        r = np.array([complex(p["rho01_re"], p["rho01_im"]) for p in qubits[qid]["points"]])
        plain, offset = detuned_decay_fits(t, r)
        misfit = 2 * offset.cost / (2 * len(t) - 6)
        f_p = stats.f.sf(max(0.0, (2 * plain.cost - 2 * offset.cost) / 2 / misfit), 2, 2 * len(t) - 6)
        rows[qid] = (f"{t[-1]:.1f}", f"{abs(r[-1]):.4f}", f"{abs(r[-1]) / abs(r[0]):.2f}",
                     f"{abs(r[-1]) / Q52_LATE_MEAN:.1f}", f"{abs(r[-1]) / SHOT_SIGMA:.1f}",
                     f"{plain.x[1]:.1f}", f_p, abs(r[0]))
    assert rows[80][:5] == ("54.1", "0.0395", "0.09", "2.1", "7.1")
    assert rows[102][:5] == ("66.0", "0.0438", "0.09", "2.4", "7.9")
    # the echo time is 1.1 to 1.2 of the records' own decay (Q80 in-sample 23.25 us, Q102 29 us)
    assert rows[102][5] == "29.2"
    t2_80, t2_102 = qubits[80]["verdict"]["T2_us"], qubits[102]["verdict"]["T2_us"]
    assert (f"{t2_80 / 23.25:.2f}", f"{t2_102 / 29.2:.2f}") == ("1.16", "1.13")
    # Q52 itself was still turning at 261 us, where its coherence stood at 1.5 times its late mean
    assert (f"{q52.t[7]:.1f}", f"{q52.a[7] / Q52_LATE_MEAN:.1f}", f"{np.degrees(np.angle(q52.r01[7])):+.1f}") == \
        ("261.0", "1.5", "+100.0")
    assert float(rows[80][3]) > 1.5 and float(rows[102][3]) > 1.5
    # a detuned decay plus a constant resolves no offset in either record (scaled by each fit's misfit)
    assert rows[80][6] > 0.1 and rows[102][6] > 0.1


def test_march_simulator_null_is_the_pre_run_synthetic():
    sim = load(MARCH_SIM_JSON)
    assert sim["experiment"] == "shadow_hunt_synthetic" and sim["seed"] == 42
    params = {q["qubit"]: (q["verdict"]["T1_us"], q["verdict"]["T2_us"]) for q in sim["qubit_results"]}
    assert params == {15: (450.0, 36), 80: (350.0, 28)}
    assert not any(q["verdict"]["quadrant_consistent"] for q in sim["qubit_results"])


def test_q80_in_sample_fits_and_their_two_baselines(march):
    qubits, _ = march
    q80 = qubits[80]
    points = q80["points"]
    T2 = q80["verdict"]["T2_us"]
    t = np.array([p["delay_us"] for p in points])
    r = np.array([complex(p["rho01_re"], p["rho01_im"]) for p in points])
    r0 = r[0]
    tl, phases = t[2:], np.unwrap(np.angle(r[2:]))
    m, phi0 = np.linalg.lstsq(np.vstack([tl, np.ones_like(tl)]).T, phases, rcond=None)[0]
    tt, rr = t[1:], r[1:]
    fixed = np.abs(rr - abs(r0) * np.exp(-tt / T2) * np.exp(1j * (m * tt + phi0))).mean()
    intercept = np.abs(rr - abs(r0) * np.exp(-tt / T2) * np.exp(1j * phases.mean())).mean()
    assert (f"{m / (2 * np.pi) * 1e3:+.2f}", f"{fixed:.4f}", f"{intercept:.4f}", f"{intercept / fixed:.1f}") == \
        ("+1.27", "0.0356", "0.0508", "1.4")

    def sse(p):
        T, dw = p
        return 1e10 if T <= 0 else np.sum(np.abs(rr - r0 * np.exp(-tt / T) * np.exp(-1j * dw * tt)) ** 2)

    free = min((minimize(sse, [T0, d0], method="Nelder-Mead") for T0 in (T2, T2 / 2.5, 15, 20)
                for d0 in (-0.01, -0.005, 0.005, 0.01)), key=lambda res: res.fun)
    free_error = np.abs(rr - r0 * np.exp(-tt / free.x[0]) * np.exp(-1j * free.x[1] * tt)).mean()
    envelope_fit = minimize_scalar(lambda T: np.sum(np.abs(rr - r0 * np.exp(-tt / T)) ** 2), bounds=(1, 200), method="bounded")
    envelope_error = np.abs(rr - r0 * np.exp(-tt / envelope_fit.x)).mean()
    calibration_error = np.abs(rr - r0.real * np.exp(-tt / T2)).mean()
    assert (f"{free.x[0]:.2f}", f"{free.x[1] / (2 * np.pi) * 1e3:.2f}", f"{free_error:.4f}") == ("23.25", "-2.58", "0.0138")
    assert (f"{envelope_error:.4f}", f"{envelope_error / free_error:.1f}") == ("0.0487", "3.5")
    assert (f"{calibration_error:.4f}", f"{calibration_error / free_error:.1f}") == ("0.0557", "4.0")
    assert free.x[0] < T2 and f"{T2:.2f}" == "27.06"
    # under rho01 ~ exp(-i dw t) a negative dw is a phase that advances: the same direction as the drift
    assert free.x[1] < 0 and m > 0


# ----------------------------------------------------------------------------------------- the cockpit
def bures(rho, sigma):
    w, v = np.linalg.eigh(rho)
    root = (v * np.sqrt(np.clip(w, 0, None))) @ v.conj().T
    fidelity = np.sum(np.sqrt(np.clip(np.linalg.eigvalsh(root @ sigma @ root), 0, None))) ** 2
    return np.sqrt(max(0.0, 2.0 * (1.0 - np.sqrt(min(1.0, fidelity)))))


def test_cockpit_readings_recompute_independently(q52, march):
    rho = [np.array(r["density_matrix_real"]) + 1j * np.array(r["density_matrix_imag"]) for r in q52.record["raw_tomography"]]
    cpsi = np.array([np.real(np.trace(p @ p)) * 2 * abs(p[0, 1]) for p in rho])
    steps = np.array([bures(rho[k], rho[k - 1]) for k in range(1, len(rho))])
    assert f"{np.corrcoef(steps, cpsi[1:])[0, 1]:.3f}" == "0.954"
    assert f"{q52.a[q52.t > 300].mean():.5f}" == "0.01828" and (q52.t > 300).sum() == 16
    qubits, _ = march
    sim = {q["qubit"]: q for q in load(MARCH_SIM_JSON)["qubit_results"]}
    hw80, sim80 = qubits[80]["points"], sim[80]["points"]
    gaps = [a["cpsi"] - b["cpsi"] for a, b in zip(hw80, sim80)]
    assert sum(g < 0 for g in gaps) == 10 and f"{np.mean(np.abs(gaps)):.4f}" == "0.0338"
    assert (f"{gaps[0]:+.4f}", f"{hw80[0]['cpsi']:.4f}", f"{sim80[0]['cpsi']:.4f}") == ("-0.0463", "0.8175", "0.8639")

    def crossing(points):
        t = [p["delay_us"] for p in points]
        c = [p["cpsi"] for p in points]
        i = next(i for i in range(1, len(c)) if c[i - 1] >= 0.25 > c[i])
        return t[i - 1] + (c[i - 1] - 0.25) / (c[i - 1] - c[i]) * (t[i] - t[i - 1])

    assert (f"{crossing(hw80):.1f}", f"{crossing(sim80):.1f}", f"{crossing(qubits[102]['points']):.1f}") == ("18.5", "23.9", "23.8")
    run3 = load(RUN3_JSON)
    assert (f"{run3['measured_crossing_us']:.2f}", run3["timestamp"][:8]) == ("15.29", "20260318")


def test_cockpit_producer_reproduces_its_committed_output(tmp_path):
    output = tmp_path / "cockpit_validation.txt"
    completed = subprocess.run(
        [sys.executable, str(ROOT / "simulations/cockpit_validation.py")], cwd=tmp_path, capture_output=True,
        text=True, encoding="utf-8", timeout=600, check=False,
        env={**os.environ, "RCPSI_REPO_ROOT": str(ROOT), "RCPSI_COCKPIT_VALIDATION_OUTPUT": str(output),
             "PYTHONIOENCODING": "utf-8", "PYTHONDONTWRITEBYTECODE": "1"})
    assert completed.returncode == 0, completed.stderr[-3000:]
    committed = (ROOT / "simulations/results/cockpit_validation.txt").read_bytes().replace(b"\r\n", b"\n")
    assert output.read_bytes().replace(b"\r\n", b"\n") == committed
    assert sorted(p.name for p in tmp_path.iterdir()) == ["cockpit_validation.txt"]


# ------------------------------------------------------------------------------------ mutation controls
def test_mutations_break_the_gates(q52, march):
    # a late sample turned by 90 degrees breaks the 17/17 direction
    record = copy.deepcopy(q52.record)
    k = int(np.flatnonzero(q52.late)[5])
    re, im = record["raw_tomography"][k]["density_matrix_real"][0][1], record["raw_tomography"][k]["density_matrix_imag"][0][1]
    record["raw_tomography"][k]["density_matrix_real"][0][1], record["raw_tomography"][k]["density_matrix_imag"][0][1] = -im, re
    assert late_direction_counts(Q52(record)) != (17, 17, 17)
    # a late tail that keeps turning at the early rate is no longer static under the drift gate
    record = copy.deepcopy(q52.record)
    for i in np.flatnonzero(q52.late):
        z = q52.r01[i] * np.exp(-1j * 0.035696 * q52.t[i])
        record["raw_tomography"][i]["density_matrix_real"][0][1] = z.real
        record["raw_tomography"][i]["density_matrix_imag"][0][1] = z.imag
    rotating = late_phase_drift(Q52(record))
    assert abs(rotating.slope) > 1.0 and rotating.pvalue < 1e-6
    # late samples turned onto the fixed point's direction no longer sit 4.4 standard errors from it
    record = copy.deepcopy(q52.record)
    r_minus = fixed_point_r_minus(q52)
    turn = np.exp(1j * (np.angle(r_minus) - np.angle(q52.r01[q52.late].mean())))
    for i in np.flatnonzero(q52.late):
        z = q52.r01[i] * turn
        record["raw_tomography"][i]["density_matrix_real"][0][1] = z.real
        record["raw_tomography"][i]["density_matrix_imag"][0][1] = z.imag
    assert abs(late_direction_distance(Q52(record), r_minus)) < 1e-6
    # a stored T2* that is not the proxy fails the proxy gate
    qubits, _ = march
    verdict = dict(qubits[80]["verdict"])
    verdict["T2_star_us"] = 11.01
    assert not stored_t2star_is_proxy(verdict)
    # a Q102 sample moved into Q80's quadrant changes the quadrant list
    q102 = copy.deepcopy(qubits[102])
    q102["points"][2]["rho01_re"], q102["points"][2]["rho01_im"] = 0.1, 0.1
    assert march_quadrants(q102) != [3, 2, 1, 1, 4, 3, 3, 2]
    # an MC value far from its estimate fails the standard-error law
    with pytest.raises(AssertionError):
        assert_mc(0.00859, 1e-5, 0.0090, 5, "mutation")
