"""Finite-grid and normalization controls for the raw frequency SFF report."""
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


def helper():
    path = Path(__file__).resolve().parents[1] / "sff_window_summary.py"
    assert path.exists(), "protected SFF window-summary helper is missing"
    spec = importlib.util.spec_from_file_location("sff_window_summary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert callable(getattr(module, "sff_frequency", None)), "protected actual SFF kernel must be import-safe"
    assert callable(getattr(module, "raw_multiset_density_scale", None)), "protected raw multiset density scale is missing"
    return module


def test_reference_tracks_mode_count_not_unit_plateau():
    h = helper()
    assert h.independent_phase_reference(8) == 0.125  # wrong normalization 1 fails
    assert h.independent_phase_reference(1) == 1.0  # zero/constant-1/8 control
    with pytest.raises(ValueError):
        h.independent_phase_reference(0)


def test_unreached_windows_are_missing_not_measured_zero():
    h = helper()
    result = h.summarize_windows(np.arange(0, 201, 50), np.zeros(5), 5810.10)
    assert result["below"] == 0.0
    assert result["intermediate"] is None
    assert result["beyond"] is None
    assert result["slope"] is None
    assert h.format_sample(result["beyond"]) == "not sampled"
    assert h.format_sample(result["below"]) == "0.0000"  # same zeros, reached window


def test_reached_windows_and_slope_use_only_their_samples():
    h = helper()
    times = np.arange(21, dtype=float)
    values = times / 50
    result = h.summarize_windows(times, values, 10.0)
    assert result["below"] == 0.0
    assert result["intermediate"] == pytest.approx(0.11, abs=1e-15)
    assert result["beyond"] == pytest.approx(0.31, abs=1e-15)
    assert result["slope"] == pytest.approx(0.2, abs=1e-15)
    # A genuine flat reached window reports a measured zero slope, not missing.
    flat = h.summarize_windows(times, np.zeros_like(times), 10.0)
    assert flat["slope"] == 0.0


def test_sparse_intermediate_bin_has_mean_but_no_fit():
    h = helper()
    result = h.summarize_windows(np.array([0., 5., 15.]), np.array([1., 2., 3.]), 10.0)
    assert result["intermediate"] == 2.0
    assert result["beyond"] == 3.0
    assert result["slope"] is None  # one point cannot establish a slope


def test_empty_below_bin_has_no_minimum():
    h = helper()
    result = h.summarize_windows(np.array([20., 30.]), np.array([4., 5.]), 10.0)
    assert result["below"] is None
    assert result["minimum"] is None  # no fallback to K[0]
    assert result["beyond"] == 4.5


def test_actual_kernel_retains_m_squared_normalization():
    # Four roots of unity: the four equally weighted discrete time samples
    # average to 1/M. A /M kernel or unit reference fails; t=0 pins K(0)=1.
    values = helper().sff_frequency(1j * np.arange(4), np.arange(4) * np.pi / 2)
    assert values[0] == 1.0
    assert np.mean(values) == pytest.approx(helper().independent_phase_reference(4), abs=1e-15)
    assert np.max(values[1:]) < 1e-28


def test_duplication_preserves_sff_but_changes_raw_density_scale():
    h = helper()
    eigenvalues = 1j * np.array([1., 2.])
    duplicated = np.repeat(eigenvalues, 2)
    times = np.linspace(0, 20, 37)
    np.testing.assert_allclose(h.sff_frequency(eigenvalues, times), h.sff_frequency(duplicated, times), atol=1e-15, rtol=0)
    assert h.raw_multiset_density_scale(eigenvalues) == pytest.approx(2 * np.pi, abs=1e-14)
    assert h.raw_multiset_density_scale(duplicated) == pytest.approx(6 * np.pi, abs=1e-14)
    assert h.raw_multiset_density_scale(duplicated) != h.raw_multiset_density_scale(eigenvalues)
    # No fictitious default scale when the nonzero absolute-frequency range vanishes.
    with pytest.raises(ValueError):
        h.raw_multiset_density_scale(1j * np.array([1., 1.]))


def test_producer_temp_output_uses_actual_helper_and_csv_windows(tmp_path):
    root = Path(__file__).resolve().parents[2]
    producer = root / "simulations/spectral_form_factor.py"
    source = producer.read_text(encoding="utf-8")
    assert "RCPSI_SFF_OUTPUT_PATH" in source, "protected temporary-output seam is missing; do not run against tracked output"
    assert "def sff_frequency(" not in source, "producer must consume the actual import-safe kernel"
    tracked = root / "simulations/results/spectral_form_factor.txt"
    before, before_mtime = tracked.read_bytes(), tracked.stat().st_mtime_ns
    destination = tmp_path / "sff.txt"
    run = subprocess.run([sys.executable, str(producer)], cwd=root, env={**os.environ, "RCPSI_SFF_OUTPUT_PATH": str(destination)},
                         capture_output=True, text=True, encoding="utf-8", timeout=60)
    assert run.returncode == 0, run.stderr
    assert tracked.read_bytes() == before and tracked.stat().st_mtime_ns == before_mtime
    report = destination.read_text(encoding="utf-8")
    assert "raw multiset density scale" in report and "multiplicity-dependent" in report
    assert "Only the raw oscillation-frequency SFF is computed" in report
    assert "no decay-weighted or connected estimator is defined" in report
    assert "K_diss" not in report and "overflow/NaN" not in report
    assert "Heisenberg" not in report and "t_H" not in report and "t_Th" not in report
    assert "first sampled exceedance" not in report, "the t=0 exceedance diagnostic is vacuous and must be absent"
    assert "2*cos" in report and "doubled and cross frequencies" in report
    phase4 = report.split("PHASE 4:", 1)[1].split("PHASE 5:", 1)[0]
    blocks = dict(re.findall(r"  N=(\d+):\n(.*?)(?=\n  N=|\Z)", phase4, re.S))
    assert set(blocks) == {"3", "4", "5", "6", "7"}
    h = helper()
    for n in range(3, 8):
        csv = np.loadtxt(root / f"simulations/results/rmt_eigenvalues_N{n}.csv", delimiter="\t", skiprows=1)
        eigenvalues = csv[:, 0] + 1j * csv[:, 1]
        scale = h.raw_multiset_density_scale(eigenvalues)
        # Independent telescoping denominator checks the actual multiset,
        # including duplicate absolute frequencies; no deduplication/unfolding.
        frequencies = np.abs(eigenvalues.imag[np.abs(eigenvalues.imag) > 1e-10])
        expected = 2 * np.pi * (len(frequencies) - 1) / np.ptp(frequencies)
        assert scale == pytest.approx(expected, rel=1e-14)
        period = 2 * np.pi / (4 * (1 - np.cos(np.pi / n)))
        end = min(3 * scale, 50 * period, 200)
        times = np.linspace(0, end, min(2000, max(500, int(end * 50))))
        windows = h.summarize_windows(times, h.sff_frequency(eigenvalues, times), scale)
        block = blocks[str(n)]
        assert f"raw scale = {scale:.2f}" in block
        assert f"reference 1/M = {1 / len(eigenvalues):.8g}" in block
        for label, key in (("Below", "below"), ("Intermediate", "intermediate"), ("Beyond", "beyond")):
            line = next(line for line in block.splitlines() if line.strip().startswith(label + " ("))
            assert f"<K> = {h.format_sample(windows[key])}" in line
        if n >= 5:
            assert windows["beyond"] is None
        if n == 7:
            assert windows["intermediate"] is None and windows["slope"] is None
            assert "slope = not sampled" in block and "not classifiable" in block


def test_zero_frequency_sff_is_constant_and_pair_square_has_doubled_frequency():
    h = helper()
    times = np.array([0., np.pi / 2, np.pi])
    zero_frequencies = np.full(6, -2.0, dtype=complex)
    np.testing.assert_array_equal(h.sff_frequency(zero_frequencies, times), np.ones(3))
    np.testing.assert_array_equal(6**2 * h.sff_frequency(zero_frequencies, times), np.full(3, 36.0))
    pair = h.sff_frequency(np.array([-1j, 1j]), times)
    np.testing.assert_allclose(pair, np.cos(times)**2, rtol=0, atol=1e-15)
    assert abs(pair[-1] - np.cos(times[-1])) > 1.9  # trace-amplitude vs SFF control


def test_producer_has_no_misdefined_decay_or_vacuous_connected_kernel():
    source = (Path(__file__).resolve().parents[1] / "spectral_form_factor.py").read_text(encoding="utf-8")
    document = (Path(__file__).resolve().parents[2] / "experiments" / "SPECTRAL_FORM_FACTOR.md").read_text(encoding="utf-8")
    for stale in ("def sff_dissipative", "def sff_connected", "K_diss", "overflow/NaN"):
        assert stale not in source + document
    assert "physical Liouvillian propagation would involve" in document
    assert "`exp(lambda*t)`" in document
    assert "requires a specified ensemble or averaging prescription" in document


def test_f2_frequency_is_block_reference_not_full_spectrum_minimum():
    root = Path(__file__).resolve().parents[2]
    source = (root / "simulations" / "spectral_form_factor.py").read_text(encoding="utf-8")
    artifact = (root / "simulations" / "results" / "spectral_form_factor.txt").read_text(encoding="utf-8")
    document = (root / "experiments" / "SPECTRAL_FORM_FACTOR.md").read_text(encoding="utf-8")
    combined = source + artifact + document

    assert "(0,1)-block k=1 reference" in combined
    assert "ω_min (slowest)" not in combined
    assert "slowest pair's trace-amplitude term" not in combined

    # The N=5 full spectrum itself supplies a stable counterexample to the
    # old global-minimum label: resolved nonzero frequencies below the F2 k=1
    # block value occur in the same committed CSV.
    data = np.loadtxt(root / "simulations" / "results" / "rmt_eigenvalues_N5.csv",
                      delimiter="\t", skiprows=1)
    resolved = np.abs(data[:, 1])
    f2_reference = 4 * (1 - np.cos(np.pi / 5))
    assert np.any((resolved > 0.1) & (resolved < 0.9 * f2_reference))
