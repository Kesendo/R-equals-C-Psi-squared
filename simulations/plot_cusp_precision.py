#!/usr/bin/env python3
"""Plot the frozen April-26 CΨ(t) samples with finite F25 comparisons.

The two figures preserve the measured scalar rows, the fitted curve, the shot-
noise guide and the contemporaneous curated T2-echo input.  The displayed
``CΨ=1/4`` is a selected reference level, not a measured Mandelbrot cusp or
phase transition.  The curated T2-echo value is explicit CLI input because its
underlying calibration snapshot is not stored with this data set.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from operator import itemgetter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


APRIL26_JSON = "data/ibm_cusp_precision_april2026/cusp_precision_ibm_kingston_20260426_115939.json"
APRIL26_SHA256 = "DAC33BA260594B7262D329C2D17E0629130D478A6C6D94927D467882EBB25211"
REQUIRED_GAMMA_T2ECHO_PER_US = 0.00165


def require_exact_json_path(json_path, expected_relative):
    repo_root = Path(__file__).parent.parent
    supplied = Path(json_path)
    resolved = supplied if supplied.is_absolute() else repo_root / supplied
    expected = repo_root / expected_relative
    if resolved.resolve() != expected.resolve():
        raise ValueError(f"expected immutable input {expected_relative}")


def load_frozen_json(relative_path, expected_sha256):
    payload = (Path(__file__).parent.parent / relative_path).read_bytes()
    actual_sha256 = hashlib.sha256(payload).hexdigest().upper()
    if actual_sha256 != expected_sha256:
        raise ValueError(f"immutable JSON digest mismatch: {relative_path}")
    parsed = json.loads(payload)
    return parsed


def cpsi_bellplus(t_us, gamma):
    f = np.exp(-4.0 * gamma * t_us)
    return f * (1.0 + f * f) / 6.0


def render_trajectory(output_path, t, cpsi, sigma, t_fine,
                      declared_curve, fitted_curve, t_cross_declared,
                      t_cross_fit, declared_gamma, gamma_fit, backend, pair):
    figure, axis = plt.subplots(figsize=(8, 5.5))
    axis.axhline(
        0.25, color="gray", linestyle=":", linewidth=1.0,
        label="selected scalar reference CΨ=1/4",
    )
    axis.axvline(
        t_cross_declared, color="C1", linestyle=":", alpha=0.5,
        label=f"declared-model reference time={t_cross_declared:.1f} us",
    )
    axis.axvline(
        t_cross_fit, color="C2", linestyle=":", alpha=0.7,
        label=f"fitted-model reference time={t_cross_fit:.2f} us",
    )
    axis.plot(
        t_fine, declared_curve, color="C1", linestyle="--", linewidth=1.2,
        alpha=0.6, label=f"F25 at declared gamma={declared_gamma:.5f}/us",
    )
    axis.plot(
        t_fine, fitted_curve, color="C2", linewidth=1.6,
        label=f"F25 at in-situ fitted gamma={gamma_fit:.5f}/us",
    )
    axis.errorbar(
        t, cpsi, yerr=sigma, fmt="o", color="black", markersize=5,
        capsize=2, elinewidth=0.8,
        label=f"saved hardware samples ({backend}, qubits {pair})",
    )
    axis.set_xlabel("t [microseconds]")
    axis.set_ylabel("CΨ")
    axis.set_title("quarter reference level CΨ=1/4")
    axis.legend(loc="upper right", fontsize=8.5)
    axis.grid(alpha=0.3)
    axis.set_xlim(left=0)
    axis.set_ylim(bottom=0)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def render_residuals(output_path, t, residual_declared, residual_fit,
                     sigma, shots, t_cross_fit, declared_gamma,
                     backend, pair):
    figure, axis = plt.subplots(figsize=(8, 4.5))
    axis.axhline(0.0, color="gray", linewidth=0.8, alpha=0.6)
    axis.axhspan(
        -sigma, sigma, color="gray", alpha=0.15,
        label=f"shot-noise guide: 0.5/sqrt({shots})={sigma:.4f}",
    )
    axis.plot(
        t, residual_declared, "o--", color="C1", markersize=5,
        label="residual against declared T2-echo comparison",
    )
    axis.plot(
        t, residual_fit, "o-", color="C2", markersize=5,
        label="residual against in-situ fit",
    )
    axis.axvline(
        t_cross_fit, color="C2", linestyle=":", alpha=0.5,
        label="fitted-model quarter-reference time",
    )
    axis.set_xlabel("t [microseconds]")
    axis.set_ylabel("measured CΨ - finite F25 model")
    axis.set_title(
        f"declared T2-echo input: {declared_gamma:.5f} /us\n"
        f"pointwise residuals on {backend}, qubits {pair}"
    )
    axis.legend(loc="best", fontsize=9)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def main(json_path, gamma_t2echo_per_us, output_dir) -> None:
    destination = Path(output_dir)
    declared_gamma = gamma_t2echo_per_us
    if gamma_t2echo_per_us != REQUIRED_GAMMA_T2ECHO_PER_US:
        raise ValueError("the frozen plot contract requires 0.00165 /us")
    require_exact_json_path(json_path, APRIL26_JSON)
    data = load_frozen_json(json_path, APRIL26_SHA256)
    cpsi_rows = data["cpsi_data"]
    t = np.array(list(map(itemgetter("t_us"), cpsi_rows)))
    cpsi = np.array(list(map(itemgetter("cpsi"), cpsi_rows)))
    gamma_fit = data["gamma_fit"]
    K_CROSS = data["K_CROSS"]
    backend = data["backend"]
    pair = data["pair"]["qubits"]
    shots = data["shots"]

    t_fine = np.linspace(0.0, t.max() * 1.05, 600)
    declared_curve = cpsi_bellplus(t_fine, declared_gamma)
    sampled_declared = cpsi_bellplus(t, declared_gamma)
    fitted_curve = cpsi_bellplus(t_fine, gamma_fit)
    sampled_fit = cpsi_bellplus(t, gamma_fit)
    residual_declared = cpsi - sampled_declared
    residual_fit = cpsi - sampled_fit
    sigma = 0.5 / math.sqrt(shots)
    t_cross_declared = K_CROSS / declared_gamma
    t_cross_fit = K_CROSS / gamma_fit

    trajectory_output = destination / "cusp_precision_trajectory.png"
    residual_output = destination / "cusp_precision_residuals.png"
    render_trajectory(
        trajectory_output, t, cpsi, sigma, t_fine,
        declared_curve, fitted_curve, t_cross_declared, t_cross_fit,
        declared_gamma, gamma_fit, backend, pair,
    )
    render_residuals(
        residual_output, t, residual_declared, residual_fit, sigma, shots,
        t_cross_fit, declared_gamma, backend, pair,
    )
    print(f"saved: {trajectory_output}")
    print(f"saved: {residual_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path")
    parser.add_argument(
        "--gamma-t2echo-per-us", dest="gamma_t2echo_per_us",
        required=True, type=float,
    )
    parser.add_argument("--output-dir", dest="output_dir", required=True)
    args = parser.parse_args()
    main(args.json_path, args.gamma_t2echo_per_us, args.output_dir)
