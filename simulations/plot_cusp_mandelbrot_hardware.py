#!/usr/bin/env python3
"""Embed the frozen April-26 hardware samples in a Mandelbrot-plane drawing.

The measured scalar samples are mapped to the real plotting coordinate
``c_plot=CΨ``.  They are not iterates of ``z_(n+1)=z_n^2+c``.  The selected
quarter reference on that real line is shown beside, but not identified with,
the cardioid's algebraic cusp at ``c=+1/4``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from operator import itemgetter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


APRIL26_JSON = "data/ibm_cusp_precision_april2026/cusp_precision_ibm_kingston_20260426_115939.json"
APRIL26_SHA256 = "DAC33BA260594B7262D329C2D17E0629130D478A6C6D94927D467882EBB25211"


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


def compute_mandelbrot(x_range, y_range, nx, ny, max_iter):
    x_values = np.linspace(x_range[0], x_range[1], nx)
    y_values = np.linspace(y_range[0], y_range[1], ny)
    x_grid, y_grid = np.meshgrid(x_values, y_values)
    c_grid = x_grid + 1j * y_grid
    escape = np.zeros_like(x_grid, dtype=int)
    z_grid = np.zeros_like(c_grid)
    active = np.ones_like(x_grid, dtype=bool)
    for iteration in range(max_iter):
        z_grid[active] = z_grid[active] ** 2 + c_grid[active]
        escaped = active & (np.abs(z_grid) > 2)
        escape[escaped] = iteration
        active = active & ~escaped
    escape[active] = max_iter
    return escape


def bellplus_readout(times, gamma):
    f = np.exp(-4.0 * gamma * times)
    return f * (1.0 + f * f) / 6.0


def period_one_cardioid(n_pts):
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return (np.exp(1j * theta) / 2.0) - (np.exp(2j * theta) / 4.0)


def make_plot(output_path, x_range, y_range, nx, ny, figsize, zoom,
              t_meas, cpsi_meas, t_theory, cpsi_theory,
              gamma_fit, backend, pair):
    escape = compute_mandelbrot(x_range, y_range, nx, ny, 200)
    figure, axis = plt.subplots(figsize=figsize)
    axis.imshow(
        np.log(escape.astype(float) + 1.0),
        extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
        origin="lower", cmap="bone", aspect="equal", interpolation="bilinear",
    )

    axis.plot(
        cpsi_theory, np.zeros_like(cpsi_theory), color="lightblue",
        linewidth=1.0, alpha=0.7,
        label=f"finite F25 model at fitted gamma={gamma_fit * 1000.0:.2f}/ms",
    )
    samples = axis.scatter(
        cpsi_meas, np.zeros_like(cpsi_meas), c=t_meas, cmap="viridis",
        s=80, zorder=5, edgecolor="black", linewidth=0.6,
        label=f"saved hardware samples ({backend}, qubits {pair})",
    )
    axis.plot(
        0.25, 0.0, "X", color="lime", markersize=18, mew=2.5,
        label="selected radial CΨ=1/4 reference", zorder=6,
    )
    axis.plot(
        cpsi_theory[0], 0.0, "o", color="red", markersize=12,
        label=f"finite model t={t_theory[0]:.0f}", zorder=6,
    )
    axis.plot(
        cpsi_theory[-1], 0.0, "s", color="cyan", markersize=10,
        label=f"finite plotted endpoint t={t_theory[-1]:.1f}", zorder=6,
    )

    cardioid = period_one_cardioid(1000)
    axis.plot(
        cardioid.real, cardioid.imag, "--", color="lime", alpha=0.4,
        linewidth=1.0, label="period-one cardioid boundary",
    )
    if zoom:
        axis.set_title(
            "quarter reference; cardioid cusp is a separate algebraic marker"
        )
    else:
        axis.set_title("hardware samples embedded as c_plot = CΨ")
    axis.set_xlabel("Re(c_plot)")
    axis.set_ylabel("Im(c_plot)")
    axis.legend(loc="upper left", fontsize=9)
    axis.grid(alpha=0.2)
    figure.colorbar(samples, ax=axis, label="t [microseconds]", shrink=0.7)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def main(json_path, output_dir) -> None:
    destination = Path(output_dir)
    require_exact_json_path(json_path, APRIL26_JSON)
    data = load_frozen_json(json_path, APRIL26_SHA256)
    cpsi_rows = data["cpsi_data"]
    t_meas = np.array(list(map(itemgetter("t_us"), cpsi_rows)))
    cpsi_meas = np.array(list(map(itemgetter("cpsi"), cpsi_rows)))
    gamma_fit = data["gamma_fit"]
    backend = data["backend"]
    pair = data["pair"]["qubits"]
    t_theory = np.linspace(0.0, 8.0, 2000)
    cpsi_theory = bellplus_readout(t_theory, gamma_fit)

    full_output = destination / "bellplus_trajectory_on_mandelbrot_hardware.png"
    zoom_output = destination / "bellplus_trajectory_on_mandelbrot_hardware_zoom.png"
    make_plot(
        full_output, (-2.2, 0.8), (-1.2, 1.2), 1200, 900, (14, 10), False,
        t_meas, cpsi_meas, t_theory, cpsi_theory, gamma_fit, backend, pair,
    )
    make_plot(
        zoom_output, (-0.10, 0.50), (-0.20, 0.20), 1200, 800, (12, 6), True,
        t_meas, cpsi_meas, t_theory, cpsi_theory, gamma_fit, backend, pair,
    )
    print(f"saved: {full_output}")
    print(f"saved: {zoom_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_path")
    parser.add_argument("--output-dir", dest="output_dir", required=True)
    args = parser.parse_args()
    main(args.json_path, args.output_dir)
