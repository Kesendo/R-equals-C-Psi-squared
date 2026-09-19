#!/usr/bin/env python3
"""Embed a finite Bell+ readout beside Mandelbrot iteration geometry.

The real scalar

    CΨ(t) = f(t) (1 + f(t)^2) / 6,  f(t) = exp(-4 gamma t),

is plotted as the chosen coordinate ``c_plot=CΨ(t)``.  This is an explicit
comparison map, not a claim that the Lindblad trajectory obeys
``z_(n+1)=z_n^2+c``.  Its selected radial value CΨ=1/4 and the cardioid cusp
c=+1/4 are different objects with the same displayed number.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


def compute_mandelbrot(x_range, y_range, nx, ny, max_iter=200):
    x = np.linspace(x_range[0], x_range[1], nx)
    y = np.linspace(y_range[0], y_range[1], ny)
    x_grid, y_grid = np.meshgrid(x, y)
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


def render_view(output_path, x_range, y_range, nx, ny,
                times, readout, zoom):
    escape = compute_mandelbrot(x_range, y_range, nx, ny)
    log_escape = np.log(escape.astype(float) + 1.0)
    figsize = (10, 8) if zoom else (14, 10)
    figure, axis = plt.subplots(figsize=figsize)
    axis.imshow(
        log_escape,
        extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
        cmap="bone", origin="lower", aspect="equal", interpolation="bilinear",
    )

    selected = ((readout >= x_range[0]) & (readout <= x_range[1]))
    shown_readout = readout[selected]
    shown_times = times[selected]
    points = np.array([shown_readout, np.zeros_like(shown_readout)]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    normalization = plt.Normalize(times.min(), times.max())
    collection = LineCollection(
        segments, cmap="plasma", norm=normalization,
        linewidths=4.0 if zoom else 3.5, zorder=5,
    )
    collection.set_array(shown_times[:-1])
    axis.add_collection(collection)
    figure.colorbar(collection, ax=axis, label="time t (gamma=1)", shrink=0.6)

    axis.plot(
        0.25, 0.0, "x", color="lime", markersize=18 if zoom else 14,
        markeredgewidth=4 if zoom else 3, zorder=10,
        label="selected CΨ=1/4 radial reference",
    )
    axis.plot(
        readout[0], 0.0, "o", color="red", markersize=10, zorder=10,
        label=f"t=0: CΨ={readout[0]:.3f}",
    )

    cardioid = period_one_cardioid(1000)
    axis.plot(
        cardioid.real, cardioid.imag, "--", color="lime", alpha=0.4,
        linewidth=1.0, label="period-one cardioid boundary",
    )
    if zoom:
        axis.set_title(
            "quarter reference CΨ=1/4; cardioid cusp is algebraic"
        )
    else:
        axis.set_title(
            "embedded comparison: c_plot = CΨ(t)\n"
            "iteration interior / escape exterior"
        )
    axis.set_xlabel("Re(c_plot)")
    axis.set_ylabel("Im(c_plot)")
    axis.legend(loc="upper left", fontsize=9)
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def main(output_dir) -> None:
    destination = Path(output_dir)
    full_output = destination / "bellplus_trajectory_on_mandelbrot.png"
    zoom_output = destination / "bellplus_trajectory_on_mandelbrot_zoom.png"
    gamma = 1.0
    times = np.linspace(0.0, 5.0, 2000)
    readout = bellplus_readout(times, gamma)
    render_view(
        full_output, (-2.2, 0.8), (-1.2, 1.2), 1200, 900,
        times, readout, False,
    )
    render_view(
        zoom_output, (-0.1, 0.5), (-0.3, 0.3), 800, 600,
        times, readout, True,
    )
    print(f"saved: {full_output}")
    print(f"saved: {zoom_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", dest="output_dir", required=True)
    args = parser.parse_args()
    main(args.output_dir)
