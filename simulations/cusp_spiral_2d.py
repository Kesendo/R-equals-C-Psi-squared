#!/usr/bin/env python3
"""Draw a radial quarter reference beside a distinct cardioid cusp.

The named two-qubit Bell+ model supplies finite ``CΨ_com`` trajectories.
They are embedded as drawing coordinates beside the independently defined
period-one cardioid.  The circle ``|CΨ_com|=1/4`` records a selected radial
crossing; only its positive-real point shares the cardioid cusp's coordinate.
The frozen April-16 and April-26 data are overlays, not Mandelbrot iterates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from operator import itemgetter
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from cpsi_complex_plane import main_cardioid, period_2_bulb, trajectory  # noqa: E402
from hardware_cpsi_cplane import extract_trajectory  # noqa: E402


APRIL16_JSON = "data/ibm_cusp_slowing_april2026/cusp_slowing_ibm_kingston_20260416_212042.json"
APRIL16_SHA256 = "7210E6F31C1211F8C1236A4DC6020E569DDCF8F0585CEFA423F7AEC70669AE36"
APRIL26_JSON = "data/ibm_cusp_precision_april2026/cusp_precision_ibm_kingston_20260426_115939.json"
APRIL26_SHA256 = "DAC33BA260594B7262D329C2D17E0629130D478A6C6D94927D467882EBB25211"
CUSP = 0.25
COLORS = ("#CC3333", "#33AACC", "#EE9922", "#AA33CC", "#33AA77")


def load_frozen_json(relative_path, expected_sha256):
    payload = (Path(__file__).parent.parent / relative_path).read_bytes()
    actual_sha256 = hashlib.sha256(payload).hexdigest().upper()
    if actual_sha256 != expected_sha256:
        raise ValueError(f"immutable JSON digest mismatch: {relative_path}")
    parsed = json.loads(payload)
    return parsed


def quarter_circle(n_pts):
    theta = np.linspace(0.0, 2.0 * np.pi, n_pts)
    return CUSP * np.exp(1j * theta)


def first_crossing(c_values):
    magnitudes = np.abs(c_values)
    inside = np.where(magnitudes <= CUSP)[0]
    if len(inside) == 0:
        return None, None
    index = int(inside[0])
    return index, float(np.degrees(np.angle(c_values[index])))


def plot_hardware_arc(axis, c_values, color, label):
    axis.plot(
        c_values.real, c_values.imag, "o-", color=color, lw=1.3,
        markersize=5, alpha=0.75, zorder=3, label=label,
    )
    axis.plot(
        c_values.real[0], c_values.imag[0], "o", mfc="white", mec=color,
        mew=1.5, markersize=10, zorder=4,
    )


def make_static(out_png, april16_data, april26_data):
    gamma = 0.05
    configs = (
        (gamma, 0.0, 0.0, 25.0),
        (gamma, 0.3, 0.0, 25.0),
        (gamma, 0.6, 0.0, 25.0),
        (gamma, 1.0, 0.0, 25.0),
        (gamma, -0.6, 0.0, 25.0),
    )
    trajectories = [trajectory(*configuration) for configuration in configs]

    pair_runs = april16_data["pair_runs"]
    pair_a = extract_trajectory(pair_runs["A_mid"])
    pair_b = extract_trajectory(pair_runs["B_high"])
    hardware_a = np.array(pair_a["cpsi_complex"])
    hardware_b = np.array(pair_b["cpsi_complex"])
    precision_rows = april26_data["cpsi_data"]
    precision_t = np.array(list(map(itemgetter("t_us"), precision_rows)))
    precision_cpsi = np.array(list(map(itemgetter("cpsi"), precision_rows)))

    figure, axes = plt.subplots(1, 3, figsize=(21, 7))
    full_axis, zoom_axis, hardware_axis = axes
    cardioid = main_cardioid()
    bulb = period_2_bulb()
    radial_circle = quarter_circle(400)
    for axis in axes:
        axis.plot(
            cardioid.real, cardioid.imag, "-", color="#999", lw=0.9,
            alpha=0.6, label="period-one cardioid",
        )
        axis.plot(bulb.real, bulb.imag, "-", color="#999", lw=0.6, alpha=0.4)
        axis.plot(
            radial_circle.real, radial_circle.imag, "--", color="red",
            lw=1.6, alpha=0.9, label="selected radial |CΨ|=1/4 circle",
        )
        axis.plot(
            0.25, 0.0, "o", color="red", markersize=7, zorder=6,
            label="cardioid cusp c=+1/4",
        )
        axis.axhline(0.0, color="gray", lw=0.3, alpha=0.5)
        axis.axvline(0.0, color="gray", lw=0.3, alpha=0.5)

    for index, model_trajectory in enumerate(trajectories):
        c_values = np.array(model_trajectory["cpsi_complex"])
        color = COLORS[index % len(COLORS)]
        omega = model_trajectory["detuning_sum"]
        winding = omega / (4.0 * model_trajectory["gamma"])
        label = f"finite model Omega={omega:+.1f}, Omega/(4 gamma)={winding:+.1f}"
        crossing_index, crossing_angle = first_crossing(c_values)
        for axis in (full_axis, zoom_axis):
            axis.plot(c_values.real, c_values.imag, "-", color=color,
                      lw=1.8, alpha=0.85, label=label)
            axis.plot(c_values.real[0], c_values.imag[0], "o",
                      color=color, markersize=7)
            if crossing_index is not None:
                axis.plot(
                    c_values.real[crossing_index], c_values.imag[crossing_index],
                    "*", color=color, markersize=15, markeredgecolor="black",
                    markeredgewidth=0.6, zorder=7,
                )
        if crossing_index is not None:
            print(f"  Omega={omega:+.1f}: finite radial crossing angle={crossing_angle:+.1f} deg")

    precision_zero = np.zeros_like(precision_cpsi)
    hardware_axis.plot(
        precision_cpsi, precision_zero, "o-", color="#2E8B57",
        markersize=6, lw=1.0, alpha=0.9, zorder=5,
        label=f"April-26 saved scalar rows ({len(precision_t)} delays)",
    )
    plot_hardware_arc(
        hardware_axis, hardware_a, "#B22222", "April-16 saved pair A density matrices"
    )
    plot_hardware_arc(
        hardware_axis, hardware_b, "#1F6FB2", "April-16 saved pair B density matrices"
    )

    full_axis.set_xlim(-1.0, 0.6)
    full_axis.set_ylim(-0.6, 0.6)
    full_axis.set_aspect("equal")
    full_axis.set_title("Finite model spirals and separate algebraic references")
    full_axis.grid(True, alpha=0.2)
    full_axis.legend(loc="lower left", fontsize=7)
    full_axis.set_xlabel("Re(c_plot)")
    full_axis.set_ylabel("Im(c_plot)")

    zoom_axis.set_xlim(-0.32, 0.40)
    zoom_axis.set_ylim(-0.32, 0.32)
    zoom_axis.set_aspect("equal")
    zoom_axis.set_title("Finite radial crossings; stars mark sampled first-inside points")
    zoom_axis.grid(True, alpha=0.2)
    zoom_axis.legend(loc="lower left", fontsize=7)
    zoom_axis.set_xlabel("Re(c_plot)")
    zoom_axis.set_ylabel("Im(c_plot)")

    hardware_axis.set_aspect("equal")
    hardware_axis.set_title(
        "selected radial reference |CΨ| = 1/4\n"
        "embedded coordinate: c_plot = CΨ_com"
    )
    hardware_axis.grid(True, alpha=0.2)
    hardware_axis.legend(loc="lower left", fontsize=6.5)
    hardware_axis.set_xlabel("Re(c_plot)")
    hardware_axis.set_ylabel("Im(c_plot)")

    figure.suptitle(
        "Finite CΨ_com paths beside the independently defined period-one cardioid",
        y=1.00,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.95])
    figure.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(figure)


def make_gif(out_gif):
    gamma, omega, phi_0, t_max = 0.05, 0.5, 0.0, 25.0
    model_trajectory = trajectory(gamma, omega, phi_0, t_max, 160)
    c_values = np.array(model_trajectory["cpsi_complex"])
    times = np.array(model_trajectory["times"])
    crossing_index, _ = first_crossing(c_values)

    figure, axis = plt.subplots(figsize=(8, 8))
    cardioid = main_cardioid()
    radial_circle = quarter_circle(400)
    axis.plot(cardioid.real, cardioid.imag, "-", color="#999", lw=0.9, alpha=0.5)
    axis.plot(
        radial_circle.real, radial_circle.imag, "--", color="red", lw=1.6,
        alpha=0.9, label="selected radial |CΨ|=1/4 circle",
    )
    axis.plot(0.25, 0.0, "o", color="red", markersize=6, zorder=6)
    axis.axhline(0.0, color="gray", lw=0.3, alpha=0.5)
    axis.axvline(0.0, color="gray", lw=0.3, alpha=0.5)
    axis.set_xlim(-0.40, 0.42)
    axis.set_ylim(-0.40, 0.40)
    axis.set_aspect("equal")
    axis.grid(True, alpha=0.2)
    axis.set_xlabel("Re(c_plot)")
    axis.set_ylabel("Im(c_plot)")
    axis.set_title(
        "c=+1/4 is the cardioid cusp; the circle is a radial readout"
    )

    line, = axis.plot([], [], "-", color="#AA33CC", lw=2.0, alpha=0.9)
    head = axis.scatter([], [], color="#AA33CC", s=70, edgecolor="black",
                        linewidths=0.6, zorder=7)
    star = axis.scatter([], [], color="gold", s=220, marker="*",
                        edgecolor="black", linewidths=0.8, zorder=8)
    readout_text = axis.text(
        0.02, 0.97, "", transform=axis.transAxes, fontsize=10, va="top",
        family="monospace", bbox=dict(facecolor="white", alpha=0.85,
                                      edgecolor="gray"),
    )
    axis.legend(loc="lower left", fontsize=9)

    def update(frame):
        line.set_data(c_values.real[:frame + 1], c_values.imag[:frame + 1])
        head.set_offsets([[c_values.real[frame], c_values.imag[frame]]])
        if crossing_index is not None and frame >= crossing_index:
            star.set_offsets([[
                c_values.real[crossing_index], c_values.imag[crossing_index]
            ]])
        magnitude = abs(c_values[frame])
        argument = np.degrees(np.angle(c_values[frame]))
        side = "inside radial reference" if magnitude <= CUSP else "outside"
        readout_text.set_text(
            f"t={times[frame]:6.2f}\n|CΨ|={magnitude:6.4f}\n"
            f"arg={argument:+6.1f} deg\n{side}"
        )
        return line, head, star, readout_text

    animation = FuncAnimation(
        figure, update, frames=160, interval=80, blit=False, repeat=True
    )
    animation.save(out_gif, writer=PillowWriter(fps=12.5))
    plt.close(figure)


def main(output_dir) -> None:
    destination = Path(output_dir)
    april16_data = load_frozen_json(APRIL16_JSON, APRIL16_SHA256)
    april26_data = load_frozen_json(APRIL26_JSON, APRIL26_SHA256)
    static_output = destination / "cusp_spiral_2d.png"
    animated_output = destination / "cusp_spiral_2d.gif"
    make_static(static_output, april16_data, april26_data)
    make_gif(animated_output)
    print(f"saved: {static_output}")
    print(f"saved: {animated_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", dest="output_dir", required=True)
    args = parser.parse_args()
    main(args.output_dir)
