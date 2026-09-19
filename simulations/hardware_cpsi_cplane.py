#!/usr/bin/env python3
"""
Hardware CΨ in the complex plane: a finite readout from the saved density
matrices of the 2026-04-16 run on ibm_kingston.

The cusp-slowing JSON stores the full 4×4 density matrix for each delay
point via rho2_real + rho2_imag. We can therefore compute complex
CΨ_com = C · (Σ ρ_{ij} off-diagonal, signed) / (d-1) without a new QPU
run. The plotted CΨ_com is a basis-chosen one-coherence readout embedded in a
copy of the Mandelbrot parameter plane. It is not an iteration orbit, and its
radial quarter reference is not the cardioid cusp.
The resolved phase change is only this sparse, basis-dependent coordinate
readout, not a standalone detuning or calibration estimate.

Plots:
  (1) c-plane trajectories for both hardware pairs with Mandelbrot cardioid
  (2) Complex CΨ evolution vs real CΨ (does the phase rotate?)
  (3) Zoom to the cusp region

Date: 2026-04-16
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

# ── Paths ─────────────────────────────────────────────────────────
APRIL16_JSON = "data/ibm_cusp_slowing_april2026/cusp_slowing_ibm_kingston_20260416_212042.json"
APRIL16_SHA256 = "7210E6F31C1211F8C1236A4DC6020E569DDCF8F0585CEFA423F7AEC70669AE36"


def load_frozen_json(relative_path, expected_sha256):
    payload = (Path(__file__).parent.parent / relative_path).read_bytes()
    actual_sha256 = hashlib.sha256(payload).hexdigest().upper()
    if actual_sha256 != expected_sha256:
        raise ValueError(f"immutable JSON digest mismatch: {relative_path}")
    parsed = json.loads(payload)
    return parsed


# ── CΨ metrics ────────────────────────────────────────────────────
def cpsi_real_from_rho(rho: np.ndarray) -> float:
    d = rho.shape[0]
    C = float(np.real(np.trace(rho @ rho)))
    diag = np.diag(np.diag(rho))
    L1 = float(np.sum(np.abs(rho - diag)))
    return C * L1 / (d - 1)


def cpsi_complex_from_rho(rho: np.ndarray) -> complex:
    d = rho.shape[0]
    C = float(np.real(np.trace(rho @ rho)))
    # Sum of upper-triangle off-diagonals (signed/phase-aware), ×2 for the
    # lower triangle counterparts. /(d−1) for the normalization.
    mask_upper = np.triu(np.ones_like(rho, dtype=bool), k=1)
    Psi_com = complex(2.0 * np.sum(rho[mask_upper]) / (d - 1))
    return C * Psi_com


# ── Mandelbrot references ─────────────────────────────────────────
def main_cardioid(n_pts: int = 400) -> np.ndarray:
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return (np.exp(1j * theta) / 2.0) - (np.exp(2j * theta) / 4.0)


def period_2_bulb(n_pts: int = 200) -> np.ndarray:
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return -1.0 + 0.25 * np.exp(1j * theta)


# ── Analysis ──────────────────────────────────────────────────────
def extract_trajectory(pair_data: dict) -> dict:
    points = [point for point in pair_data["trajectory"]
              if "rho2_real" in point]
    matrices = [
        np.array(point["rho2_real"], dtype=complex)
        + 1j * np.array(point["rho2_imag"], dtype=complex)
        for point in points
    ]
    return {"t_us": [point["t_us"] for point in points],
            "cpsi_complex": [cpsi_complex_from_rho(rho) for rho in matrices],
            "cpsi_real": [cpsi_real_from_rho(rho) for rho in matrices],
            "pair": pair_data.get("pair", {})}


def plot_all(pair_a: dict, pair_b: dict, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    ax_full, ax_zoom, ax_phase = axes

    c_card = main_cardioid()
    c_bulb = period_2_bulb()

    for ax in (ax_full, ax_zoom):
        ax.plot(c_card.real, c_card.imag, "-", color="#888", linewidth=0.9,
                alpha=0.6, label="Mandelbrot cardioid")
        ax.plot(c_bulb.real, c_bulb.imag, "-", color="#888", linewidth=0.5,
                alpha=0.4, label="period-2 bulb")
        ax.plot(0.25, 0, "o", color="red", markersize=10, zorder=5,
                label="cardioid cusp c=+1/4 (algebraic reference)")
        ax.axhline(0, color="gray", linewidth=0.3, alpha=0.4)
        ax.axvline(0, color="gray", linewidth=0.3, alpha=0.4)

    colors = {"A_mid": "#CC5533", "B_high": "#3355CC"}

    for label, tr in [("A_mid", pair_a), ("B_high", pair_b)]:
        c_vals = np.array(tr["cpsi_complex"])
        col = colors[label]
        pair = tr["pair"]
        qubits = pair.get("qubits", ["?"])
        label_str = f"{label}: qubits {qubits}"
        for ax in (ax_full, ax_zoom):
            ax.plot(c_vals.real, c_vals.imag, "o-", color=col,
                    linewidth=2, markersize=7, alpha=0.9, label=label_str)
            ax.plot(c_vals.real[0], c_vals.imag[0], "s", color=col,
                    markersize=11, markerfacecolor="white", markeredgewidth=2,
                    zorder=6, label=f"{label}: t=0")
            ax.plot(c_vals.real[-1], c_vals.imag[-1], "x", color=col,
                    markersize=12, markeredgewidth=3, zorder=6)

    ax_full.set_xlim(-0.05, 0.40)
    ax_full.set_ylim(-0.15, 0.15)
    ax_full.set_aspect("equal")
    ax_full.set_title("Hardware c-plane trajectories\n(Mandelbrot cardioid context)")
    ax_full.grid(True, alpha=0.2)
    ax_full.legend(loc="upper left", fontsize=7)
    ax_full.set_xlabel("Re(CΨ_com)")
    ax_full.set_ylabel("Im(CΨ_com)")

    # Zoom around the shared numerical value, without identifying the objects.
    ax_zoom.set_xlim(0.15, 0.30)
    ax_zoom.set_ylim(-0.04, 0.04)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title("Radial hardware samples beside the algebraic cusp")
    ax_zoom.grid(True, alpha=0.2)
    ax_zoom.set_xlabel("Re(CΨ_com)")
    ax_zoom.set_ylabel("Im(CΨ_com)")
    ax_zoom.legend(loc="upper right", fontsize=7)

    # Third panel: phase and magnitude of CΨ_com vs t
    for label, tr in [("A_mid", pair_a), ("B_high", pair_b)]:
        c_vals = np.array(tr["cpsi_complex"])
        ts = np.array(tr["t_us"])
        col = colors[label]
        ax_phase.plot(ts, np.abs(c_vals), "o-", color=col,
                      linewidth=1.5, markersize=6,
                      label=f"{label}: |CΨ_com|")
        # Small arg(CΨ_com) trace, rescaled onto the same axis via twin
    ax_phase_twin = ax_phase.twinx()
    for label, tr in [("A_mid", pair_a), ("B_high", pair_b)]:
        c_vals = np.array(tr["cpsi_complex"])
        ts = np.array(tr["t_us"])
        col = colors[label]
        args = np.degrees(np.angle(c_vals))
        ax_phase_twin.plot(ts, args, "--", color=col, alpha=0.6, linewidth=1,
                            label=f"{label}: arg(CΨ_com)")
    ax_phase.axhline(0.25, color="red", linestyle=":", linewidth=1.2,
                     alpha=0.7, label="|CΨ| = 1/4 (radial reference)")
    ax_phase.set_xlabel("delay t (μs)")
    ax_phase.set_ylabel("|CΨ_com|", color="#444")
    ax_phase_twin.set_ylabel("arg(CΨ_com) [°]", color="#888")
    ax_phase.set_title(
        "|CΨ_com|(t) and arg(CΨ_com)(t)\n"
        "not standalone detuning/calibration evidence"
    )
    ax_phase.grid(True, alpha=0.2)
    ax_phase.legend(loc="upper right", fontsize=7)
    ax_phase_twin.legend(loc="center right", fontsize=7)

    fig.suptitle(
        "saved April-16 density matrices\n"
        "radial quarter reference; cardioid cusp is separate\n"
        "embedded coordinate c_plot = CΨ_com",
        y=1.00,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_png}")


# ── Main ──────────────────────────────────────────────────────────
def main(output_path) -> None:
    destination = Path(output_path)
    print("=" * 70)
    print("  Hardware CΨ in the complex plane (c-plane extension)")
    print("=" * 70)

    data = load_frozen_json(APRIL16_JSON, APRIL16_SHA256)
    pair_runs = data["pair_runs"]
    pair_a = extract_trajectory(pair_runs["A_mid"])
    pair_b = extract_trajectory(pair_runs["B_high"])

    print(f"\n  Pair A (mid-T2): {len(pair_a['t_us'])} delay points")
    print(f"\n  Pair B (high-T2): {len(pair_b['t_us'])} delay points")

    plot_all(pair_a, pair_b, destination)

    print()
    print("Interpretation:")
    print("  The small resolved phase changes belong to this sparse,")
    print("  basis-dependent coordinate readout.")
    print("  They do not by themselves estimate residual detuning or calibration.")
    print("  The red marker is the separate algebraic cusp c=+1/4.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", dest="output_path", required=True)
    args = parser.parse_args()
    main(args.output_path)
