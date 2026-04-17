#!/usr/bin/env python3
"""
Hardware CΨ in the complex plane: extracts complex CΨ_com from the
saved density matrices of the 2026-04-16 cusp-slowing run on ibm_kingston.

The cusp-slowing JSON stores the full 4×4 density matrix for each delay
point via rho2_real + rho2_imag. We can therefore compute complex
CΨ_com = C · (Σ ρ_{ij} off-diagonal, signed) / (d-1) without a new QPU
run: the 2D trajectory is already in the data if Kingston's natural
detuning on qubits 124-125 / 14-15 produced any phase rotation.

Plots:
  (1) c-plane trajectories for both hardware pairs with Mandelbrot cardioid
  (2) Complex CΨ evolution vs real CΨ (does the phase rotate?)
  (3) Zoom to the cusp region

Date: 2026-04-16
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ── Paths ─────────────────────────────────────────────────────────
# Public repo: reads hardware JSON from data/ibm_cusp_slowing_april2026/
# (self-contained, no dependency on the private run-script location).
REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / "data" / "ibm_cusp_slowing_april2026"
PUBLIC_RESULTS = Path(__file__).parent / "results"
PUBLIC_RESULTS.mkdir(exist_ok=True)
OUT = PUBLIC_RESULTS / "hardware_cpsi_cplane.png"


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
    traj = pair_data.get("trajectory", [])
    ts, cpsi_c, cpsi_r = [], [], []
    for p in traj:
        if "rho2_real" not in p:
            continue
        rho = (np.array(p["rho2_real"], dtype=complex)
               + 1j * np.array(p["rho2_imag"], dtype=complex))
        ts.append(p["t_us"])
        cpsi_c.append(cpsi_complex_from_rho(rho))
        cpsi_r.append(cpsi_real_from_rho(rho))
    return {"t_us": ts, "cpsi_complex": cpsi_c, "cpsi_real": cpsi_r,
            "pair": pair_data.get("pair", {})}


def plot_all(pair_a: dict, pair_b: dict, freeze_sim: dict | None,
             out_png: Path) -> None:
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
                label="cusp c = 1/4")
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

    # Zoom to cusp region
    ax_zoom.set_xlim(0.15, 0.30)
    ax_zoom.set_ylim(-0.04, 0.04)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title("Zoom to the cusp c = 1/4 (fold crossing)")
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
                     alpha=0.7, label="|CΨ| = 1/4 (fold)")
    ax_phase.set_xlabel("delay t (μs)")
    ax_phase.set_ylabel("|CΨ_com|", color="#444")
    ax_phase_twin.set_ylabel("arg(CΨ_com) [°]", color="#888")
    ax_phase.set_title("|CΨ_com|(t) and arg(CΨ_com)(t)\n"
                       "(phase ≈ constant ⇒ Kingston has little Z-drift)")
    ax_phase.grid(True, alpha=0.2)
    ax_phase.legend(loc="upper right", fontsize=7)
    ax_phase_twin.legend(loc="center right", fontsize=7)

    fig.suptitle(
        "Hardware CΨ in the complex plane, ibm_kingston Bell⁺ pairs\n"
        "Real-axis trajectory = the 1D case of BOUNDARY_NAVIGATION; "
        "deviation from the real axis is Kingston Z-detuning",
        y=1.00,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_png}")


# ── Main ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("  Hardware CΨ in the complex plane (c-plane extension)")
    print("=" * 70)

    jsons = sorted(DATA_DIR.glob("cusp_slowing_*.json"))
    if not jsons:
        print(f"ERROR: no cusp_slowing_*.json found in {DATA_DIR}")
        sys.exit(1)
    json_path = jsons[-1]
    print(f"\n  input: {json_path.relative_to(REPO_ROOT)}")

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    pair_runs = data.get("pair_runs", {})
    pair_a = extract_trajectory(pair_runs.get("A_mid", {}))
    pair_b = extract_trajectory(pair_runs.get("B_high", {}))

    print(f"\n  Pair A (mid-T2): {len(pair_a['t_us'])} delay points")
    print(f"    {'t (μs)':>8}  {'|CΨ_com|':>10}  {'arg (deg)':>10}  "
          f"{'Re':>10}  {'Im':>10}")
    for t, c in zip(pair_a["t_us"], pair_a["cpsi_complex"]):
        print(f"    {t:>8.2f}  {abs(c):>10.4f}  "
              f"{np.degrees(np.angle(c)):>+10.2f}  "
              f"{c.real:>+10.4f}  {c.imag:>+10.4f}")

    print(f"\n  Pair B (high-T2): {len(pair_b['t_us'])} delay points")
    print(f"    {'t (μs)':>8}  {'|CΨ_com|':>10}  {'arg (deg)':>10}  "
          f"{'Re':>10}  {'Im':>10}")
    for t, c in zip(pair_b["t_us"], pair_b["cpsi_complex"]):
        print(f"    {t:>8.2f}  {abs(c):>10.4f}  "
              f"{np.degrees(np.angle(c)):>+10.2f}  "
              f"{c.real:>+10.4f}  {c.imag:>+10.4f}")

    plot_all(pair_a, pair_b, None, OUT)

    print()
    print("Interpretation:")
    print("  If all arg(CΨ_com) are ≈ 0, Kingston's rotating frame was well-")
    print("  calibrated: Bell+ stays on the real axis, 1D trajectory.")
    print("  If arg drifts linearly with t, there's residual Z-detuning and")
    print("  the trajectory is a 2D spiral in the c-plane.")
    print("  Either way, the cusp at c = 1/4 sits at the red marker.")
