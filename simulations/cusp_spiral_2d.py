#!/usr/bin/env python3
"""The 2D spiral into the cusp circle: the graphical representation.

The interior axis read in the complex plane. With a common Z-drift Ω under
the dephasing, the Bell+ coherence CΨ_com becomes complex and winds inward.
The cusp ¼, a point on the real line, is the circle |CΨ| = ¼ here: every
spiral crosses it (the radial magnitude law is Ω-independent), and only the
crossing angle φ₀ − Ω·t_cross is free, the one IBM Kingston steered.

Produces:
  simulations/results/cusp_spiral_2d/cusp_spiral_2d.png   (static, annotated)
  simulations/results/cusp_spiral_2d/cusp_spiral_2d.gif   (one spiral winding in)

The static figure overlays the real Kingston spirals (best-effort: skipped
cleanly if the hardware JSON is absent). Date: 2026-06-03.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from cpsi_complex_plane import trajectory, main_cardioid, period_2_bulb  # noqa: E402

CUSP = 0.25
COLORS = ["#CC3333", "#33AACC", "#EE9922", "#AA33CC", "#33AA77"]


def quarter_circle(n_pts: int = 400) -> np.ndarray:
    """The cusp circle |CΨ| = 1/4: the 1D cusp point seen edge-on."""
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return CUSP * np.exp(1j * theta)


def first_crossing(c_vals: np.ndarray):
    """First index where |CΨ_com| <= 1/4 (the spiral entering the circle), and
    the crossing angle in degrees. Returns (idx, angle_deg) or (None, None)."""
    mags = np.abs(c_vals)
    inside = np.where(mags <= CUSP)[0]
    if len(inside) == 0:
        return None, None
    idx = int(inside[0])
    return idx, float(np.degrees(np.angle(c_vals[idx])))


def load_hardware():
    """Best-effort: the two real Kingston spirals. Returns a list of
    (label, c_vals) or [] if the data/loader is unavailable."""
    try:
        import json
        from hardware_cpsi_cplane import extract_trajectory
        data_dir = Path(__file__).parent.parent / "data" / "ibm_cusp_slowing_april2026"
        jsons = sorted(data_dir.glob("cusp_slowing_*.json"))
        if not jsons:
            return []
        with open(jsons[-1], "r", encoding="utf-8") as f:
            data = json.load(f)
        pair_runs = data.get("pair_runs", {})
        out = []
        for key, label in [("A_mid", "Kingston A (clockwise)"),
                           ("B_high", "Kingston B (counter-cw)")]:
            tr = extract_trajectory(pair_runs.get(key, {}))
            c = np.array(tr["cpsi_complex"])
            if len(c) > 0:
                out.append((label, c))
        return out
    except Exception as e:  # noqa: BLE001  best-effort overlay
        print(f"  (hardware overlay skipped: {e})")
        return []


def make_static(out_png: Path) -> None:
    # Idealized spirals: fixed γ, a fan of Ω (Ω=0 is the real-axis 1D baseline).
    gamma = 0.05
    configs = [
        {"gamma": gamma, "detuning_sum": 0.0, "phi_0": 0.0, "t_max": 25.0},
        {"gamma": gamma, "detuning_sum": 0.3, "phi_0": 0.0, "t_max": 25.0},
        {"gamma": gamma, "detuning_sum": 0.6, "phi_0": 0.0, "t_max": 25.0},
        {"gamma": gamma, "detuning_sum": 1.0, "phi_0": 0.0, "t_max": 25.0},
    ]
    trajs = [trajectory(**c) for c in configs]
    hardware = load_hardware()

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    ax_full, ax_zoom = axes

    c_card = main_cardioid()
    c_bulb = period_2_bulb()
    c_circ = quarter_circle()
    for ax in axes:
        ax.plot(c_card.real, c_card.imag, "-", color="#999", lw=0.9, alpha=0.6,
                label="Mandelbrot cardioid")
        ax.plot(c_bulb.real, c_bulb.imag, "-", color="#999", lw=0.6, alpha=0.4)
        # The new element: the cusp as a CIRCLE.
        ax.plot(c_circ.real, c_circ.imag, "--", color="red", lw=1.6, alpha=0.9,
                label="cusp circle |CΨ| = 1/4")
        ax.plot(0.25, 0, "o", color="red", markersize=7, zorder=6,
                label="the 1D cusp point")
        ax.axhline(0, color="gray", lw=0.3, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.3, alpha=0.5)

    for i, tr in enumerate(trajs):
        c = np.array(tr["cpsi_complex"])
        col = COLORS[i % len(COLORS)]
        omega = tr["detuning_sum"]
        wind = omega / (4 * tr["gamma"])
        lbl = f"Ω={omega:.1f} (Ω/4γ={wind:.1f})" if omega > 0 else "Ω=0 (the 1D real-axis line)"
        idx, ang = first_crossing(c)
        for ax in axes:
            ax.plot(c.real, c.imag, "-", color=col, lw=1.8, alpha=0.85, label=lbl)
            ax.plot(c.real[0], c.imag[0], "o", color=col, markersize=7)
            if idx is not None:
                ax.plot(c.real[idx], c.imag[idx], "*", color=col, markersize=15,
                        markeredgecolor="black", markeredgewidth=0.6, zorder=7)
        if idx is not None:
            print(f"  Ω={omega:.1f}: crosses |CΨ|=1/4 at angle {ang:+.1f}°")

    for j, (label, c) in enumerate(hardware):
        ax_zoom.plot(c.real, c.imag, "s-", color="#222", lw=1.3, markersize=5,
                     alpha=0.7, label=label if j == 0 else None)

    ax_full.set_xlim(-1.0, 0.6)
    ax_full.set_ylim(-0.6, 0.6)
    ax_full.set_aspect("equal")
    ax_full.set_title("The c-plane: every spiral crosses the same ¼-circle")
    ax_full.grid(True, alpha=0.2)
    ax_full.legend(loc="lower left", fontsize=7)
    ax_full.set_xlabel("Re(CΨ_com)")
    ax_full.set_ylabel("Im(CΨ_com)")

    ax_zoom.set_xlim(-0.32, 0.40)
    ax_zoom.set_ylim(-0.32, 0.32)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title("Zoom: the cusp is a circle, the angle is the free thing\n"
                      "(stars = crossings; squares = real Kingston spirals)")
    ax_zoom.grid(True, alpha=0.2)
    ax_zoom.legend(loc="lower left", fontsize=7)
    ax_zoom.set_xlabel("Re(CΨ_com)")
    ax_zoom.set_ylabel("Im(CΨ_com)")

    fig.suptitle(
        "The interior axis in 2D: the cusp ¼ is a circle, every spiral crosses it,\n"
        "only the crossing angle is free (the one IBM Kingston steered)", y=1.00)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_png}")


def make_gif(out_gif: Path) -> None:
    gamma, omega, phi0, t_max = 0.05, 0.5, 0.0, 25.0
    tr = trajectory(gamma=gamma, detuning_sum=omega, phi_0=phi0, t_max=t_max, n_steps=160)
    c = np.array(tr["cpsi_complex"])
    times = np.array(tr["times"])
    idx_cross, _ = first_crossing(c)

    fig, ax = plt.subplots(figsize=(8, 8))
    c_card = main_cardioid()
    c_circ = quarter_circle()
    ax.plot(c_card.real, c_card.imag, "-", color="#999", lw=0.9, alpha=0.5)
    ax.plot(c_circ.real, c_circ.imag, "--", color="red", lw=1.6, alpha=0.9,
            label="cusp circle |CΨ| = 1/4")
    ax.plot(0.25, 0, "o", color="red", markersize=6, zorder=6)
    ax.axhline(0, color="gray", lw=0.3, alpha=0.5)
    ax.axvline(0, color="gray", lw=0.3, alpha=0.5)
    ax.set_xlim(-0.40, 0.42)
    ax.set_ylim(-0.40, 0.40)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("Re(CΨ_com)")
    ax.set_ylabel("Im(CΨ_com)")
    ax.set_title(f"A spiral winding into the cusp circle (γ={gamma}, Ω={omega})")

    line, = ax.plot([], [], "-", color="#AA33CC", lw=2.0, alpha=0.9)
    head = ax.scatter([], [], color="#AA33CC", s=70, edgecolor="black",
                      linewidths=0.6, zorder=7)
    star = ax.scatter([], [], color="gold", s=220, marker="*",
                      edgecolor="black", linewidths=0.8, zorder=8)
    text = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=10,
                   va="top", family="monospace",
                   bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"))
    ax.legend(loc="lower left", fontsize=9)

    def update(frame: int):
        line.set_data(c.real[:frame + 1], c.imag[:frame + 1])
        head.set_offsets([[c.real[frame], c.imag[frame]]])
        if idx_cross is not None and frame >= idx_cross:
            star.set_offsets([[c.real[idx_cross], c.imag[idx_cross]]])
        mag = abs(c[frame])
        arg = np.degrees(np.angle(c[frame]))
        text.set_text(f"t   = {times[frame]:6.2f}\n"
                      f"|CΨ| = {mag:6.4f}\n"
                      f"arg = {arg:+6.1f}°\n"
                      f"{'INSIDE ¼' if mag <= CUSP else 'outside'}")
        return line, head, star, text

    print(f"  rendering GIF ({len(c)} frames)...")
    anim = FuncAnimation(fig, update, frames=len(c), interval=80, blit=False, repeat=True)
    anim.save(out_gif, writer=PillowWriter(fps=12))
    plt.close()
    print(f"  saved: {out_gif}")


if __name__ == "__main__":
    out_dir = Path(__file__).parent / "results" / "cusp_spiral_2d"
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("  The 2D spiral into the cusp circle")
    print("=" * 70)
    make_static(out_dir / "cusp_spiral_2d.png")
    make_gif(out_dir / "cusp_spiral_2d.gif")
    print("\nThe cusp ¼ is a circle; every spiral crosses it; only the angle is free.")
