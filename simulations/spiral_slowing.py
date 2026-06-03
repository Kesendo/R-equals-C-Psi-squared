#!/usr/bin/env python3
"""The spiral slows, but the carrier does not: the slowing is ours.

A reading of one thing the eye catches in cusp_spiral_2d: the spiral visibly slows as it
winds into the cusp. It really does, but not because of the rotation, and not because of
a brake at ¼. Split the coherence into its two factors,

    CΨ = purity · coherence,   purity = Tr(ρ²) = ½(1+f²),   coherence = f/3,   f = e^(−4γt),

and the story is clean: the bare coherence (the off-diagonal, the carrier, the eigenvalue
−γ₀) descends at a perfectly constant rate −4γ. It never slows. What slows is CΨ, because
it folds in the purity, which collapses fast early (pure Bell → mixed floor ½) then
flattens. So CΨ's log-rate halves, −8γ → −4γ, and in the linear-radius picture the main
figure draws, that flattening looks even more dramatic. The rotation only bends this
constant-rate fall into a logarithmic spiral that crowds geometrically toward the center.

Same lesson as the interior axis: the slowing is the observable's and the lens's, not the
carrier's. The carrier's own clock runs steady.

Produces: simulations/results/cusp_spiral_2d/spiral_slowing.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CUSP = 0.25


def main() -> None:
    gamma, omega, phi0 = 0.05, 0.4, 0.0
    wind = omega / (4 * gamma)

    # Fine grid for curves; equal-time dots for the "crowding" panel.
    t = np.linspace(0.0, 30.0, 600001)
    f = np.exp(-4 * gamma * t)
    mag = f * (1 + f**2) / 6          # |CΨ|, the observable
    coh = f / 3.0                     # the bare coherence (the carrier), exp at rate 4γ
    c = mag * np.exp(1j * (phi0 - omega * t))

    dlog_cpsi = -np.gradient(np.log(mag), t) / gamma    # descent rate of CΨ, in units of γ
    dlog_coh = -np.gradient(np.log(coh), t) / gamma     # descent rate of the carrier (≡ 4)
    i_cross = int(np.argmin(np.abs(mag - CUSP)))
    t_cross = t[i_cross]

    print("=" * 80)
    print(f"  spiral slowing   γ={gamma}  Ω={omega}  winding Ω/4γ={wind:.1f}")
    print("=" * 80)
    print(f"  CΨ descent rate:   {dlog_cpsi[0]:.2f}γ at start  ->  {dlog_cpsi[-1]:.2f}γ late")
    print(f"  carrier rate:      {dlog_coh[0]:.2f}γ throughout (constant)")
    print(f"  per-turn shrink:   e^(-2π/{wind:.0f}) = {np.exp(-2*np.pi/wind):.4f}")
    print(f"  crosses ¼ at t={t_cross:.2f}, |CΨ| still falling at finite pace there")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19, 6.2))

    # ── Panel 1: the spiral with equal-time dots (the slowing you see) ──
    circ = CUSP * np.exp(1j * np.linspace(0, 2 * np.pi, 400))
    ax1.plot(circ.real, circ.imag, "--", color="red", lw=1.4, alpha=0.9,
             label="cusp circle |CΨ| = 1/4")
    ax1.plot(0.0, 0.0, "+", color="gray", markersize=9)
    ax1.plot(c.real, c.imag, "-", color="#AA33CC", lw=1.3, alpha=0.6, zorder=3)
    n_dots = 46
    td = np.linspace(0.0, t[-1], n_dots)
    fd = np.exp(-4 * gamma * td)
    magd = fd * (1 + fd**2) / 6
    cd = magd * np.exp(1j * (phi0 - omega * td))
    ax1.scatter(cd.real, cd.imag, c=td, cmap="viridis", s=42, edgecolor="black",
                linewidths=0.4, zorder=5)
    ax1.plot(c.real[i_cross], c.imag[i_cross], "*", color="gold", markeredgecolor="black",
             markeredgewidth=0.7, markersize=20, zorder=6, label="crosses ¼")
    ax1.set_xlim(-0.36, 0.40)
    ax1.set_ylim(-0.38, 0.38)
    ax1.set_aspect("equal")
    ax1.grid(True, alpha=0.2)
    ax1.set_title("Equal time-steps, crowding inward:\nthe spiral slows as it winds in")
    ax1.set_xlabel("Re(CΨ_com)")
    ax1.set_ylabel("Im(CΨ_com)")
    ax1.legend(loc="lower left", fontsize=8)

    # ── Panel 2: linear-looking decay is a log straight line for the carrier ──
    ax2.set_yscale("log")
    ax2.plot(t, mag, "-", color="#AA33CC", lw=2.0, label="|CΨ| (the observable): bends")
    ax2.plot(t, coh, "-", color="#2E8B57", lw=2.0,
             label="coherence f/3 (the carrier): straight = constant rate")
    ax2.axhline(CUSP, color="red", ls=":", lw=1.0, alpha=0.7)
    ax2.plot(t_cross, CUSP, "*", color="gold", markeredgecolor="black",
             markeredgewidth=0.6, markersize=16, zorder=6)
    ax2.annotate("|CΨ| = 1/4", (t_cross, CUSP), textcoords="offset points",
                 xytext=(8, 6), fontsize=8, color="red")
    ax2.annotate("carrier: constant −4γ\n(a straight line in log)", (22, coh[int(22/t[-1]*len(t))]),
                 textcoords="offset points", xytext=(-150, 20), fontsize=8, color="#2E8B57")
    ax2.set_xlim(0, t[-1])
    ax2.grid(True, alpha=0.2, which="both")
    ax2.set_title("On a log axis the carrier is dead straight;\nCΨ bends (steeper early), so it 'slows'")
    ax2.set_xlabel("t")
    ax2.set_ylabel("magnitude (log)")
    ax2.legend(loc="lower left", fontsize=8)

    # ── Panel 3: the descent rate, and the gap that is the depurification ──
    ax3.plot(t, dlog_cpsi, "-", color="#AA33CC", lw=2.0, label="CΨ descent rate: 8γ → 4γ (slows)")
    ax3.plot(t, dlog_coh, "-", color="#2E8B57", lw=2.0, label="carrier rate: 4γ (never slows)")
    ax3.fill_between(t, dlog_coh, dlog_cpsi, color="#AA33CC", alpha=0.15)
    ax3.annotate("the gap = the purity collapsing\n(this is the slowing, and it is ours)",
                 (3.0, 6.0), fontsize=8, color="#7722AA")
    ax3.axvline(t_cross, color="red", ls=":", lw=1.0, alpha=0.7)
    ax3.annotate("¼", (t_cross, 7.6), textcoords="offset points", xytext=(4, 0),
                 fontsize=9, color="red")
    ax3.set_xlim(0, t[-1])
    ax3.set_ylim(3.0, 8.6)
    ax3.grid(True, alpha=0.2)
    ax3.set_title("The slowing lives in the gap:\nthe depurification, not the carrier")
    ax3.set_xlabel("t")
    ax3.set_ylabel("descent rate −d ln(·)/dt  [units of γ]")
    ax3.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "Does the spiral slow? Yes, but the carrier does not. The carrier (the coherence, the\n"
        "eigenvalue −γ₀) falls at a constant rate; CΨ only seems to slow, through the purity and the "
        "linear lens. The rotation just winds it into a spiral. The slowing is ours.",
        y=1.02, fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.94])

    out_dir = Path(__file__).parent / "results" / "cusp_spiral_2d"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "spiral_slowing.png"
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
