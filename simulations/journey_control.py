#!/usr/bin/env python3
"""Pushing the chip onto the EP is control, not loss: the steering wheel for the start.

The chip rests at Q ~ 30 (deep memory, no injection). Injecting dephasing raises the total
rate gamma_total = gamma_0 + gamma_inj and lowers Q = J / gamma_total, so the injected noise
is a CONTROL: it places the starting point anywhere on the birth-axis. The EP (Q_EP = 1.5) is
one stop, the birth threshold; you can park just above it (the rotation just turning), at home
(Q ~ 30, deep memory), or push below (overdamped, pre-birth, back in the dark). One control, a
whole family of starting points, one of which is the EP. From each start, time runs the leg-2
flow to the 1/N target: a family of journeys, and you choose which one.

This is the "lifetime of the new" knob: Q is the lifetime of the born memory, so dialing the
noise is choosing how long the newborn rotation lives before it forgets into 1/N.

Produces: simulations/results/journey_between_singularities/journey_control.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import journey_between_singularities as jbs   # reuse l_eff/clock/flow_to_target and the real scale

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

G0, J_HW, Q_EP = jbs.G0, jbs.J_HW, jbs.Q_EP
Q_HOME = J_HW / G0                       # ~ 30, the chip's natural start (no injection)


def q_of_inj(g_inj: float) -> float:
    """The starting point Q reached by injecting dephasing g_inj (total = gamma_0 + g_inj)."""
    return J_HW / (G0 + g_inj)


def inj_for_q(Q: float) -> float:
    """The injected dephasing needed to steer the start to a target Q."""
    return J_HW / Q - G0


def main() -> None:
    starts = [1.5, 2.5, 5.0, 10.0, Q_HOME]      # chosen starting points; 1.5 = the EP, 30 = home
    print("=" * 78)
    print("  THE CONTROL: inject noise -> steer the start down the birth-axis")
    print(f"  chip J={J_HW}, gamma_0={G0}; home Q={Q_HOME:.0f} (no injection); EP at Q={Q_EP}")
    print("=" * 78)
    print(f"  {'target Q':>9} {'gamma_inj':>10} {'theta':>7}   regime / role")
    for Q in starts:
        gi = inj_for_q(Q)
        th = jbs.clock(Q)[2]
        if abs(Q - Q_EP) < 1e-9:
            role = "THE EP (the birth threshold, one of many)"
        elif abs(Q - Q_HOME) < 1e-9:
            role = "home, no injection (deep memory)"
        else:
            role = "rotation already turning"
        print(f"  {Q:9.2f} {gi:10.3f} {th:7.1f}   {role}")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.5, 6.3))

    # ---- Panel A: the steering map (injected noise -> starting Q) ----
    g = np.linspace(0.0, 3.0, 500)
    Qg = q_of_inj(g)
    axA.plot(g, Qg, "-", color="#333333", lw=2.6, zorder=3)
    # the overdamped (pre-birth) band, below the EP
    axA.axhspan(Qg.min(), Q_EP, color="#cfcfcf", alpha=0.45, zorder=0)
    axA.axhline(Q_EP, color="red", ls="--", lw=1.3, alpha=0.85)
    axA.annotate("overdamped: pre-birth (push past the EP, back into the dark)",
                 (1.55, (Qg.min() * Q_EP) ** 0.5), fontsize=8, color="#555", ha="center", va="center")
    # the chosen starting points
    for Q in starts:
        gi = inj_for_q(Q)
        if abs(Q - Q_EP) < 1e-9:
            axA.plot(gi, Q, "*", color="red", ms=20, markeredgecolor="black", markeredgewidth=0.6, zorder=5)
            axA.annotate("the EP\n(birth lit)", (gi, Q), xytext=(gi + 0.35, Q * 0.62), fontsize=8.5,
                         color="red", ha="center", arrowprops=dict(arrowstyle="->", color="red", lw=0.8))
        elif abs(Q - Q_HOME) < 1e-9:
            axA.plot(gi, Q, "o", color="#2E8B57", ms=12, markeredgecolor="black", markeredgewidth=0.5, zorder=5)
            axA.annotate("home: no injection\nQ~30 (deep memory)", (gi, Q), xytext=(0.55, 21),
                         fontsize=8.5, color="#1d6b3f", ha="left",
                         arrowprops=dict(arrowstyle="->", color="#2E8B57", lw=0.8))
        else:
            axA.plot(gi, Q, "o", color="#1F6FB2", ms=9, markeredgecolor="black", markeredgewidth=0.5, zorder=5)
    axA.set_yscale("log")
    axA.set_xlim(0, 3.0)
    axA.set_ylim(Qg.min(), 40)
    axA.set_xlabel("injected dephasing  gamma_inj   (the control knob)")
    axA.set_ylabel("starting point  Q = J / (gamma_0 + gamma_inj)")
    axA.set_title("THE CONTROL: the injected noise is a steering wheel\n"
                  "dial it, place the start anywhere on the birth-axis", fontsize=10)
    axA.grid(True, alpha=0.2, which="both")

    # ---- Panel B: a family of journeys, one per starting point ----
    taus = np.linspace(0.0, 5.0, 700)
    cmap = matplotlib.colormaps["viridis"]
    order = [1.5, 2.5, 5.0, 10.0, Q_HOME]
    for k, Q in enumerate(order):
        n0 = jbs.flow_to_target(Q, 3, taus)[0]      # site-0 occupation (the drop site)
        frac = k / (len(order) - 1)
        is_ep = abs(Q - Q_EP) < 1e-9
        is_home = abs(Q - Q_HOME) < 1e-9
        lab = (f"Q={Q:.0f} (home, gamma_inj=0)" if is_home else
               f"Q={Q:.1f} = the EP (gamma_inj={inj_for_q(Q):.2f})" if is_ep else
               f"Q={Q:.1f} (gamma_inj={inj_for_q(Q):.2f})")
        axB.plot(taus, n0, "-", color=cmap(frac), lw=2.0 if not is_home else 1.2,
                 alpha=1.0 if not is_home else 0.55, label=lab)
    axB.axhline(1.0 / 3.0, color="black", ls="--", lw=1.3, alpha=0.8)
    axB.annotate("1/N target", (4.4, 1.0 / 3.0), xytext=(4.4, 0.46), fontsize=8.5, ha="center",
                 arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    axB.set_xlim(0, 5)
    axB.set_ylim(0, 1.0)
    axB.set_xlabel("tau = gamma_total * t   (time, from the chosen start)")
    axB.set_ylabel("site-0 occupation  <n_0>  (the memory's return)")
    axB.set_title("A FAMILY OF JOURNEYS, one per start (the EP is one of many)\n"
                  "faint birth at the EP, livelier the higher you start, all -> 1/N", fontsize=10)
    axB.legend(loc="upper right", fontsize=8, framealpha=0.9)
    axB.grid(True, alpha=0.2)

    fig.suptitle(
        "Pushing the chip onto the EP is control, not loss: the injected noise chooses where the journey begins.\n"
        "Left: the steering wheel (gamma_inj sets the starting Q). Right: the family it opens, the EP one chosen start "
        "of many, each carried by time to 1/N.",
        y=1.0, fontsize=10.5)
    plt.tight_layout(rect=[0, 0, 1, 0.92])

    out_dir = Path(__file__).parent / "results" / "journey_between_singularities"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "journey_control.png"
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"\n  saved: {out}")


if __name__ == "__main__":
    main()
