#!/usr/bin/env python3
"""Illustrative dephasing control for the single-excitation population flow.

At fixed coupling, increasing the canonical Lindblad rate lowers Q=J/gamma. This script shows
that mathematical control relation and propagates population trajectories for several chosen Q.
It is not a calibration of the Kingston runner and makes no critical-damping or EP claim.

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

Q_BASE = 40.0


def q_of_inj(relative_increase: float) -> float:
    """Canonical Q after increasing gamma by a dimensionless fraction."""
    return Q_BASE / (1.0 + relative_increase)


def inj_for_q(Q: float) -> float:
    """Relative canonical-rate increase needed to steer from Q_BASE to Q."""
    return Q_BASE / Q - 1.0


def control_domain(q_values) -> tuple[float, float]:
    """Plot domain containing every selected canonical-Q marker."""
    injections = [inj_for_q(Q) for Q in q_values]
    return min(injections), max(injections)


def main() -> None:
    starts = [3.0, 5.0, 10.0, 20.0, Q_BASE]
    min_injection, max_injection = control_domain(starts)
    print("=" * 78)
    print("  CANONICAL CONTROL: increase gamma -> lower Q")
    print(f"  illustrative base Q={Q_BASE:.0f}; no hardware or spectral-transition calibration")
    print("=" * 78)
    print(f"  {'target Q':>9} {'gamma_inj':>10} {'theta':>7}   regime / role")
    for Q in starts:
        gi = inj_for_q(Q)
        th = jbs.clock(Q)[2]
        if abs(Q - Q_BASE) < 1e-9:
            role = "baseline, no added dephasing"
        else:
            role = "rotation already turning"
        print(f"  {Q:9.2f} {gi:10.3f} {th:7.1f}   {role}")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.5, 6.3))

    # ---- Panel A: the steering map (injected noise -> starting Q) ----
    g = np.linspace(min_injection, max_injection, 500)
    Qg = q_of_inj(g)
    axA.plot(g, Qg, "-", color="#333333", lw=2.6, zorder=3)
    # the chosen starting points
    for Q in starts:
        gi = inj_for_q(Q)
        if abs(Q - Q_BASE) < 1e-9:
            axA.plot(gi, Q, "o", color="#2E8B57", ms=12, markeredgecolor="black", markeredgewidth=0.5, zorder=5)
            axA.annotate("baseline: no added dephasing", (gi, Q), xytext=(0.55, 28),
                         fontsize=8.5, color="#1d6b3f", ha="left",
                         arrowprops=dict(arrowstyle="->", color="#2E8B57", lw=0.8))
        else:
            axA.plot(gi, Q, "o", color="#1F6FB2", ms=9, markeredgecolor="black", markeredgewidth=0.5, zorder=5)
    axA.set_yscale("log")
    axA.set_xlim(min_injection, max_injection)
    axA.set_ylim(min(starts), Q_BASE)
    axA.set_xlabel("relative increase of canonical Lindblad gamma")
    axA.set_ylabel("canonical Q = Q_base / (1 + relative increase)")
    axA.set_title("ILLUSTRATIVE CONTROL MAP\nno hardware or EP calibration", fontsize=10)
    axA.grid(True, alpha=0.2, which="both")

    # ---- Panel B: a family of journeys, one per starting point ----
    taus = np.linspace(0.0, 5.0, 700)
    cmap = matplotlib.colormaps["viridis"]
    order = starts
    for k, Q in enumerate(order):
        n0 = jbs.flow_to_target(Q, 3, taus)[0]      # site-0 occupation (the drop site)
        frac = k / (len(order) - 1)
        is_base = abs(Q - Q_BASE) < 1e-9
        lab = f"Q={Q:.0f} (relative increase={inj_for_q(Q):.2f})"
        axB.plot(taus, n0, "-", color=cmap(frac), lw=1.2 if is_base else 2.0,
                 alpha=0.55 if is_base else 1.0, label=lab)
    axB.axhline(1.0 / 3.0, color="black", ls="--", lw=1.3, alpha=0.8)
    axB.annotate("1/N target", (4.4, 1.0 / 3.0), xytext=(4.4, 0.46), fontsize=8.5, ha="center",
                 arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    axB.set_xlim(0, 5)
    axB.set_ylim(0, 1.0)
    axB.set_xlabel("tau = gamma_per_site * t   (uniform per-site rate)")
    axB.set_ylabel("site-0 occupation  <n_0>  (population return)")
    axB.set_title("POPULATION TRAJECTORIES AT CHOSEN Q\nall approach the 1/N fixed point", fontsize=10)
    axB.legend(loc="upper right", fontsize=8, framealpha=0.9)
    axB.grid(True, alpha=0.2)

    fig.suptitle(
        "Canonical dephasing control and the single-excitation population flow.\n"
        "Illustrative model only; the Kingston population handover has separate rate-book provenance.",
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
