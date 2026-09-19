#!/usr/bin/env python3
"""The exact two-qubit decay family: how the scalar-quarter crossing depends on the start.

Sweep the partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ (Bell+ is α=π/4).
Each member is an exact two-exponential decay curve (verified inline against the Lindblad
evolution to machine precision):

    CΨ(α,t) = w₀·e^(−4γt) + w₁·e^(−12γt),   w₀ = s(1−s²/2)/3,  w₁ = s³/6,  s = sin2α.

The scaling laws this draws:
  - s is the pure-state concurrence of the initial state; CΨ(0)=s/3 is a distinct linear
    readout equal to exactly one third of it;
  - a genuine temporal downward crossing occurs iff γ>0 and s>3/4; s=3/4 only touches
    ¼ at t=0, and at γ=0 the curve is constant and does not cross;
  - for every s>0, the cubic 12γ term is present with relative weight s²/2 of the start;
    that share grows quadratically and becomes appreciable toward Bell+ (s=1), the 50/50 member. At s=0 the total
    start vanishes, and s²/2 is the continuous value of the shape parameter rather than a ratio;
  - every nonzero member (0 < s ≤ 1) has the late-time 4γ term; at the endpoint s=0,
    w₀=w₁=0 and the curve is identically zero, so there is no late-time term.

Produces: simulations/results/approach_family/approach_family.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cpsi_complex_plane import site_op, liouvillian, cpsi_real, Z  # noqa: E402

SCALAR_QUARTER = 0.25


def weights(s):
    """The closed-form two-exponential weights for s = sin2α."""
    return s * (1 - 0.5 * s**2) / 3.0, s**3 / 6.0


def cpsi_closed(s, gamma, t):
    w0, w1 = weights(s)
    return w0 * np.exp(-4 * gamma * t) + w1 * np.exp(-12 * gamma * t)


def has_downward_crossing(s, gamma):
    """A sign-changing temporal passage through ¼, excluding a t=0 touch."""
    return gamma > 0.0 and s > 0.75


def cpsi_evolved(alpha, gamma, times):
    """The real CΨ(t) from the actual Lindblad evolution, for the inline check."""
    psi = np.zeros(4, dtype=complex)
    psi[0], psi[3] = np.cos(alpha), np.sin(alpha)
    rho0 = np.outer(psi, psi.conj())
    L = liouvillian(np.zeros((4, 4), dtype=complex),
                    [np.sqrt(gamma) * site_op(Z, 0, 2), np.sqrt(gamma) * site_op(Z, 1, 2)])
    rv0 = rho0.flatten(order="F")
    return np.array([cpsi_real((expm(L * t) @ rv0).reshape(4, 4, order="F")) for t in times])


def main(output_path) -> None:
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    quarter_level = SCALAR_QUARTER
    gamma = 0.05
    s_crit = 0.75                                   # at γ>0, crosses ¼ iff s > 3/4
    t = np.linspace(0, 26, 1400)

    crossing_cases = (
        (1.00, gamma, True),
        (0.80, gamma, True),
        (0.75, gamma, False),
        (0.50, gamma, False),
        (1.00, 0.0, False),
    )
    for s, gamma_case, expected in crossing_cases:
        assert has_downward_crossing(s, gamma_case) is expected
    # Mutation control: the old s-only condition misclassifies the constant γ=0 curve.
    assert (1.00 > s_crit) and not has_downward_crossing(1.00, 0.0)
    assert np.all(cpsi_closed(1.00, 0.0, t) == 1.0 / 3.0)

    # Inline verification: closed form vs Lindblad (machine precision).
    tcheck = np.linspace(0, 20, 40)
    max_err = max(np.max(np.abs(cpsi_closed(np.sin(2 * a), gamma, tcheck)
                                - cpsi_evolved(a, gamma, tcheck)))
                  for a in (np.pi / 6, np.pi / 4, np.pi / 3))
    print("=" * 78)
    print(f"  approach family   γ={gamma}   closed-form vs Lindblad max error = {max_err:.1e}")
    print("=" * 78)
    assert max_err < 1e-9, f"closed form does not match Lindblad: {max_err}"

    s_family = np.array([0.30, 0.45, 0.60, 0.75, 0.85, 0.92, 1.00])
    cmap = matplotlib.colormaps["viridis"]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19, 6.2))

    # ── Panel 1: the family of shapes (linear) ──
    for s in s_family:
        col = cmap((s - 0.25) / 0.78)
        y = cpsi_closed(s, gamma, t)
        special = abs(s - 1.0) < 1e-9 or abs(s - s_crit) < 1e-9
        ax1.plot(t, y, "-", color=col, lw=2.4 if special else 1.6,
                 alpha=0.95 if special else 0.8,
                 label=(f"s={s:.2f}" + (" (Bell+)" if abs(s-1) < 1e-9 else
                        " (threshold, starts at ¼)" if abs(s-s_crit) < 1e-9 else "")))
    ax1.axhline(
        quarter_level, color="red", ls="--", lw=1.3, alpha=0.8,
        label="scalar quarter level CΨ=1/4; not a cardioid cusp",
    )
    ax1.set_xlim(0, t[-1])
    ax1.set_ylim(0, 0.35)
    ax1.grid(True, alpha=0.2)
    ax1.set_title("The family at γ>0: each initial concurrence s gives a shape\n"
                  "for γ>0, downward crossing iff s>3/4; s=3/4 only touches ¼ at t=0")
    ax1.set_xlabel("t")
    ax1.set_ylabel("CΨ(t)")
    ax1.legend(loc="upper right", fontsize=7.5)

    # ── Panel 2: the nonzero members' shared late-time exponent (log) ──
    for s in s_family:
        col = cmap((s - 0.25) / 0.78)
        ax2.plot(t, cpsi_closed(s, gamma, t), "-", color=col, lw=1.8, alpha=0.85)
    # the e^(-4γt) slope, anchored at the Bell+ coefficient 1/6
    ax2.plot(t, (1 / 6) * np.exp(-4 * gamma * t), "--", color="black", lw=1.4,
             alpha=0.8, label="late-time exponent e^(−4γt)")
    ax2.axhline(quarter_level, color="red", ls=":", lw=1.0, alpha=0.6)
    ax2.set_yscale("log")
    ax2.set_xlim(0, t[-1])
    ax2.set_ylim(1e-3, 0.4)
    ax2.grid(True, alpha=0.2, which="both")
    ax2.set_title("Late-time every nonzero member (0<s≤1) runs parallel to the 4γ term;\n"
                  "the weights set the shape and the cubic 12γ term fades first")
    ax2.set_xlabel("t")
    ax2.set_ylabel("CΨ(t)  (log)")
    ax2.legend(loc="lower left", fontsize=8)

    # ── Panel 3: the shape parameter (the weights / cubic-term fraction) ──
    ss = np.linspace(0, 1, 400)
    w0s, w1s = weights(ss)
    cubic_frac = np.divide(w1s, w0s + w1s, out=np.zeros_like(ss), where=(w0s + w1s) > 0)
    ax3.plot(ss, w0s, "-", color="#2E8B57", lw=2.0, label="carrier weight w₀ = s(1−s²/2)/3")
    ax3.plot(ss, w1s, "-", color="#AA33CC", lw=2.0, label="cubic f³ weight w₁ = s³/6")
    ax3.plot(ss, cubic_frac, "-", color="#CC7711", lw=2.0,
             label="cubic-term fraction (s>0) = s²/2; continuous value 0 at s=0")
    ax3.axvline(s_crit, color="red", ls="--", lw=1.0, alpha=0.7)
    ax3.annotate("s=3/4\n(starts at ¼)", (s_crit, 0.45), fontsize=8, color="red", ha="center")
    ax3.axvline(1.0, color="gray", ls=":", lw=1.0, alpha=0.7)
    ax3.annotate("Bell+\n(50/50)", (1.0, 0.55), fontsize=8, color="gray", ha="center")
    ax3.set_xlim(0, 1.0)
    ax3.set_ylim(0, 0.62)
    ax3.grid(True, alpha=0.2)
    ax3.set_title("The shape parameter: the cubic term has relative share s²/2 for s>0;\n"
                  "the cubic f³ contribution belongs to CΨ; it is not a Liouvillian mode")
    ax3.set_xlabel("s = sin 2α  (initial pure-state concurrence)")
    ax3.set_ylabel("weight")
    ax3.legend(loc="upper left", fontsize=8)

    fig.suptitle(
        "Exact two-qubit decay family: |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩, CΨ(α,t) = w₀e^(−4γt) + w₁e^(−12γt).\n"
        "s is the initial pure-state concurrence; CΨ(0)=s/3 is one third of it. A temporal downward crossing occurs "
        "iff γ>0 and s>3/4; every nonzero member has the 4γ term.\n"
        "s=0 is identically zero; no late-time term",
        y=1.02, fontsize=10.5)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    fig.savefig(destination, dpi=170, bbox_inches="tight")
    plt.close()
    print("  s is the initial pure-state concurrence; CΨ(0)=s/3 is one third of it; "
          "temporal downward crossing iff γ>0 and s>3/4; cubic-term fraction for s>0 is s²/2 "
          "(Bell+ s=1: 1/2; continuous value 0 at s=0)")
    print(f"  saved: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", dest="output_path", required=True)
    args = parser.parse_args()
    main(args.output_path)
