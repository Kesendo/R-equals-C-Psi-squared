#!/usr/bin/env python3
"""The family of approach shapes: how the approach to the cusp ¼ depends on the start.

Sweep the partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ (Bell+ is α=π/4).
Each member is a two-exponential approach to ¼, closed-form (verified inline against the
Lindblad evolution to machine precision):

    CΨ(α,t) = w₀·e^(−4γt) + w₁·e^(−12γt),   w₀ = s(1−s²/2)/3,  w₁ = s³/6,  s = sin2α.

The scaling laws this draws:
  - the start CΨ(0) = s/3: the starting height IS the entanglement;
  - it crosses ¼ only if s > 3/4 (a threshold in entanglement); s = 3/4 starts exactly on ¼;
  - the fast mode (12γ) carries a fraction s²/2 of the start, growing quadratically: only strong
    entanglement excites it; Bell+ (s=1) is the 50/50 member;
  - every member shares the carrier 4γ (the slowest mode, the eigenvalue −γ₀) and collapses onto
    it at late time. The shape is the early harmonic transient; the carrier is universal. This is the
    "the slowing is ours" reading (spiral_slowing.py) across the whole family.

Produces: simulations/results/approach_family/approach_family.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from cpsi_complex_plane import site_op, liouvillian, cpsi_real, Z  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CUSP = 0.25


def weights(s):
    """The closed-form two-exponential weights for s = sin2α."""
    return s * (1 - 0.5 * s**2) / 3.0, s**3 / 6.0


def cpsi_closed(s, gamma, t):
    w0, w1 = weights(s)
    return w0 * np.exp(-4 * gamma * t) + w1 * np.exp(-12 * gamma * t)


def cpsi_evolved(alpha, gamma, times):
    """The real CΨ(t) from the actual Lindblad evolution, for the inline check."""
    psi = np.zeros(4, dtype=complex)
    psi[0], psi[3] = np.cos(alpha), np.sin(alpha)
    rho0 = np.outer(psi, psi.conj())
    L = liouvillian(np.zeros((4, 4), dtype=complex),
                    [np.sqrt(gamma) * site_op(Z, 0, 2), np.sqrt(gamma) * site_op(Z, 1, 2)])
    rv0 = rho0.flatten(order="F")
    return np.array([cpsi_real((expm(L * t) @ rv0).reshape(4, 4, order="F")) for t in times])


def main() -> None:
    gamma = 0.05
    s_crit = 0.75                                   # crosses ¼ iff s > 3/4
    t = np.linspace(0, 26, 1400)

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
    ax1.axhline(CUSP, color="red", ls="--", lw=1.3, alpha=0.8, label="the cusp ¼")
    ax1.set_xlim(0, t[-1])
    ax1.set_ylim(0, 0.35)
    ax1.grid(True, alpha=0.2)
    ax1.set_title("The family: each entanglement s a shape\n(s>3/4 crosses ¼; s<3/4 starts below it)")
    ax1.set_xlabel("t")
    ax1.set_ylabel("CΨ(t)")
    ax1.legend(loc="upper right", fontsize=7.5)

    # ── Panel 2: the carrier collapse (log) ──
    for s in s_family:
        col = cmap((s - 0.25) / 0.78)
        ax2.plot(t, cpsi_closed(s, gamma, t), "-", color=col, lw=1.8, alpha=0.85)
    # the bare carrier slope e^(-4γt), anchored at the Bell+ carrier weight 1/6
    ax2.plot(t, (1 / 6) * np.exp(-4 * gamma * t), "--", color="black", lw=1.4,
             alpha=0.8, label="carrier slope e^(−4γt)")
    ax2.axhline(CUSP, color="red", ls=":", lw=1.0, alpha=0.6)
    ax2.set_yscale("log")
    ax2.set_xlim(0, t[-1])
    ax2.set_ylim(1e-3, 0.4)
    ax2.grid(True, alpha=0.2, which="both")
    ax2.set_title("Late-time every member runs parallel to the carrier 4γ;\n"
                  "the shape lives only in the early harmonic (steeper start)")
    ax2.set_xlabel("t")
    ax2.set_ylabel("CΨ(t)  (log)")
    ax2.legend(loc="lower left", fontsize=8)

    # ── Panel 3: the shape parameter (the weights / harmonic fraction) ──
    ss = np.linspace(0, 1, 400)
    w0s, w1s = weights(ss)
    harm_frac = np.divide(w1s, w0s + w1s, out=np.zeros_like(ss), where=(w0s + w1s) > 0)
    ax3.plot(ss, w0s, "-", color="#2E8B57", lw=2.0, label="carrier weight w₀ = s(1−s²/2)/3")
    ax3.plot(ss, w1s, "-", color="#AA33CC", lw=2.0, label="harmonic weight w₁ = s³/6")
    ax3.plot(ss, harm_frac, "-", color="#CC7711", lw=2.0,
             label="harmonic fraction at t=0 = s²/2")
    ax3.axvline(s_crit, color="red", ls="--", lw=1.0, alpha=0.7)
    ax3.annotate("s=3/4\n(starts at ¼)", (s_crit, 0.45), fontsize=8, color="red", ha="center")
    ax3.axvline(1.0, color="gray", ls=":", lw=1.0, alpha=0.7)
    ax3.annotate("Bell+\n(50/50)", (1.0, 0.55), fontsize=8, color="gray", ha="center")
    ax3.set_xlim(0, 1.0)
    ax3.set_ylim(0, 0.62)
    ax3.grid(True, alpha=0.2)
    ax3.set_title("The shape parameter: the fast mode (12γ) carries s²/2,\n"
                  "growing quadratically; only strong entanglement excites it")
    ax3.set_xlabel("s = sin 2α  (the entanglement)")
    ax3.set_ylabel("weight")
    ax3.legend(loc="upper left", fontsize=8)

    fig.suptitle(
        "The family of approach shapes: |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩, CΨ(α,t) = w₀e^(−4γt) + w₁e^(−12γt).\n"
        "The start height is the entanglement (CΨ(0)=s/3), the crossing is a threshold (s>3/4), the shape "
        "is the early harmonic (s³); every member shares the carrier 4γ.",
        y=1.02, fontsize=10.5)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    out_dir = Path(__file__).parent / "results" / "approach_family"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "approach_family.png"
    plt.savefig(out, dpi=170, bbox_inches="tight")
    plt.close()
    print(f"  CΨ(0)=s/3; crosses ¼ iff s>3/4; harmonic fraction at t=0 = s²/2 (Bell+ s=1: 1/2)")
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
