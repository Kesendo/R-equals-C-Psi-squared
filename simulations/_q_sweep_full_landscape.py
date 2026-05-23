"""Full Q-landscape sweep: J ∈ [1e-10, 20] at γ₀ = 0.05, N = 5 chain.

Spans 11+ decades of Q to show the three regimes of the slow mode:

  1. Noise/kernel floor (Q < ~1e-5):
     Lebensader mode below machine precision; eigensolver returns it as kernel.
     Slowest non-kernel mode = F50 weight-1 mode, |Re(λ)| = 2γ₀ = 0.1 exactly.

  2. Q² Lebensader regime (~1e-5 < Q < ~7):
     Lebensader mode resolvable, ⟨n_XY⟩ ≈ 0.55·Q²/N², |Re(λ)| ≈ 2γ₀·0.55·Q²/N².

  3. F50 saturation ceiling (Q > ~7 at N=5):
     Lebensader hits F50 ceiling, ⟨n_XY⟩ = 1 exactly, |Re(λ)| = 2γ₀ = 0.1 exactly.

The "canonical Q-band" [0.2, 2.5] sits inside region 2, well below saturation.

Output: simulations/results/q_scaling_visual/q_sweep_full_landscape_N5.png
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, "simulations")

from _q_scaling_as_gamma_distribution import (  # noqa: E402
    all_pauli_matrices_flat,
    all_pauli_strings_n_xy,
    compute_eigenmode_n_xy,
)


def main() -> None:
    N = 5
    gamma = 0.05
    n_points = 50
    J_min, J_max = 1e-10, 20.0

    print("=" * 100)
    print(f"Full Q-landscape sweep: γ₀={gamma}, N={N} chain")
    print(f"  J range: [{J_min:.0e}, {J_max}]  ({n_points} log-spaced points)")
    print(f"  Q range: [{J_min/gamma:.2e}, {J_max/gamma:.1f}]")
    print("=" * 100)

    J_values = np.logspace(math.log10(J_min), math.log10(J_max), n_points)
    Q_values = J_values / gamma

    pauli_mats = all_pauli_matrices_flat(N)
    n_xy_pauli = all_pauli_strings_n_xy(N).astype(float)

    slow_re = np.zeros(n_points)
    slow_nxy = np.zeros(n_points)

    t_total = time.time()
    for i, J in enumerate(J_values):
        t0 = time.time()
        re_abs, n_xy = compute_eigenmode_n_xy(
            N, J=J, gamma=gamma,
            pauli_mats=pauli_mats, n_xy_pauli=n_xy_pauli)
        non_kernel = re_abs > 1e-13
        if np.any(non_kernel):
            idx = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
            slow_re[i] = re_abs[idx]
            slow_nxy[i] = n_xy[idx]
        else:
            slow_re[i] = np.nan
            slow_nxy[i] = np.nan
        if i % 5 == 0 or i == n_points - 1:
            print(f"  [{i+1:2d}/{n_points}] J={J:.3e} (Q={Q_values[i]:.3e}): "
                  f"|Re(λ)|_min={slow_re[i]:.3e}, ⟨n_XY⟩_slow={slow_nxy[i]:.6f}  "
                  f"({time.time()-t0:.1f}s)")
    print(f"  Total: {time.time()-t_total:.1f}s")

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    Q_ref = np.logspace(math.log10(Q_values[0] * 0.5),
                        math.log10(Q_values[-1] * 2), 400)

    ax = axes[0]
    ax.axvspan(0.2, 2.5, color="lightgreen", alpha=0.25,
               label="canonical Q-band [0.2, 2.5]")
    ax.loglog(Q_values, slow_nxy, 'o-', color="black", markersize=6,
              label="⟨n_XY⟩_slow (measured)")
    ax.loglog(Q_ref, 0.55 * Q_ref**2 / N**2, 'r--', lw=2,
              label="0.55·Q²/N² (Lebensader prediction)")
    ax.axhline(1.0, color="blue", ls=":", lw=2,
               label="F50 ceiling (⟨n_XY⟩ = 1)")
    ax.set_xlabel("Q = J / γ₀")
    ax.set_ylabel("⟨n_XY⟩ of slowest non-kernel mode")
    ax.set_title(f"Slow-mode light content over 11 decades of Q\n"
                 f"(N={N} chain, γ₀={gamma})")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_xlim(Q_values[0] * 0.5, Q_values[-1] * 2)
    ax.set_ylim(1e-15, 5)

    ax = axes[1]
    ax.axvspan(0.2, 2.5, color="lightgreen", alpha=0.25,
               label="canonical Q-band [0.2, 2.5]")
    ax.loglog(Q_values, slow_re, 'o-', color="black", markersize=6,
              label="|Re(λ)|_min (measured)")
    ax.loglog(Q_ref, 2 * gamma * 0.55 * Q_ref**2 / N**2, 'r--', lw=2,
              label="2γ₀·0.55·Q²/N² (Lebensader)")
    ax.axhline(2 * gamma, color="blue", ls=":", lw=2,
               label=f"F50 floor (2γ₀ = {2*gamma})")
    ax.set_xlabel("Q = J / γ₀")
    ax.set_ylabel("|Re(λ)| of slowest non-kernel mode")
    ax.set_title("Three regimes:\n"
                 "noise floor (left) → Q² rising → F50 ceiling (right)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_xlim(Q_values[0] * 0.5, Q_values[-1] * 2)
    ax.set_ylim(1e-16, 1)

    plt.tight_layout()
    out_dir = Path("simulations/results/q_scaling_visual")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "q_sweep_full_landscape_N5.png"
    plt.savefig(out_file, dpi=120, bbox_inches="tight")
    print(f"\nPlot saved: {out_file}")

    print()
    print(f"Strategic samples across the landscape (γ₀={gamma}, N={N}):")
    print(f"{'Q':>14}  {'J':>14}  {'|Re(λ)|_min':>14}  {'⟨n_XY⟩_slow':>14}  {'regime':>22}")

    sample_indices = list(range(0, n_points, max(1, n_points // 12)))
    if (n_points - 1) not in sample_indices:
        sample_indices.append(n_points - 1)

    for i in sample_indices:
        Q = Q_values[i]
        J = J_values[i]
        nxy = slow_nxy[i]
        if np.isnan(nxy) or nxy < 1e-13:
            regime = "noise/kernel"
        elif nxy < 1e-3:
            regime = "Q² below band"
        elif 0.2 <= Q <= 2.5:
            regime = "canonical Q-band"
        elif nxy < 0.99:
            regime = "transition"
        else:
            regime = "F50 saturation"
        print(f"{Q:14.3e}  {J:14.3e}  {slow_re[i]:14.6e}  "
              f"{nxy:14.6e}  {regime:>22}")

    print()
    print("Reading:")
    print(f"  Lower edge of resolvable Lebensader: Q ≈ 1e-5 (set by float precision).")
    print(f"  Upper edge (F50 ceiling): Q² ≈ N²/0.55 → Q ≈ {math.sqrt(N**2/0.55):.2f} for N={N}.")
    print(f"  The Lebensader band is geometrically bounded: between numerical floor")
    print(f"  and the F50 saturation, ~7 decades wide. Canonical [0.2, 2.5] sits inside.")


if __name__ == "__main__":
    main()
