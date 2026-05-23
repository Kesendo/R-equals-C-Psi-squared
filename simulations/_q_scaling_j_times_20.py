"""Re-run the Q-scaling table with J literally multiplied by 20.

Original Q anchors {0.5, 1.0, 1.5, √3, 2.0, 2.5} at γ₀=0.05 give
J in {0.025 .. 0.125}. Multiply by 20: J in {0.5 .. 2.5}, which at
γ₀=0.05 means Q in {10 .. 50} — far outside the small-Q (Lebensader)
regime.

Question: does the Absorption Theorem still hold (yes, by structure)?
Does the slow-mode ⟨n_XY⟩ still scale as Q²·0.55/N² (no, saturation)?
What does the angle look like when sin²(angle) ≈ 1?
"""
from __future__ import annotations

import math
import sys

import numpy as np

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
    print("=" * 100)
    print("Q-scaling table with J × 20 (γ₀ = 0.05 fixed, N = 5 chain)")
    print("=" * 100)
    print("Original J × 20 → Q × 20. Same machinery, just steeper angles.")
    print()

    pauli_mats = all_pauli_matrices_flat(N)
    n_xy_pauli = all_pauli_strings_n_xy(N).astype(float)

    sqrt3 = math.sqrt(3.0)
    J_orig = [0.025, 0.050, 0.075, 0.05 * sqrt3, 0.100, 0.125]
    q_orig_labels = ["0.5", "1.0", "1.5", "√3", "2.0", "2.5"]

    print(f"{'Q_orig':>8} {'J_orig':>8} {'J×20':>8} {'Q×20':>8} "
          f"{'|Re(λ)|_min':>14} {'⟨n_XY⟩_slow':>14} {'⟨n_XY⟩/Q²':>14} "
          f"{'⟨n_XY⟩/N':>10}")
    for j_orig, lbl in zip(J_orig, q_orig_labels):
        j_scaled = j_orig * 20.0
        re_abs, n_xy = compute_eigenmode_n_xy(
            N, J=j_scaled, gamma=gamma,
            pauli_mats=pauli_mats, n_xy_pauli=n_xy_pauli)
        non_kernel = re_abs > 1e-9
        idx = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
        slow_nxy = n_xy[idx]
        slow_re = re_abs[idx]
        q_scaled = j_scaled / gamma
        ratio_q2 = slow_nxy / (q_scaled * q_scaled) if q_scaled > 0 else float("nan")
        ratio_N = slow_nxy / N
        print(f"{lbl:>8} {j_orig:8.4f} {j_scaled:8.4f} {q_scaled:8.2f} "
              f"{slow_re:14.6e} {slow_nxy:14.6f} {ratio_q2:14.8f} "
              f"{ratio_N:10.4f}")

    print()
    print("Reading:")
    print("  γ₀ unchanged.  AT line slope 2γ₀ = 0.1 unchanged.")
    print("  ⟨n_XY⟩/Q² no longer constant — small-angle Q² law has saturated.")
    print("  ⟨n_XY⟩/N shows how close the slow mode is to bulk light content.")
    print("  N/2 = 2.5 = bulk-average light content (uniform Pauli mixture).")


if __name__ == "__main__":
    main()
