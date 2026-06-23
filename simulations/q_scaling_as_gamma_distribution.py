"""Q-scaling visualised as a γ-distribution sweep at fixed γ₀ = 0.05.

What does scaling Q = J/γ₀ actually change? Not γ₀ (substrate clock, fixed).
Not the Pauli basis (orthogonal labels, fixed). What changes is *which*
operators are Liouvillian eigenmodes, and therefore *how much γ each absorbs*.

Absorption Theorem (F2):
    |Re(λ_mode)| = 2γ₀ · ⟨n_XY⟩_mode.

γ₀ is the slope of the line (2γ₀ = 0.1). ⟨n_XY⟩ is the mode's light content.
J does not tune the line; J redistributes the dots along the line.

Three panels:
  1. Scatter (⟨n_XY⟩, |Re(λ)|) for all 1024 modes across 6 Q values, AT line.
  2. Histogram of ⟨n_XY⟩ at each Q (operator-population shift).
  3. Slow-mode ⟨n_XY⟩(Q) with chain prediction 0.55·Q²/N² overlaid.

Output: simulations/results/q_scaling_visual/q_scaling_at_fixed_gamma_N5.png
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

from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

from f1_chain_heisenberg_small_n_anchor import (  # noqa: E402
    chain_bonds,
    heisenberg_graph_h,
)


_PAULI_SINGLE = np.array([
    [[1, 0], [0, 1]],
    [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]],
    [[1, 0], [0, -1]],
], dtype=complex)


def pauli_string_matrix(p: int, N: int) -> np.ndarray:
    digits = []
    q = p
    for _ in range(N):
        digits.append(q & 3)
        q >>= 2
    M = _PAULI_SINGLE[digits[-1]]
    for d in reversed(digits[:-1]):
        M = np.kron(M, _PAULI_SINGLE[d])
    return M


def all_pauli_strings_n_xy(N: int) -> np.ndarray:
    n_xy = np.zeros(4 ** N, dtype=np.int32)
    for p in range(4 ** N):
        q = p
        c = 0
        for _ in range(N):
            d = q & 3
            if d == 1 or d == 2:
                c += 1
            q >>= 2
        n_xy[p] = c
    return n_xy


def all_pauli_matrices_flat(N: int) -> np.ndarray:
    d = 1 << N
    out = np.zeros((4 ** N, d * d), dtype=complex)
    for p in range(4 ** N):
        out[p] = pauli_string_matrix(p, N).reshape(-1)
    return out


def compute_eigenmode_n_xy(N: int, J: float, gamma: float,
                            pauli_mats: np.ndarray,
                            n_xy_pauli: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = 1 << N
    bonds = chain_bonds(N)
    H = heisenberg_graph_h(N, bonds, J=J)
    L = lindbladian_z_dephasing(H, [gamma] * N)

    evals, evecs = np.linalg.eig(L)
    re_abs = np.abs(evals.real)

    norm_factor = 1.0 / d
    amps = (pauli_mats.conj() @ evecs) * math.sqrt(norm_factor)
    overlaps = np.abs(amps) ** 2  # (4^N, d²)
    totals = overlaps.sum(axis=0)
    safe = np.where(totals > 1e-12, totals, 1.0)
    n_xy_per_mode = (overlaps * n_xy_pauli[:, None]).sum(axis=0) / safe
    n_xy_per_mode = np.where(totals > 1e-12, n_xy_per_mode, np.nan)

    return re_abs, n_xy_per_mode


def main() -> None:
    print("=" * 80)
    print("Q-scaling as γ-distribution: N=5 chain, γ₀=0.05 fixed")
    print("=" * 80)
    print("Reading:")
    print("  γ₀ fixed = substrate clock. Same on every site, every Q.")
    print("  J = Q·γ₀ varies = Hamiltonian rotates the operator basis.")
    print("  |Re(λ_mode)| = 2γ₀ · ⟨n_XY⟩_mode  (Absorption Theorem).")
    print("  γ₀ is the line's slope. Q changes WHERE on the line each mode sits.")
    print()

    N = 5
    gamma_substrate = 0.05
    sqrt3 = math.sqrt(3.0)
    q_values = [0.5, 1.0, 1.5, sqrt3, 2.0, 2.5]
    q_labels = ["0.5", "1.0", "1.5", "√3", "2.0", "2.5"]

    print(f"  Precomputing Pauli basis matrices (4^{N}={4**N} strings)...",
          end=" ", flush=True)
    t0 = time.time()
    pauli_mats = all_pauli_matrices_flat(N)
    n_xy_pauli = all_pauli_strings_n_xy(N).astype(float)
    print(f"{time.time()-t0:.2f}s")

    all_data = {}
    for Q in q_values:
        J = Q * gamma_substrate
        print(f"  Q={Q:.4f} (J={J:.4f}, γ₀={gamma_substrate})...",
              end=" ", flush=True)
        t0 = time.time()
        re_abs, n_xy = compute_eigenmode_n_xy(
            N, J=J, gamma=gamma_substrate,
            pauli_mats=pauli_mats, n_xy_pauli=n_xy_pauli)
        all_data[Q] = (re_abs, n_xy)
        print(f"{time.time()-t0:.2f}s ({len(re_abs)} eigenmodes)")

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.95, len(q_values)))

    ax = axes[0]
    for (Q, (re_abs, n_xy)), color, label in zip(all_data.items(), colors, q_labels):
        ax.scatter(n_xy, re_abs, s=14, alpha=0.45, color=color, label=f"Q={label}")
    x_line = np.linspace(0, N, 100)
    ax.plot(x_line, 2 * gamma_substrate * x_line, 'k--', lw=2,
            label="AT line: 2γ₀·⟨n_XY⟩")
    ax.set_xlabel("⟨n_XY⟩  (light content of eigenmode)")
    ax.set_ylabel("|Re(λ)|  (decay rate)")
    ax.set_title("All 1024 eigenmodes ride on the Absorption-Theorem line\n"
                 f"(N=5 chain, γ₀={gamma_substrate} fixed)\n"
                 "Q changes WHERE on the line each mode sits.")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    bins = np.linspace(0, N, 41)
    for (Q, (re_abs, n_xy)), color, label in zip(all_data.items(), colors, q_labels):
        valid = ~np.isnan(n_xy)
        ax.hist(n_xy[valid], bins=bins, alpha=0.4, color=color,
                label=f"Q={label}", density=False)
    ax.set_xlabel("⟨n_XY⟩")
    ax.set_ylabel("# eigenmodes per bin")
    ax.set_title("Operator-population shift with Q\n"
                 "low Q: clusters at integer ⟨n_XY⟩ (Pauli basis)\n"
                 "high Q: spreads (H-eigenbasis mixing)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    q_arr = np.array(q_values)
    slow_n_xy = []
    slow_re = []
    for Q in q_values:
        re_abs, n_xy = all_data[Q]
        non_kernel = re_abs > 1e-9
        valid_re = np.where(non_kernel, re_abs, np.inf)
        idx_min = int(np.argmin(valid_re))
        slow_n_xy.append(n_xy[idx_min])
        slow_re.append(re_abs[idx_min])
    slow_n_xy_arr = np.array(slow_n_xy)
    slow_re_arr = np.array(slow_re)

    ax.scatter(q_arr, slow_n_xy_arr, s=80, color="black",
               zorder=5, label="slow-mode ⟨n_XY⟩ (measured)")
    q_dense = np.linspace(0, q_arr.max() * 1.05, 200)
    pred = 0.55 * q_dense * q_dense / (N * N)
    ax.plot(q_dense, pred, 'r--', lw=2, label="chain prediction: 0.55·Q²/N²")
    ax.set_xlabel("Q = J / γ₀")
    ax.set_ylabel("⟨n_XY⟩ of slowest mode")
    ax.set_title("Slow-mode admixture grows as Q²\n"
                 "γ₀ fixed; J shifts how much γ the slow mode absorbs.\n"
                 "decay rate = 2γ₀·⟨n_XY⟩, scales as Q².")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_dir = Path("simulations/results/q_scaling_visual")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "q_scaling_at_fixed_gamma_N5.png"
    plt.savefig(out_file, dpi=120, bbox_inches="tight")
    print(f"\nPlot saved: {out_file}")

    print()
    print("Summary: slow-mode trace across Q (γ₀=0.05 throughout)")
    print(f"{'Q':>6}  {'J':>8}  {'|Re(λ)|_min':>14}  {'⟨n_XY⟩_slow':>14}  "
          f"{'0.55·Q²/N²':>14}  {'2γ₀·⟨n_XY⟩':>14}")
    for Q, label in zip(q_values, q_labels):
        J = Q * gamma_substrate
        re_abs, n_xy = all_data[Q]
        non_kernel = re_abs > 1e-9
        idx_min = int(np.argmin(np.where(non_kernel, re_abs, np.inf)))
        slow_nxy = n_xy[idx_min]
        slow_re_val = re_abs[idx_min]
        pred_n_xy = 0.55 * Q * Q / (N * N)
        at_pred = 2 * gamma_substrate * slow_nxy
        print(f"{label:>6}  {J:8.4f}  {slow_re_val:14.6e}  {slow_nxy:14.6f}  "
              f"{pred_n_xy:14.6f}  {at_pred:14.6e}")

    print()
    print("Reading: γ₀ stayed at 0.05 throughout. The Absorption-Theorem line")
    print("(slope = 2γ₀ = 0.1) didn't move. What moved is which operators")
    print("are eigenmodes, i.e. how much γ each absorbs. J redistributes γ.")


if __name__ == "__main__":
    main()
