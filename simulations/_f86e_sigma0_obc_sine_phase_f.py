"""F86e σ_0 Phase F: characterise the CHANGE (switch), not the state.

Per the contract/refrigerator reframing:
the shared contract (N-universal structure) carries zero information; the
only information is the change. So instead of fitting σ_0(∞), characterise
the switch  Δσ(N) = σ_0(N+2) − σ_0(N).

Phase F:
  1. σ_0(N) for N = 5..19 (extends Phase A's 5..15).
  2. First switches Δσ(N), step 2, odd and even N separate.
  3. log(Δσ) vs N  → exponential decay rate (smooth-periodic-integrand
     trapezoidal rule converges exponentially, NOT power-law).
  4. Δσ(N)·2^(N/2) → prefactor convergence.
  5. Switch ratio Δσ(N+2)/Δσ(N) → geometric ratio.
  6. Second difference Δ²σ (the change of the change).
  7. Geometric-tail extrapolation σ_0(∞) = σ_0(N) + Δσ(N)·r/(1−r).

Output: simulations/results/f86e_sigma0_obc_sine_phase_f/
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
import framework as fw  # noqa: E402


def hd_subspace_projector(N: int, n: int, hd_value: int) -> np.ndarray:
    P_n = fw.popcount_states(N, n)
    P_np1 = fw.popcount_states(N, n + 1)
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    cols = []
    for p in P_n:
        for q in P_np1:
            if bin(p ^ q).count("1") == hd_value:
                idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
                v = np.zeros(len(P_n) * Mnp1, dtype=complex)
                v[idx] = 1.0
                cols.append(v)
    return (np.column_stack(cols) if cols else
            np.zeros((len(P_n) * Mnp1, 0), dtype=complex))


def sigma_0(N: int, n: int = 1, gamma_0: float = 0.05) -> float:
    _, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    M_H_total = sum(M_H_per_bond)
    P_HD1 = hd_subspace_projector(N, n, 1)
    P_HD3 = hd_subspace_projector(N, n, 3)
    V_inter = P_HD1.conj().T @ M_H_total @ P_HD3
    return float(np.linalg.svd(V_inter, compute_uv=False)[0])


def main() -> None:
    gamma_0 = 0.05
    n = 1
    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_f")
    out_dir.mkdir(parents=True, exist_ok=True)

    N_values = list(range(5, 20))
    print("=" * 100)
    print("F86e Phase F: characterising the CHANGE, the switch Δσ(N) = σ_0(N+2) − σ_0(N)")
    print("=" * 100)

    sig = {}
    for N in N_values:
        t0 = time.time()
        sig[N] = sigma_0(N, n, gamma_0)
        print(f"  N={N:2d}: σ_0 = {sig[N]:.10f}  ({time.time()-t0:.1f}s)")

    # First switches, step 2
    print()
    print("First switches Δσ(N) = σ_0(N+2) − σ_0(N):")
    print(f"{'N':>4}  {'parity':>6}  {'Δσ(N)':>14}  {'Δσ·2^(N/2)':>13}  "
          f"{'ratio Δσ(N)/Δσ(N-2)':>20}  {'log Δσ':>10}")
    switches = {}
    for N in N_values:
        if N + 2 in sig:
            d = sig[N + 2] - sig[N]
            switches[N] = d
    prev_by_parity = {}
    rows = []
    for N in sorted(switches):
        d = switches[N]
        parity = "odd" if N % 2 == 1 else "even"
        scaled = d * 2 ** (N / 2)
        ratio = (d / prev_by_parity[parity]
                 if parity in prev_by_parity else float("nan"))
        prev_by_parity[parity] = d
        logd = math.log(abs(d)) if d != 0 else float("nan")
        print(f"{N:>4}  {parity:>6}  {d:>14.8e}  {scaled:>13.6f}  "
              f"{ratio:>20.6f}  {logd:>10.4f}")
        rows.append((N, parity, d, scaled, ratio, logd))

    # log-linear fit of switch decay rate, per parity
    print()
    print("Exponential decay rate (log Δσ ~ −rate·N):")
    for parity in ("odd", "even"):
        pts = [(N, math.log(abs(switches[N])))
               for N in sorted(switches) if (N % 2 == 1) == (parity == "odd")]
        if len(pts) >= 2:
            Ns = np.array([p[0] for p in pts], dtype=float)
            logs = np.array([p[1] for p in pts])
            slope, intercept = np.polyfit(Ns, logs, 1)
            # use last two points for the asymptotic rate
            rate_last = (logs[-1] - logs[-2]) / (Ns[-1] - Ns[-2])
            print(f"  {parity:>5}: global slope {slope:.6f}, "
                  f"last-pair slope {rate_last:.6f}  "
                  f"(ln2/2 = {math.log(2)/2:.6f}, ln(1/√2) = {math.log(1/math.sqrt(2)):.6f})")

    # Second difference: change of the change
    print()
    print("Second difference Δ²σ(N) = Δσ(N+2) − Δσ(N)  (change of the change):")
    for N in sorted(switches):
        if N + 2 in switches:
            dd = switches[N + 2] - switches[N]
            print(f"  N={N:2d}: Δ²σ = {dd:+.8e}  "
                  f"ratio Δσ(N+2)/Δσ(N) = {switches[N+2]/switches[N]:.6f}")

    # Geometric-tail extrapolation: if Δσ(N) ~ C·r^N, then
    #   σ_0(∞) = σ_0(N) + Δσ(N)·r²/(1−r²)   (step-2 series)
    print()
    print("Geometric-tail extrapolation of σ_0(∞):")
    print(f"{'from N':>7}  {'r (step-2)':>12}  {'σ_0(∞) est':>14}")
    odd_sorted = sorted([N for N in switches if N % 2 == 1])
    for N in odd_sorted:
        if N - 2 in switches and N % 2 == 1:
            r2 = switches[N] / switches[N - 2]   # step-2 geometric ratio
            # remaining tail from σ_0(N+2): Δσ(N+2)+Δσ(N+4)+... ≈ Δσ(N)·r2²/(1−r2)
            tail = switches[N] * r2 / (1 - r2)
            est = sig[N + 2] + tail
            print(f"{N:>7}  {r2:>12.6f}  {est:>14.8f}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    ax = axes[0]
    odd_N = [N for N in sorted(switches) if N % 2 == 1]
    even_N = [N for N in sorted(switches) if N % 2 == 0]
    ax.semilogy(odd_N, [switches[N] for N in odd_N], 'o-', label="odd N")
    ax.semilogy(even_N, [switches[N] for N in even_N], 's-', label="even N")
    ax.set_xlabel("N")
    ax.set_ylabel("Δσ(N)  (log scale)")
    ax.set_title("Switch Δσ(N) = σ_0(N+2)−σ_0(N)\n"
                 "straight line on log = exponential decay")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(odd_N, [switches[N] * 2 ** (N / 2) for N in odd_N], 'o-',
            label="odd N")
    ax.plot(even_N, [switches[N] * 2 ** (N / 2) for N in even_N], 's-',
            label="even N")
    ax.set_xlabel("N")
    ax.set_ylabel("Δσ(N) · 2^(N/2)")
    ax.set_title("Prefactor: Δσ·2^(N/2)\n"
                 "flat = pure geometric 2^(−N/2); drift = corrections")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ratios = [(N, switches[N + 2] / switches[N])
              for N in sorted(switches) if N + 2 in switches]
    ax.plot([r[0] for r in ratios], [r[1] for r in ratios], 'o-')
    ax.axhline(0.5, color="red", ls="--", label="1/2")
    ax.set_xlabel("N")
    ax.set_ylabel("Δσ(N+2) / Δσ(N)")
    ax.set_title("Switch ratio → geometric limit")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = out_dir / "phase_f_switches.png"
    plt.savefig(plot_file, dpi=120, bbox_inches="tight")
    print(f"\nPlot saved: {plot_file}")

    np.savez(out_dir / "switches.npz",
             N_values=np.array(N_values),
             sigma_0=np.array([sig[N] for N in N_values]),
             switch_N=np.array(sorted(switches)),
             switch_val=np.array([switches[N] for N in sorted(switches)]))
    print(f"Data saved: {out_dir / 'switches.npz'}")


if __name__ == "__main__":
    main()
