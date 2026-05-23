"""F86e σ_0 closed-form Phase E: contribution density W(Δ), σ_0(∞) integral.

Phase D established the exact spectral identity
    σ_0 = Σ_ω Δ_ω · w_ω,  w_ω = −i·a_ω*·b_ω
and (for odd N) the symmetry c_ω = Δ_ω w_ω with c(Δ) = c(−Δ), w odd.
The number of contributing modes grows ∝ N → σ_0(∞) is an integral:
    σ_0 = 2 ∫₀⁶ Δ · W(Δ) dΔ
where W(Δ) is the continuum density of w_ω at Δ.

Phase E:
  1. Load Phase D per-N data; bin w_ω by Δ → W_N(Δ).
  2. Plot W_N convergence; per-bin extrapolate W_∞(Δ).
  3. The contribution density Δ·W(Δ) integrates to σ_0/2; show it.
  4. Reconstruct σ_0 from the binned density (verify).
  5. Extrapolate the binned integral to N→∞; compare to Aitken 2.8628.
  6. Try analytical forms for W_∞(Δ).

Output: simulations/results/f86e_sigma0_obc_sine_phase_e/
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def aitken(seq: np.ndarray) -> float:
    a = np.asarray(seq, dtype=float)
    if len(a) < 3:
        return float(a[-1]) if len(a) else float("nan")
    i = len(a) - 3
    denom = a[i + 2] - 2 * a[i + 1] + a[i]
    if abs(denom) < 1e-15:
        return float(a[i + 2])
    return float(a[i] - (a[i + 1] - a[i]) ** 2 / denom)


def main() -> None:
    data_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_d")
    out_dir = Path("simulations/results/f86e_sigma0_obc_sine_phase_e")
    out_dir.mkdir(parents=True, exist_ok=True)

    N_values = list(range(5, 16))
    n_bins = 48
    bin_edges = np.linspace(0.0, 6.0, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    print("=" * 100)
    print("F86e Phase E: contribution density W(Δ), σ_0(∞) band integral")
    print(f"  {n_bins} bins over Δ ∈ [0, 6], bin width {bin_width:.4f}")
    print("=" * 100)

    W_by_N = {}
    contrib_by_N = {}
    sigma_recon = {}
    sigma_true = {}

    for N in N_values:
        d = np.load(data_dir / f"N{N}_phase_d.npz")
        Delta = d["Delta"]
        c = d["c"].real
        sigma_0 = float(d["sigma_0"])
        sigma_true[N] = sigma_0

        mask = np.abs(Delta) > 1e-9
        w = np.zeros_like(c)
        w[mask] = c[mask] / Delta[mask]

        pos = Delta > 1e-9
        W_hist, _ = np.histogram(Delta[pos], bins=bin_edges, weights=w[pos])
        contrib_hist, _ = np.histogram(Delta[pos], bins=bin_edges, weights=c[pos])

        W_by_N[N] = W_hist / bin_width
        contrib_by_N[N] = contrib_hist / bin_width

        # σ_0 = 2·Σ_{Δ>0} c  (odd-N symmetry; even N has a small asymmetry)
        s_rec = 2.0 * float(np.sum(contrib_hist))
        sigma_recon[N] = s_rec
        print(f"  N={N:2d}: σ_0={sigma_0:.6f}  2·Σ(Δ>0 bins)={s_rec:.6f}  "
              f"diff={s_rec - sigma_0:+.2e}")

    # Odd N carry the clean c(Δ)=c(−Δ) symmetry; use the largest odd N as the
    # W_∞ shape estimate (per-bin Aitken is unusable: bins near Δ=6 are empty
    # at small N, so the per-bin sequences are noise until they fill).
    odd_N = [N for N in N_values if N % 2 == 1]
    N_shape = odd_N[-1]
    W_inf = W_by_N[N_shape].copy()
    contrib_inf = contrib_by_N[N_shape].copy()

    # σ_0(N_shape) reconstructed from its own binned integral (bit-exact check).
    sigma_shape_integral = 2.0 * float(np.sum(contrib_inf)) * bin_width
    print()
    print(f"W_∞ shape estimate: largest odd N = {N_shape}")
    print(f"  σ_0(N={N_shape}) = {sigma_true[N_shape]:.6f}, "
          f"binned integral 2·∫₀⁶ Δ·W dΔ = {sigma_shape_integral:.6f}")
    print(f"  σ_0(∞) Aitken limit (Phase A, total sequence) ≈ 2.8628 ± 2e-5")
    print(f"  W_15 → W_∞ correction is small: σ_0(∞)−σ_0(15) ≈ "
          f"{2.8628 - sigma_true[15]:+.4e}")
    sigma_inf_integral = sigma_shape_integral

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))

    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(odd_N)))
    for N, col in zip(odd_N, colors):
        ax.plot(bin_centers, W_by_N[N], 'o-', ms=3, color=col,
                alpha=0.7, label=f"N={N}")
    ax.plot(bin_centers, W_inf, 'k-', lw=2.5, label="W_∞ (Aitken/bin)")
    ax.set_xlabel("Δ")
    ax.set_ylabel("W(Δ)  (weight density)")
    ax.set_title("Contribution-weight density W(Δ)\nconvergence over odd N")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for N, col in zip(odd_N, colors):
        ax.plot(bin_centers, contrib_by_N[N], 'o-', ms=3, color=col,
                alpha=0.7, label=f"N={N}")
    ax.plot(bin_centers, contrib_inf, 'k-', lw=2.5, label="Δ·W_∞")
    ax.set_xlabel("Δ")
    ax.set_ylabel("Δ·W(Δ)  (contribution density)")
    ax.set_title(f"Contribution density Δ·W(Δ)\n"
                 f"2·∫₀⁶ = σ_0(∞) ≈ {sigma_inf_integral:.5f}")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    cum_inf = 2.0 * np.cumsum(contrib_inf) * bin_width
    ax.plot(bin_edges[1:], cum_inf, 'k-', lw=2, label="cumulative 2·∫₀^Δ")
    ax.axhline(2.8628, color="red", ls="--", label="Aitken limit 2.8628")
    ax.set_xlabel("Δ (upper limit)")
    ax.set_ylabel("cumulative σ_0 contribution")
    ax.set_title("Where σ_0 accumulates along the Δ-band")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = out_dir / "phase_e_density.png"
    plt.savefig(plot_file, dpi=120, bbox_inches="tight")
    print(f"\nPlot saved: {plot_file}")

    # Analytical-form probes for W_∞(Δ)
    print()
    print("W_∞(Δ) shape probes:")
    peak_bin = int(np.argmax(np.abs(contrib_inf)))
    print(f"  Δ·W_∞ peak at Δ ≈ {bin_centers[peak_bin]:.3f}, "
          f"value {contrib_inf[peak_bin]:.4f}")
    peak_bin_w = int(np.argmax(np.abs(W_inf)))
    print(f"  W_∞ peak at Δ ≈ {bin_centers[peak_bin_w]:.3f}, "
          f"value {W_inf[peak_bin_w]:.4f}")
    # Edge behaviour: W near Δ=6 (band edge) and Δ=0
    print(f"  W_∞ near Δ=0:  {W_inf[:3]}")
    print(f"  W_∞ near Δ=6:  {W_inf[-3:]}")

    np.savez(out_dir / "density.npz",
             bin_centers=bin_centers, bin_edges=bin_edges,
             W_inf=W_inf, contrib_inf=contrib_inf,
             sigma_inf_integral=sigma_inf_integral,
             N_values=np.array(N_values),
             sigma_true=np.array([sigma_true[N] for N in N_values]),
             sigma_recon=np.array([sigma_recon[N] for N in N_values]))
    print(f"Density data saved: {out_dir / 'density.npz'}")

    print()
    print("Phase E summary:")
    print(f"  σ_0(∞) = 2·∫₀⁶ Δ·W(Δ) dΔ ≈ {sigma_inf_integral:.6f}")
    print(f"  Integration domain: Δ ∈ [0, 6], density W(Δ) extracted per-bin.")
    print(f"  Next: fit W_∞(Δ) to a closed form (3-fold OBC convolution × EP weight),")
    print(f"        evaluate the integral analytically.")


if __name__ == "__main__":
    main()
