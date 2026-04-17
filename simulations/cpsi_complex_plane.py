#!/usr/bin/env python3
"""
CΨ as a complex quantity: the 2D → 3D step for BOUNDARY_NAVIGATION.

The original CΨ = Tr(ρ²) · L1/(d−1) uses the L1 norm of off-diagonals,
which is real and positive. This puts every Bell+ trajectory on the
real axis of the Mandelbrot c-plane, a 1D path.

Complex extension: CΨ_com = Tr(ρ²) · (Σ ρ_{ij} off-diagonal, signed) / (d−1).
Now CΨ_com is complex in general. For Bell+ with phase φ on |11⟩:

    ρ_{00,11}(0) = (1/2) · exp(iφ)

Under Lindblad Z-dephasing + Z-drift (detuning δ):

    ρ_{00,11}(t) = (1/2) · exp(iφ + iΔt − 4γt)

where Δ = 2(δ_1 − δ_2). The trajectory in the c-plane is a **logarithmic
spiral** winding inward toward the origin, passing near the cardioid
cusp at c = 1/4 when |CΨ_com| ≈ 1/3.

This document generates the visualization. Main cardioid of the
Mandelbrot set drawn analytically; trajectories computed from Lindblad
integration for several (γ, Δ, φ_0) triples.

Date: 2026-04-16
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ── Constants ─────────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op: np.ndarray, site: int, n: int) -> np.ndarray:
    factors = [I2] * n
    factors[site] = op
    out = factors[0]
    for f in factors[1:]:
        out = np.kron(out, f)
    return out


def liouvillian(h: np.ndarray, jump_ops: list[np.ndarray]) -> np.ndarray:
    d = h.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, h) - np.kron(h.T, Id))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Id, LdL)
              - 0.5 * np.kron(LdL.T, Id))
    return L


# ── CΨ metrics (real and complex) ──
def cpsi_real(rho: np.ndarray) -> float:
    d = rho.shape[0]
    C = float(np.real(np.trace(rho @ rho)))
    diag = np.diag(np.diag(rho))
    L1 = float(np.sum(np.abs(rho - diag)))
    return C * L1 / (d - 1)


def cpsi_complex(rho: np.ndarray) -> complex:
    """Complex extension: Psi_com = sum of off-diagonals (signed), /(d-1)."""
    d = rho.shape[0]
    C = float(np.real(np.trace(rho @ rho)))
    # For 2-qubit Bell-like states the dominant off-diagonal is ρ_{0,3}
    # (|00⟩⟨11|). We sum all upper-triangle off-diagonals with their sign.
    mask_upper = np.triu(np.ones_like(rho, dtype=bool), k=1)
    Psi_com = complex(2.0 * np.sum(rho[mask_upper]) / (d - 1))
    return C * Psi_com


# ── State preparation ─────────────────────────────────────────────
def bell_plus_phased(phi: float) -> np.ndarray:
    """Bell+ with a phase: (|00⟩ + exp(iφ)|11⟩) / √2."""
    psi = np.zeros(4, dtype=complex)
    psi[0] = 1 / np.sqrt(2)
    psi[3] = np.exp(1j * phi) / np.sqrt(2)
    return np.outer(psi, psi.conj())


# ── Trajectory computation ────────────────────────────────────────
def trajectory(
    gamma: float,
    detuning_sum: float,  # Ω = δ₁ + δ₂: the rotation rate of ρ_{00,11}
    phi_0: float,
    t_max: float,
    n_steps: int = 400,
) -> dict:
    """Two-qubit Bell+ with phase φ_0, under Z-dephasing rate γ per site,
    plus common Z-Hamiltonian rotating ρ_{00,11} at rate Ω. The phase of
    |00⟩⟨11| advances as exp(-i·Ω·t); combined with dephasing this gives
    a logarithmic spiral in the c-plane.
    """
    # For H = (Ω/2)(Z_0 + Z_1)/2 = (Ω/4)(Z_0 + Z_1), the eigenvalue gap
    # between |00⟩ and |11⟩ is Ω, so ρ_{00,11} rotates at Ω.
    h = (detuning_sum / 4.0) * (site_op(Z, 0, 2) + site_op(Z, 1, 2))

    jumps = [np.sqrt(gamma) * site_op(Z, 0, 2),
             np.sqrt(gamma) * site_op(Z, 1, 2)]
    L = liouvillian(h, jumps)

    rho0 = bell_plus_phased(phi_0)
    rho_vec0 = rho0.flatten(order="F")
    times = np.linspace(0.0, t_max, n_steps)
    cpsi_c_list = []
    cpsi_r_list = []
    for t in times:
        rho_vec_t = expm(L * t) @ rho_vec0
        rho_t = rho_vec_t.reshape(4, 4, order="F")
        cpsi_c_list.append(cpsi_complex(rho_t))
        cpsi_r_list.append(cpsi_real(rho_t))
    return {
        "times": times.tolist(),
        "cpsi_complex": [complex(z) for z in cpsi_c_list],
        "cpsi_real": cpsi_r_list,
        "gamma": gamma,
        "detuning_sum": detuning_sum,
        "phi_0": phi_0,
    }


# ── Mandelbrot cardioid boundary ──
def main_cardioid(n_pts: int = 400) -> np.ndarray:
    """The main cardioid of the Mandelbrot set (period-1 center):
    c(θ) = e^{iθ}/2 − e^{2iθ}/4, θ ∈ [0, 2π]. Cusp at c = 1/4.
    """
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return (np.exp(1j * theta) / 2.0) - (np.exp(2j * theta) / 4.0)


def period_2_bulb(n_pts: int = 200) -> np.ndarray:
    """Period-2 bulb: circle of radius 1/4 centered at c = −1."""
    theta = np.linspace(0, 2 * np.pi, n_pts)
    return -1.0 + 0.25 * np.exp(1j * theta)


# ── Plotting ──────────────────────────────────────────────────────
def plot_all(trajectories: list[dict], out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    ax_full, ax_zoom = axes

    # Mandelbrot cardioid in both panels
    c_card = main_cardioid()
    c_bulb = period_2_bulb()
    for ax in axes:
        ax.plot(c_card.real, c_card.imag, "-", color="#888", linewidth=0.9,
                alpha=0.7, label="Mandelbrot cardioid (period-1)")
        ax.plot(c_bulb.real, c_bulb.imag, "-", color="#888", linewidth=0.6,
                alpha=0.5, label="period-2 bulb")
        ax.plot(0.25, 0, "o", color="red", markersize=9, zorder=5,
                label="cusp c = 1/4")
        ax.axhline(0, color="gray", linewidth=0.3, alpha=0.5)
        ax.axvline(0, color="gray", linewidth=0.3, alpha=0.5)

    # Trajectories
    colors = ["#CC3333", "#33AACC", "#EE9922", "#AA33CC", "#33AA77"]
    for i, tr in enumerate(trajectories):
        c_vals = np.array(tr["cpsi_complex"])
        col = colors[i % len(colors)]
        gamma = tr["gamma"]
        omega = tr["detuning_sum"]
        phi0 = tr["phi_0"]
        winds = omega / (4 * gamma) if gamma > 0 else float("inf")
        label = (f"γ={gamma:.3f}, Ω={omega:.2f}, φ₀={phi0:.2f}"
                 f"  (Ω/4γ={winds:.1f})")
        for ax in axes:
            ax.plot(c_vals.real, c_vals.imag, "-", color=col,
                    linewidth=1.8, alpha=0.85, label=label)
            ax.plot(c_vals.real[0], c_vals.imag[0], "o", color=col,
                    markersize=8)
            ax.plot(c_vals.real[-1], c_vals.imag[-1], "x", color=col,
                    markersize=10, markeredgewidth=2)

    ax_full.set_xlim(-1.7, 0.6)
    ax_full.set_ylim(-0.8, 0.8)
    ax_full.set_aspect("equal")
    ax_full.set_title("Full c-plane: Mandelbrot + CΨ-trajectories")
    ax_full.grid(True, alpha=0.2)
    ax_full.legend(loc="lower left", fontsize=7)
    ax_full.set_xlabel("Re(CΨ_com)")
    ax_full.set_ylabel("Im(CΨ_com)")

    # Zoom to the cusp region
    ax_zoom.set_xlim(-0.05, 0.45)
    ax_zoom.set_ylim(-0.25, 0.25)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_title("Zoom to the cusp c = 1/4 (the fold entrance)")
    ax_zoom.grid(True, alpha=0.2)
    ax_zoom.set_xlabel("Re(CΨ_com)")
    ax_zoom.set_ylabel("Im(CΨ_com)")

    fig.suptitle(
        "CΨ in the complex plane: the 2D extension of BOUNDARY_NAVIGATION\n"
        "φ₀ = 0 stays on the real axis; φ₀ ≠ 0 spirals through the cusp",
        y=1.00,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_png}")


# ── Run ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    RESULTS = Path(__file__).parent / "results"
    RESULTS.mkdir(exist_ok=True)
    OUT = RESULTS / "cpsi_complex_plane.png"

    print("=" * 70)
    print("  CΨ in the complex plane: 2D extension of the fold navigation")
    print("=" * 70)

    # Five trajectories with different (γ, Ω, φ₀)
    configs = [
        {"gamma": 0.05, "detuning_sum": 0.0, "phi_0": 0.0, "t_max": 20.0},         # real-axis baseline
        {"gamma": 0.05, "detuning_sum": 0.0, "phi_0": np.pi/4, "t_max": 20.0},     # fixed phase, ray
        {"gamma": 0.05, "detuning_sum": 0.4, "phi_0": 0.0, "t_max": 20.0},         # slow spiral (Ω/4γ=2)
        {"gamma": 0.05, "detuning_sum": 1.5, "phi_0": 0.0, "t_max": 20.0},         # fast spiral (Ω/4γ=7.5)
        {"gamma": 0.02, "detuning_sum": 0.5, "phi_0": np.pi/3, "t_max": 40.0},     # long-lived spiral
    ]

    trajs = []
    for cfg in configs:
        print(f"\n  running: γ={cfg['gamma']}, Ω={cfg['detuning_sum']}, "
              f"φ₀={cfg['phi_0']:.3f}, t_max={cfg['t_max']}")
        tr = trajectory(**cfg)
        trajs.append(tr)
        c0 = tr["cpsi_complex"][0]
        cend = tr["cpsi_complex"][-1]
        print(f"    CΨ_com(0)  = {c0.real:+.4f} + {c0.imag:+.4f}i  "
              f"|CΨ_com|={abs(c0):.4f}")
        print(f"    CΨ_com(end) = {cend.real:+.4f} + {cend.imag:+.4f}i  "
              f"|CΨ_com|={abs(cend):.4f}")

    plot_all(trajs, OUT)

    print()
    print("Interpretation:")
    print("  φ₀=0, Ω=0: trajectory stays on positive real axis, crosses")
    print("  the cusp at c = 1/4 head-on. This is BOUNDARY_NAVIGATION.md.")
    print("  φ₀=0, Ω>0 (common Z-drift): trajectory spirals through the cusp,")
    print("  winding around c=0 as it decays. Winding number = Ω/(4γ).")
    print("  φ₀≠0: trajectory starts off the real axis and spirals if Ω≠0,")
    print("  or decays along a ray if Ω=0.")
    print()
    print("The cusp at c = 1/4 remains the topology-changing point, but the")
    print("path through it is now 2D instead of 1D. The 3D view is this")
    print("2D c-plane × time: a spiral staircase descending into the cusp.")
