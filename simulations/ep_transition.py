#!/usr/bin/env python3
"""The exceptional point, in detail between the anchors: the birth of the rotation.

The F86a slow-pair effective Liouvillian (k=1, the slowest channel pair, γ₀=1):

    L_eff(Q) = [[ -2,        i·J·g_eff ],
                [ i·J·g_eff, -6        ]],   J = Q·γ₀,   x = Q·g_eff/2 = Q/Q_EP

has eigenvalues λ_± = -4γ₀ ± √(4γ₀² − J²g_eff²): two real below the EP, a defective double root
−4γ₀ AT the EP (Q_EP = 2/g_eff), a complex conjugate pair above it. Read through the clock we
built (Takt = the radial decay −Re λ, Rotation = the angle arctan(|Im λ|/|Re λ|) = the F95 angle),
the EP is exactly where the Rotation hand lifts off the Takt axis. This draws the four readings
between the anchors Q=0, Q_EP, Q_peak:

  1. the eigenvalue pair in the complex plane: two reals converging, coalescing at −4γ₀, splitting
     vertically into the oscillating pair (the √-branch-point, the rotation born);
  2. the clock hands vs Q: the decay pins at 4γ₀ at the EP, the Rotation angle lifts off 0;
  3. the defectiveness vs Q: the eigenvectors coalesce (overlap → 1) and the Petermann factor
     spikes (→ ∞) at the EP, the fragile-bridge pinch;
  4. the hardware: IBM Kingston (ep_onset_may2026, job d8drjbfd0j8c73f4mobg) swept Q and the memory
     revival stayed at the 1/N floor until Q crossed ~Q_EP, then lifted off.

Console sibling (clock + defectiveness, text): simulations/f86_ep_through_the_clock.py.
Produces: simulations/results/ep_transition/ep_transition.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

G0 = 1.0
G_EFF = 4.0 / 3.0           # Q_EP = 2/g_eff = 1.5 (the c=2 peak orbit, matching the hardware)
Q_EP = 2.0 / G_EFF
X_PEAK = 2.196910329331     # C2BareDoubledPtfClosedForm: the resonance peak in x = Q/Q_EP units
Q_PEAK = X_PEAK * Q_EP


def l_eff(Q: float, k: int = 1) -> np.ndarray:
    J = Q * G0
    return np.array([[-2.0 * G0 * (2 * k - 1), 1j * J * G_EFF],
                     [1j * J * G_EFF, -2.0 * G0 * (2 * k + 1)]], dtype=complex)


def clock(Q: float):
    """(decay, omega, theta_deg) of the slowest mode, read through the clock."""
    ev = np.linalg.eigvals(l_eff(Q))
    decay = -ev.real
    omega = np.abs(ev.imag)
    i = int(np.argmin(decay))            # the slowest (longest-lived) mode
    return decay[i], omega[i], np.degrees(np.arctan2(omega[i], decay[i]))


def defectiveness(Q: float):
    """(eigenvector overlap |cos|, Petermann factor) of L_eff(Q)."""
    w, V = np.linalg.eig(l_eff(Q))
    v0, v1 = V[:, 0], V[:, 1]
    cosang = abs(np.vdot(v0, v1)) / (np.linalg.norm(v0) * np.linalg.norm(v1))
    Vinv = np.linalg.inv(V)
    petermann = max(float(np.linalg.norm(Vinv[n, :]) ** 2) for n in range(2))
    return min(cosang, 1.0), petermann


def main() -> None:
    print("=" * 74)
    print(f"  EP transition   γ₀={G0}  g_eff={G_EFF:.4f}  Q_EP={Q_EP:.3f}  Q_peak={Q_PEAK:.3f}")
    print("=" * 74)

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    ax_ev, ax_clock, ax_def, ax_hw = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    cmap = matplotlib.colormaps["viridis"]

    # ── Panel 1: the eigenvalue pair in the complex plane (coalescence + lift-off) ──
    Qs = np.linspace(0.0, 2.5 * Q_EP, 400)
    lam_hi = np.empty(len(Qs), dtype=complex)
    lam_lo = np.empty(len(Qs), dtype=complex)
    for j, Q in enumerate(Qs):
        ev = np.linalg.eigvals(l_eff(Q))
        ev = ev[np.argsort(ev.imag)]      # lo = most negative Im, hi = most positive Im
        lam_lo[j], lam_hi[j] = ev[0], ev[1]
    for lam in (lam_hi, lam_lo):
        ax_ev.scatter(lam.real, lam.imag, c=Qs, cmap=cmap, s=10, zorder=4)
    ax_ev.plot(-4.0, 0.0, "*", color="red", markersize=20, markeredgecolor="black",
               markeredgewidth=0.7, zorder=6, label=f"the EP: λ=−4γ₀ (Q_EP={Q_EP:.2f})")
    ax_ev.axhline(0, color="gray", lw=0.4, alpha=0.5)
    ax_ev.set_xlim(-6.5, -1.5)
    ax_ev.set_ylim(-3.2, 3.2)
    ax_ev.set_title("The eigenvalue pair: two reals coalesce at −4γ₀,\nthen split vertically (the rotation born)")
    ax_ev.set_xlabel("Re(λ)  (−decay)")
    ax_ev.set_ylabel("Im(λ)  (±ω, the oscillation)")
    ax_ev.grid(True, alpha=0.2)
    ax_ev.legend(loc="upper left", fontsize=9)
    sc = ax_ev.scatter([], [], c=[], cmap=cmap, vmin=0, vmax=2.5 * Q_EP)
    cb = fig.colorbar(sc, ax=ax_ev, fraction=0.046, pad=0.02)
    cb.set_label("Q", fontsize=9)

    # ── Panel 2: the clock hands vs Q (decay pins, rotation lifts off) ──
    Qc = np.linspace(0.01, 2.5 * Q_EP, 400)
    decay = np.array([clock(Q)[0] for Q in Qc])
    theta = np.array([clock(Q)[2] for Q in Qc])
    ax_clock.plot(Qc, decay, "-", color="#2E8B57", lw=2.2, label="Takt: decay −Re(λ) (pins at 4γ₀)")
    ax_clock.axhline(4.0, color="#2E8B57", ls=":", lw=1.0, alpha=0.6)
    ax_clock.set_xlabel("Q")
    ax_clock.set_ylabel("decay −Re(λ)  [γ₀]", color="#2E8B57")
    ax_clock.set_ylim(0, 6.5)
    ax_t = ax_clock.twinx()
    ax_t.plot(Qc, theta, "-", color="#AA33CC", lw=2.2, label="Rotation: angle θ=arctan(ω/gap) (the F95 angle)")
    ax_t.set_ylabel("Rotation angle θ  [deg]", color="#AA33CC")
    ax_t.set_ylim(0, 70)
    for ax in (ax_clock,):
        ax.axvline(Q_EP, color="red", ls="--", lw=1.2, alpha=0.8)
        ax.axvline(Q_PEAK, color="orange", ls="--", lw=1.0, alpha=0.7)
    ax_clock.annotate("Q_EP\n(rotation\nlifts off)", (Q_EP, 5.6), fontsize=8, color="red", ha="center")
    ax_clock.annotate("Q_peak", (Q_PEAK, 5.9), fontsize=8, color="orange", ha="center")
    ax_clock.set_title("The clock hands: below the EP pure Takt (θ=0),\nat the EP the decay pins at 4γ₀, above it θ lifts off")
    l1, lab1 = ax_clock.get_legend_handles_labels()
    l2, lab2 = ax_t.get_legend_handles_labels()
    ax_clock.legend(l1 + l2, lab1 + lab2, loc="center right", fontsize=8)

    # ── Panel 3: the defectiveness vs Q (the fragile-bridge pinch) ──
    Qd = np.linspace(0.05, 2.5 * Q_EP, 600)
    Qd = Qd[np.abs(Qd - Q_EP) > 1e-3]      # skip exactly the EP (eig is singular there)
    overlap = np.array([defectiveness(Q)[0] for Q in Qd])
    peter = np.array([defectiveness(Q)[1] for Q in Qd])
    ax_def.plot(Qd, overlap, "-", color="#CC3333", lw=2.2, label="eigenvector overlap |⟨v₊|v₋⟩| → 1")
    ax_def.set_xlabel("Q")
    ax_def.set_ylabel("eigenvector overlap |cos|", color="#CC3333")
    ax_def.set_ylim(0, 1.05)
    ax_p = ax_def.twinx()
    ax_p.plot(Qd, peter, "-", color="#333333", lw=1.6, alpha=0.8, label="Petermann factor K → ∞")
    ax_p.set_yscale("log")
    ax_p.set_ylabel("Petermann factor K (log)", color="#333333")
    ax_def.axvline(Q_EP, color="red", ls="--", lw=1.2, alpha=0.8)
    ax_def.annotate("Q_EP: the eigenvectors coalesce\n(defective, Jordan block)", (Q_EP, 0.5),
                    fontsize=8, color="red", ha="center")
    ax_def.set_title("The defectiveness: the eigenvectors collapse to one and the\nPetermann sensitivity spikes at the EP (the fragile-bridge pinch)")
    l1, lab1 = ax_def.get_legend_handles_labels()
    l2, lab2 = ax_p.get_legend_handles_labels()
    ax_def.legend(l1 + l2, lab1 + lab2, loc="upper right", fontsize=8)

    # ── Panel 4: the hardware (Kingston EP onset, Part B) ──
    hw_Q = np.array([0.5, 1.0, 1.5, 2.5, 5.0, 20.0])
    hw_rev = np.array([0.30, 0.36, 0.34, 0.49, 0.56, 0.70])   # data/ibm_ep_onset_may2026 README, job d8drjbfd0j8c73f4mobg
    ax_hw.axhline(1.0 / 3.0, color="gray", ls=":", lw=1.2, alpha=0.7, label="1/N equipartition floor")
    ax_hw.plot(hw_Q, hw_rev, "o-", color="#1F6FB2", lw=1.6, markersize=9, markeredgecolor="black",
               markeredgewidth=0.5, label="IBM Kingston revival (max ⟨n₀⟩)")
    ax_hw.axvline(Q_EP, color="red", ls="--", lw=1.2, alpha=0.8)
    ax_hw.annotate("Q_EP ≈ 1.5", (Q_EP, 0.66), fontsize=9, color="red", ha="center")
    ax_hw.annotate("overdamped\n(forgotten)", (0.7, 0.30), fontsize=8, color="#555", ha="center")
    ax_hw.annotate("memory lifts off →", (4.0, 0.45), fontsize=8, color="#1F6FB2")
    ax_hw.set_xscale("log")
    ax_hw.set_xlim(0.4, 25)
    ax_hw.set_ylim(0.25, 0.75)
    ax_hw.set_title("The hardware: the memory switches on as Q crosses the EP\n(IBM Kingston, 2026-05-31, the rotation born on a real chip)")
    ax_hw.set_xlabel("Q = J/γ₀  (log)")
    ax_hw.set_ylabel("revival (memory return)")
    ax_hw.grid(True, alpha=0.2, which="both")
    ax_hw.legend(loc="upper left", fontsize=8)

    fig.suptitle(
        "The exceptional point, in detail between the anchors (Q=0, Q_EP, Q_peak): the birth of the rotation.\n"
        "Two real decay channels converge, coalesce defectively at the EP (the Takt pins at 4γ₀, the eigenvectors "
        "collapse), and split into an oscillating pair, the Rotation hand lifting off, confirmed on IBM Kingston.",
        y=0.99, fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    out_dir = Path(__file__).parent / "results" / "ep_transition"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "ep_transition.png"
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()

    # console summary of the anchors
    for Q in (0.5 * Q_EP, Q_EP, 1.5 * Q_EP, Q_PEAK):
        d, w, th = clock(Q)
        print(f"  Q={Q:5.2f} (Q/Q_EP={Q/Q_EP:4.2f}): decay={d:5.2f}γ₀  ω={w:5.2f}  θ={th:5.1f}°")
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()
