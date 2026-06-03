#!/usr/bin/env python3
"""The journey between the two singularities, closed: born in Q, dies in tau.

THE_FLOW_BETWEEN_TWO_SINGULARITIES says the flow is bracketed by two singularities
of DIFFERENT kinds, and that difference is the whole point:

  - the EP (the BIRTH): parameter space, at Q_EP = 2/g_eff = 1.5, DEFECTIVE (a Jordan
    block, the two slow eigenvectors coalesce, the Petermann sensitivity diverges). You
    reach it by tuning the coupling: gamma_0 is the universal carrier and stays constant,
    so raising Q = J/gamma_0 means raising J. Below it: overdamped, no memory (the revival
    sits on the 1/N floor). Above it: the rotation is born, the memory sloshes back.

  - the target (the DEATH): state space, the 1/N equipartitioned fixed point, SIMPLE (the
    lambda = 0 kernel of L). You reach it by letting TIME run at a fixed Q above the EP: the
    reborn memory sloshes, fades, and the excitation spreads to 1/N per site.

So the journey is not a single 1D path. It is born in Q (parameter space, the defective
EP) and dies in tau (state space, the simple 1/N kernel). This script draws both legs on
one canvas, with both hardware halves (IBM Kingston, 2026-05-31): Part B (revival vs Q, the
birth) on the left, Part A (the population flow to 1/N, the death) on the right.

The closed forms sit at the two endpoints; the middle of the second leg is the un-closed-form
loop. We close the ARC (both legs, both singularities, both hardware halves in one view); the
deep reading of the middle stays parked, per the experiment doc.

Produces: simulations/results/journey_between_singularities/journey.png
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

# ----------------------------------------------------------------------------------------
# Leg 1: the birth (parameter space, the Q-axis). The slow-pair effective Liouvillian.
# ----------------------------------------------------------------------------------------
G0 = 1.0
G_EFF = 4.0 / 3.0          # Q_EP = 2/g_eff = 1.5 (the c=2 peak orbit, matching the hardware)
Q_EP = 2.0 / G_EFF
X_PEAK = 2.196910329331    # C2BareDoubledPtfClosedForm resonance peak, in x = Q/Q_EP units
Q_PEAK = X_PEAK * Q_EP


def l_eff(Q: float) -> np.ndarray:
    J = Q * G0
    return np.array([[-2.0 * G0, 1j * J * G_EFF],
                     [1j * J * G_EFF, -6.0 * G0]], dtype=complex)


def clock(Q: float):
    """(decay -Re lambda, omega |Im lambda|, theta deg, eigenvector overlap) of the slow mode."""
    w, V = np.linalg.eig(l_eff(Q))
    decay = -w.real
    omega = np.abs(w.imag)
    i = int(np.argmin(decay))                 # the slowest (longest-lived) mode
    v0, v1 = V[:, 0], V[:, 1]
    overlap = abs(np.vdot(v0, v1)) / (np.linalg.norm(v0) * np.linalg.norm(v1))
    return decay[i], omega[i], np.degrees(np.arctan2(omega[i], decay[i])), min(overlap, 1.0)


# Part B hardware: the revival (max <n_0> for t >= 2us) vs Q (job d8drjbfd0j8c73f4mobg)
HW_Q = np.array([0.5, 1.0, 1.5, 2.5, 5.0, 20.0])
HW_REV = np.array([0.30, 0.36, 0.34, 0.49, 0.56, 0.70])

# ----------------------------------------------------------------------------------------
# Leg 2: the death (state space, the tau-axis at a fixed Q above the EP). N=3 single
# excitation, XY chain under Z-dephasing, the flow to 1/N. (gamma = 1, J = Q, tau = gamma*t.)
# ----------------------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op: np.ndarray, l: int, N: int) -> np.ndarray:
    mats = [op if k == l else I2 for k in range(N)]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def flow_to_target(Q: float, N: int, taus: np.ndarray):
    """<n_site>(tau) for the single excitation on site 0, relaxing to 1/N. Returns array [site][tau]."""
    d = 1 << N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        H += Q * (site_op(X, b, N) @ site_op(X, b + 1, N) + site_op(Y, b, N) @ site_op(Y, b + 1, N))
    Id = np.eye(d, dtype=complex)
    # column-stacking vec: vec(A rho B) = (B^T kron A) vec(rho)
    Lsup = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = site_op(Z, l, N)
        Lsup += np.kron(Zl.T, Zl) - np.kron(Id, Id)
    # initial state: excitation on site 0 (leftmost = MSB), index 2^(N-1)
    idx = 1 << (N - 1)
    rho0 = np.zeros((d, d), dtype=complex)
    rho0[idx, idx] = 1.0
    vec0 = rho0.reshape(-1, order="F")
    w, V = np.linalg.eig(Lsup)
    c = np.linalg.solve(V, vec0)
    n_ops = [(Id - site_op(Z, l, N)) * 0.5 for l in range(N)]
    out = np.zeros((N, len(taus)))
    for j, tau in enumerate(taus):
        vt = V @ (np.exp(w * tau) * c)
        rho = vt.reshape(d, d, order="F")
        for l in range(N):
            out[l, j] = float(np.real(np.trace(n_ops[l] @ rho)))
    return out


def main() -> None:
    print("=" * 80)
    print(f"  THE JOURNEY  gamma0={G0}  g_eff={G_EFF:.4f}  Q_EP={Q_EP:.3f}  Q_peak={Q_PEAK:.3f}")
    print("=" * 80)

    # ---- Leg 1: the birth, the full Q-axis to 20 (the stretch the EP scout did not reach) ----
    print("\n  LEG 1 - the birth (parameter space, Q-axis), the closed-form clock:")
    print(f"  {'Q':>6} {'x=Q/Q_EP':>9} {'decay':>7} {'omega':>7} {'theta':>7} {'overlap':>8}   regime")
    Qs_table = [0.3, 0.5, 0.75, 1.0, 1.5, 2.5, Q_PEAK, 5.0, 10.0, 20.0]
    for Q in Qs_table:
        d, w, th, ov = clock(Q)
        x = Q / Q_EP
        if Q < Q_EP - 1e-9:
            reg = "overdamped (theta=0, pre-birth)"
        elif abs(Q - Q_EP) < 1e-6:
            reg = "THE EP (defective pinch, overlap=1)"
        else:
            reg = "rotation born (theta lifts, overlap falls)"
        print(f"  {Q:6.2f} {x:9.3f} {d:7.3f} {w:7.3f} {th:7.1f} {ov:8.3f}   {reg}")

    print("\n  Part B hardware revival (IBM Kingston, max <n_0> vs Q):")
    for q, r in zip(HW_Q, HW_REV):
        print(f"  Q={q:5.1f}  revival={r:.2f}   {'<= 1/N floor (forgotten)' if r < 0.40 else 'memory present'}")

    # ---- Leg 2: the death, the flow to 1/N at a representative Q above the EP ----
    N = 3
    Q_death = 5.0
    taus = np.linspace(0.0, 4.0, 240)
    flow = flow_to_target(Q_death, N, taus)
    print(f"\n  LEG 2 - the death (state space, tau-axis), N={N}, Q={Q_death} (above the EP):")
    print(f"  {'tau':>6} " + " ".join(f"<n{l}>" for l in range(N)))
    for tau in (0.0, 0.2, 0.5, 1.0, 2.0, 4.0):
        j = int(np.argmin(np.abs(taus - tau)))
        print(f"  {taus[j]:6.2f} " + " ".join(f"{flow[l, j]:5.3f}" for l in range(N)))
    print(f"  target 1/N = {1.0 / N:.4f} (every site converges here: the equipartitioned death)")

    # ====================================================================================
    # The figure: the two legs on one canvas
    # ====================================================================================
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(16, 6.5))

    # ---- Panel A (left): the birth, parameter space, the Q-axis ----
    Qc = np.logspace(np.log10(0.3), np.log10(20.0), 400)
    theta = np.array([clock(Q)[2] for Q in Qc])
    overlap = np.array([clock(Q)[3] for Q in Qc])

    axL.axhline(1.0 / 3.0, color="gray", ls=":", lw=1.2, alpha=0.7)
    axL.plot(HW_Q, HW_REV, "o-", color="#1F6FB2", lw=1.8, markersize=10, markeredgecolor="black",
             markeredgewidth=0.5, label="IBM Kingston revival (memory return)", zorder=5)
    axL.annotate("1/N floor\n(forgotten)", (0.34, 0.335), fontsize=8, color="#555", ha="center", va="bottom")
    axL.set_xscale("log")
    axL.set_xlim(0.3, 20)
    axL.set_ylim(0.0, 1.0)
    axL.set_xlabel("Q = J / gamma_0   (gamma_0 const -> raising Q raises the coupling J)")
    axL.set_ylabel("revival  (memory return)", color="#1F6FB2")
    axL.axvline(Q_EP, color="red", ls="--", lw=1.3, alpha=0.85)
    axL.axvline(Q_PEAK, color="orange", ls="--", lw=1.0, alpha=0.7)
    axL.annotate("Q_EP = 1.5\nthe defective pinch\n(the rotation born)", (1.5, 0.10), fontsize=8.5,
                 color="red", ha="center")
    axL.annotate("Q_peak", (Q_PEAK, 0.04), fontsize=8, color="orange", ha="center")

    axLt = axL.twinx()
    axLt.plot(Qc, theta, "-", color="#AA33CC", lw=2.2, label="Rotation angle theta (the F95 angle)")
    axLt.plot(Qc, overlap * 90.0, "-", color="#CC3333", lw=1.8, alpha=0.8,
              label="eigenvector overlap min(x,1/x)  [x90, peaks =1 at the EP]")
    axLt.set_ylabel("theta [deg]   /   overlap [x90]", color="#7733AA")
    axLt.set_ylim(0, 95)
    axLt.annotate("theta -> 90 deg\n(pure rotation)", (13.0, 66), fontsize=8, color="#AA33CC", ha="center")

    lA, labA = axL.get_legend_handles_labels()
    lAt, labAt = axLt.get_legend_handles_labels()
    axL.legend(lA + lAt, labA + labAt, loc="center left", fontsize=8, framealpha=0.9)
    axL.set_title("THE BIRTH  -  parameter space (tune Q through the EP)\n"
                  "below: overdamped, no memory.  at Q_EP: the defective pinch.  above: the rotation born.",
                  fontsize=10)

    # ---- Panel B (right): the death, state space, the tau-axis ----
    colors = ["#1F6FB2", "#2E8B57", "#CC7722"]
    labels = ["site 0 (edge, the drop)", "site 1 (bulk)", "site 2 (edge)"]
    for l in range(N):
        axR.plot(taus, flow[l], "-", color=colors[l], lw=2.2, label=labels[l])
    axR.axhline(1.0 / N, color="black", ls="--", lw=1.4, alpha=0.8)
    axR.annotate("1/N target  (the simple lambda=0 kernel:\nthe equipartitioned death, the memory forgotten)",
                 (2.0, 1.0 / N), xytext=(1.7, 0.52), fontsize=8.5, ha="center",
                 arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    axR.set_xlim(0, 4)
    axR.set_ylim(0, 1.0)
    axR.set_xlabel("tau = gamma_0 * t   (time, at fixed Q above the EP)")
    axR.set_ylabel("per-site occupation  <n_site>")
    axR.set_title(f"THE DEATH  -  state space (let time run, Q={Q_death:.0f} above the EP)\n"
                  "the reborn memory sloshes, fades, and spreads to 1/N: the flow into the future.",
                  fontsize=10)
    axR.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    axR.grid(True, alpha=0.2)

    fig.suptitle(
        "The journey between the two singularities, closed: born in Q, dies in tau.\n"
        "Left: the BIRTH in parameter space, the rotation switching on as the coupling crosses the defective EP "
        "(IBM Kingston, Part B).\n"
        "Right: the DEATH in state space, that reborn memory relaxing over time to the simple 1/N kernel "
        "(the flow watched on Kingston, Part A).",
        y=1.0, fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.92])

    out_dir = Path(__file__).parent / "results" / "journey_between_singularities"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "journey.png"
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"\n  saved: {out}")


if __name__ == "__main__":
    main()
