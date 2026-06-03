#!/usr/bin/env python3
"""The journey between the two singularities, closed: born in Q, dies in tau.

THE_FLOW_BETWEEN_TWO_SINGULARITIES says the flow is bracketed by two singularities
of DIFFERENT kinds, and that difference is the whole point:

  - the EP (the BIRTH): parameter space, at Q_EP = 2/g_eff = 1.5, DEFECTIVE (a Jordan
    block, the two slow eigenvectors coalesce, the Petermann sensitivity diverges). Below
    it: overdamped, no memory (the revival sits on the 1/N floor). Above it: the rotation
    is born, the memory sloshes back.

  - the target (the DEATH): state space, the 1/N equipartitioned fixed point, SIMPLE (the
    lambda = 0 kernel of L). Reached by letting TIME run at a fixed Q above the EP: the
    reborn memory sloshes, fades, and the excitation spreads to 1/N per site.

The carrier gamma_0 is REAL: ~0.05 on the chip, not the toy 1.0. That matters here. Q_EP =
2/g_eff = 1.5 is gamma_0-invariant (J = Q*gamma_0 cancels it), and so is the whole
dimensionless clock (theta, the eigenvector overlap, the flow vs tau): gamma_0 is a pure
SCALE, the journey's SHAPE does not depend on it. But the physical coupling does: J_EP =
Q_EP*gamma_0 = 0.075, while the chip runs at J = 1.5 rad/us, i.e. Q = J/gamma_0 ~ 30, deep
above the EP. Setting gamma_0 = 1 hides this by forcing J = Q, which makes J_EP = 1.5
coincide with the chip's J = 1.5 and falsely reads as "the chip sits on the EP". It does not:
the chip lives in the deep-memory regime, and the EP is reached only by INJECTING noise
(Part B pushed Q down from ~30 toward 1.5). That is the important thing the toy value buried.

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
# The REAL scale. gamma_0 is the chip's carrier (~0.05), not the toy 1.0.
# ----------------------------------------------------------------------------------------
G0 = 0.05                       # the hardware carrier Gamma_0 ~ 0.05 (the real value)
G_EFF = 4.0 / 3.0               # Q_EP = 2/g_eff = 1.5, gamma_0-invariant
Q_EP = 2.0 / G_EFF
X_PEAK = 2.196910329331         # C2BareDoubledPtfClosedForm resonance peak, in x = Q/Q_EP units
Q_PEAK = X_PEAK * Q_EP
J_HW = 1.5                      # the hardware coupling (rad/us), fixed; Part B injected gamma to scan Q
Q_HOME = J_HW / G0              # the chip's natural operating point ~ 30 (deep memory, no injection)


def l_eff(Q: float, gamma0: float = G0) -> np.ndarray:
    """The slow-pair effective Liouvillian at coupling Q (J = Q*gamma0, threaded honestly)."""
    J = Q * gamma0
    return np.array([[-2.0 * gamma0, 1j * J * G_EFF],
                     [1j * J * G_EFF, -6.0 * gamma0]], dtype=complex)


def clock(Q: float, gamma0: float = G0):
    """(decay -Re lambda, omega |Im lambda|, theta deg, eigenvector overlap) of the slow mode."""
    w, V = np.linalg.eig(l_eff(Q, gamma0))
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
# Leg 2: the death (state space). N=3 single excitation, XY chain under Z-dephasing,
# the flow to 1/N. Read dimensionlessly (gamma_0-invariant): L'(Q) = -iQ[H_unit,.] +
# Sum_l(Z_l rho Z_l - rho), tau = gamma_0 * t. (Any gamma_0 gives the same curve vs tau.)
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
    print("=" * 84)
    print(f"  THE JOURNEY  gamma0={G0}  g_eff={G_EFF:.4f}  Q_EP={Q_EP:.3f}  Q_peak={Q_PEAK:.3f}")
    print(f"  J_EP = Q_EP*gamma0 = {Q_EP * G0:.4f}   |   chip: J={J_HW} -> Q_home = J/gamma0 = {Q_HOME:.1f}")
    print("=" * 84)

    # ---- gamma_0-invariance check: the dimensionless clock is identical at 1.0 and 0.05 ----
    print("\n  gamma_0-INVARIANCE of the clock (theta, overlap) vs Q -- the carrier is a pure scale:")
    print(f"  {'Q':>6} | {'theta(g0=1.0)':>13} {'theta(g0=0.05)':>14} | {'overlap(1.0)':>12} {'overlap(0.05)':>13}")
    for Q in (0.75, 1.5, 2.5, 5.0, 20.0):
        _, _, th1, ov1 = clock(Q, 1.0)
        _, _, th2, ov2 = clock(Q, 0.05)
        print(f"  {Q:6.2f} | {th1:13.4f} {th2:14.4f} | {ov1:12.5f} {ov2:13.5f}")

    # ---- Leg 1: the birth, the closed-form clock, with PHYSICAL J alongside Q ----
    print("\n  LEG 1 - the birth (parameter space). Q, x=Q/Q_EP, physical J=Q*gamma0, the clock:")
    print(f"  {'Q':>6} {'x':>6} {'J':>7} {'decay':>7} {'omega':>7} {'theta':>7} {'overlap':>8}   regime")
    for Q in [0.3, 0.5, 0.75, 1.0, 1.5, 2.5, Q_PEAK, 5.0, 10.0, 20.0, Q_HOME]:
        d, w, th, ov = clock(Q)
        x = Q / Q_EP
        J = Q * G0
        if Q < Q_EP - 1e-9:
            reg = "overdamped (pre-birth)"
        elif abs(Q - Q_EP) < 1e-6:
            reg = "THE EP (defective pinch)"
        elif abs(Q - Q_HOME) < 1e-6:
            reg = "<- the chip's natural home (deep memory)"
        else:
            reg = "rotation born"
        print(f"  {Q:6.2f} {x:6.2f} {J:7.4f} {d:7.4f} {w:7.4f} {th:7.1f} {ov:8.3f}   {reg}")

    # ---- Leg 2: the death, the flow to 1/N at a representative Q above the EP ----
    N = 3
    Q_death = 5.0
    taus = np.linspace(0.0, 4.0, 240)
    flow = flow_to_target(Q_death, N, taus)
    print(f"\n  LEG 2 - the death (state space, tau-axis), N={N}, Q={Q_death} (above the EP):")
    for tau in (0.0, 0.2, 0.5, 1.0, 2.0, 4.0):
        j = int(np.argmin(np.abs(taus - tau)))
        print(f"  tau={taus[j]:5.2f}  " + " ".join(f"<n{l}>={flow[l, j]:5.3f}" for l in range(N)))
    print(f"  target 1/N = {1.0 / N:.4f}  (the chip's home Q~{Q_HOME:.0f} sloshes far longer, same arc)")

    # ====================================================================================
    # The figure: the two legs on one canvas
    # ====================================================================================
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(16, 6.6))

    # ---- Panel A (left): the birth, parameter space, the Q-axis (to the chip's home ~30) ----
    Qc = np.logspace(np.log10(0.3), np.log10(40.0), 500)
    theta = np.array([clock(Q)[2] for Q in Qc])
    overlap = np.array([clock(Q)[3] for Q in Qc])

    axL.axhline(1.0 / 3.0, color="gray", ls=":", lw=1.2, alpha=0.7)
    axL.plot(HW_Q, HW_REV, "o-", color="#1F6FB2", lw=1.8, markersize=9, markeredgecolor="black",
             markeredgewidth=0.5, label="IBM Kingston revival (Part B: inject noise, scan Q)", zorder=5)
    axL.annotate("1/N floor\n(forgotten)", (0.34, 0.335), fontsize=8, color="#555", ha="center", va="bottom")
    axL.set_xscale("log")
    axL.set_xlim(0.3, 40)
    axL.set_ylim(0.0, 1.0)
    axL.set_xlabel("Q = J / gamma_0   (the chip rests at Q~30; the EP is reached by ADDING noise, lowering Q)")
    axL.set_ylabel("revival  (memory return)", color="#1F6FB2")
    axL.axvline(Q_EP, color="red", ls="--", lw=1.3, alpha=0.85)
    axL.axvline(Q_PEAK, color="orange", ls="--", lw=1.0, alpha=0.6)
    axL.axvline(Q_HOME, color="#2E8B57", ls="-.", lw=1.6, alpha=0.85)
    axL.annotate("Q_EP = 1.5\n(J_EP = 0.075)\nthe defective pinch", (1.5, 0.085), fontsize=8.5,
                 color="red", ha="center")
    axL.annotate("Q_peak", (Q_PEAK, 0.50), fontsize=8, color="orange", ha="center", rotation=90)
    axL.annotate("the chip's home\nJ=1.5, gamma_0~0.05\nQ~30 (deep memory)", (30, 0.27), fontsize=8.5,
                 color="#1d6b3f", ha="center")

    # physical coupling J on a secondary top axis: J_EP=0.075 vs the chip's J=1.5 both visible
    secax = axL.secondary_xaxis("top", functions=(lambda q: q * G0, lambda j: j / G0))
    secax.set_xlabel("physical coupling  J = Q * gamma_0   [rad/us]   (gamma_0 ~ 0.05, the chip's carrier)")

    axLt = axL.twinx()
    axLt.plot(Qc, theta, "-", color="#AA33CC", lw=2.2, label="Rotation angle theta (the F95 angle)")
    axLt.plot(Qc, overlap * 90.0, "-", color="#CC3333", lw=1.8, alpha=0.8,
              label="eigenvector overlap min(x,1/x)  [x90, peaks =1 at the EP]")
    axLt.set_ylabel("theta [deg]   /   overlap [x90]", color="#7733AA")
    axLt.set_ylim(0, 95)
    axLt.annotate("theta -> 90 deg\n(pure rotation)", (24.0, 70), fontsize=8, color="#AA33CC", ha="center")

    lA, labA = axL.get_legend_handles_labels()
    lAt, labAt = axLt.get_legend_handles_labels()
    axL.legend(lA + lAt, labA + labAt, loc="center left", fontsize=7.8, framealpha=0.9)
    axL.set_title("THE BIRTH  -  parameter space.  The dimensionless clock is gamma_0-invariant\n"
                  "(a pure scale), but the chip rests at Q~30: the EP is a NOISE-degraded regime.",
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
        "The journey between the two singularities, closed: born in Q, dies in tau   (gamma_0 ~ 0.05, the real carrier).\n"
        "Left: the BIRTH in parameter space; the clock's SHAPE is gamma_0-invariant, but the chip rests at Q~30 and the "
        "EP (Q=1.5, J=0.075) is reached only by injecting noise (Kingston Part B).\n"
        "Right: the DEATH in state space; the reborn memory relaxing over time to the simple 1/N kernel (Kingston Part A).",
        y=1.0, fontsize=10.5)
    plt.tight_layout(rect=[0, 0, 1, 0.90])

    out_dir = Path(__file__).parent / "results" / "journey_between_singularities"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "journey.png"
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"\n  saved: {out}")


if __name__ == "__main__":
    main()
