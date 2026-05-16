"""Phase 1 of the F86 Γ_pair derivation attempt:

Scan HWHM_left/Q_peak for a GENERIC 2×2 cluster-pair sub-block over its parameter
space (a, b, X, Δδ). The bare-doubled-PTF closed form 0.671535 (today's
C2BareDoubledPtfClosedForm) is the special case (a=-2, b=-6, Δδ=0).

The cluster-pair sub-block in JW basis at the (i, j)-pair:

    L_sub(J) = [[ γa + iJ·Δδ/2,  +iJ·X     ],
                [ +iJ·X,         γb − iJ·Δδ/2 ]]

with γ = γ₀ = 1, t_peak = 1/(4γ₀) = 1/4. K_b(J) = 2·Re⟨ρ(t_peak)|dρ/dJ⟩ with
probe ρ_0 = (1, 0) and V = ∂L/∂J. Direct numerical Duhamel evaluator (same
machinery as today's brute test), parabolic Q_peak refinement, leftward
half-max search for HWHM_left.

Goal: see if HWHM_left/Q_peak is universal at 0.671535 (bare floor) across
generic (a, b, X, Δδ), or if it depends on the parameters. If parameter-
dependent, identify the dependence pattern.
"""

from __future__ import annotations

import sys
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import math
import numpy as np
from scipy.optimize import brentq


def k_b_cluster_pair(J, a, b, X, dDelta, gamma=1.0, t_peak=None):
    """K_b(J) for the 2x2 cluster-pair sub-block.

    L_sub(J) = [[γa + iJ·Δδ/2, +iJ·X], [+iJ·X, γb − iJ·Δδ/2]]
    V = ∂L/∂J = [[+i·Δδ/2, +i·X], [+i·X, −i·Δδ/2]]
    """
    if t_peak is None:
        t_peak = 1.0 / (4.0 * gamma)
    L = np.array([
        [gamma * a + 1j * J * dDelta / 2, +1j * J * X],
        [+1j * J * X, gamma * b - 1j * J * dDelta / 2],
    ], dtype=complex)
    Vb = np.array([
        [+1j * dDelta / 2, +1j * X],
        [+1j * X, -1j * dDelta / 2],
    ], dtype=complex)
    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    rho0 = np.array([1.0, 0.0], dtype=complex)
    expLam = np.exp(evals * t_peak)
    rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)
    n = 2
    I_mat = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            diff = evals[j] - evals[i]
            if abs(diff) < 1e-10:
                I_mat[i, j] = t_peak * expLam[i]
            else:
                I_mat[i, j] = (expLam[j] - expLam[i]) / diff
    V_eig = Rinv @ Vb @ R
    c0 = Rinv @ rho0
    drho_eig = np.zeros(n, dtype=complex)
    for i in range(n):
        for j in range(n):
            drho_eig[i] += V_eig[i, j] * I_mat[i, j] * c0[j]
    drho = R @ drho_eig
    return 2.0 * np.real(np.vdot(rho_t, drho))


def compute_q_peak_hwhm(a, b, X, dDelta, gamma=1.0, q_max=8.0, n_grid=4000):
    """Numerical Q_peak, |K_b|_peak, HWHM_left via grid scan + parabolic refinement."""
    Q_grid = np.linspace(0.01, q_max, n_grid)
    K_vals = np.array([k_b_cluster_pair(Q, a, b, X, dDelta, gamma) for Q in Q_grid])
    K_abs = np.abs(K_vals)
    i_max = int(np.argmax(K_abs))
    if i_max == 0 or i_max == n_grid - 1:
        return None  # peak at boundary
    y0, y1, y2 = K_abs[i_max - 1], K_abs[i_max], K_abs[i_max + 1]
    denom = y0 - 2 * y1 + y2
    delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-15 else 0.0
    dQ = Q_grid[1] - Q_grid[0]
    Q_peak = Q_grid[i_max] + delta * dQ

    # Brent-refine K_max
    def neg_K_abs(q):
        return -abs(k_b_cluster_pair(q, a, b, X, dDelta, gamma))
    try:
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(neg_K_abs, bracket=(Q_grid[max(0, i_max-3)], Q_peak, Q_grid[min(n_grid-1, i_max+3)]), method="brent", options={"xtol": 1e-10})
        Q_peak_refined = res.x
        K_max = -res.fun
    except Exception:
        Q_peak_refined = Q_peak
        K_max = K_abs[i_max]

    half = K_max / 2
    # Search LEFTWARD from Q_peak: |K| = K_max at peak, decreases to <half going left
    hwhm_left = None
    for j in range(i_max, 0, -1):
        if K_abs[j - 1] < half <= K_abs[j]:
            try:
                Q_half = brentq(
                    lambda q: abs(k_b_cluster_pair(q, a, b, X, dDelta, gamma)) - half,
                    Q_grid[j - 1], Q_grid[j], xtol=1e-12,
                )
                hwhm_left = Q_peak_refined - Q_half
                break
            except ValueError:
                continue
    if hwhm_left is None:
        return None
    return {
        "Q_peak": Q_peak_refined,
        "K_max": K_max,
        "HWHM_left": hwhm_left,
        "ratio": hwhm_left / Q_peak_refined,
    }


def main():
    # Sanity-check: bare-doubled-PTF parameters (a=-2, b=-6, X=2, Δδ=0).
    # Today's closed form gives ratio = 0.671535517861.
    print("=" * 72)
    print("Sanity: bare-doubled-PTF (a=-2, b=-6, X=2, Δδ=0)")
    print("=" * 72)
    res = compute_q_peak_hwhm(a=-2, b=-6, X=2, dDelta=0)
    if res:
        print(f"  Q_peak    = {res['Q_peak']:.6f}  (expected ≈ 2.197 in units where Q_EP=1)")
        print(f"  HWHM_left = {res['HWHM_left']:.6f}")
        print(f"  ratio     = {res['ratio']:.6f}  (expected 0.671535)")
    print()

    # Phase 1A: vary Δδ with (a, b, X) fixed at bare-doubled-PTF defaults
    print("=" * 72)
    print("Phase 1A: vary Δδ ∈ {0, 0.5, 1, 2, 4, 8}, fix (a=-2, b=-6, X=2)")
    print("=" * 72)
    print(f"{'Δδ':>8} | {'Q_peak':>10} | {'HWHM_left':>10} | {'ratio':>10}")
    print("-" * 50)
    for dDelta in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
        res = compute_q_peak_hwhm(a=-2, b=-6, X=2, dDelta=dDelta)
        if res:
            print(f"{dDelta:>8.2f} | {res['Q_peak']:>10.4f} | {res['HWHM_left']:>10.4f} | {res['ratio']:>10.6f}")
        else:
            print(f"{dDelta:>8.2f} | (no peak found)")
    print()

    # Phase 1B: vary (a-b) with average fixed at -4 (= bare-doubled-PTF average), X=2, Δδ=0
    print("=" * 72)
    print("Phase 1B: vary (a-b)/2 ∈ {0.5, 1, 2, 4, 6}, fix mean (a+b)/2 = -4, X=2, Δδ=0")
    print("=" * 72)
    print(f"{'(a-b)/2':>10} | {'a':>6} | {'b':>6} | {'Q_peak':>10} | {'HWHM':>8} | {'ratio':>10}")
    print("-" * 65)
    for halfDiff in [0.5, 1.0, 2.0, 4.0, 6.0]:
        a = -4 + halfDiff
        b = -4 - halfDiff
        res = compute_q_peak_hwhm(a=a, b=b, X=2, dDelta=0)
        if res:
            print(f"{halfDiff:>10.2f} | {a:>6.2f} | {b:>6.2f} | {res['Q_peak']:>10.4f} | {res['HWHM_left']:>8.4f} | {res['ratio']:>10.6f}")
        else:
            print(f"{halfDiff:>10.2f} | {a:>6.2f} | {b:>6.2f} | (no peak)")
    print()

    # Phase 1C: vary (a+b) with diff fixed at 4 (= bare), X=2, Δδ=0
    print("=" * 72)
    print("Phase 1C: vary (a+b)/2 ∈ {-1, -2, -4, -8, -16}, fix (a-b)/2 = 2, X=2, Δδ=0")
    print("=" * 72)
    print(f"{'(a+b)/2':>10} | {'a':>6} | {'b':>6} | {'Q_peak':>10} | {'HWHM':>8} | {'ratio':>10}")
    print("-" * 65)
    for halfSum in [-1.0, -2.0, -4.0, -8.0, -16.0]:
        a = halfSum + 2
        b = halfSum - 2
        res = compute_q_peak_hwhm(a=a, b=b, X=2, dDelta=0)
        if res:
            print(f"{halfSum:>10.2f} | {a:>6.2f} | {b:>6.2f} | {res['Q_peak']:>10.4f} | {res['HWHM_left']:>8.4f} | {res['ratio']:>10.6f}")
        else:
            print(f"{halfSum:>10.2f} | {a:>6.2f} | {b:>6.2f} | (no peak)")
    print()

    # Phase 1D: vary X with (a, b, Δδ) fixed
    print("=" * 72)
    print("Phase 1D: vary X ∈ {0.5, 1, 2, 4, 8}, fix (a=-2, b=-6, Δδ=0)")
    print("=" * 72)
    print(f"{'X':>6} | {'Q_peak':>10} | {'HWHM':>8} | {'ratio':>10}")
    print("-" * 45)
    for X in [0.5, 1.0, 2.0, 4.0, 8.0]:
        res = compute_q_peak_hwhm(a=-2, b=-6, X=X, dDelta=0, q_max=20)
        if res:
            print(f"{X:>6.2f} | {res['Q_peak']:>10.4f} | {res['HWHM_left']:>8.4f} | {res['ratio']:>10.6f}")
        else:
            print(f"{X:>6.2f} | (no peak)")
    print()

    # Phase 1E: combined sweep (representative JW cluster-pair parameters at N=5)
    print("=" * 72)
    print("Phase 1E: combined sweep — what does the ratio do at non-trivial Δδ?")
    print("=" * 72)
    print(f"{'a':>6} | {'b':>6} | {'X':>6} | {'Δδ':>6} | {'Q_peak':>8} | {'ratio':>10}")
    print("-" * 60)
    cases = [
        (-2, -6, 2, 0, "bare-doubled-PTF"),
        (-2, -6, 2, 1, "small cluster gap"),
        (-2, -6, 2, 4, "large cluster gap = X"),
        (-3, -5, 2, 2, "tighter (a-b), Δδ=X"),
        (-4, -4, 2, 2, "degenerate a=b, Δδ=X"),
        (-4, -4, 1, 2, "degenerate, X<Δδ"),
        (-2, -6, 4, 2, "X dominant"),
    ]
    for (a, b, X, dDelta, label) in cases:
        res = compute_q_peak_hwhm(a=a, b=b, X=X, dDelta=dDelta, q_max=20)
        if res:
            print(f"{a:>6.2f} | {b:>6.2f} | {X:>6.2f} | {dDelta:>6.2f} | {res['Q_peak']:>8.4f} | {res['ratio']:>10.6f}  ({label})")
        else:
            print(f"{a:>6.2f} | {b:>6.2f} | {X:>6.2f} | {dDelta:>6.2f} | (no peak)  ({label})")


if __name__ == "__main__":
    main()
