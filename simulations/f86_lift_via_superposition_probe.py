"""Phase 4 of the F86 Γ_pair derivation attempt:

Phase 1 found that single-cluster-pair K_b curves give HWHM_left/Q_peak ≈ 0.671535
universally (across (a+b)/2, X, mild Δδ dependence). The empirical Interior 0.7506
and Endpoint 0.7728 must therefore come from SUPERPOSITION over multiple cluster
pairs — each contributing a Lorentzian-like K_b with its own Q_EP_pair.

This probe constructs simple multi-pair sums and measures the resulting HWHM_left/
Q_peak ratio, checking whether realistic JW-cluster-pair distributions can lift
from 0.6715 to 0.7506/0.7728.

Model: K_b_total(Q) = Σ_p w_p · K_b_bare(Q/Q_EP_p) using today's bare closed form
where each pair has its own Q_EP_p and weight w_p (Frobenius² of the cross-block).
"""

from __future__ import annotations

import sys
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def k_b_bare(x: float) -> float:
    """Bare-doubled-PTF K_b at dimensionless x = Q/Q_EP. Both regimes, stable across EP.
    Reproduces C2BareDoubledPtfClosedForm.EvaluateKb (γ₀=1, normalised scale)."""
    eMinus2 = math.exp(-2.0)
    small = 0.05
    if x > 1.0:
        xi = math.sqrt(x ** 2 - 1)
        if xi < small:
            xi2 = xi ** 2
            bracket = -5 * xi2 ** 2 / 12 + 11 * xi2 ** 3 / 360
        else:
            bracket = (xi ** 2 + 2) * math.cos(xi) - 2
        return eMinus2 * x * bracket / xi ** 4
    if x < 1.0:
        mu = math.sqrt(1 - x ** 2)
        if mu < small:
            mu2 = mu ** 2
            bracket = -5 * mu2 ** 2 / 12 - 11 * mu2 ** 3 / 360
        else:
            bracket = (2 - mu ** 2) * math.cosh(mu) - 2
        return eMinus2 * x * bracket / mu ** 4
    return -5 * eMinus2 / 12  # EP


def k_b_total(Q: float, pair_specs):
    """K_b_total(Q) = Σ_p w_p · K_b_bare(Q/Q_EP_p). pair_specs = [(Q_EP_p, w_p), ...]"""
    total = 0.0
    for Q_EP, w in pair_specs:
        total += w * k_b_bare(Q / Q_EP)
    return total


def scan_hwhm_ratio(pair_specs, q_max=8.0, n_grid=4000):
    Q_grid = np.linspace(0.05, q_max, n_grid)
    K_vals = np.array([k_b_total(Q, pair_specs) for Q in Q_grid])
    K_abs = np.abs(K_vals)
    i_max = int(np.argmax(K_abs))
    if i_max == 0 or i_max == n_grid - 1:
        return None

    # Parabolic refine
    y0, y1, y2 = K_abs[i_max - 1], K_abs[i_max], K_abs[i_max + 1]
    denom = y0 - 2 * y1 + y2
    delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-15 else 0.0
    dQ = Q_grid[1] - Q_grid[0]
    Q_peak = Q_grid[i_max] + delta * dQ

    # Brent refine
    try:
        res = minimize_scalar(
            lambda q: -abs(k_b_total(q, pair_specs)),
            bracket=(Q_grid[max(0, i_max - 3)], Q_peak, Q_grid[min(n_grid - 1, i_max + 3)]),
            method="brent", options={"xtol": 1e-10},
        )
        Q_peak = res.x
        K_max = -res.fun
    except Exception:
        K_max = K_abs[i_max]

    half = K_max / 2
    hwhm_left = None
    for j in range(i_max, 0, -1):
        if K_abs[j - 1] < half <= K_abs[j]:
            try:
                Q_half = brentq(
                    lambda q: abs(k_b_total(q, pair_specs)) - half,
                    Q_grid[j - 1], Q_grid[j], xtol=1e-12,
                )
                hwhm_left = Q_peak - Q_half
                break
            except Exception:
                continue
    if hwhm_left is None:
        return None
    return {"Q_peak": Q_peak, "K_max": K_max, "HWHM_left": hwhm_left, "ratio": hwhm_left / Q_peak}


def main():
    print("=" * 72)
    print("Phase 4 superposition probe — can multi-pair sums lift HWHM ratio above 0.6715?")
    print("=" * 72)
    print()

    # Sanity: single pair, should give 0.671535
    print("Sanity: single pair at Q_EP=1, w=1")
    res = scan_hwhm_ratio([(1.0, 1.0)])
    print(f"  Q_peak = {res['Q_peak']:.4f}, HWHM = {res['HWHM_left']:.4f}, ratio = {res['ratio']:.6f}")
    print()

    # Two pairs at different Q_EP
    print("Two pairs at different Q_EP, equal weight w=0.5 each:")
    print(f"{'Q_EP_1':>8} {'Q_EP_2':>8} | {'Q_peak':>8} {'HWHM':>8} {'ratio':>10}")
    print("-" * 50)
    for Q1, Q2 in [(1.0, 1.0), (1.0, 1.2), (1.0, 1.5), (1.0, 2.0), (1.0, 3.0), (0.5, 1.5), (0.5, 2.0)]:
        res = scan_hwhm_ratio([(Q1, 0.5), (Q2, 0.5)])
        if res:
            print(f"{Q1:>8.2f} {Q2:>8.2f} | {res['Q_peak']:>8.4f} {res['HWHM_left']:>8.4f} {res['ratio']:>10.6f}")
    print()

    # Three pairs spread
    print("Three pairs spread across Q_EP, equal weight w=1/3:")
    print(f"{'Q_EPs':>20} | {'Q_peak':>8} {'HWHM':>8} {'ratio':>10}")
    print("-" * 55)
    for triplet in [[1.0, 1.0, 1.0], [0.8, 1.0, 1.2], [0.5, 1.0, 1.5], [0.5, 1.0, 2.0], [1.0, 1.5, 2.0], [0.5, 1.5, 3.0]]:
        specs = [(q, 1.0/3) for q in triplet]
        res = scan_hwhm_ratio(specs, q_max=10)
        if res:
            print(f"{str(triplet):>20} | {res['Q_peak']:>8.4f} {res['HWHM_left']:>8.4f} {res['ratio']:>10.6f}")
    print()

    # Unequal weights (one dominant, others as small contributions)
    print("Dominant pair + small contribution (Endpoint-like):")
    print(f"{'Q_EP_dom':>10} {'Q_EP_aux':>10} {'w_aux':>8} | {'ratio':>10}")
    print("-" * 50)
    for w_aux in [0.0, 0.05, 0.1, 0.2, 0.4, 0.6]:
        specs = [(1.0, 1.0 - w_aux), (0.5, w_aux)]
        res = scan_hwhm_ratio(specs)
        if res:
            print(f"{1.0:>10.2f} {0.5:>10.2f} {w_aux:>8.2f} | {res['ratio']:>10.6f}")
        else:
            print(f"{1.0:>10.2f} {0.5:>10.2f} {w_aux:>8.2f} | (no peak)")
    print()

    # Search: target ratio 0.7506 (Interior) or 0.7728 (Endpoint)
    # by tuning a 2-pair configuration
    print("Search for ratio = 0.7506 (Interior target) via 2-pair sum w_1 + w_2 = 1:")
    targets = [0.7506, 0.7728]
    for target in targets:
        best = None
        for Q1 in np.linspace(0.3, 3.0, 20):
            for Q2 in np.linspace(0.3, 3.0, 20):
                if abs(Q1 - Q2) < 0.01:
                    continue
                for w1 in np.linspace(0.1, 0.9, 9):
                    w2 = 1.0 - w1
                    specs = [(Q1, w1), (Q2, w2)]
                    res = scan_hwhm_ratio(specs, q_max=8)
                    if res:
                        err = abs(res['ratio'] - target)
                        if best is None or err < best[0]:
                            best = (err, Q1, Q2, w1, w2, res['ratio'])
        if best:
            err, Q1, Q2, w1, w2, ratio = best
            print(f"  target {target}: best ratio = {ratio:.6f} at Q1={Q1:.2f}, Q2={Q2:.2f}, w1={w1:.2f}, w2={w2:.2f}, err={err:.4f}")
        else:
            print(f"  target {target}: no match found")


if __name__ == "__main__":
    main()
