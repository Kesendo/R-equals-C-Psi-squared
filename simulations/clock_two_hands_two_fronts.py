#!/usr/bin/env python3
"""The clock's two hands read the two live-open F-fronts.

The two open fronts that sit on this week's clock are NOT the same system, but the clock (a Takt
hand for the radial decay, a Rotation hand for the angular ω) is the shared instrument , each
front rides one hand:

  Rotation hand -> F86b3 (the universal resonance shape). The 2-level EP rotation depends on
    x = Q/Q_EP alone (the proof's eigenvector mixing τ² = (Q−Q_EP)/(Q+Q_EP) = (x−1)/(x+1)), so the
    clock angle θ(x) is IDENTICAL for different g_eff (different Q_EP) , it collapses. That is why
    the resonance shape is universal in Q/Q_EP. (The bare-shape constants x_peak = 2.1969 and the
    HWHM_left/Q_peak floor 0.6715 are derived on C2BareDoubledPtfClosedForm; the +0.08/0.10 lift to
    the empirical 0.756/0.770 rides g_eff, the blocked residue, not the clock.)

  Takt hand -> F87 (the soft/hard break). At γ=0 the Takt is stopped, L = −i[H,·] is symmetric
    about −σ=0, everything is soft. The break of a hard pair (odd hopping cycle) grows FIRST-ORDER
    in the γ-tick (residual/γ -> const); a soft pair stays ~0 at every γ.

This is a seeing, not a solve: the clock unifies the VIEW of the two fronts, not their solutions
(F86b3's lift is g_eff-blocked; F87's converse is the set-level statement we sharpened per-block).
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain


# ---- Rotation hand -> F86b3: the resonance shape rides Q/Q_EP --------------------------------

def slow_mode_theta(Q, g_eff, g0=1.0, k=1):
    """The clock angle θ = arctan(|Im λ|/|Re λ|) of the slowest mode of the F86a 2-level L_eff(k)."""
    J = Q * g0
    L = np.array([[-2 * g0 * (2 * k - 1), 1j * J * g_eff],
                  [1j * J * g_eff,        -2 * g0 * (2 * k + 1)]], dtype=complex)
    ev = np.linalg.eigvals(L)
    decay, omega = -ev.real, np.abs(ev.imag)
    i = int(np.argmin(decay))
    return np.degrees(np.arctan2(omega[i], decay[i]))


def rotation_hand_f86b3():
    print("=" * 70)
    print("ROTATION hand  ->  F86b3: the resonance shape rides Q/Q_EP")
    print("  θ(x) read for two g_eff (two Q_EP); plotted vs x = Q/Q_EP it COLLAPSES.")
    print("=" * 70)
    g1, g2 = 4 / 3, 0.8                                  # Q_EP = 1.5 and 2.5
    print(f"  {'x=Q/Q_EP':>9}  {'θ (g_eff=4/3)':>14}  {'θ (g_eff=0.8)':>14}  {'τ²=(x−1)/(x+1)':>15}")
    for x in [0.5, 0.9, 1.0, 1.2, 1.5, 2.0, 2.197, 3.0]:
        t1 = slow_mode_theta(x * 2.0 / g1, g1)           # Q = x · Q_EP, Q_EP = 2/g_eff
        t2 = slow_mode_theta(x * 2.0 / g2, g2)
        tau2 = (x - 1.0) / (x + 1.0)
        mark = "  <- x_peak (bare-shape peak)" if abs(x - 2.197) < 1e-2 else ("  <- EP" if abs(x - 1.0) < 1e-9 else "")
        print(f"  {x:9.3f}  {t1:13.3f}°  {t2:13.3f}°  {tau2:15.4f}{mark}")
    print("  => θ(x) identical for both g_eff: the EP rotation, hence the shape, is universal in")
    print("     x = Q/Q_EP. Derived bare-shape constants: x_peak = 2.1969, HWHM_left/Q_peak floor")
    print("     = 0.6715 (C2BareDoubledPtfClosedForm); the lift to 0.756/0.770 rides g_eff (blocked).")


# ---- Takt hand -> F87: the soft/hard break rides the γ-tick ----------------------------------

def pair_residual(pair, gamma, N=4):
    """Optimal-transport asymmetry of Spec(L) about −σ (the palindrome-pairing residual)."""
    if gamma == 0.0:
        H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
        L = lindbladian_pauli_dephasing(H, [0.0] * N, dephase_letter='Z')
    else:
        H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
        L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    ev = np.linalg.eigvals(L)
    sigma = N * gamma
    tgt = -ev - 2 * sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def takt_hand_f87():
    print("\n" + "=" * 70)
    print("TAKT hand  ->  F87: the soft/hard break rides the γ-tick")
    print("  soft pair XXZ+ZXX (bipartite) vs hard pair XXZ+XZX (odd cycle), N=4.")
    print("=" * 70)
    soft = [('X', 'X', 'Z'), ('Z', 'X', 'X')]
    hard = [('X', 'X', 'Z'), ('X', 'Z', 'X')]
    print(f"  {'γ':>8}  {'soft residual':>14}  {'hard residual':>14}  {'hard/γ':>10}")
    for gamma in [0.0, 1e-4, 1e-3, 1e-2, 1e-1]:
        rs = pair_residual(soft, gamma)
        rh = pair_residual(hard, gamma)
        hg = f"{rh / gamma:10.4f}" if gamma > 0 else f"{'(—)':>10}"
        print(f"  {gamma:8.0e}  {rs:14.3e}  {rh:14.3e}  {hg}")
    print("  => γ=0: Takt stopped, both soft (residual ~0). hard/γ -> const ≈ 0.256: the hard break")
    print("     is first-order in the γ-tick. The soft pair stays ~0 at every γ.")


def main():
    rotation_hand_f86b3()
    takt_hand_f87()
    print("\n" + "=" * 70)
    print("Two fronts, two hands, one clock. Rotation reads the F86b3 shape (universal in Q/Q_EP);")
    print("Takt reads the F87 break (first-order in γ). The clock unifies the view, not the")
    print("solutions , a seeing.")
    print("=" * 70)


if __name__ == "__main__":
    main()
