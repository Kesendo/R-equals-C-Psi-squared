#!/usr/bin/env python3
"""F86a's exceptional point, read through the clock (Takt / Rotation).

F86a (PROOF_F86A_EP_MECHANISM): inside each coherence block the two slowest adjacent rate
channels form a same-sign-imaginary 2x2 effective Liouvillian

    L_eff(k) = [[ -2γ₀(2k-1),  iJ·g_eff ],
                [  iJ·g_eff,   -2γ₀(2k+1) ]]

with eigenvalues λ_±(k) = -4γ₀k ± sqrt(4γ₀² - J²·g_eff²). The discriminant vanishes at the
exceptional point J·g_eff = 2γ₀, i.e. Q_EP = 2/g_eff.

Read that 2-level through the clock we built on MirrorSystem: a Takt hand (the radial decay
-Re λ) and a Rotation hand (the angular |Im λ|, angle θ = arctan(|Im|/|Re|)). The point: the EP
is exactly where the Rotation hand lifts off the Takt axis. Below Q_EP the pair is real (θ=0,
pure Takt, overdamped, the slow mode is the longer-lived one); at Q_EP it coalesces; above Q_EP
the decay is PINNED at -4γ₀k (the Takt holds) and only the Rotation opens (θ grows). Reading
F86a's Q_EP at the real g_eff values puts the F86 peak Q's (1.5, 2.5) on the clock as
Rotation-onset points.
"""
import numpy as np

GAMMA0 = 1.0
K = 1  # the slowest channel pair (sets t_peak = 1/(4γ₀))


def slow_mode_clock(Q, g_eff, k=K, g0=GAMMA0):
    """Read the slowest mode of L_eff(k) through the clock: (decay, omega, theta_deg)."""
    J = Q * g0
    L = np.array([[-2 * g0 * (2 * k - 1), 1j * J * g_eff],
                  [1j * J * g_eff,        -2 * g0 * (2 * k + 1)]], dtype=complex)
    ev = np.linalg.eigvals(L)
    decay = -ev.real            # Takt: radial decay rate
    omega = np.abs(ev.imag)     # Rotation: angular frequency
    i = int(np.argmin(decay))   # the slowest (longest-lived) mode, as the Rotation voice reads
    theta = np.degrees(np.arctan2(omega[i], decay[i]))
    return decay[i], omega[i], theta


def main():
    for g_eff in [4 / 3, 0.8]:                  # Q_EP = 1.5 (c=2 peak) and 2.5 (Endpoint orbit)
        q_ep = 2.0 / g_eff
        print(f"\ng_eff = {g_eff:.4f}  ->  Q_EP = 2/g_eff = {q_ep:.4f}   (k={K}, γ₀={GAMMA0})")
        print(f"  {'Q':>7}  {'decay(-Re)':>10}  {'ω(|Im|)':>8}  {'θ deg':>6}")
        for Q in [0.5 * q_ep, 0.9 * q_ep, q_ep, 1.1 * q_ep, 1.5 * q_ep, 2.0 * q_ep]:
            decay, omega, theta = slow_mode_clock(Q, g_eff)
            if abs(Q - q_ep) < 1e-9:
                mark = "  <- EP (modes coalesce, θ still 0)"
            elif Q > q_ep:
                mark = "  <- Rotation on (Takt pinned at 4γ₀)"
            else:
                mark = ""
            print(f"  {Q:7.3f}  {decay:10.3f}  {omega:8.3f}  {theta:6.1f}{mark}")


if __name__ == "__main__":
    main()
