"""_takt_absorption_overlay.py - the Absorption Theorem laid over the takt.

Tom, 2026-05-29: "Wir haben die Formel für den Takt gehabt, Absorption Theorem, leg
das drüber." The Takt is γ₀, the carrier. The Absorption Theorem is its formula:
Re(λ) = −2γ₀·popcount(i XOR j). Lay it over the basic Q = J/γ₀ exchange.

The exchange |10> <-> |01> lives in the 2-level single-excitation block (the XX+YY
hop and the Z-dephasing both conserve excitation number, so nothing leaks to |00>,|11>;
the 2-level reduction is EXACT). Under H = J·(XX+YY)/2 the coherent drive is a σ_x
rotation at ω₀ = 2J. The coherence ρ[10,01] sits at Hamming distance
popcount(10 XOR 01) = 2, so the Absorption Theorem damps it at

    Γ = 2·γ₀·popcount = 4·γ₀.

The population difference s_z = ρ_11 − ρ_10 then obeys the damped oscillator

    s_z'' + Γ·s_z' + ω₀²·s_z = 0,   ω₀ = 2J,  Γ = 4γ₀,

with s_z(0) = −1, s_z'(0) = 0, and transfer T(t) = (1 + s_z)/2. Critical damping is
Γ = 2ω₀, i.e. 4γ₀ = 4J, i.e. J = γ₀, Q = 1. So the Absorption Theorem PREDICTS the
swing-death at Q = 1 with no free parameter: the popcount-2 damping rate IS the takt's
threshold. We check the closed form against the full 2-qubit Lindblad bit-exactly.

Run: python simulations/_takt_absorption_overlay.py
"""
from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np
from scipy.linalg import expm

GAMMA_0 = 0.05
TOL = 1e-9
_ok = []

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def popcount(x):
    return bin(x).count("1")


def lindbladian(H, c_list, gamma):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c in c_list:
        cdc = c.conj().T @ c
        L = L + gamma * (np.kron(c, c.conj()) - 0.5 * (np.kron(cdc, Id) + np.kron(Id, cdc.T)))
    return L


def exact_T(J, gamma, ts):
    """Full 2-qubit Lindblad transfer T(t) = P(qubit 1 excited) from |10>."""
    H = J * (np.kron(X, X) + np.kron(Y, Y)) / 2.0
    L = lindbladian(H, [np.kron(Z, I2), np.kron(I2, Z)], gamma)
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[2, 2] = 1.0  # |10><10|
    P1 = np.kron(I2, (I2 - Z) / 2.0)  # qubit 1 excited
    v0 = rho0.reshape(-1, order="F")
    return np.array([float(np.trace(P1 @ (expm(L * t) @ v0).reshape(4, 4, order="F")).real)
                     for t in ts])


def absorption_T(J, gamma, ts):
    """Closed form from the Absorption Theorem: damped oscillator with
    ω₀ = 2J and Γ = 2·γ·popcount(10 XOR 01) = 4γ. No free parameter."""
    w0 = 2.0 * J
    Gamma = 2.0 * gamma * popcount(0b10 ^ 0b01)   # = 4γ, popcount = 2
    out = []
    for t in ts:
        if J > gamma:                              # underdamped, Q > 1
            wd = math.sqrt(w0 * w0 - (Gamma / 2.0) ** 2)
            sz = -math.exp(-Gamma * t / 2.0) * (math.cos(wd * t) + (Gamma / (2.0 * wd)) * math.sin(wd * t))
        elif abs(J - gamma) < 1e-12:               # critical, Q = 1
            sz = -math.exp(-Gamma * t / 2.0) * (1.0 + (Gamma / 2.0) * t)
        else:                                      # overdamped, Q < 1
            r = math.sqrt((Gamma / 2.0) ** 2 - w0 * w0)
            sz = -math.exp(-Gamma * t / 2.0) * (math.cosh(r * t) + (Gamma / (2.0 * r)) * math.sinh(r * t))
        out.append((1.0 + sz) / 2.0)
    return np.array(out)


def main():
    print("=" * 78)
    print(f"THE ABSORPTION THEOREM LAID OVER THE TAKT  (γ₀ = {GAMMA_0})")
    print("=" * 78)
    print("  exchange |10>↔|01>: coherence at Hamming distance popcount(10⊕01) = 2")
    print(f"  Absorption damping Γ = 2γ₀·popcount = 4γ₀ = {4*GAMMA_0};  drive ω₀ = 2J")
    print(f"  critical damping Γ = 2ω₀  =>  4γ₀ = 4J  =>  J = γ₀ = {GAMMA_0}, Q = 1\n")

    ts = np.linspace(0.0, 80.0, 200)
    print(f"  {'J':>7} {'Q':>5} {'regime':>12}  closed-form vs exact Lindblad")
    for J in (0.0125, 0.025, 0.05, 0.10, 0.20, 0.40):
        Texact = exact_T(J, GAMMA_0, ts)
        Tabs = absorption_T(J, GAMMA_0, ts)
        dev = float(np.max(np.abs(Texact - Tabs)))
        Q = J / GAMMA_0
        regime = "overdamped" if J < GAMMA_0 - 1e-12 else ("critical" if abs(J - GAMMA_0) < 1e-12 else "underdamped")
        report(f"J={J:<6} Q={Q:<4.2f} {regime:>12}: Absorption closed form = full Lindblad",
               dev < TOL, f"   max|Δ| = {dev:.2e}")

    # The threshold the Absorption Theorem predicts, read off max T (overshoot past 1/2)
    print("\n  Absorption-predicted max transfer (overshoot past ½ ⟺ underdamped ⟺ Q>1):")
    for J in (0.025, 0.05, 0.10, 0.20):
        Tabs = absorption_T(J, GAMMA_0, np.linspace(0, 200, 400))
        print(f"    J={J:<6} Q={J/GAMMA_0:<4.2f}  max T = {float(np.max(Tabs)):.4f}  "
              f"{'swings (overshoots ½)' if np.max(Tabs) > 0.5 + 1e-6 else 'no overshoot'}")

    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 78)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 78)
    print("""
The Absorption Theorem IS the formula for the takt:
  The exchange coherence sits at popcount 2, so the carrier γ₀ damps it at Γ = 4γ₀.
  Against the coherent drive ω₀ = 2J that is a damped oscillator whose critical point
  Γ = 2ω₀ falls exactly at J = γ₀, Q = 1 — the swing-death the hardware showed. No fit,
  no free parameter: the popcount-2 absorption rate sets the threshold. The beat's
  formula sits exactly on the line where the dance falls into step.
""")


if __name__ == "__main__":
    main()
