"""Verify (or falsify) the path-2 closed form by comparing numerical evaluation
of the sympy-derived expression against the bond-isolate (2)-topology CSV at N=7.

If they match: the H_B-evolved + dephasing-mask approach was exact for this case.
If they don't: we need full Liouvillian eigendecomposition (exact).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV = REPO / "simulations" / "results" / "bond_isolate" / "N7_b0-1_J0.0750_gamma0.0500_probe-coherence.csv"


def closed_form_path2(N: int, J: float, gamma: float, t: np.ndarray) -> np.ndarray:
    """Path-2 topology S(t) from sympy derivation.

    Note: this is the approximate "elementwise dephasing on H_B-evolved" version.
    Verify against bond-isolate full RK4 CSV; if match → exact.
    """
    sqrt2 = np.sqrt(2)
    cos2 = np.cos(2 * sqrt2 * J * t)
    cos4 = np.cos(4 * sqrt2 * J * t)
    cos6 = np.cos(6 * sqrt2 * J * t)
    cos8 = np.cos(8 * sqrt2 * J * t)

    # Block contribution: dispersion 2√2 J base frequency, multiple harmonics
    # Approximate form from sympy expansion (real part of complex exp combination)
    # S_block * exp(4γt)  ≈  block_const + Σ harmonic terms
    # I'll just numerically reconstruct S_block from the formula structure

    # From sympy output (rearranged):
    # Block part contains polynomial of N times exp(±k·sqrt(2)·iJ t) terms / N²√(N-1)
    # Pulling real part:
    #   block * exp(4γt) = (1/(16 N² √(N-1))) · [bare_envelope * exp(0) +
    #                       harmonic_real(N, 4Jt√2)]
    # Rather than recompute symbolically, do:

    bare = (N - 3) * (N - 1) / (N**2)
    # Heuristic: try the exp(-4γt) envelope with harmonic terms scaled by N
    # Re-derive concrete form by substituting into the sympy output
    # For now just compute and compare structure:

    # Apply the exact formula extracted from sympy (real-part of complex exp combo).
    # Pre-factor at exp(-4γt): (16(N-3)(N-1)^{3/2} + harmonic_imag_real_combo(N, J·t)) / (16 N² √(N-1))
    # The harmonic part from sympy was an unhandled complex polynomial; let me eval numerically.

    # CHEATING here: I'll just evaluate the sympy expression directly via lambdify pattern.
    # But simpler: import sympy and eval.
    import sympy as sp
    Ns, Js, gs, ts = sp.symbols("N J gamma t", positive=True, real=True)

    # Reconstruct the sympy formula from the printed output
    # S_total_str_factored = (... long expression in Ns, Js, gs, ts ...)
    # For simplicity, I'll re-import the python script and call its computation.
    # That's a hack; let me just compute S_total directly here.

    sqrt2_s = sp.sqrt(2)
    block_complex = (
        48 * Ns**2 * sp.exp(8 * sqrt2_s * sp.I * Js * ts)
        + 4 * Ns * sp.exp(14 * sqrt2_s * sp.I * Js * ts)
        + 92 * Ns * sp.exp(10 * sqrt2_s * sp.I * Js * ts)
        - 288 * Ns * sp.exp(8 * sqrt2_s * sp.I * Js * ts)
        + 92 * Ns * sp.exp(6 * sqrt2_s * sp.I * Js * ts)
        + 4 * Ns * sp.exp(2 * sqrt2_s * sp.I * Js * ts)
        + sp.exp(16 * sqrt2_s * sp.I * Js * ts)
        - 12 * sp.exp(14 * sqrt2_s * sp.I * Js * ts)
        + 12 * sp.exp(12 * sqrt2_s * sp.I * Js * ts)
        - 276 * sp.exp(10 * sqrt2_s * sp.I * Js * ts)
        + 598 * sp.exp(8 * sqrt2_s * sp.I * Js * ts)
        - 276 * sp.exp(6 * sqrt2_s * sp.I * Js * ts)
        + 12 * sp.exp(4 * sqrt2_s * sp.I * Js * ts)
        - 12 * sp.exp(2 * sqrt2_s * sp.I * Js * ts)
        + 1
    )
    # Bare term in symbolic form (real)
    bare_complex = 16 * (Ns - 3) * (Ns - 1)**(sp.Rational(3, 2)) * sp.exp(8 * sqrt2_s * sp.I * Js * ts)
    # Common factor
    overall_envelope = sp.exp(-4 * gs * ts) * sp.exp(-8 * sqrt2_s * sp.I * Js * ts)
    # Conjugate factor for the block_complex part to make it real
    # The original printed expression has conjugate(1/sqrt(N-1)) — for N real positive, conjugate is just 1/sqrt(N-1)
    norm = 1 / (16 * Ns**2 * sp.sqrt(Ns - 1))

    S_total_sym = norm * overall_envelope * (bare_complex + block_complex / sp.sqrt(Ns - 1))

    # Take real part
    S_real = sp.re(S_total_sym)

    # Evaluate at t values
    S_func = sp.lambdify((Ns, Js, gs, ts), S_real, modules="numpy")
    return S_func(N, J, gamma, t)


def main() -> None:
    if not CSV.exists():
        print(f"CSV not found: {CSV}")
        return

    data = np.loadtxt(CSV, delimiter=",", skiprows=1)
    t = data[:, 0]
    S_csv = data[:, -1]

    print("# Path-2 topology N=7 J=0.075 γ=0.05: closed-form vs bond-isolate CSV")
    print()
    print(f"# bond-isolate CSV: t in [{t[0]}, {t[-1]}], {len(t)} samples")

    S_pred = closed_form_path2(7, 0.075, 0.05, t)

    diff = S_pred - S_csv
    print(f"# max absolute diff: {np.max(np.abs(diff)):.4e}")
    print(f"# mean abs diff:     {np.mean(np.abs(diff)):.4e}")
    print()
    print("# Spot-check at sample times:")
    print("|   t  |  S_csv (bond-isolate RK4)  |  S_pred (sympy approx)  |  diff |")
    for i in [0, 30, 50, 100, 150, 200, 300]:
        print(f"| {t[i]:5.2f} | {S_csv[i]:.6f} | {S_pred[i]:.6f} | {S_pred[i] - S_csv[i]:+.4e} |")


if __name__ == "__main__":
    main()
