"""F89 path-2: derive closed-form N-scaling for the 4 populated mode amplitudes.

Path-2 has 4 distinct populated mode-groups at J/γ=1.5:
  (rate 2γ, freq 0)
  (rate 2γ, freq ±2√2 J)
  (rate 3.04γ, freq 0)
  (rate 3.48γ, freq ±5.45 J)

Each mode-group has a per-N amplitude. Numerically we have N=5, 7, 11 data;
the goal is to fit closed-form rational expressions in N for each amplitude,
analog of the all-isolated formula's (N-1)/N + 4m(N-2)·(...)/(N²(N-1))
combinatorial prefactors.

Approach:
  1. Compute amplitudes numerically at N ∈ {5, 6, 7, 8, 9, 10, 11, 12, 15, 20}.
  2. For each mode-group, fit ratio amplitude / (1/(N²(N-1))) and check if it
     reduces to a polynomial in N (likely candidates: 1, (N-2), (N-2)², N_E²).
  3. Verify fitted closed form against new N values.

Anchor: per F89 partial-trace prefactor is 1/(2√(N·C(N,2))) = 1/√(2·N²(N-1)),
so squared per-coherence amplitude scales as 1/(N²(N-1)). N_E = N - 3 for
path-2 enters as a prefactor (term-2 of ρ_block(0) carries N_E factor).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

from _f89_pathk_lib import (
    build_block_L,
    compute_rho_block_0,
    per_site_reduction_matrix,
)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
N_BLOCK = 3


def numerical_amplitude_per_mode(N: int, J: float, gamma: float):
    """Return (eigvals, sigs) where sigs[k] = Σ_l |a[l, k]|² is the total
    population of L_super eigenmode k, summed over block sites."""
    L = build_block_L(J, gamma, N_BLOCK)
    rho = compute_rho_block_0(N_BLOCK, N)
    vec = rho.flatten(order="F")
    eigvals, R = np.linalg.eig(L)
    c = np.linalg.solve(R, vec)
    w = per_site_reduction_matrix(N_BLOCK)
    M = w @ R
    a = M * c[None, :]
    sigs = np.sum(np.abs(a) ** 2, axis=0)
    return eigvals, sigs


def grouped_amplitudes(N: int, J: float, gamma: float, threshold: float = 1e-12):
    """Group amplitudes by (rate Γ/γ, |freq|/J)."""
    eigvals, sigs = numerical_amplitude_per_mode(N, J, gamma)
    rates = -eigvals.real / gamma
    freqs = np.abs(eigvals.imag) / J
    contributing = np.where(sigs > threshold)[0]
    groups = {}
    for k in contributing:
        key = (round(rates[k], 4), round(freqs[k], 4))
        groups[key] = groups.get(key, 0.0) + sigs[k]
    return groups


def fit_rational_in_n(ns: list[int], values: list[float], j_gamma_factor: float = 1.0):
    """Try to fit values[i] = polynomial(N=ns[i]) / (N²(N-1)) for some integer poly.

    Returns the polynomial coefficients (low to high degree) if the fit is
    clean rational with integer/rational numerators, else None.

    j_gamma_factor: optional rational scale to factor out (e.g. 1/2, √2 etc.).
    """
    # Multiply values by N²(N-1)/j_gamma_factor → should be a polynomial in N
    n_arr = np.array(ns, dtype=float)
    v_arr = np.array(values, dtype=float)
    poly_target = v_arr * (n_arr * n_arr * (n_arr - 1)) / j_gamma_factor

    # Try increasing polynomial degree in N
    for degree in range(0, 5):
        coefs = np.polynomial.polynomial.polyfit(n_arr, poly_target, degree)
        residual = poly_target - np.polynomial.polynomial.polyval(n_arr, coefs)
        rel_err = np.max(np.abs(residual)) / max(1e-15, np.max(np.abs(poly_target)))
        if rel_err < 1e-9:
            return coefs, degree, rel_err
    return None, None, None


def main() -> None:
    J, gamma = 0.075, 0.05
    print(f"# F89 path-2 amplitude N-scaling, J={J}, γ={gamma}\n")

    Ns = [5, 6, 7, 8, 9, 10, 11, 12, 15, 20]
    print(f"## Computing amplitudes at N ∈ {Ns}")

    # Collect per-N amplitudes for each mode-group
    all_groups: dict[tuple, list[float]] = {}
    for N in Ns:
        groups = grouped_amplitudes(N, J, gamma)
        for key, val in groups.items():
            all_groups.setdefault(key, []).append((N, val))

    # Filter: only keep groups present at all N values
    consistent = {k: vals for k, vals in all_groups.items() if len(vals) == len(Ns)}
    print(f"# {len(consistent)} mode-groups consistent across all N values:\n")

    print("| (rate Γ/γ, |freq|/J) | N=5 | N=7 | N=11 | N=20 |")
    print("|---|---|---|---|---|")
    for key in sorted(consistent.keys()):
        vals_dict = dict(consistent[key])
        v5 = vals_dict.get(5, 0.0)
        v7 = vals_dict.get(7, 0.0)
        v11 = vals_dict.get(11, 0.0)
        v20 = vals_dict.get(20, 0.0)
        print(f"| ({key[0]}, {key[1]}) | {v5:.4e} | {v7:.4e} | {v11:.4e} | {v20:.4e} |")
    print()

    # Try fitting each group as rational in N: amplitude · N²(N-1) = polynomial(N)
    print("## Attempting closed-form fit: amplitude(N) · N²·(N−1) = polynomial(N)?\n")
    print("| Mode | Polynomial fit (low→high degree) | Degree | rel err |")
    print("|---|---|---|---|")
    for key in sorted(consistent.keys()):
        ns = [n for n, _ in consistent[key]]
        vals = [v for _, v in consistent[key]]
        coefs, deg, err = fit_rational_in_n(ns, vals)
        if coefs is not None:
            # Try to round to rationals
            rounded = [sp.nsimplify(c, rational=True, tolerance=1e-6) for c in coefs]
            poly_str = " + ".join(f"({c})·N^{i}" for i, c in enumerate(rounded) if c != 0)
            print(f"| ({key[0]}, {key[1]}) | {poly_str} | {deg} | {err:.2e} |")
        else:
            print(f"| ({key[0]}, {key[1]}) | (no fit ≤ degree 4) | – | – |")
    print()

    # Investigate: maybe scale by N_E = N-3 helps for some modes
    print("## Alternative fit: amplitude(N) · N²(N−1) / (N−3)² = polynomial(N)?")
    print("# (rationale: dominant (vac,SE) mode scales as N_E²/(N²(N-1)) per Parseval analysis)\n")
    print("| Mode | Polynomial fit | Degree | rel err |")
    print("|---|---|---|---|")
    for key in sorted(consistent.keys()):
        ns = [n for n, _ in consistent[key]]
        vals = [v for _, v in consistent[key]]
        n_arr = np.array(ns, dtype=float)
        v_arr = np.array(vals, dtype=float)
        # Divide by (N-3)²
        adjusted = v_arr * (n_arr * n_arr * (n_arr - 1)) / ((n_arr - 3) ** 2)
        for degree in range(0, 4):
            coefs = np.polynomial.polynomial.polyfit(n_arr, adjusted, degree)
            residual = adjusted - np.polynomial.polynomial.polyval(n_arr, coefs)
            rel_err = np.max(np.abs(residual)) / max(1e-15, np.max(np.abs(adjusted)))
            if rel_err < 1e-9:
                rounded = [sp.nsimplify(c, rational=True, tolerance=1e-6) for c in coefs]
                poly_str = " + ".join(f"({c})·N^{i}" for i, c in enumerate(rounded) if c != 0)
                print(f"| ({key[0]}, {key[1]}) | {poly_str} | {degree} | {rel_err:.2e} |")
                break
        else:
            print(f"| ({key[0]}, {key[1]}) | (no fit ≤ degree 3) | – | – |")
    print()

    # Extract q-dependence: redo at multiple q values, see if amplitude / N-scaling
    # gives a clean function of q = J/γ
    print("## q-dependence scan: amplitude(N=11, q) · N²(N-1) / N_E^p")
    print("# (N=11 fixed; varying q to extract J/γ scaling of the per-mode coefficient)\n")
    qs = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
    N_fixed = 11
    print("| Mode key (rate Γ/γ, |freq|/J at q=1.5 — note rates shift with q) | " +
          " | ".join(f"q={q}" for q in qs) + " |")
    print("|---|" + "|".join(["---"] * len(qs)) + "|")
    for q in qs:
        gamma_q = 0.05
        J_q = q * gamma_q
        groups_q = grouped_amplitudes(N_fixed, J_q, gamma_q)
        # Show amplitude · N²(N-1) for each q
    # Print per-q mode survey separately (modes at different q have different rates/freqs)
    for q in qs:
        gamma_q = 0.05
        J_q = q * gamma_q
        groups_q = grouped_amplitudes(N_fixed, J_q, gamma_q)
        print(f"\n### q = {q} (J = {J_q}, γ = {gamma_q}, N = {N_fixed})")
        print("| (rate Γ/γ, |freq|/J) | amplitude · N²(N-1) | / (N-3)² |")
        print("|---|---|---|")
        for key in sorted(groups_q.keys()):
            amp = groups_q[key]
            scaled = amp * N_fixed * N_fixed * (N_fixed - 1)
            scaled_NE2 = scaled / ((N_fixed - 3) ** 2)
            print(f"| ({key[0]}, {key[1]}) | {scaled:.6f} | {scaled_NE2:.6f} |")


if __name__ == "__main__":
    main()
