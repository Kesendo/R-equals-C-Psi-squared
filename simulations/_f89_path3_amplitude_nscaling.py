"""F89 path-3: derive closed-form N-scaling for the populated mode amplitudes.

Path-3 has 10 distinct populated mode-groups at q = J/γ = 1.5 (per
`_f89_path3_multi_exp_derive.py`). Goal: for each mode-group, find a
closed-form rational A(N) = polynomial(N) / [N²(N−1)] using numerical
amplitudes at multiple N values, analog of path-2's
A_Bloch(N) = 3·(N−3)²/(2·N²(N−1)).

The 4 AT-locked modes (rates 2γ and 6γ, freqs J·(−1±√5)) are q-INDEPENDENT
in rate, so their amplitudes should admit clean rational forms in N.
The 8 octic-derived modes have q-dependent rates; at fixed q their
amplitudes vs N may still be rational, but rates shift with q so any
closed form in (N, q) requires the cubic-octic substitution.

Per F89, the partial-trace prefactor is 1/(2√(N·C(N,2))) = 1/√(2·N²(N−1)),
so squared per-coherence amplitudes scale as 1/(N²(N−1)). Path-3's N_E
factor (where N_E = N − 4 for path-3) enters from the term-2 of
ρ_block(0) carrying N_E weight.
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
N_BLOCK = 4


def numerical_amplitude_per_mode(N: int, J: float, gamma: float):
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
    eigvals, sigs = numerical_amplitude_per_mode(N, J, gamma)
    rates = -eigvals.real / gamma
    freqs = np.abs(eigvals.imag) / J
    contributing = np.where(sigs > threshold)[0]
    groups: dict[tuple, float] = {}
    for k in contributing:
        key = (round(rates[k], 4), round(freqs[k], 4))
        groups[key] = groups.get(key, 0.0) + sigs[k]
    return groups


def fit_polynomial_in_n(ns: list[int], values: list[float], divisor_fn) -> tuple:
    """Fit values[i] · divisor_fn(N) = polynomial(N) of integer degree.

    Returns (coefs_low_to_high, degree, rel_err) of best-fitting polynomial
    up to degree 4, or (None, None, None) if no clean fit.
    """
    n_arr = np.array(ns, dtype=float)
    v_arr = np.array(values, dtype=float)
    target = v_arr * np.array([divisor_fn(N) for N in ns], dtype=float)
    for degree in range(0, 5):
        coefs = np.polynomial.polynomial.polyfit(n_arr, target, degree)
        residual = target - np.polynomial.polynomial.polyval(n_arr, coefs)
        rel_err = np.max(np.abs(residual)) / max(1e-15, np.max(np.abs(target)))
        if rel_err < 1e-9:
            return coefs, degree, rel_err
    return None, None, None


def main() -> None:
    J, gamma = 0.075, 0.05
    print(f"# F89 path-3 amplitude N-scaling, J={J}, γ={gamma} (q=1.5)\n")

    # Path-3 needs N ≥ 5 (4-qubit block + ≥1 bare site for N_E ≥ 1)
    Ns = [5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20]
    print(f"## Computing amplitudes at N ∈ {Ns}")

    all_groups: dict[tuple, list[tuple[int, float]]] = {}
    for N in Ns:
        groups = grouped_amplitudes(N, J, gamma)
        for key, val in groups.items():
            all_groups.setdefault(key, []).append((N, val))

    consistent = {k: v for k, v in all_groups.items() if len(v) == len(Ns)}
    print(f"# {len(consistent)} mode-groups consistent across all N values:\n")

    print("| (rate Γ/γ, |freq|/J) | N=5 | N=7 | N=11 | N=20 |")
    print("|---|---|---|---|---|")
    for key in sorted(consistent.keys()):
        d = dict(consistent[key])
        print(f"| ({key[0]}, {key[1]}) | {d.get(5, 0):.4e} | {d.get(7, 0):.4e} | "
              f"{d.get(11, 0):.4e} | {d.get(20, 0):.4e} |")
    print()

    # Try multiple divisor patterns analogous to path-2's findings
    divisor_patterns = [
        ("N²(N−1)", lambda N: N * N * (N - 1)),
        ("N²(N−1)/(N−4)²", lambda N: N * N * (N - 1) / max((N - 4) ** 2, 1)),
        ("N²(N−1)/(N−4)", lambda N: N * N * (N - 1) / max(N - 4, 1)),
        ("N²(N−1)·N²", lambda N: N * N * (N - 1) * N * N),
    ]

    for label, fn in divisor_patterns:
        print(f"\n## Fit attempt: amplitude(N) · [{label}] = polynomial(N)")
        print("| Mode | Polynomial fit (low→high degree) | Degree | rel err |")
        print("|---|---|---|---|")
        for key in sorted(consistent.keys()):
            ns = [n for n, _ in consistent[key]]
            vals = [v for _, v in consistent[key]]
            coefs, deg, err = fit_polynomial_in_n(ns, vals, fn)
            if coefs is not None:
                rounded = [sp.nsimplify(c, rational=True, tolerance=1e-6) for c in coefs]
                poly_str = " + ".join(
                    f"({c})·N^{i}" for i, c in enumerate(rounded) if c != 0
                ) or "0"
                print(f"| ({key[0]}, {key[1]}) | {poly_str} | {deg} | {err:.2e} |")
            else:
                print(f"| ({key[0]}, {key[1]}) | (no fit ≤ degree 4) | – | – |")

    # Hamming-complement pair structure: pairs (Γ, 8γ−Γ) at fixed |freq|
    print("\n## Hamming-complement pair structure (Γ_a + Γ_b = 8γ at fixed |freq|/J)")
    pairs = [
        ((3.3488, 1.206), (4.6512, 1.206)),
        ((3.5989, 2.93), (4.4011, 2.93)),
        ((3.777, 5.178), (4.223, 5.178)),
        ((4.0, 0.5944), (4.0, 7.5024)),  # both at rate 4γ but different freq
    ]
    print("| Pair (Γ_a, |ω|/J) ↔ (Γ_b, |ω|/J) | sum Γ/γ | A_a (N=11) | A_b (N=11) | A_a/A_b |")
    print("|---|---|---|---|---|")
    for ka, kb in pairs:
        a = dict(consistent.get(ka, [])).get(11, 0.0)
        b = dict(consistent.get(kb, [])).get(11, 0.0)
        ratio = a / b if b > 1e-15 else float("inf")
        print(f"| ({ka[0]}, {ka[1]}) ↔ ({kb[0]}, {kb[1]}) | {ka[0] + kb[0]:.4f} | "
              f"{a:.3e} | {b:.3e} | {ratio:.2f} |")


if __name__ == "__main__":
    main()
