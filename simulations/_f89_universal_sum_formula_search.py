"""F89: search for a universal closed form for Sum F_a · N²(N−1) as f(N_block).

Empirical values across paths:
  N_block=4 (path-3): 22/3
  N_block=5 (path-4): 25/2
  N_block=6 (path-5): 483/25
  N_block=7 (path-6): 256/9

Test various candidate forms:
  - polynomial in N_block
  - rational P(N_block) / Q(N_block)
  - involving N_block-related products like N_block(N_block-1)
"""

from __future__ import annotations

import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def main() -> None:
    print("# F89 universal sum-formula search\n")

    # Empirical sums (Sum F_a · N²(N-1) at q=1.5, q-independent)
    data = {
        4: sp.Rational(22, 3),
        5: sp.Rational(25, 2),
        6: sp.Rational(483, 25),
        7: sp.Rational(256, 9),
    }
    n_blocks = sorted(data.keys())
    sums = [data[n] for n in n_blocks]

    print("## Empirical data:")
    print("| N_block | Sum F_a · N²(N-1) | Float |")
    print("|---|---|---|")
    for n, s in zip(n_blocks, sums):
        print(f"| {n} | {s} | {float(s):.6f} |")

    # Try polynomial fits
    print("\n## Polynomial fits in N_block:")
    n_arr = np.array(n_blocks, dtype=float)
    s_arr = np.array([float(s) for s in sums])
    for deg in range(1, 5):
        coefs = np.polynomial.polynomial.polyfit(n_arr, s_arr, deg)
        residual = s_arr - np.polynomial.polynomial.polyval(n_arr, coefs)
        rel_err = np.max(np.abs(residual)) / np.max(np.abs(s_arr))
        print(f"# Degree {deg}: rel_err = {rel_err:.2e}, coefs = {coefs}")
        if rel_err < 1e-9:
            print(f"#   → Polynomial fit found at degree {deg}")

    # Try rational P/Q
    print("\n## Rational candidate forms:")
    candidates = [
        ("(N(N-1))²/((N-2)·something)", None),
        ("N²(N-1)/(N+1)", lambda N: N * N * (N - 1) / (N + 1)),
        ("(N-2)(N-1)N", lambda N: (N - 2) * (N - 1) * N),
        ("N²(N-2)", lambda N: N * N * (N - 2)),
        ("N(N-1)(N-2)", lambda N: N * (N - 1) * (N - 2)),
        ("(N-1)²(N-2)", lambda N: (N - 1) ** 2 * (N - 2)),
        ("(N-1)·something", None),
    ]
    for label, fn in candidates:
        if fn is None:
            continue
        ratios = [float(s) / fn(n) for n, s in zip(n_blocks, sums)]
        print(f"# Sum / [{label}]: {[f'{r:.6f}' for r in ratios]}")
        if all(abs(r - ratios[0]) < 1e-9 for r in ratios):
            print(f"#   → CONSTANT! Sum = {ratios[0]:.10f} · {label}")

    # Compare with N_overlap_pairs (= N_block(N_block-1) = number of overlap basis pairs)
    print("\n## Compare with N_overlap_pairs = N_block(N_block-1):")
    print("| N_block | Sum | N_overlap = N_b(N_b-1) | Sum/N_overlap | nsimplify |")
    print("|---|---|---|---|---|")
    for n, s in zip(n_blocks, sums):
        n_ov = n * (n - 1)
        ratio = s / n_ov
        rat_simp = sp.nsimplify(ratio, rational=True)
        print(f"| {n} | {s} | {n_ov} | {ratio} = {float(ratio):.6f} | {rat_simp} |")

    # Try (Sum) · (denom of sum) / something
    print("\n## Numerator/denominator pattern:")
    print("| N_block | Sum | Numer | Denom |")
    print("|---|---|---|---|")
    for n, s in zip(n_blocks, sums):
        print(f"| {n} | {s} | {s.p} | {s.q} |")

    # Numerators: 22, 25, 483, 256
    # Denominators: 3, 2, 25, 9
    # Try if Sum · (some N_block-dependent quantity) is rational with simpler structure

    # Check if Sum · (N_block+1) gives simpler:
    print("\n## Sum · (N_block + 1):")
    for n, s in zip(n_blocks, sums):
        v = s * (n + 1)
        print(f"# N_block={n}: Sum · {n+1} = {v} = {float(v):.6f}")

    print("\n## Sum · (N_block - 1):")
    for n, s in zip(n_blocks, sums):
        v = s * (n - 1)
        print(f"# N_block={n}: Sum · {n-1} = {v} = {float(v):.6f}")

    # Try guessing via cyclotomic structure: maybe Sum involves φ(N_block+1) somehow
    from sympy import totient
    print("\n## Sum vs φ(N_block+1) and Φ_{N_block+1}(1):")
    print("| N_block | N_block+1 | φ(N_block+1) | Sum | Sum/φ |")
    print("|---|---|---|---|---|")
    for n, s in zip(n_blocks, sums):
        phi = totient(n + 1)
        print(f"| {n} | {n+1} | {phi} | {s} | {s/phi} = {float(s/phi):.6f} |")


if __name__ == "__main__":
    main()
