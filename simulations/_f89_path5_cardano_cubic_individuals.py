"""F89 path-5: identify individual F_a amplitudes as Cardano-cubic algebraic forms.

Empirical sigs · N²(N-1) for path-5 F_a modes:
  ω/J = +2.494 (= +4cos(2π/7)): 16.5745
  ω/J = -0.890 (= +4cos(4π/7)): 2.6525
  ω/J = -3.604 (= +4cos(6π/7)): 0.0930

The three frequencies y_n = 4cos(2πn/7) for n=1, 2, 3 are roots of
y³ + 2y² − 8y − 8 = 0 (Cardano-cubic, irreducible over Q).

Hypothesis: each amplitude has form a + b·y_n + c·y_n² for rational
(a, b, c) shared across the three modes. Solve linear system to extract.
"""

from __future__ import annotations

import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def main() -> None:
    print("# F89 path-5 F_a amplitudes as Cardano-cubic forms\n")

    # Empirical (sigs · N²(N-1)) values
    y = [4 * np.cos(2 * np.pi * n / 7) for n in [1, 2, 3]]
    amps = [16.57449911, 2.65249486, 0.09300603]

    print("## Empirical:")
    print("| n | y_n = 4cos(2πn/7) | sigs · N²(N-1) |")
    print("|---|---|---|")
    for n, yi, ai in zip([1, 2, 3], y, amps):
        print(f"| {n} | {yi:+.10f} | {ai:.10f} |")

    # Solve a + b·y + c·y² = amp for each n; 3 unknowns, 3 equations
    A = np.array([[1, yi, yi ** 2] for yi in y])
    b = np.array(amps)
    coefs = np.linalg.solve(A, b)
    print(f"\n## Solving a + b·y + c·y² = amp:")
    print(f"# a = {coefs[0]:.10f}")
    print(f"# b = {coefs[1]:.10f}")
    print(f"# c = {coefs[2]:.10f}")

    # Verify
    for n, yi, ai in zip([1, 2, 3], y, amps):
        pred = coefs[0] + coefs[1] * yi + coefs[2] * yi ** 2
        print(f"# n={n}: predicted = {pred:.10f}, actual = {ai:.10f}, diff = {pred - ai:+.2e}")

    # Recognise as rationals
    print(f"\n## Symbolic recognition (nsimplify):")
    for label, val in [("a", coefs[0]), ("b", coefs[1]), ("c", coefs[2])]:
        rat = sp.nsimplify(val, rational=True, tolerance=1e-9)
        print(f"# {label} = {rat} ≈ {float(rat):.10f}  (numerical {val:.10f})")

    # Verify with exact symbolic y values
    print("\n## Verification with symbolic y_n:")
    a_sym = sp.nsimplify(coefs[0], rational=True, tolerance=1e-9)
    b_sym = sp.nsimplify(coefs[1], rational=True, tolerance=1e-9)
    c_sym = sp.nsimplify(coefs[2], rational=True, tolerance=1e-9)
    for n in [1, 2, 3]:
        y_sym = 4 * sp.cos(2 * sp.pi * n / 7)
        amp_sym = a_sym + b_sym * y_sym + c_sym * y_sym ** 2
        amp_simplified = sp.simplify(amp_sym)
        amp_float = float(amp_simplified)
        actual = amps[n - 1]
        print(f"# n={n}: amp = {amp_simplified} = {amp_float:.10f}  (empirical {actual:.10f}, diff {amp_float - actual:+.2e})")

    print(f"\n## Final closed forms:")
    print(f"# y_n = 4cos(2πn/7) for n=1, 2, 3  (roots of y³ + 2y² − 8y − 8 = 0)")
    print(f"# sigs[F_a:n](N) = (a + b·y_n + c·y_n²) / [N²(N-1)] with:")
    print(f"#   a = {a_sym}")
    print(f"#   b = {b_sym}")
    print(f"#   c = {c_sym}")

    # Sum
    print(f"\n## Sum check (should = 483/25):")
    sum_pred = sum(a_sym + b_sym * y[i] + c_sym * y[i] ** 2 for i in range(3))
    print(f"# Σ amp = {sum_pred:.10f}, expected 483/25 = {483/25:.10f}")
    # Symbolic sum
    sum_sym = 3 * a_sym + b_sym * sum(4 * sp.cos(2 * sp.pi * n / 7) for n in [1, 2, 3]) + c_sym * sum((4 * sp.cos(2 * sp.pi * n / 7)) ** 2 for n in [1, 2, 3])
    sum_simplified = sp.nsimplify(sp.simplify(sum_sym), rational=True)
    print(f"# Symbolic Σ = {sum_simplified}")


if __name__ == "__main__":
    main()
