#!/usr/bin/env python3
"""F50 K_3 N=3 anomaly: testing the class-sum-scalar conjecture (and finding it falsified).

After identifying the K_3 N=3 weight-1 anomaly as living in the 2-dim standard
irrep of S_3 = Aut(K_3), a natural conjecture is:

  CONJECTURE: For graph G with bond set forming a single conjugacy class of
  Aut(G), the dimension excess `dim(ker[H_G, ·]|_{weight-1}) − 2N` equals
  Σ_{ρ ≠ trivial} [mult(ρ, weight-1-as-Aut(G)-module)] × dim(ρ) where the sum
  is over irreps ρ for which the bond-class-sum acts as scalar 0 under
  Schur's lemma (i.e., χ_ρ(t) = 0 on the transposition class).

For K_3 N=3: standard rep (dim 2, χ on transposition = 0) has mult 2 in
weight-1 → predicted excess = 2×2 / 2 = 2. Matches empirical (2 extras).

For K_4 N=4: (2,2) irrep (dim 2, χ on transposition = 0) has mult 4 in
weight-1 (verified via Frobenius reciprocity on each (a, c) orbit's stabilizer).
Predicted excess = 4×2 = 8. **Empirical excess = 0.** Conjecture falsified.

This script documents the test and the falsification: the class-sum scalar
condition is necessary (Schur's lemma) but NOT sufficient for ker-contribution.
The mechanism that selects K_3 N=3 as the unique anomaly remains open.

The likely resolution: the class-sum scalar under Schur acts via
LEFT-multiplication on the group algebra. But `[H, A] = 0` is a MATRIX
commutator condition on the operator space, not a left-multiplication
condition. Even if class-sum · v = 0 (left mult, vanishing on (2,2)-isotypic
sub-rep), this does not imply `[class-sum, A] = 0` for operators A.

The right framework involves the operator-space rep structure under
conjugation by Aut(G), combined with the matrix-commutator condition with
H. This requires deeper rep-theory than the naive class-sum check.

References:
- PROOF_WEIGHT1_DEGENERACY § Appendix (2026-05-17): the K_3 N=3 anomaly + Step-5 derivation gap
- experiments/WEIGHT2_KERNEL.md: the Trivial/Alternating/Mixed framework that today's K_3 finding extends to weight-1
- f50_topology_anomaly_sweep.py: the empirical topology sweep that established the anomaly
"""
from __future__ import annotations

import sys
from itertools import product, permutations

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# S_3 character table (rows: irreps, cols: classes e | (ab) | (abc))
# Class sizes: 1 | 3 | 2
S3_CHARS = {
    "trivial": [1, 1, 1],
    "sign":    [1, -1, 1],
    "standard (2-dim)": [2, 0, -1],
}
S3_CLASS_SIZES = [1, 3, 2]

# S_4 character table (rows: irreps, cols: classes e | (ab) | (ab)(cd) | (abc) | (abcd))
# Class sizes: 1 | 6 | 3 | 8 | 6
S4_CHARS = {
    "trivial":           [1, 1, 1, 1, 1],
    "sign":              [1, -1, 1, 1, -1],
    "(3,1) standard":    [3, 1, -1, 0, -1],
    "(2,2) two-row":     [2, 0, 2, -1, 0],
    "(2,1,1) skew-std":  [3, -1, -1, 0, 1],
}
S4_CLASS_SIZES = [1, 6, 3, 8, 6]


def class_sum_scalar(chi_irrep, dim_irrep, class_idx, class_size):
    """Scalar by which the class sum acts on an irrep (Schur's lemma):

        λ_ρ(class) = (class size) · χ_ρ(class) / dim(ρ)
    """
    return class_size * chi_irrep[class_idx] / dim_irrep


def predict_class_sum_zero_irreps(chars, class_sizes, transposition_idx):
    """For each irrep, compute the class-sum scalar at the transposition class.
    Return list of (name, scalar) pairs, highlighting those with scalar 0."""
    results = []
    for name, chi in chars.items():
        dim = chi[0]  # identity class
        lam = class_sum_scalar(chi, dim, transposition_idx, class_sizes[transposition_idx])
        results.append((name, dim, lam))
    return results


def report():
    print("=" * 70)
    print("F50 K_3 N=3 anomaly: class-sum scalar conjecture test")
    print("=" * 70)
    print()
    print("CONJECTURE: irreps with class-sum scalar = 0 on transpositions contribute to ker[H_KN, ·]|_{w=1}.")
    print()
    print("K_3 = full S_3 transposition class. Class size = 3.")
    print(f"  {'irrep':>20} {'dim':>4} {'class-sum scalar':>20}")
    for name, dim, lam in predict_class_sum_zero_irreps(S3_CHARS, S3_CLASS_SIZES, transposition_idx=1):
        marker = " ← ZERO" if abs(lam) < 1e-9 else ""
        print(f"  {name:>20} {dim:>4} {lam:>20.4f}{marker}")
    print("  Standard (dim 2) annihilated → predicted excess from standard = 2 (mult) × 2 (dim) / 2 = 2")
    print("  Empirical excess at K_3 N=3: 2 ✓ MATCHES")
    print()
    print("K_4 = full S_4 transposition class. Class size = 6.")
    print(f"  {'irrep':>20} {'dim':>4} {'class-sum scalar':>20}")
    for name, dim, lam in predict_class_sum_zero_irreps(S4_CHARS, S4_CLASS_SIZES, transposition_idx=1):
        marker = " ← ZERO" if abs(lam) < 1e-9 else ""
        print(f"  {name:>20} {dim:>4} {lam:>20.4f}{marker}")
    print()
    print("  Compute mult of (2,2) in weight-1 sector via Frobenius reciprocity:")
    print("    - (a, c=0) orbit: stab = S_3 ⊂ S_4, Σ_{h ∈ S_3} χ_(2,2)(h) = 2 + 3·0 + 2·(-1) = 0 → mult = 0")
    print("    - (a, c=1) orbit: stab = S_2 = {e, (k,l)}, Σ = 2 + 0 = 2 → mult = 1")
    print("    - (a, c=2) orbit: stab = S_2, mult = 1")
    print("    - (a, c=3) orbit: stab = S_3, mult = 0")
    print("  Total mult((2,2), weight-1 at N=4) = 2 a-types · (0 + 1 + 1 + 0) = 4")
    print("  Predicted excess from (2,2) = 4 × 2 = 8")
    print("  Empirical excess at K_4 N=4: 0 ✗ CONJECTURE FALSIFIED")
    print()
    print("=" * 70)
    print("Lesson: class-sum scalar = 0 (Schur's lemma on left-multiplication)")
    print("is NOT sufficient for [H, A] = 0 (matrix commutator) on the operator space.")
    print()
    print("The two actions of Aut(G) on the operator space differ:")
    print("  - Conjugation: A → g A g^{-1}.  Decomposes operator space into irrep-isotypic parts.")
    print("  - Left mult:    A → g · A.       Different decomposition.")
    print()
    print("[H, A] = 0 lives in the matrix-commutator framework, not pure left-mult or pure conj.")
    print("The right structural predictor for K_3 N=3 anomaly remains open.")


if __name__ == "__main__":
    report()
