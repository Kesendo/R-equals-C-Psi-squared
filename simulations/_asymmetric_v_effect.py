#!/usr/bin/env python3
"""First framework-based calculation: asymmetric V-Effect emergent exchange.

Generalize the symmetric result δE_GS = -(3/8) α²/J (EXCHANGE_FROM_V_EFFECT)
to two atoms with different intra-pair couplings J_A and J_B.

Setup: N=4 chain.
  Pair A: qubits {0, 1}, intra-Heisenberg with strength J_A
  Pair B: qubits {2, 3}, intra-Heisenberg with strength J_B
  Bridge: Heisenberg on bond (1, 2) with strength α (V-Effect bridge)

Prediction (second-order PT):
  Both pairs must flip singlet→triplet to be reachable by V.
  Pair A gap: 4 J_A. Pair B gap: 4 J_B. Total gap: 4(J_A + J_B).
  Matrix element norm: ⟨(σ_1·σ_2)²⟩_singlet-singlet = 3 (Pauli identity,
  since ⟨σ_1·σ_2⟩ = 0 between two independent singlets).
  Therefore:  δE_GS^(2) = -3 α² / (4 (J_A + J_B))

Symmetric limit J_A = J_B = J: δE = -3α²/(8J), matches EXCHANGE_FROM_V_EFFECT.

This is the FIRST calculation that builds on framework.py primitives.
"""
import math
import sys

import numpy as np

import framework as fw

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


def build_asymmetric_h(N, J_A, J_B, alpha):
    """N=4 chain with pair-A on (0,1), pair-B on (2,3), bridge on (1,2).

    Each part is a Heisenberg coupling J σ·σ = J(XX + YY + ZZ).
    """
    H = (
        fw._build_bilinear(N, [(0, 1)], [('X', 'X', J_A), ('Y', 'Y', J_A), ('Z', 'Z', J_A)])
        + fw._build_bilinear(N, [(2, 3)], [('X', 'X', J_B), ('Y', 'Y', J_B), ('Z', 'Z', J_B)])
        + fw._build_bilinear(N, [(1, 2)], [('X', 'X', alpha), ('Y', 'Y', alpha), ('Z', 'Z', alpha)])
    )
    return H


def predict_delta_e_gs(J_A, J_B, alpha):
    """Second-order PT prediction: δE_GS = -3α² / (4(J_A + J_B))."""
    return -3.0 * alpha ** 2 / (4.0 * (J_A + J_B))


def main():
    N = 4
    J_pairs = [
        (1.0, 1.0),    # symmetric (sanity, matches -3/(8J))
        (0.5, 1.0),    # mild asymmetry
        (1.0, 2.0),    # 2x asymmetry
        (0.5, 2.0),    # 4x asymmetry
        (1.0, 5.0),    # 5x asymmetry
        (0.3, 3.0),    # 10x asymmetry
        (2.0, 3.0),    # both heavy
    ]
    alphas = [0.025, 0.05, 0.10, 0.20]

    print("=" * 90)
    print("Asymmetric V-Effect emergent exchange (first framework-based calculation)")
    print(f"N = {N}, pair A on (0,1) with J_A, pair B on (2,3) with J_B, bridge on (1,2) with α")
    print("Prediction: δE_GS = -3α² / (4(J_A + J_B))")
    print("=" * 90)
    print(f"\n{'J_A':>5s} {'J_B':>5s} {'α':>6s} {'E_0(0)':>10s} {'E_0(α)':>10s} "
          f"{'δE':>11s} {'δE/α²':>10s} {'predicted':>12s} {'rel err':>10s}")

    summary = []
    for J_A, J_B in J_pairs:
        H_0 = build_asymmetric_h(N, J_A, J_B, 0.0)
        E0_zero = float(np.linalg.eigvalsh(H_0).min())
        # Sanity check: at α=0, GS energy = -3J_A - 3J_B
        expected_E0 = -3 * (J_A + J_B)
        assert abs(E0_zero - expected_E0) < 1e-10, \
            f"GS at α=0 should be -3(J_A+J_B), got {E0_zero} vs {expected_E0}"

        for alpha in alphas:
            H = build_asymmetric_h(N, J_A, J_B, alpha)
            E0 = float(np.linalg.eigvalsh(H).min())
            delta_E = E0 - E0_zero
            scaled = delta_E / alpha ** 2 if alpha > 0 else 0
            predicted = predict_delta_e_gs(J_A, J_B, alpha)
            rel_err = abs(delta_E - predicted) / abs(predicted) if predicted != 0 else 0
            print(f"{J_A:5.2f} {J_B:5.2f} {alpha:6.3f} {E0_zero:10.4f} {E0:10.4f} "
                  f"{delta_E:11.6f} {scaled:10.5f} {predicted:12.6f} {rel_err:10.4%}")
            if alpha == 0.025:
                summary.append((J_A, J_B, scaled, -3.0 / (4 * (J_A + J_B))))
        print()

    print("\nSmall-α limit (α=0.025): δE/α² approaches the asymptotic prediction")
    print(f"\n{'J_A':>5s} {'J_B':>5s} {'δE/α² numeric':>16s} {'-3/(4(J_A+J_B))':>18s} {'rel err':>10s}")
    for J_A, J_B, num, pred in summary:
        rel = abs(num - pred) / abs(pred) if pred != 0 else 0
        print(f"{J_A:5.2f} {J_B:5.2f} {num:16.6f} {pred:18.6f} {rel:10.4%}")

    print()
    print("=" * 90)
    print("Reading: at small α, δE/α² → -3/(4(J_A + J_B)) regardless of asymmetry.")
    print("Symmetric limit J_A=J_B=J recovers -3/(8J) (EXCHANGE_FROM_V_EFFECT result).")
    print("Higher-order corrections grow with α and asymmetry; the leading-order")
    print("prediction is exact in the perturbative regime.")
    print("=" * 90)

    # Verify framework.py's v_effect_emergent_exchange in the symmetric case
    print()
    print("Cross-check: framework.v_effect_emergent_exchange uses the symmetric formula.")
    fw_predict = fw.v_effect_emergent_exchange(0.10, J_intra=1.0)
    here_predict = predict_delta_e_gs(1.0, 1.0, 0.10)
    print(f"  framework symmetric formula at α=0.10, J=1: {fw_predict:.6f}")
    print(f"  this script's symmetric prediction:        {here_predict:.6f}")
    assert abs(fw_predict - here_predict) < 1e-12, "Framework and this script should agree at symmetric point"
    print("  ✓ agreement at symmetric point.")


if __name__ == "__main__":
    main()
