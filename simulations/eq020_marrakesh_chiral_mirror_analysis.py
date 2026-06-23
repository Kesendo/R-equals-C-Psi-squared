#!/usr/bin/env python3
"""Re-analyze existing Marrakesh K-partnership data for the chiral-mirror Bloch
prediction.

Existing data: data/ibm_k_partnership_april2026/ — 5 bonding-mode receivers
(k=1..5) at N=5 on Marrakesh path [48, 49, 50, 51, 58], t = 0.8, 3 Trotter
steps, XX-only. Each receiver has tomographed ρ_AB (4×4) on (qubit 0, qubit 4).

Today's chiral mirror prediction:
  ⟨X_a⟩(ψ_{N+1-k}) = +(−1)^a · ⟨X_a⟩(ψ_k)
  ⟨Y_a⟩(ψ_{N+1-k}) = −(−1)^a · ⟨Y_a⟩(ψ_k)
  ⟨Z_a⟩(ψ_{N+1-k}) = +⟨Z_a⟩(ψ_k)
  P_a (purity) identical

For sites a ∈ {0, 4} (the pair endpoints in this data), N=5:
  K_1 sign at site 0: (−1)^0 = +1
  K_1 sign at site 4: (−1)^4 = +1
So at the K_1-paired states (k=1, k=5) and (k=2, k=4):
  ⟨X_0⟩, ⟨X_4⟩ should match (sign +1)
  ⟨Y_0⟩, ⟨Y_4⟩ should sign-flip (−1 from energy reversal)
  ⟨Z_0⟩, ⟨Z_4⟩ should match
"""
from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

DATA_DIR = Path(__file__).parent.parent / "data" / "ibm_k_partnership_april2026"

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def load_rho_ab(receiver_data):
    """Reconstruct ρ_AB (4×4) from the JSON real/imag parts."""
    re = np.array(receiver_data["rho_2q_real"])
    im = np.array(receiver_data["rho_2q_imag"])
    return re + 1j * im


def reduced_singles(rho_ab):
    """Return (ρ_A, ρ_B) for the 2-qubit ρ_AB.

    Convention: ρ_AB on (A, B) = (qubit 0 of pair, qubit 1 of pair).
    """
    rho_4 = rho_ab.reshape(2, 2, 2, 2)
    rho_A = np.trace(rho_4, axis1=1, axis2=3)
    rho_B = np.trace(rho_4, axis1=0, axis2=2)
    return rho_A, rho_B


def bloch(rho_1q):
    """⟨X⟩, ⟨Y⟩, ⟨Z⟩ from 2×2 ρ."""
    return (
        float(np.real(np.trace(rho_1q @ X))),
        float(np.real(np.trace(rho_1q @ Y))),
        float(np.real(np.trace(rho_1q @ Z))),
    )


def purity(rho_1q):
    return float(np.real(np.trace(rho_1q @ rho_1q)))


def analyze(label, fname):
    print(f"\n{'=' * 80}")
    print(f"{label}: {fname.name}")
    print(f"{'=' * 80}")
    with open(fname) as f:
        d = json.load(f)
    print(f"  Backend: {d['backend']}, path: {d['path']}, N={d['n']}, "
          f"t={d['t_evol']}, xx_only={d['xx_only']}")
    print(f"  Shots: {d['shots']}")
    print()

    # Extract Bloch components per receiver, for sites 0 (endpoint A) and 4 (endpoint B)
    print(f"  Per-receiver single-qubit Bloch components (sites 0 and 4):")
    print(f"  {'k':<10} {'⟨X_0⟩':>9} {'⟨Y_0⟩':>9} {'⟨Z_0⟩':>9} {'P_0':>8}   "
          f"{'⟨X_4⟩':>9} {'⟨Y_4⟩':>9} {'⟨Z_4⟩':>9} {'P_4':>8}")
    blochs = {}
    for k in [1, 2, 3, 4, 5]:
        rho_ab = load_rho_ab(d["receivers"][f"bonding:{k}"])
        rho_A, rho_B = reduced_singles(rho_ab)
        bA = bloch(rho_A)
        bB = bloch(rho_B)
        pA = purity(rho_A)
        pB = purity(rho_B)
        blochs[k] = {"A": bA, "B": bB, "P_A": pA, "P_B": pB}
        print(f"  bonding:{k:<3} {bA[0]:>+9.4f} {bA[1]:>+9.4f} {bA[2]:>+9.4f} "
              f"{pA:>8.4f}   {bB[0]:>+9.4f} {bB[1]:>+9.4f} {bB[2]:>+9.4f} {pB:>8.4f}")

    # Apply chiral mirror tests for K_1 pairs (1, 5) and (2, 4)
    print()
    print("  Chiral mirror predictions:")
    print(f"    Sites 0 and 4: K_1 sign factor (−1)^a is +1 for both (a=0, 4)")
    print()
    print(f"    Predicted:")
    print(f"      ⟨X_a⟩(k=N+1−k) = +⟨X_a⟩(k)")
    print(f"      ⟨Y_a⟩(k=N+1−k) = −⟨Y_a⟩(k)   ← sharp Y sign flip")
    print(f"      ⟨Z_a⟩(k=N+1−k) = +⟨Z_a⟩(k)")
    print(f"      P_a identical")
    print()
    print(f"    Test: |observed − predicted|, expressed as relative error to mean:")
    print(f"  {'pair':<14} {'site':<6} "
          f"{'X-mirror':>11} {'Y-mirror':>11} {'Z-mirror':>11} {'P-mirror':>11}")
    for k_pair in [(1, 5), (2, 4)]:
        k_a, k_b = k_pair
        for site_label, key in [("0", "A"), ("4", "B")]:
            xa, ya, za = blochs[k_a][key]
            xb, yb, zb = blochs[k_b][key]
            pa = blochs[k_a][f"P_{key}"]
            pb = blochs[k_b][f"P_{key}"]
            # Relative deviation from prediction:
            #   X: x_b − (+x_a),
            #   Y: y_b − (−y_a) = y_b + y_a,
            #   Z: z_b − z_a,
            #   P: p_b − p_a.
            def reldev(o, p):
                base = max(abs(o), abs(p), 1e-3)
                return abs(o - p) / base
            print(f"  ({k_a},{k_b}){'':>6} site {site_label:<3} "
                  f"{reldev(xb, +xa):>11.2%} "
                  f"{reldev(yb, -ya):>11.2%} "
                  f"{reldev(zb, +za):>11.2%} "
                  f"{reldev(pb, +pa):>11.2%}")

    return blochs


def main():
    print("=" * 80)
    print("EQ-020: chiral mirror Bloch test on EXISTING Marrakesh data")
    print("=" * 80)
    print()
    print("Setup: bonding-mode receivers k=1..5 at N=5, XX-only Hamiltonian,")
    print("3 Trotter steps to t=0.8. ρ_AB tomography on (qubit 0, qubit 4).")
    print("Path on Marrakesh: [48, 49, 50, 51, 58].")
    print()
    print("Today's chiral mirror prediction (sharper than the original")
    print("|·|²-observable identity in PROOF_K_PARTNERSHIP):")
    print("  ⟨X_a⟩(k') = +(−1)^a · ⟨X_a⟩(k)")
    print("  ⟨Y_a⟩(k') = −(−1)^a · ⟨Y_a⟩(k)   for K_1 pair k' = N+1−k")
    print("  ⟨Z_a⟩(k') = +⟨Z_a⟩(k)")

    aer_blochs = analyze(
        "Aer reference (Marrakesh noise model)",
        DATA_DIR / "k_partnership_marrakesh_aer_20260425_140311.json",
    )
    hw_blochs = analyze(
        "Hardware (live Marrakesh)",
        DATA_DIR / "k_partnership_marrakesh_20260425_140913.json",
    )


if __name__ == "__main__":
    main()
