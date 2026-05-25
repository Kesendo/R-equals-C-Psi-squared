"""Polarity {−1/2, 0, +1/2} demo via F81 Pi-decomposition (Schicht 1).

Computes M = Π·L·Π⁻¹ + L + 2σ·I and decomposes it into:
    M_sym  = (M + Π·M·Π⁻¹) / 2   (Π-symmetric, "on the 0-axis" component)
    M_anti = (M − Π·M·Π⁻¹) / 2   (Π-antisymmetric, the +1/2 vs −1/2 polarized residue)

Per F81: M_anti = L_{H_odd} when the only dissipator is Z-dephasing.
Frobenius-orthogonal: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².

For each Hamiltonian family, reports the polarity ratio
    sym_share  = ‖M_sym‖² / ‖M‖²    (= 0-axis fraction)
    anti_share = ‖M_anti‖² / ‖M‖²   (= polarized fraction)

Tom 2026-05-25: demo to decide whether to type a PolarityCoordinates primitive.
"""

import sys
sys.path.insert(0, 'simulations')
import framework as fw

import numpy as np

N = 3
gamma = 0.05

# Test cases: each is (label, terms-list).
# terms format per pi_decompose_M docstring: list of (a, b) Pauli letter tuples.
# The framework's letter type is the C# PauliLetter mirrored; framework uses strings.
cases = [
    ("Heisenberg (F1 truly, Pi^2-Z-even closed under canonical Pi)",
     [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]),
    ("XX only (Pi^2-Z-even, truly subset)",
     [('X', 'X')]),
    ("YZ + ZY (Pi^2-Z-EVEN non-truly; F108 anomaly under canonical Pi)",
     [('Y', 'Z'), ('Z', 'Y')]),
    ("XY (Pi^2-Z-ODD pure)",
     [('X', 'Y')]),
    ("XY + YX (Pi^2-Z-ODD pair)",
     [('X', 'Y'), ('Y', 'X')]),
    ("Heisenberg + XY perturbation (mixed even+odd)",
     [('X', 'X'), ('Y', 'Y'), ('Z', 'Z'), ('X', 'Y')]),
]

print(f"{'='*88}")
print(f"Polarity demo at N={N}, gamma_Z={gamma} (uniform Z-dephasing)")
print(f"{'='*88}\n")

chain = fw.ChainSystem(N=N, gamma_0=gamma)

print(f"{'Case':<60} {'||M||^2':>10} {'sym%':>8} {'anti%':>8}")
print(f"{'-'*60} {'-'*10} {'-'*8} {'-'*8}")
for label, terms in cases:
    result = fw.pi_decompose_M(chain, terms, gamma_z=gamma)
    norm_sq_M = result['norm_sq']['M']
    norm_sq_sym = result['norm_sq']['M_sym']
    norm_sq_anti = result['norm_sq']['M_anti']
    sym_share = norm_sq_sym / norm_sq_M if norm_sq_M > 1e-15 else 0.0
    anti_share = norm_sq_anti / norm_sq_M if norm_sq_M > 1e-15 else 0.0
    print(f"{label:<60} {norm_sq_M:>10.4f} {100*sym_share:>7.2f}% {100*anti_share:>7.2f}%")

print()
print("Reading:")
print("  - sym%  high (-> 100%)  : polarity concentrated on 0-axis (Pi-symmetric component)")
print("  - anti% high (-> 100%)  : polarity in the +/- 1/2 residue (Pi-antisymmetric)")
print("  - ||M||^2 = 0            : M vanishes; full F1 palindrome (no polarity residue)")
print()
print("F81 50/50 prediction: pure Pi^2-Z-odd H gives sym = anti = 50% at any N (Step 8).")
print("F1 prediction: Heisenberg-truly gives ||M||^2 = 0 (no residue at all).")
print("F108 motivation: YZ+ZY is Pi^2-Z-EVEN non-truly, canonical Pi does NOT palindromize")
print("                 it (||M||^2 > 0), but Pi_5bilinear does (the F108 closure).")
