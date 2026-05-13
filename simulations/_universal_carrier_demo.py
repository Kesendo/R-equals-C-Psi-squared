"""Universal Carrier Engine Demo (2026-05-12).

Reads existing rmt eigenvalue export at γ=0.05 (engine: dotnet run -- rmt),
applies the Universal Carrier reading (every rate is a₀·γ·n_XY where
a₀=2 is the Pi2 dyadic ladder anchor and n_XY is dimensionless), then
verifies γ-scaling by computing the N=3 Liouvillian directly at γ=0.10
and showing the doubled grid.

Universal Carrier prediction: rates scale as γ; dimensionless n_XY values
(γ-invariant) are the same. Same Pi2 anchor a₀=2 in both cases.
"""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw

# 1. Read engine data at default γ=0.05
csv_path = Path(__file__).parent / "results" / "rmt_eigenvalues_N3.csv"
data = np.loadtxt(csv_path, delimiter="\t", skiprows=1)
re_05 = np.unique(np.round(data[:, 0], decimals=6))
re_05 = re_05[re_05 < 0]  # drop kernel (Re=0); focus on non-trivial decay rates

# 2. Apply Universal Carrier reading: n_XY = -Re / (2γ)
gamma_05 = 0.05
absorption_quantum_05 = 2 * gamma_05  # = a₀ · γ = 0.10
n_xy_values = -re_05 / absorption_quantum_05

print("═" * 78)
print("Universal Carrier Engine Demo: N=3 chain, Heisenberg + Z-dephasing")
print("═" * 78)
print()
print(f"Engine output: {csv_path.name} ({len(data)} eigenvalues at γ=0.05)")
print(f"Unique non-trivial Re(λ) values: {len(re_05)}")
print(f"Absorption quantum = a₀·γ = 2·{gamma_05} = {absorption_quantum_05}")
print()
print(f"  {'Re(λ) at γ=0.05':>18}  {'n_XY = -Re/(2γ)':>18}  {'reading':>30}")
print(f"  {'-'*18}  {'-'*18}  {'-'*30}")
for re, nxy in zip(re_05, n_xy_values):
    if abs(nxy - round(nxy)) < 0.001:
        reading = f"a₀·γ·{int(round(nxy))}  pure-weight rung"
    elif abs(nxy - 4/3) < 0.001:
        reading = "a₀·γ·4/3  F33 mixed (8γ/3)"
    elif abs(nxy - 5/3) < 0.001:
        reading = "a₀·γ·5/3  F33 mixed (10γ/3)"
    else:
        reading = f"a₀·γ·{nxy:.4f}"
    print(f"  {re:>18.6f}  {nxy:>18.6f}  {reading:>30}")
print()

# 3. Predict γ=0.10 grid (doubled)
gamma_10 = 0.10
absorption_quantum_10 = 2 * gamma_10  # = 0.20
print("─" * 78)
print(f"Universal Carrier prediction at γ=0.10 (carrier-rate doubled):")
print(f"  Absorption quantum = a₀·γ = 2·{gamma_10} = {absorption_quantum_10}")
print(f"  Dimensionless n_XY values: γ-invariant (same as above)")
print(f"  Re(λ) values: scaled by 2 (= 0.10/0.05 ratio)")
print()

# 4. Verify by computing N=3 Liouvillian at γ=0.10
print("─" * 78)
print("Verification: build N=3 chain Liouvillian at γ=0.10 in Python framework...")

chain = fw.ChainSystem(N=3, J=1.0, gamma_0=gamma_10)
L = chain.L
eigs = np.linalg.eigvals(L)
re_10 = np.real(eigs)
re_10_unique = np.unique(np.round(re_10, decimals=6))
re_10_unique = re_10_unique[re_10_unique < -1e-10]

print(f"Computed {len(eigs)} eigenvalues, {len(re_10_unique)} unique non-trivial Re values")
print()
print(f"  {'Re(λ) at γ=0.05':>18}  {'×2 prediction':>18}  {'Re(λ) at γ=0.10':>18}  {'match':>8}")
print(f"  {'-'*18}  {'-'*18}  {'-'*18}  {'-'*8}")
for re_05_val in re_05:
    pred = 2 * re_05_val
    # find closest match in computed γ=0.10 grid
    closest = re_10_unique[np.argmin(np.abs(re_10_unique - pred))]
    match = "✓" if abs(closest - pred) < 1e-4 else "✗"
    print(f"  {re_05_val:>18.6f}  {pred:>18.6f}  {closest:>18.6f}  {match:>8}")

print()
print("═" * 78)
print("Conclusion (Universal Carrier reading):")
print("  • a₀ = 2 (Pi2 dyadic ladder, n=0 anchor) is γ-invariant")
print("  • n_XY values are γ-invariant (Pauli-letter exposure counts)")
print("  • Re(λ) = a₀·γ·n_XY scales linearly with γ")
print("  • What's universal: a₀ and the n_XY structure (the framework's")
print("    own structural integers). What scales: γ alone (the carrier).")
print("═" * 78)
