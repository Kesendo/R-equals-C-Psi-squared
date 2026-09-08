"""Universal Carrier Engine Demo (2026-05-12).

Reads existing rmt eigenvalue export at γ=0.05 (engine: dotnet run -- rmt),
applies the Universal Carrier reading (every rate is a₀·γ·n_XY where
a₀=2 is the Pi2 dyadic ladder anchor and n_XY is dimensionless), then
verifies the exact joint Lindblad scaling by computing the N=3 Liouvillian
at (J,γ)=(2,0.10), keeping Q=J/gamma=20 fixed, and showing the doubled grid.

Universal Carrier prediction at fixed Q: rates scale as γ; dimensionless
n_XY values are unchanged. Same Pi2 anchor a₀=2 in both cases. A γ-only
sweep at fixed J would change Q and is not this scaling identity.
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

# 3. Predict the fixed-Q (J,γ) -> 2(J,γ) grid
gamma_10 = 0.10
J_10 = 2.0
absorption_quantum_10 = 2 * gamma_10  # = 0.20
print("─" * 78)
print(f"Universal Carrier prediction at γ=0.10 (fixed Q=J/gamma=20):")
print(f"  Absorption quantum = a₀·γ = 2·{gamma_10} = {absorption_quantum_10}")
print(f"  Dimensionless n_XY values: fixed-Q invariant (same as above)")
print(f"  Re(λ) values: scaled by 2 under joint J,γ scaling")
print()

# 4. Verify with the two directly built N=3 generators. The committed
# γ=0.05 CSV above is rounded to six decimals, so it is descriptive input,
# not the numerical oracle for an exact scaling identity.
print("─" * 78)
print("Verification: build N=3 generators at (J,γ)=(1,0.05) and (2,0.10)...")

base_chain = fw.ChainSystem(N=3, J=1.0, gamma_0=gamma_05)
chain = fw.ChainSystem(N=3, J=2.0, gamma_0=gamma_10)
generator_residual = np.max(np.abs(chain.L - 2 * base_chain.L))
if generator_residual > 1e-12:
    raise AssertionError(f"fixed-Q generator scaling failed: {generator_residual:.3e}")
base_eigs = np.linalg.eigvals(base_chain.L)
eigs = np.linalg.eigvals(chain.L)
re_05_direct = np.real(base_eigs)
re_10 = np.real(eigs)
re_05_direct_unique = np.unique(np.round(re_05_direct, decimals=10))
re_05_direct_unique = re_05_direct_unique[re_05_direct_unique < -1e-10]
re_10_unique = np.unique(np.round(re_10, decimals=10))
re_10_unique = re_10_unique[re_10_unique < -1e-10]

print(f"Generator identity max|L(2J,2γ)-2L(J,γ)| = {generator_residual:.3e}")
print(f"Computed {len(eigs)} eigenvalues, {len(re_10_unique)} unique non-trivial Re values")
print()
print(f"  {'direct Re(λ), J=1':>18}  {'×2 prediction':>18}  {'Re(λ) at γ=0.10':>18}  {'match':>8}")
print(f"  {'-'*18}  {'-'*18}  {'-'*18}  {'-'*8}")
for re_05_val in re_05_direct_unique:
    pred = 2 * re_05_val
    # Eigenvalues inherit the exact generator scale; use direct producer
    # values on both sides of this displayed check.
    closest = re_10_unique[np.argmin(np.abs(re_10_unique - pred))]
    match = "✓" if abs(closest - pred) < 1e-8 else "✗"
    print(f"  {re_05_val:>18.6f}  {pred:>18.6f}  {closest:>18.6f}  {match:>8}")

print()
print("═" * 78)
print("Conclusion (Universal Carrier reading):")
print("  • a₀ = 2 (Pi2 dyadic ladder, n=0 anchor) is γ-invariant")
print("  • n_XY values are invariant at fixed Q (Pauli-letter expectations)")
print("  • Re(λ) = −a₀·γ·n_XY scales linearly under joint J,γ scaling")
print("  • What's universal: a₀ and the n_XY structure (the framework's")
print("    own structural data). A γ-only sweep at fixed J is not this identity.")
print("═" * 78)
