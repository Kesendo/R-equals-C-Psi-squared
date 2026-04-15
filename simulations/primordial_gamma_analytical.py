"""
Primordial Gamma: Analytical Derivation of g(r) = <n_XY>_B
============================================================

The function g(r) = gamma_eff / gamma_B for the 3-qubit chain S-M-B
is the <n_XY>_B of the slowest S-coherence eigenmode of L.

Strategy: the XX+YY Hamiltonian conserves total excitation number.
In the single-excitation sector (which contains the S-coherence modes),
the Hamiltonian is a 3x3 tridiagonal matrix. Its eigenvectors give
the site amplitudes, and |amplitude_B|^2 = <n_XY>_B for that mode.

For general couplings J_SM and J_MB, the single-excitation Hamiltonian:

    H_1 = [[0,     J_SM,  0   ],
            [J_SM,  0,     J_MB],
            [0,     J_MB,  0   ]]

Eigenvalues: 0, +-sqrt(J_SM^2 + J_MB^2)
Eigenvectors analytically known.

The slowest S-coherence mode is the one with the smallest <n_XY>_B
AND nonzero S-site amplitude. This determines g(r).

Date: 2026-04-15
Authors: Tom and Claude (Code)
"""

import numpy as np
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# =========================================================================
# Analytical solution of the single-excitation sector
# =========================================================================
print("=" * 80)
print("Analytical derivation of g(r)")
print("=" * 80)

print("""
3-qubit chain S-M-B, XX+YY coupling:
  H_1 = [[0,     J_SM,  0   ],
          [J_SM,  0,     J_MB],
          [0,     J_MB,  0   ]]

Characteristic polynomial: -E(E^2 - J_SM^2 - J_MB^2) = 0

Eigenvalues:
  E_0 = 0
  E_+ = +sqrt(J_SM^2 + J_MB^2)
  E_- = -sqrt(J_SM^2 + J_MB^2)

Eigenvectors (unnormalized):
  |0>: (-J_MB, 0, J_SM)          [zero mode, antisymmetric]
  |+>: (J_SM, E_+, J_MB)         [bonding]
  |->: (J_SM, E_-, J_MB)         [antibonding]

Normalized |+> (= |->  up to sign of M-component):
  a_S = J_SM / sqrt(2(J_SM^2 + J_MB^2))
  a_M = +-1/sqrt(2)
  a_B = J_MB / sqrt(2(J_SM^2 + J_MB^2))

Site-B weight of bonding/antibonding modes:
  |a_B|^2 = J_MB^2 / (2(J_SM^2 + J_MB^2))

With r = J_SM/J_MB:
  |a_B|^2 = 1 / (2(r^2 + 1))

Zero-mode site-B weight:
  |a_B^(0)|^2 = J_SM^2 / (J_SM^2 + J_MB^2) = r^2 / (r^2 + 1)
""")

def g_analytical(r):
    """Analytical g(r) = |a_B|^2 for the bonding/antibonding modes."""
    return 1.0 / (2.0 * (r**2 + 1))

def g_zero_mode(r):
    """<n_XY>_B for the zero mode."""
    return r**2 / (r**2 + 1)

print("Analytical predictions:")
print(f"  g_bonding(r) = 1 / (2(r^2 + 1))")
print(f"  g_zero(r)    = r^2 / (r^2 + 1)")
print()

# =========================================================================
# Which mode is the slowest S-coherence mode?
# =========================================================================
print("=" * 80)
print("Which mode is slowest for S-coherence?")
print("=" * 80)

print("""
The Absorption Theorem: Re(lambda) = -2*gamma_B * <n_XY>_B.
Slowest = smallest <n_XY>_B (among modes with S-coherence content).

Three single-excitation modes:
  Zero mode:  <n_XY>_B = r^2/(r^2+1).     S-amplitude: J_MB/sqrt(J_SM^2+J_MB^2)
  Bonding:    <n_XY>_B = 1/(2(r^2+1)).     S-amplitude: J_SM/sqrt(2(J_SM^2+J_MB^2))
  Antibonding: same <n_XY>_B as bonding.   Same S-amplitude.

For r < 1 (J_SM < J_MB, strong bath coupling):
  g_bonding = 1/(2(r^2+1)) ~ 1/2 for r << 1
  g_zero    = r^2/(r^2+1)  ~ r^2 for r << 1
  Zero mode has SMALLER <n_XY>_B -> zero mode is the slowest.

For r > 1 (J_SM > J_MB, weak bath coupling):
  g_bonding = 1/(2(r^2+1)) ~ 1/(2r^2) for r >> 1
  g_zero    = r^2/(r^2+1)  ~ 1 for r >> 1
  Bonding mode has SMALLER <n_XY>_B -> bonding mode is the slowest.

For r = 1 (equal coupling):
  g_bonding = 1/4
  g_zero    = 1/2
  Bonding mode is slower (1/4 < 1/2).

Crossover: g_bonding = g_zero when 1/(2(r^2+1)) = r^2/(r^2+1)
  -> 1/2 = r^2 -> r = 1/sqrt(2) ~ 0.707

So:
  g(r) = r^2/(r^2+1)      for r < 1/sqrt(2)   [zero mode dominates]
  g(r) = 1/(2(r^2+1))     for r > 1/sqrt(2)   [bonding mode dominates]
""")

def g_combined(r):
    """The actual g(r): minimum of zero-mode and bonding-mode <n_XY>_B,
    restricted to modes with S-coherence content."""
    g_bond = 1.0 / (2.0 * (r**2 + 1))
    g_zero = r**2 / (r**2 + 1)
    # But: zero mode has S-amplitude proportional to J_MB, which vanishes
    # as J_MB -> 0 (r -> inf). Need to check if S-coherence content is nonzero.
    # Zero mode: S-amplitude = J_MB/sqrt(J_SM^2+J_MB^2) = 1/sqrt(r^2+1)
    # Bonding:   S-amplitude = J_SM/sqrt(2(J_SM^2+J_MB^2)) = r/sqrt(2(r^2+1))
    # Both are nonzero for finite r. So pick the smaller <n_XY>_B.
    return min(g_bond, g_zero)

# =========================================================================
# Verify against numerical data
# =========================================================================
print("=" * 80)
print("Verification: analytical vs numerical (V2 data)")
print("=" * 80)

# V2 data (gamma_B = 0.01 column, good-cavity regime)
numerical_data = [
    (100.0, 0.000050),
    (33.33, 0.000449),
    (10.0,  0.004949),
    (3.33,  0.041272),
    (1.0,   0.249987),
    (0.33,  0.100000),
    (0.1,   0.009901),
]

print(f"\n{'r':>10} {'g_numerical':>12} {'g_analytical':>13} {'g_zero':>10} {'g_bonding':>10} {'which':>10} {'error':>10}")
print("-" * 80)

max_err = 0
for r, g_num in numerical_data:
    g_z = g_zero_mode(r)
    g_b = g_analytical(r)
    g_pred = g_combined(r)
    which = "zero" if g_z < g_b else "bonding"
    err = abs(g_num - g_pred) / g_num if g_num > 0 else 0
    max_err = max(max_err, err)
    print(f"{r:10.4f} {g_num:12.6f} {g_pred:13.6f} {g_z:10.6f} {g_b:10.6f} {which:>10} {err:10.4e}")

print(f"\nMax relative error: {max_err:.4e}")

# =========================================================================
# The complete formula
# =========================================================================
print("\n" + "=" * 80)
print("THE FORMULA")
print("=" * 80)

r_cross = 1.0 / np.sqrt(2)
print(f"""
For a 3-qubit chain S-M-B with XX+YY coupling (J_SM, J_MB) and
Z-dephasing only on B at rate gamma_B, in the good-cavity regime
(gamma_B << J_MB):

  gamma_eff = gamma_B * g(r),    r = J_SM / J_MB

where:
               {{ r^2 / (r^2 + 1)       for r < 1/sqrt(2) = {r_cross:.4f}
  g(r) =      {{
               {{ 1 / (2(r^2 + 1))      for r >= 1/sqrt(2)

Asymptotics:
  r << 1: g ~ r^2 = (J_SM/J_MB)^2          [perturbative S->MB leakage]
  r = 1:  g = 1/4                            [equal coupling]
  r >> 1: g ~ 1/(2r^2) = J_MB^2/(2*J_SM^2)  [perturbative B->SM coupling]

Crossover at r = 1/sqrt(2) where both branches give g = 1/3.

Physical meaning:
  g(r) is the fraction of XY-Pauli weight that the slowest S-coherence
  eigenmode carries on the dissipative site B. Below the crossover,
  the slowest mode is the antisymmetric (zero-energy) mode. Above it,
  the slowest mode is the bonding/antibonding mode.

Verification: max relative error vs full Liouvillian numerics = {max_err:.2e}.
""")

# Save
results_dir = Path("simulations/results/primordial_gamma")
with open(results_dir / 'analytical_g_r.txt', 'w', encoding='utf-8') as f:
    f.write("Analytical derivation of g(r) for 3-qubit chain\n")
    f.write("=" * 80 + "\n\n")
    f.write("g(r) = min(r^2/(r^2+1), 1/(2(r^2+1)))\n")
    f.write(f"Crossover at r = 1/sqrt(2) = {r_cross:.6f}\n")
    f.write(f"g(1) = 1/4 exactly\n\n")
    f.write("Verification:\n")
    for r, g_num in numerical_data:
        g_pred = g_combined(r)
        err = abs(g_num - g_pred) / g_num if g_num > 0 else 0
        f.write(f"  r={r:.4f}: numerical={g_num:.6f}, analytical={g_pred:.6f}, err={err:.2e}\n")
    f.write(f"\nMax relative error: {max_err:.4e}\n")

print(f"Results saved to {results_dir / 'analytical_g_r.txt'}")
