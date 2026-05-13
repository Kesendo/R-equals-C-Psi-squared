"""Falsification witness: γ₀ = 1/8 has no operational signature beyond
substrate-invariance (2026-05-13).

Speculation tested: γ₀ = 1/8 = a_4 on the Pi2 dyadic ladder would put
2γ₀ = 1/4 = a_3 at the fold-boundary anchor. Numerologically attractive
but structurally suspect (dimensional mismatch: 2γ₀ is a rate, 1/4 is
CΨ density; substrate-invariance says any γ₀ choice gives same
dimensionless physics).

Test (N=3 chain, J=1.0, γ₀ ∈ {0.05, 0.10, 0.125=1/8}):
- 9 distinct non-trivial Re values at all three γ₀; no novel degeneracy
- Pure-weight rungs (n_XY ∈ {1,2,3}, multiplicities 14,14,4) γ-invariant
- F33 mixed Q-dependent in the same pattern at all three γ₀
- 14 modes at Re = -1/4 at γ₀ = 1/8 = n_XY=1 cluster at γ-scaled position

Result: γ₀ = 1/8 is operationally indistinguishable from γ₀ = 0.05 or 0.10.
The 2γ₀ = 1/4 match is pure number coincidence across different physical
units. Substrate-invariance per UniversalCarrierClaim is confirmed.

See AbsorptionTheoremClaim, UniversalCarrierClaim.DefaultGammaZero, and
ValidateAgainstPythonStepFTests for the substrate-invariance authority.
"""
import sys
import numpy as np
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw

print("=" * 78)
print("Falsification witness: gamma_0 = 1/8 has no operational signature")
print("=" * 78)
print()

results = {}
for gamma in [0.05, 0.10, 0.125]:
    chain = fw.ChainSystem(N=3, J=1.0, gamma_0=gamma)
    eigs = np.linalg.eigvals(chain.L)
    re = np.real(eigs)
    re_dim = re / gamma
    buckets_dim = Counter(round(r, 5) for r in re_dim if r < -1e-8)
    results[gamma] = (re, buckets_dim)

print("Dimensionless rates (Re / gamma) at three gamma_0 values:")
print(f"  Q = J/gamma:  {[round(1.0/g, 3) for g in [0.05, 0.10, 0.125]]}")
print()

all_keys = set()
for _, buckets in results.values():
    all_keys.update(buckets.keys())

print(f"  {'Re/gamma (dimensionless)':>26}  {'g=0.05':>8}  {'g=0.10':>8}  {'g=0.125':>8}  {'g-invariant?':>14}")
print(f"  {'-'*26}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*14}")
for key in sorted(all_keys):
    counts = [results[g][1].get(key, 0) for g in [0.05, 0.10, 0.125]]
    invariant = "YES" if counts[0] == counts[1] == counts[2] else "no (Q-dependent)"
    print(f"  {key:>26.5f}  {counts[0]:>8}  {counts[1]:>8}  {counts[2]:>8}  {invariant:>14}")

print()
print("-" * 78)
print("Specific check at gamma_0 = 1/8: how many eigenvalues at Re = -1/4 exactly?")
re_125 = results[0.125][0]
exact_quarter = sum(1 for r in re_125 if abs(r - (-0.25)) < 1e-10)
re_05 = results[0.05][0]
at_2gamma_05 = sum(1 for r in re_05 if abs(r - (-0.10)) < 1e-10)
print(f"  Modes with Re = exactly -1/4 at gamma_0 = 0.125:  {exact_quarter}")
print(f"  Modes with Re = exactly -2*gamma at gamma_0 = 0.05: {at_2gamma_05}")
match = "yes (same n_XY=1 cluster, just gamma-scaled)" if exact_quarter == at_2gamma_05 else "no (different multiplicity!)"
print(f"  Same multiplicity? {match}")

print()
print("-" * 78)
print("Distinct non-trivial Re values per gamma_0 (looking for novel degeneracy):")
for gamma in [0.05, 0.10, 0.125]:
    re_g = results[gamma][0]
    distinct = len(set(round(r, 8) for r in re_g if r < -1e-8))
    print(f"  gamma_0 = {gamma:>5}:  {distinct} distinct non-trivial Re values")

print()
print("=" * 78)
print("Result: gamma_0 = 1/8 has NO operational signature beyond substrate-invariance.")
print()
print("Pure-weight rungs are gamma-invariant (Re/gamma identical across all gamma_0);")
print("F33-mixed rates are Q-dependent (same pattern at all gamma_0); distinct count")
print("is 9 at all gamma_0; no novel degeneracy at gamma_0 = 1/8. The numerical match")
print("2*gamma = 1/4 at gamma_0 = 1/8 is a pure number coincidence across different")
print("physical units (rate vs CPsi density).")
print()
print("Hypothesis 'gamma_0 = 1/8 as on-ladder convention' is FALSIFIED.")
print("=" * 78)
