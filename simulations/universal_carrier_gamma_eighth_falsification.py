"""Falsification witness: γ₀ = 1/8 has no operational signature beyond
substrate-invariance (2026-05-13).

# What this script falsifies

During the Universal Carrier session (2026-05-12, see commits dd27b6f /
aec0772 / f15fc8f), Claude (assistant) proposed the speculation:

    "If γ₀ = 1/8 = a_4 on the Pi2 dyadic ladder, then 2γ₀ = 1/4 = a_3
     would coincide with the fold-boundary anchor. That would be the
     'on-ladder' choice for γ₀."

The speculation was numerologically attractive (everything on the Pi2
ladder) but had three a-priori problems:

1. Dimensional mismatch: 2γ₀ is a rate (1/time), 1/4 is a CΨ density
   (dimensionless). Same number, different physical things.
2. Substrate-invariance (per F86KnowledgeBase γ-tests + UniversalCarrierClaim):
   γ₀ is a code/substrate convention; dimensionless physics (Q = J/γ₀,
   K = γ·t, n_XY-grid) is γ-invariant. No γ₀ choice should be physically
   privileged.
3. Pattern-matching urge: the speculation was driven by "wouldn't it be
   neat if everything aligned" rather than by any framework constraint.

This script tests whether γ₀ = 1/8 produces ANY operational signature that
other γ₀ values do not.

# What was tested

For N=3 chain, J=1.0, three γ₀ values {0.05, 0.10, 0.125 = 1/8}:
- Liouvillian eigenvalue spectrum (full 4³ = 64 modes)
- Dimensionless rate Re/γ per cluster
- Multiplicity preservation across γ₀ values
- Distinct-Re count (looking for novel degeneracies at γ₀ = 1/8)
- Direct count of modes at Re = -1/4 specifically at γ₀ = 1/8

# Conclusion

γ₀ = 1/8 is operationally indistinguishable from γ₀ = 0.05 or 0.10.

- Pure-weight rungs (n_XY ∈ {1, 2, 3}, multiplicity 14, 14, 4) show
  identical Re/γ values across all three γ₀: γ-invariant per Absorption
  Theorem.
- F33-mixed clusters shift slightly with Q = J/γ across all three γ₀
  (J=1 fixed, varying γ varies Q), in the same Q-dependent pattern; no
  γ₀ = 1/8 specific anomaly.
- 9 distinct non-trivial Re values at all three γ₀: no novel degeneracy
  at γ₀ = 1/8.
- The "14 modes at Re = -1/4" at γ₀ = 1/8 is just the n_XY = 1 cluster
  at its γ-scaled position; identical to "14 modes at -0.10" at γ₀ =
  0.05 or "14 modes at -0.20" at γ₀ = 0.10.

The 2γ₀ = 1/4 numerical match at γ₀ = 1/8 is a pure number coincidence
across different physical units (rate vs CΨ density). It produces no
observable signature.

The hypothesis "γ₀ = 1/8 as on-ladder convention" is empty; substrate-
invariance (per UniversalCarrierClaim docstring) is the right reading.

# Why this script is committed (rather than archived as WIP)

Per CLAUDE.md "Negative results matter: document failures with the same
rigor as successes." This file is the falsification witness. Future readers
of the project (including future-Claude) are protected from re-spinning the
same pattern-matching speculation by being able to re-run this and confirm
the empty result.

# Anchors

- Speculation source: conversation 2026-05-12 (memory layer; see
  user-local Claude memory project_algebra_handshake_at_quarter.md
  open-territory section, where γ₀ = 1/8 was listed and is now removed
  per this falsification).
- Substrate-invariance authority: compute/RCPsiSquared.Core/Symmetry/
  UniversalCarrierClaim.cs DefaultGammaZero const + ExtraChild
  "substrate invariance" + ValidateAgainstPythonStepFTests covering
  γ₀ ∈ {0.025, 0.05, 0.10}.
- Sibling demos (positive witnesses for substrate-invariance):
  simulations/_universal_carrier_demo.py and
  simulations/_universal_carrier_zoom.py.
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
