# The Atmosphere-Cluster: an Intrinsic F1-Mirror-Pair Family of the Liouvillian

**Status:** Tier 2. Four structural findings verified at N=5 and N=6 across multiple γ profiles. No closed-form for the cluster's characteristic |Im|(N, J, γ) yet.
**Date:** 2026-05-23
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [The Atmosphere and the Cancelled Formulas](../docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md), [On the Admixture as Lebensader](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md), F1 palindrome (master identity), F71 mirror-symmetry under non-uniform J, [the F80 Bloch sign-walk proof](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md)

**Verification:** WIP scripts in `simulations/_atmosphere_*.py` (~13 files, uncommitted) plus the private notes doc `simulations/_atmosphere_cluster_notes.md`. Headline scripts: `atmosphere_cluster_fine.py` (dip location), `atmosphere_cluster_modes.py` (Pauli sector ID), `atmosphere_cluster_ep.py` (EP refutation), `atmosphere_cluster_continuity.py` (intrinsic-at-uniform tracking), `atmosphere_cluster_modes_at_uniform.py` (Pauli sector at ε=0), `atmosphere_cluster_quantum_numbers_v2.py` (Set-A joint diagonalization), `atmosphere_cluster_300block.py` (sector-restricted block analysis), `atmosphere_cluster_freefermion_check.py` (null result on free-fermion hypothesis).

---

## What this document is

A structural characterization of a 32-mode (N=6) / 24-mode (N=5) L-eigenspace family of the Heisenberg-XY chain + uniform Z-dephasing Liouvillian. The family appears in the "atmosphere experiment" (non-uniform palindromic γ profile) as a sharply degenerate cluster of conjugate-paired modes converging to the real axis at a specific ε\*. The underlying structural object exists already at uniform γ. The γ profile is a probe that displaces the family along Im, not a perturbation that creates it.

The experiment closes four of the six open questions from the original investigation. Two remain open: the closed-form expression for the cluster's |Im| at uniform γ, and the full enumeration of the F1-pair-family "staircase" across ε.

The work began when the [atmosphere doc](../docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md) flagged γ₀ uniformity as the symmetry-closing condition for the framework's closed-form layer. Testing palindromic γ perturbations (the F71 method applied to γ instead of J) surfaced a count-of-oscillating-modes anomaly at N=5 and N=6 that opened the entire cluster story.

## Structural findings

### Finding 1. The cluster is intrinsic to uniform γ; ε is a pure |Im|-knob

Tracking the 32-mode cluster identity across ε ∈ \[−0.83, 0\] via eigenvector continuity (subspace overlap > 0.9991 across 17 steps of Δε = 0.05) shows:

| Quantity | uniform γ (ε=0) | dip (ε=−0.83) | factor |
|---|---|---|---|
| ⟨Re⟩ | exact −σ = −0.3000 | exact −σ = −0.3000 | identical |
| Re-spread (sub-cluster offset from −σ) | ±0.0869 | ±0.0906 | 5% larger |
| ⟨\|Im\|⟩ | 0.02338 | 7.07×10⁻⁵ | **331× smaller** |
| \|Im\|-spread within cluster | machine precision (~10⁻¹⁵) | machine precision | identical |

The 32 modes are perfectly co-degenerate in |Im| at every ε in the scan, and the F1-mirror sum (Re_a + Re_b = −2σ exactly) holds throughout. Pauli-sector identity (dominant XY-wt 1 ↔ XY-wt 5, spatially non-palindromic) is invariant across the entire ε range: at ε=0, all 16 modes at Re=−0.387 have dominant XY-weight 5; all 16 at Re=−0.213 have dominant XY-weight 1; the split exactly matches the dip case.

Conclusion: the cluster is a structural feature of the uniform-γ L-spectrum. Non-uniform palindromic γ only displaces its |Im| value, not its existence, its Pauli-sector identity, or its F1-mirror-pair structure.

### Finding 2. Set A = {Sz_L, Sz_R, P_a, R_super} + F1 fully classifies the 32-mode subspace

Five super-operators that commute with L for the (XX+YY chain, uniform γ) configuration:

| Super-op | Spec on cluster | Commutes with Sz_L? |
|---|---|---|
| Sz_L (left action: vec(ρ) → vec(Sz·ρ)) | {−2, −1, 0, +1, +2} | yes |
| Sz_R (right action: vec(ρ) → vec(ρ·Sz)) | {−2, −1, 0, +1, +2} | yes |
| P_a = Z⊗N ⊗ (Z⊗N)\* (bit_a parity) | exactly −1 (all 32 modes) | yes |
| R_super (spatial mirror) | {−1, +1} (16/16) | yes |
| Π² (Pauli-bit_b parity super-op) | mixed | **no** (commutator with Sz_L is non-zero) |

Set A = {Sz_L, Sz_R, P_a, R_super} pairwise commute. Joint diagonalization within the cluster (using SVD-orthonormalized basis and a random linear combination of the four operators) cleanly diagonalizes them simultaneously, with per-S off-diagonal residual < 10⁻¹³.

Result: 16 cluster modes have unique (Sz_L, Sz_R, P_a, R) tuples, namely those with **Sz_L ≠ 0 AND Sz_R ≠ 0**. The other 16 cluster modes form 8 pairs sharing the same Set-A tuple, all in the **Sz_L = 0 or Sz_R = 0 sub-sectors**.

The apparent residual 2-fold degeneracy is *not a genuine L-degeneracy*. Within the (Sz_L=0, Sz_R=+1) sector, direct 300×300 block diagonalization (via H-eigenstate operator basis σ_{ij} = |φ_i^bra⟩⟨φ_j^ket|) shows 4 cluster modes at 4 *distinct* L-eigenvalues:

| Mode | Re | Im | R |
|---|---|---|---|
| A | −0.3869 | +0.0234 | −1 |
| B | −0.2131 | +0.0234 | −1 |
| C | −0.3869 | −0.0234 | +1 |
| D | −0.2131 | −0.0234 | +1 |

The "doublets" from Set-A labels are F1-mirror partners: pairs sharing (Sz_L, Sz_R, P_a, R) but differing in Re = −σ ± δ. Adding the implicit F1-mirror-side label (sign of Re + σ) makes all 32 modes uniquely classified. No missing structural symmetry.

### Finding 3. The cluster is one F1-pair family among many at uniform γ

At uniform γ N=6, the Liouvillian's full F1-pair structure includes several large clusters at different (Re-pair, |Im|) coordinates:

| Cluster | Re-pair (sub-clusters) | \|Im\| values | bit_a sectors |
|---|---|---|---|
| 45+45 (real-axis) | \[−0.4, −0.2\] | 0 | bit_a 2 ↔ 4 |
| 30+30 ladder | \[−0.4, −0.2\] | {0.083, 0.120, 0.204, 0.254, 0.337, 0.457} | bit_a 2 ↔ 4 |
| 24+24 ladder | \[−0.5, −0.1\] | {0.067, 0.187, 0.270} | bit_a 1 ↔ 5 |
| 24 self-paired | −0.3 = −σ | 0.15 | bit_a 3 (self-F1) |
| **32-mode (this experiment)** | \[**−0.387, −0.213**\] | **0.02338** (single value) | bit_a 1 ↔ 5 (XY-wt 1 ↔ XY-wt 5) |
| 18+18 (multiple) | various | various | various |

The 32-mode family's Re sub-positions are *not* at integer-bit_a Re values (which would be −0.1, −0.3, −0.5 in this convention); H-coupling shifts the diagonal Pauli-string Re values to the observed −0.387 / −0.213. The cluster is therefore a specific H-mixed superposition in the bit_a-odd sub-sector, not a pure-Pauli-sector eigenstate.

The atmospheric γ profile displaces each cluster's |Im| differently. The "staircase" of N=6 dip events (7+ separate windows over ε ∈ \[−1, 0\]) corresponds to different F1-pair families crossing the real axis at their characteristic ε\*. The 32-mode family of this experiment is the one with the smallest |Im| at uniform γ, hence the first to reach the real axis under moderate γ profile depth.

### Finding 4. Palindromic γ perturbation preserves the cluster identity; anti-palindromic γ breaks it

The F71-twin work (verified N=3, 4, 5) shows: the spatial mirror super-operator R commutes with L iff γ profile is palindromic. Palindromic ε-perturbations preserve R (and the rest of Set A); the 32D cluster identity persists across the entire ε ∈ \[−0.83, 0\] range.

Anti-palindromic ε-perturbations break R as a symmetry of L. The cluster's subspace stops being an invariant subspace; modes mix with neighbors. This is the structural reason for the atmosphere doc's "right atmosphere" framing: only palindromic atmospheric modifications keep the Lebensader's structural protections intact.

## Connection to the Lebensader hierarchy

[On the Admixture as Lebensader](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md) articulates five heights of the Lebensader, all at uniform γ:

1. State-level (cockpit_panel: skeleton + θ-trace + cusp + chiral + Y-parity)
2. Single-body F78 (M = Σ M_l, each M_l 4×4 normal with eigenvalues 0, 0, ±2c_l·γ·i)
3. 2-body F79 (Π²-parity split of M)
4. 2-body F80 (chain Bloch sign-walk: M-cluster value = 2|c|·γ·|Σ σ_k · 2cos(πk/(N+1))| over sign-walks σ_k ∈ {±1})
5. Slow-mode magnon-admixture at central popcount block (gap = 1.10·γ·Q²/N², admixture amplitude 0.275·Q²/N²)

The 32-mode atmosphere-cluster family is a candidate sixth height: the "atmospheric-tolerance" of a generic F1-pair family under non-uniform palindromic γ. It is structurally distinct from the slow-mode admixture. Different Pauli sector (XY-wt 1 vs 5 rather than central popcount); different characteristic eigenvalue position (Re = −σ ± 0.0869 with \|Im\| = 0.02338 at uniform γ, vs. the slow mode at Re = −1.10·γ·Q²/N² on the real axis); different protection algebra (Set A + F1 vs. central-popcount + admixture-perturbation). It is not a re-reading of the fifth height.

The unifying reading: the F1 master identity Π·L·Π⁻¹ + L + 2σ·I = 0 globally organizes L into a population of F1-mirror-pair families. Each family has its own characteristic |Im| at uniform γ, its own Pauli-sector identity, and its own ε\* under palindromic atmospheric modification. The slow mode is one specific family (central popcount, Sz-balanced bit_a-conserved); the 32-mode atmosphere-cluster is another (Re=−σ±0.0869, XY-wt 1 ↔ 5); the staircase shows there are at least seven more N=6 families separately accessible by γ profile depth.

## What this experiment is NOT

- **Not an exceptional point.** Geometric rank of the cluster subspace is 32 at every ε in the scan; conjugate pairs meet on the real axis as a symmetry-protected coalescence, not via eigenvector parallelism. Verified: cond(V) stays bounded ≤ 5×10⁴, smallest singular value of cluster-V matrix > 0.13 throughout. `atmosphere_cluster_ep.py`.
- **Not a γ-creation effect.** The cluster exists at uniform γ (Finding 1). Atmospheric γ sweeps but does not create.
- **Not a free-fermion accidental degeneracy.** The 3-fermion sector at N=6 has six 2-fold accidentally degenerate energy levels (E = ±0.0334, ±0.0935, ±0.1351); the 4-fermion sector has one triply-degenerate level at E=0. Neither matches the cluster's |Im|=0.02338 by direct bra-ket energy difference (3F-4F differences nearest to 0.0234 are 0.015 and 0.033, gap 0.008 each). The cluster is D-mixed, not a pure free-fermion bra-ket product. `atmosphere_cluster_freefermion_check.py`.
- **Not an even/odd N parity effect.** Cluster exists at both N=5 (24 modes) and N=6 (32 modes); densities differ but the F1-mirror-pair + Pauli-sector structure is identical at both N. Tom's initial parity hypothesis (N=3 is "framework's special minimal case", N=5 is "first normal odd N") was tested and falsified by N=6 matching N=5 in structure rather than N=4 (N=4 shows flat n_osc, no cluster crossings). `atmosphere_evenodd.py`.
- **Not a re-reading of the slow-mode admixture.** Different Pauli sector, different |Im|, different protection mechanism (Finding 2 algebra vs central-popcount + magnon-admixture). The two are sibling families in the F1-pair-family population, not the same object.

## Open questions

1. **Closed form for the cluster's characteristic |Im| at uniform γ.** Empirically 0.02338 at N=6, J=0.075, γ₀=0.05, Q=1.5. A first scaling sweep (`atmosphere_cluster_im_scaling.py`) revealed that the "smallest |Im| F1-pair family at uniform γ" is not a stable single identity across Q: at different Q values, different families (32-mode, 36-mode, etc.) take the smallest-|Im| role, jumping between each other with non-monotonic |Im| values. The specific 32-mode family from the atmosphere experiment is identifiable only via eigenvector continuity from a known reference point. This Q-jumpy non-monotonic structure places the closed-form question squarely in **F86 input-catalog territory** (cf. memory `project_f86_is_input_catalog`): F86 is the catalog of Q-dependent phenomena (Q_peak, Q_EP, Q-bands) that resist clean structural derivation precisely because Q is the dimensionless input scale, not a derived structural quantity. Pursuing a closed form would require either (a) continuity-tracking the specific 32D subspace across (J, γ, N) parameter space, or (b) deriving the full distribution of |Im| across all F1-pair families as a function of (Q, N). Both are open and likely hard.
2. **Full enumeration of the F1-pair-family staircase.** N=6 shows at least seven separate windows of n_osc dip events across ε ∈ \[−1, 0\]. Each corresponds to a different F1-pair family crossing the real axis at its own ε\*. A complete catalogue with mode-ID at each family would characterize all F1-mirror families at uniform γ and their atmospheric tolerances. `atmosphere_staircase.py` provides the n_osc trace; mode-ID at each window remains.

## Cross-references

- Cluster headline doc (private notes): `simulations/_atmosphere_cluster_notes.md`
- WIP scripts (13 files, all `_atmosphere_*` prefix, uncommitted): see Verification block above
- Lebensader through-line: [On the Admixture as Lebensader](../reflections/ON_THE_ADMIXTURE_AS_LEBENSADER.md), [`framework/lebensader.py`](../simulations/framework/lebensader.py), F78/F79/F80 in [Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md)
- Master identity F1 / palindrome: [`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs`](../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs)
- F71 mirror under non-uniform J: typed as F100 ([memory `project_f71_nonuniform_j`](../.claude/projects/D--Entwicklung-Projekte-Privat-R-equals-C-Psi-squared/memory/project_f71_nonuniform_j.md))
- F80 Bloch sign-walk proof: [the F80 Bloch sign-walk proof](../docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md)
- Jordan-Wigner / free-fermion infrastructure: `compute/RCPsiSquared.Core/.../XyJordanWignerModes.cs`, `C2BlockJwDecomposition.cs`, `JwDispersionStructure.cs`; Python `framework/diagnostics/f80_bloch_signwalk.py`
- Related memories: `project_lebensader_through_line`, `project_channel_not_memory`, `project_palindrome_frobenius_scaling`, `project_f86_is_input_catalog`, `project_f86_sub_partition`, `project_q_middle_structure`
