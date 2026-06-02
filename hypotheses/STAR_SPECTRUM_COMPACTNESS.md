# Star Spectrum Compactness

**Tier 3 reading, partially resolved 2026-05-19.** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e). The star topology produces an unusually compact Liouvillian spectrum at N=8, both in the imaginary range (exactly ±σ = ±N/2 when J=2γ) and in the distinct-eigenvalue count (~30× fewer unique values than the chain at the same N, γ, J). The hub geometry is the obvious suspect.

**Resolution status (2026-05-19):**
- **Reading 1 (MaxImag = σ is a hub-induced spectral cap):** **RESOLVED** via cross-N data + SU(2)/Schur-Weyl derivation in [`STAR_CONFOCAL_LIMIT.md`](STAR_CONFOCAL_LIMIT.md). Verified at N=3,4,5,6,8 (5 anchors). The cap is `max |Im(λ)| ≤ J·N/2`, saturated when `J = 2γ` giving `Im/σ = 1.0` exactly. Star is the unique tested topology that saturates this bound.
- **Reading 2 (S_(N−1) representation theory accounts for 30× degeneracy ratio):** **SHARPENED 2026-06-02** ([probe](../simulations/star_degeneracy_repn.py)). The leaf-permutation S_(N−1) gives a clean upper bound on the distinct count. The Liouville space is `C^4` (hub, S_(N−1)-trivial) `⊗ (C^4)^{⊗(N−1)}` (leaves, S_(N−1) permuting factors), so by Schur-Weyl L block-diagonalises into blocks of size `M_λ = 4·dim S^λ(C^4)`, and `distinct ≤ Σ_{λ⊢(N−1)} M_λ = 4·Σ_{λ⊢(N−1)} dim S^λ(C^4)` (= 4400 at N=8; verified as a bound at N=3,4,5,6,8). But the exact 2 275 lies **below** this cap (≈0.52×) and is **not** a clean S_7 count: adding the U(1) `(Sz_left, Sz_right)` grading refines the blocks without lowering the sum, and the residual reduction is **cross-block value-coincidence from the non-normal dephasing spectrum** (the Re-decay ladder), not a representation-theoretic effect. The tempting `distinct = |Re|·|Im|` factorisation is refuted: `Im(λ_L) ≠ H-gaps` because L is non-normal, and `distinct/(|Re|·|Im|)` falls 0.24→0.01 across N=3..6. So S_(N−1) is necessary and dominant for the compactness but not sufficient for the exact count; the clean derivable object is the `Σ M_λ` upper bound.
- **Reading 3 (integer bound at other N):** **RESOLVED** ([`STAR_CONFOCAL_LIMIT`](STAR_CONFOCAL_LIMIT.md)). `MaxImag = J·N/2` for all (N, J, γ), verified at 29 anchors (a 24-point Q-sweep plus 5 Marrakesh-convention), typed as [`StarImMaxBoundClaim`](../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs) Tier 1 derived. The `MaxImag = N/2` reading at γ=0.5 was the J=2γ specialization, not a separate rule.

## Observed compactness at N=8

From the four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files (commit 89f725e). All four configurations share Heisenberg XXX (XX+YY+ZZ) at J=1.0 with uniform Z-dephasing γ=0.5. Spectrum size is 65 536 (d²=4^N for N=8) in every case.

| Topology | MaxImag | DistinctBinnedEigenvalueCount | Compactness |
|---|---|---|---|
| chain N=8       | 5.1249 | 35 135 | 54% |
| ring N=8        | 5.6511 | 21 845 | 33% |
| star N=8        | **4.0000** | **2 275** | **3.5%** |
| K_4 + disjoint 4-chain N=8 | 5.3660 | 5 812 | 9% |

The star is the outlier on both axes. Two oddities:

1. **MaxImag = 4.000000000000002 ≈ 4 exactly = N/2.** The numerical value differs from N/2 by 2·10⁻¹⁵, i.e. one bit at machine precision. Chain, ring, K_4 + disjoint all give MaxImag in the 5.1–5.7 range, none of which look like simple integers. The star alone hits an integer bound, suggesting a structural reason for the exact ±N/2 cap rather than a numerical accident.
2. **2 275 / 65 536 = 3.5% distinct binned eigenvalues.** The chain (54%) and ring (33%) are well within the "no large degeneracy" regime; K_4 + disjoint (9%) has some component-additive degeneracy from the disconnected structure. The star at 3.5% is ~30× more degenerate than the chain at the same N, γ, J. The hub-and-leaves geometry creates massive symmetry-driven degeneracy that the other topologies do not have.

## Why this is not just "more symmetry → more degeneracy"

All four topologies have non-trivial automorphism groups (ring: dihedral D_N; chain: Z_2 reflection; K_4: S_4 × identity-on-4-chain; star: S_(N−1) permuting leaves). The star's permutation symmetry S_7 on the 7 leaves is the largest among the four, so larger degeneracy is expected, but the order-of-magnitude separation (30× fewer distinct values than chain, 10× fewer than ring) suggests more than the automorphism count alone.

The exact ±N/2 imaginary bound is the more striking observation: it is a structural integer, not a continuous function of (J, γ, topology). With Z-dephasing at rate γ=0.5 the σ shift is Nγ = 4 and the F1 palindrome predicts the spectrum is symmetric around Re = −σ = −4. The N=4 imaginary bound matches σ exactly. Possible reading: in the hub geometry the slowest oscillating modes are pinned at Im(λ) = ±σ = ±Nγ, while in chain/ring/disconnected geometries they exceed σ.

## Open structural questions

Three readings the compactness data is consistent with:

1. **MaxImag = ±N/2 is a hub-induced spectral cap that requires σ = Nγ exactly.** Testable by computing star at the same N=8 with a different γ (e.g. γ=0.25, predicted MaxImag = 2 if reading is correct; γ=1.0, predicted MaxImag = 8). If MaxImag tracks σ = Nγ tightly for the star and not for the other topologies, the reading is confirmed.

2. **The S_(N−1) symmetry of the star reduces the effective spectrum to permutation-irrep sectors.** Testable by counting the number of S_(N−1)-irreps of the Liouvillian and checking whether 2 275 distinct binned values matches the predicted count. The 30× ratio is suspicious; there could be a clean S_(N−1) representation-theoretic accounting.

3. **The integer bound at MaxImag = N/2 generalises.** At other N (N=4, 5, 6, 7 star) does MaxImag = N/2 exactly? At N=4 with γ=0.5, σ = Nγ = 2: does the star MaxImag = 2? At N=5 with γ=0.5, σ = 2.5: does it cap there? The N=8 data point is one anchor; the rule could be (a) MaxImag = σ for any star at any (N, γ), (b) MaxImag = N/2 only when σ = N/2 (i.e. γ=0.5 specifically), or (c) a coincidence at N=8.

## Typed-claim status

**Readings 1 and 3 are now typed.** The imaginary cap `MaxImag = J·N/2` is [`StarImMaxBoundClaim`](../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs) Tier 1 derived (proof [`PROOF_STAR_OPTICAL_CONFOCAL_SATURATION`](../docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md)). Note the historical phrasing `MaxImag = σ = Nγ` was the J=2γ coincidence; the universal statement is `J·N/2`, γ-independent.

**The compactness (Reading 2) is not yet typed.** The clean part, the `Σ M_λ = 4·Σ_{λ⊢(N−1)} dim S^λ(C^4)` upper bound on the distinct-eigenvalue count, could land as a typed bound. The exact distinct count (2 275 at N=8) has no closed form via the routes tried: it is non-normal dephasing value-coincidence structure (Reading 2 above), not a representation-theoretic multiplicity. An exact-count derivation would need explicit joint-block diagonalisation plus the dephasing Re-ladder coincidence count, with uncertain payoff.

## Cross-references

- Anchor data: `simulations/results/f1_n8_n9_metrics/star_N{3..6}_python.json` + `star_N8.json` (`MaxImag`, `DistinctBinnedEigenvalueCount`); same fields in the chain/ring/K_4+disjoint sister files for contrast.
- **Sharpening of Reading 1 (the MaxImag = σ saturation):** [`STAR_CONFOCAL_LIMIT.md`](STAR_CONFOCAL_LIMIT.md) (Tier 1 candidate, 2026-05-19): full cross-N table + SU(2)/Schur-Weyl derivation + connection to OPTICAL_CAVITY_ANALYSIS point-focus limit.
- Sister Tier 3 reading from the same sweep: [F1_DISSIPATION_GAP_PATTERN](F1_DISSIPATION_GAP_PATTERN.md) (extended 2026-05-19 with cross-topology cross-N gap × N² ≈ 2.20 finding for chain).
- Companion typed claim from the same sweep (closed form that promoted to Tier 1 derived): [F4KernelDimensionByComponentsClaim](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs).
- Parent framework that today's findings sharpen: [OPTICAL_CAVITY_ANALYSIS](../experiments/OPTICAL_CAVITY_ANALYSIS.md) (April 2026, the Fabry-Perot resonator structure of qubit chain + Z-deph).
- F1 verification record that produced the data: [F1GeneralTopologyVerifiedClaim](../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs).
- Related: F1 palindrome identity [F1PalindromeIdentity](../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs) (the spectrum is palindromic around Re = −σ; the ±N/2 imaginary cap is orthogonal to the F1 identity's reflection axis, so the two readings live in different parts of the spectrum and could be checked independently).
