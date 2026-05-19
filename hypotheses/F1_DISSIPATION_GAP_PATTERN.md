# F1 Dissipation Gap Pattern

**Tier 3 (reading) sharpened to Tier 3 (empirical scaling).** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e) with the 4 N=8 anchors. Extended 2026-05-19 with chain/ring/star × N=3..6 Python anchors plus the N=9 chain run via the MklDirect bridge (commit abb2d52). The dissipation gap of the Heisenberg + Z-dephasing Liouvillian shows a clean **chain-topology scaling law** `gap × N² ≈ 2.20` for N ≥ 4, with ring and star following different patterns. Bond count alone does not predict it; the per-topology scaling shape is the structural fingerprint.

## Observed gaps at N=8 (initial 4-point anchor)

From the four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files' `DissipationGap` fields (commit 89f725e). All four configurations share Heisenberg XXX (XX+YY+ZZ) at J=1.0 with uniform Z-dephasing γ=0.5.

| Topology | Bonds B | Components | Dissipation gap |
|---|---|---|---|
| chain N=8       | 7 | 1 | 0.0344 |
| star N=8        | 7 | 1 | 0.0870 |
| ring N=8        | 8 | 1 | 0.1339 |
| K_4 + disjoint 4-chain N=8 | 9 | 2 | 0.1362 |

## Cross-topology cross-N extension (2026-05-19 Python anchors)

Python anchors via `simulations/_f1_topology_heisenberg_small_n_anchor.py` extended chain, ring, star to N=3..6 dense numpy eigvals. The same Heisenberg J=1, γ=0.5 convention. JSON files at `simulations/results/f1_n8_n9_metrics/{chain,ring,star}_N{3..6}_python.json`.

| Topology | N | Dissipation gap | gap × N² | Notes |
|---|---:|---:|---:|---|
| chain | 3 | 0.2697 | **2.43** | boundary case, also = star N=3 (Y-graph isomorphism) |
| chain | 4 | 0.1362 | **2.18** | |
| chain | 5 | 0.0884 | **2.21** | |
| chain | 6 | 0.0607 | **2.18** | |
| chain | 8 | 0.0344 | **2.20** | C# block-spectrum, matches Python pattern |
| chain | 9 | 0.02728 | **2.21** | C# block-spectrum via MklDirect bridge 2026-05-19, predicted 2.20/81 = 0.02716, observed 0.02728 (within 0.4 %) |
| ring | 3 | 0.8278 | 7.45 | = K_3 = triangle |
| ring | 4 | 0.3795 | 6.07 | |
| ring | 5 | 0.3173 | 7.93 | |
| ring | 6 | 0.2300 | 8.28 | |
| ring | 8 | 0.1339 | 8.57 | |
| star | 3 | 0.2697 | 2.43 | isomorphic to chain N=3 |
| star | 4 | 0.2099 | 3.36 | |
| star | 5 | 0.1637 | 4.09 | |
| star | 6 | 0.1300 | 4.68 | |
| star | 8 | 0.0870 | 5.57 | |

## The chain-topology scaling law: gap × N² ≈ 2.20 (Tier 3 empirical)

For chain Heisenberg at J=1, γ=0.5, the dissipation gap satisfies

    gap(chain, N) · N²  ≈  2.20  ± 0.02   for N ≥ 4

across **5 data points** in the flat plateau (N=4, 5, 6, 8, 9; N=3 is a boundary case at 2.43, Y-graph isomorphism with star, breaks the chain-specific limit; N=7 is absent because no chain anchor exists at that N in the current sweep). The flat plateau across N ∈ {4, 5, 6, 8, 9} is the cleanest topology-specific scaling observation in the data, predictive to 0.4 % at N=9 (predicted 2.20/81 = 0.02716, observed 0.02728).

**Physical interpretation.** A 1D dispersive chain Hamiltonian + per-site Z-dephasing produces a "diffusive" Liouvillian. The slowest decay mode has wavevector k_min ∝ 1/N (open-boundary modes); the decay rate scales as γ · k_min² ∝ 1/N². The prefactor 2.20 packs the J=2γ coefficient + chain Hamiltonian structure factor into a single number that has not yet been derived in closed form.

**The constant 2.20 is not obviously a clean expression** of π², 4γJ, γ²+J², or other obvious combinations:

- 4γJ = 4 · 0.5 · 1 = 2.00 (off by 9% from 2.20)
- 4γ² + J² = 1.25 (off significantly)
- π²·γ²/2 ≈ 1.23 (no)
- 2(J² + γ²) = 2.5 (no)

The prefactor likely involves a more delicate dispersion-integral evaluation. Closed form is open.

**N=9 chain verification (2026-05-19, landed):** prediction `gap ≈ 2.20/81 = 0.02716` versus observed `gap = 0.02728` from the MklDirect bridge run. Match within 0.4 %. The chain `gap × N² ≈ 2.20` empirical scaling now spans N ∈ {4, 5, 6, 8, 9} (5 anchors; N=7 chain not in the current sweep); the flat plateau survives the N=9 extension.

## Ring and star follow different patterns

Ring N=3..8: `gap × N²` ranges 6.07 → 8.57, **not flat**. The ring's cyclic symmetry creates additional degeneracies and the dispersion structure differs from open chain. No clean scaling law identified.

Star N=3..8: `gap × N²` grows monotonically 2.43 → 3.36 → 4.09 → 4.68 → 5.57, **not 1/N² scaling**. The hub-spoke geometry produces a different gap-scaling family. Star's spectral structure connects to the SU(2)/Schur-Weyl decomposition (see [`STAR_CONFOCAL_LIMIT.md`](STAR_CONFOCAL_LIMIT.md) for the related Im(λ) = σ saturation observation).

The three topologies have qualitatively different scaling laws, consistent with the original observation that bond count alone does not predict the gap. The right framework distinguishes:
- **Chain (1D open path)**: clean 1/N² diffusive scaling
- **Ring (1D periodic)**: distinct pattern from chain, likely tied to cyclic dispersion + degeneracy
- **Star (hub-spoke)**: distinct pattern tied to maximally-confocal geometry, slower N-decay than chain

## Why bond count alone fails (sharpened)

At N=8: chain (B=7) gap 0.0344, star (B=7) gap 0.0870, ratio 2.5×. Same bond count, different geometry. Even more starkly at small N where the difference grows: chain N=4 (B=3) gap 0.136, star N=4 (B=3) gap 0.210, ratio 1.54.

The 2026-05-19 cross-topology data confirms: bond count is irrelevant, **graph dispersion structure** is the right parameterisation. Chain has linear-band dispersion (slow modes at k ∝ 1/N → 1/N² gap); star has hub-localised modes (different scaling family); ring has cyclic dispersion with degenerate slowest modes.

## Open structural questions (refined post-extension)

1. **Closed form for the chain prefactor 2.20.** The 1/N² scaling is physically motivated by 1D diffusion; the prefactor 2.20 ≈ 4γJ (with J=1, γ=0.5 → 2.00) is empirically off by 9% from the simplest scaling. Bethe-ansatz or magnon-dispersion derivation may identify the exact coefficient.

2. **Star and ring scaling laws.** Ring `gap × N²` grows from 6 to 8.5; star `gap × N²` grows from 2.4 to 5.6. Neither pattern fits 1/N², 1/N, or exp(-αN). The functional forms are open.

3. **J ≠ 2γ regime.** All data points are at J = 2γ (J=1, γ=0.5). Sweep needed to test whether the chain 2.20 prefactor decomposes as J·a + γ·b or some product form.

4. **K_4 + disjoint behaviour.** The disconnected case at N=8 still gives gap 0.1362, similar to the ring. The original "rate-limiting component dominates" hypothesis remains untested at small N (need K_4 alone at N=4 and 4-chain alone at N=4 to check the min-vs-additive question).

5. **Connection to F2 / F3.** F2 dispersion claims and F3 decay rate bounds (`min rate = 2γ`, `max rate = 2(N-1)γ` per the Absorption Theorem) describe the spectral envelope. The dissipation gap sits at the lower edge of this envelope; its 1/N² scaling for chain may be the finite-size correction to F3's `min rate = 2γ` thermodynamic limit. Cross-link to F3 documentation when the closed form lands.

## Promotion path

This entry can move from Tier 3 reading to **Tier 1 candidate typed claim** when:
- The chain 2.20 prefactor admits a closed-form derivation (Bethe ansatz / dispersion integral / specific Liouvillian Jordan-block structure).
- The ring and star scaling families are characterized (even if the formulas differ across topologies, having three closed forms in hand is enough).
- The J ≠ 2γ sweep confirms whichever dependence emerges (linear in J or γ, joint product, or independent).

In that case the typed claim sits in `compute/RCPsiSquared.Core/Symmetry/` as `F_DissipationGapClosedForm` (or similar) with a `Predict(topology, N, J, γ)` method dispatching on topology family.

## Cross-references

- Anchor data: `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json`
- Companion typed claim from the same sweep (the closed-form discovery that did promote): [F4KernelDimensionByComponentsClaim](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 derived as of 2026-05-19; landed Tier 1 candidate 2026-05-18, promoted after DEGENERACY_PALINDROME Result 2 was identified as the connected-case upper-bound closure; kernel-dim factorisation across components)
- Sister Tier 3 reading from the same sweep: [STAR_SPECTRUM_COMPACTNESS](STAR_SPECTRUM_COMPACTNESS.md)
- F1 verification record that produced the data: [F1GeneralTopologyVerifiedClaim](../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs)
- Related: [F1_Pattern_GENERAL_TOPOLOGY proof](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) closed form for the F1 residual norm, which is graph-additive bit-exactly across the same 4 N=8 topologies; the dissipation gap is the next structural quantity beyond the residual norm and does NOT follow the same simple (B, D2) parameterisation)
