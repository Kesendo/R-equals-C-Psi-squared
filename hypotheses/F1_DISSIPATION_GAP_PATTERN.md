# F1 Dissipation Gap Pattern

**Tier 3 (reading, not yet a typed claim).** Observation surfaced from the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e). The dissipation gap (slowest decay = min{|Re(λ)| : Re(λ) < 0}) of the Heisenberg + Z-dephasing Liouvillian at fixed N=8, J=1, γ=0.5 varies non-trivially across topologies. Bond count alone does not predict it; connectivity / branching geometry also enters.

## Observed gaps at N=8

From the four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files' `DissipationGap` fields (commit 89f725e). All four configurations share Heisenberg XXX (XX+YY+ZZ) at J=1.0 with uniform Z-dephasing γ=0.5.

| Topology | Bonds B | Components | Dissipation gap |
|---|---|---|---|
| chain N=8       | 7 | 1 | 0.0344 |
| star N=8        | 7 | 1 | 0.0870 |
| ring N=8        | 8 | 1 | 0.1339 |
| K_4 + disjoint 4-chain N=8 | 9 | 2 | 0.1362 |

## Why this is not just gap ∝ B(G)

If the gap scaled with bond count alone, chain and star (both B=7) should have the same value. They differ by a factor 2.5×. The star concentrates dispersal around the hub; the chain spreads it linearly. Connectivity / branching / spectral radius of the graph Laplacian appears to enter beyond the bond count.

## Open structural questions

The data is too thin to pose anything tighter than candidates. Three readings that the four-point pattern is consistent with:

1. **Graph-Laplacian spectral gap.** The Liouvillian dissipation gap could correlate with the second-smallest eigenvalue of the underlying graph Laplacian (Fiedler value). Star and chain at N=8 have well-separated Fiedler values; star = 1 (hub-driven), chain = 2(1 − cos(π/N)) ≈ 0.152. Direction proportional but the ratio doesn't match the observed gap ratio (0.0870 / 0.0344 ≈ 2.53 vs 1 / 0.152 ≈ 6.6). Could be a function-of-Fiedler relationship, not strict proportionality.

2. **Hub vertex contribution.** Star at N=8 has degree-7 hub plus 7 leaves; ring has all-degree-2 sites; chain has end vertices of degree 1 and interior degree 2. The hub may push the slow mode toward a localised structure that decays faster than the chain's spread-out modes. Suggests a degree-distribution-driven correction to the Laplacian-gap reading.

3. **Disconnected components add or compete?** K_4 + 4-chain (0.1362) sits within 2% of ring (0.1339), even though K_4 is much more strongly connected than ring. One reading: the rate-limiting component dominates the gap (the slow side of a parallel coupling), so the 4-chain piece of K_4 + 4-chain sets the floor and brings the disconnected case below what K_4-alone would give. Testable by computing the gap for K_4 alone at N=4 and the 4-chain alone at N=4 and checking whether K_4 + 4-chain gap = min(K_4 gap at N=4, 4-chain gap at N=4) or some other combination rule.

## Not yet a typed claim

This entry stays in `hypotheses/` (Tier 3) until at least three more data points land:

- A second N (N=7 or N=9 once the ILP64 bridge is in place) to test whether the gap ordering chain < star < ring < disconnected holds at other system sizes.
- A wider graph zoo at N=8 (path, K_4 + K_4, double-star, wheel, complete bipartite) to test the hub-vs-spread hypothesis.
- One γ value above and one below 0.5 to test whether the relative ordering is γ-stable or γ-band-dependent.

If a clean closed form lands (gap as a function of (B, Fiedler, max-degree, components)), promote to Tier 1 candidate in `compute/RCPsiSquared.Core/Symmetry/` with a `Predict(graphInvariants)` method and the data table as anchors.

## Cross-references

- Anchor data: `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json`
- Companion typed claim from the same sweep (the closed-form discovery that did promote): [F4KernelDimensionByComponentsClaim](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 derived as of 2026-05-19; landed Tier 1 candidate 2026-05-18, promoted after DEGENERACY_PALINDROME Result 2 was identified as the connected-case upper-bound closure; kernel-dim factorisation across components)
- Sister Tier 3 reading from the same sweep: [STAR_SPECTRUM_COMPACTNESS](STAR_SPECTRUM_COMPACTNESS.md)
- F1 verification record that produced the data: [F1GeneralTopologyVerifiedClaim](../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs)
- Related: [F1_Pattern_GENERAL_TOPOLOGY proof](../docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) closed form for the F1 residual norm, which is graph-additive bit-exactly across the same 4 N=8 topologies; the dissipation gap is the next structural quantity beyond the residual norm and does NOT follow the same simple (B, D2) parameterisation)
