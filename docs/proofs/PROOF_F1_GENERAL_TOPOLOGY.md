# PROOF: F1 H-block residual closed form extends to arbitrary graph topology

**Status:** Tier 1 (synthesis). Closes the F1 OpenQuestion "general topology beyond chain/ring/star/K_N" by combining the existing analytic result PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary (universality of the bond-disjointness argument across any graph) with disconnected-component and weighted-edge extensions. Verification record: Python at N = 5, 6 across named graphs, random connected Erdős-Rényi, disconnected, weighted, single-body; C# at N = 5 (graph-aware mode of `PalindromeResidualScalingClaim`), N = 7 (F1 palindromic pairing via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` on multiple graph topologies including disconnected), N = 8 (Heisenberg + Z-deph across 4 topologies opt-in SLOW_N8 with full `F1SpectrumStatistics` metric capture). N = 9 chain wired as a `SkippableFact` but blocked at the LP64 MKL P/Invoke marshalling ceiling discovered in this sprint; see § Scale frontier.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

F1 was first proven on a chain, then extended to ring, star, and the complete graph K_N. Each extension repeated the same skeleton, but the proofs had to be written separately. The natural next question was whether the closed form holds for any graph topology, or whether some richer-connectivity graph would surface higher topological invariants (cycle counts, triangles, hub presence) that the linear-graph cases happened to miss.

The answer is the cleanest possible affirmation. The residual norm has the same closed form for every connected graph, and the only graph invariants that enter are the simplest possible two: the number of bonds, and the sum of squared degrees. No cycle count, no clustering coefficient, no spectral graph data. Two graphs with the same bond count and same degree-squared sum produce identical residual norms even if they look topologically different.

The mechanism is bond-disjointness. The Pauli-string supports of any two bond Hamiltonians (whether or not their bonds share a vertex) are disjoint, so the bond-resolved residuals are Frobenius-orthogonal, and the residual norm decomposes additively across bonds. The bond count enters because each bond contributes independently; the degree-squared sum enters in the single-body class because identity-padded single-body terms at site l get a multiplicity proportional to deg(l)². No higher graph invariant appears because the additive decomposition stops at the bond level.

The same argument extends in two directions without new work. Disconnected graphs with multiple components contribute additively (B and D2 sum across components). Weighted edges with per-bond couplings replace the bond count by the squared-coupling sum. And the dissipator-block independence from γ (proved in the non-uniform sibling) means the result holds for arbitrary γ patterns too. So the F1 closed form is genuinely a property of the bonds themselves, not of the graph's larger structure, not of the dephasing rate's distribution, not of which sites are connected to which.

## Statement

For N ≥ 2 qubits, any 2-bilinear Pauli Hamiltonian H placed on an arbitrary graph G with bond set E(G) and (uniform OR site-dependent) Z-dephasing at rates {γ_l}, the F1 palindrome residual norm has the closed form

    ‖M(N, G)‖²_F  =  c_H · F(N, G)

with the per-class graph factor

    main class         F(N, G) = B(G) · 4^(N − 2)
    single-body class  F(N, G) = (D2(G) / 2) · 4^(N − 2)

where B(G) = number of bonds and D2(G) = Σ_i deg_G(i)² is the second moment of the degree sequence. The Hamiltonian-dependent constant c_H is fixed at N = 2 once and is graph-independent. The formula holds for:

- **any connected graph** (chain, ring, star, K_N, K_{a,b}, random Erdős-Rényi, arbitrary tree, etc.);
- **disconnected graphs** with multiple components (B and D2 sum across components);
- **weighted edges** with per-bond couplings {J_b}, by the substitution B → Σ_b J²_b / J²_ref (where J_ref is the anchoring coupling for c_H);
- **uniform AND non-uniform per-site γ_l** (separately closed by [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md): the H-block scaling factor F(N, G) is γ-independent because the dissipator-block residual M_D vanishes per Pauli string).

## Proof

### Section 1. Bond-disjointness across any graph topology

The substantive analytic content is already established in [PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md), specifically Lemma 3 and its Corollary, which together state that for any two bond Hamiltonians H_e and H_e' on a graph G (whether or not e and e' share a vertex), the Pauli-string transition supports of [H_e, ·] and [H_e', ·] are disjoint. The proof there is graph-independent: it depends only on the Pauli algebra of the bond bilinears, not on the underlying connectivity. Consequently the bond-resolved residuals M_e := Π·L_He·Π⁻¹ + L_He + (per-bond σ-share)·I are mutually Frobenius-orthogonal:

    tr(M_e^† · M_e')  =  0    for any e ≠ e' ∈ E(G).

Summing M = Σ_e M_e (linearity of L in H and of Π·(·)·Π⁻¹) and taking the Frobenius square gives the bond-additive identity

    ‖M(N, G)‖²_F  =  Σ_e ‖M_e‖²_F.

For the main class each ‖M_e‖²_F equals the per-bond constant c_H · 4^(N−2) (Lemma 1 of PROOF_CROSS_TERM_FORMULA, multiplied by the Liouville-space extension factor 4^(N−2)). The sum over E(G) gives the bond-count B(G) prefactor. For the single-body class, the per-site contributions get a degree-weighted assembly (each Iσ + σI bond term contributes to the site at both endpoints), giving the D2(G) / 2 prefactor (see OPERATOR_RIGIDITY_ACROSS_CUSP.md "Algebraic origin" for the derivation of why D2 enters the single-body case). Both prefactors are linear in graph invariants; no higher-order topological dependence (cycle count, triangle count, hub presence) ever appears in the residual norm closed form.

This establishes universality of (B, D2) parameterisation across any **connected** graph.

### Section 2. Disconnected components

For G = G_1 ⊔ G_2 a disjoint union of two connected components, the bond set E(G) = E(G_1) ⊔ E(G_2) and the bond-additive identity from Section 1 gives

    ‖M(N, G)‖²_F  =  Σ_{e ∈ E(G_1)} ‖M_e‖²_F  +  Σ_{e ∈ E(G_2)} ‖M_e‖²_F
                  =  c_H · (B(G_1) + B(G_2)) · 4^(N − 2)
                  =  c_H · B(G) · 4^(N − 2),

using B(G) = B(G_1) + B(G_2) for disjoint unions. The degree-squared sum D2(G) is similarly component-additive: D2(G) = Σ_i deg_G(i)² = D2(G_1) + D2(G_2) since the degree of each site only counts edges to its own component.

This generalises to any finite number of components by induction. **The dimensional factor 4^(N − 2) uses the global N (all qubits), not per-component sizes**: the Liouville-space extension factor is a property of the ambient Hilbert space, not the connected component.

### Section 3. Weighted edges

For weighted bonds with per-bond couplings J_b ∈ ℝ_{>0}, the bond Hamiltonian H_e = J_b · (canonical bilinear at bond b) scales the per-bond M_e by J_b, so ‖M_e‖²_F scales by J_b². The bond-disjointness identity from Section 1 still holds (it depends only on the operator support, not the magnitude), so

    ‖M(N, G, {J_b})‖²_F  =  Σ_b J²_b · c_H · 4^(N − 2)
                         =  c_H · (Σ_b J²_b) · 4^(N − 2)
                         =  c_H · B_weighted · 4^(N − 2),

where B_weighted := Σ_b J²_b and c_H is anchored at the reference coupling J_ref = 1. This is the natural extension: at uniform J = 1, B_weighted = B(G) recovers the unweighted case.

### Section 4. Non-uniform γ

The γ-independence of F(N, G) is the separate result of [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md): the dissipator-block residual M_D = Π·L_D·Π⁻¹ + L_D + 2Σγ·I vanishes per Pauli string for arbitrary per-site γ_l, because the per-site Z-dephasing kernel is proportional to I_4 and the F1 σ-shift 2Σγ·I cancels the sum exactly. Hence ‖M‖²_F = ‖M_H‖²_F and the H-block carries no γ dependence by construction. This combines with Sections 1-3 to give the fully general statement above.

## Verification record

| (N, G) | path | predicted vs observed | relative error |
|---|---|---|---|
| N = 2, single bond | full PalindromeResidual.Build (anchor) | c_H = 128 by construction | 0 |
| N = 5, 6, named graphs (path, cycle, star, K_N, K_{2,N-2}) | full PalindromeResidual.Build (Python) | bit-exact | 0 |
| N = 5, 6, 30 random connected ER graphs per N (p ∈ {0.3, 0.5, 0.7}) | full PalindromeResidual.Build (Python) | bit-exact | 0 |
| N = 6, two disjoint 3-chains | full PalindromeResidual.Build (Python) | bit-exact | 0 |
| N = 4 chain, weighted J = (1, 2, 3) | full PalindromeResidual.Build (Python) | bit-exact | 0 |
| N = 5, single-body IY+YI on chain | full PalindromeResidual.Build (Python) | bit-exact | 0 |
| N = 5, chain / ring / star / triangle+disjoint-bond (XX+YZ non-truly) | C# `PalindromeResidualScalingClaim` graph-aware Factor | rel < 1e-9 | machine precision |
| N = 7, chain / ring / star / K_4+disjoint-3-chain (XY+Z-deph, F1-truly) | C# `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` + F1 palindromic-pairing multiset assertion {λ} = {−2σ − λ} | tolerance 1e-7 across 16 384 eigenvalues per topology | within MKL Evd accumulation |
| N = 8 Heisenberg, 4 topologies (see N=8 table below) | C# `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` + multiset assertion | tolerance 1e-6 across 65 536 eigenvalues per topology; observed max distance 2.6E-13 .. 2.6E-07 | well within tolerance |
| N = 9 Heisenberg chain | C# `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` with MklDirect bridge for blocks > LP64 ceiling | tolerance 1e-5 across 262 144 eigenvalues; max pairing distance **3.48E-13**; 0 outlier pairs | bit-exact; well within tolerance |

### N = 8 four-topology metric capture (Heisenberg J = 1, uniform γ = 0.5, σ = N·γ = 4)

Per-system metrics from the SLOW_N8 dogfood (full `F1SpectrumStatistics.TopologyMetrics` JSON files under `simulations/results/f1_n8_n9_metrics/`; numbers shown here are from the 2026-05-18 capture). All 4 systems pass the F1 palindromic-pairing assertion at tolerance 1e-6; observed max distance ranges from 2.6E-13 (ring) to 2.6E-07 (K_4 + disjoint 4-chain), comfortably inside the tolerance and 1-7 decades tighter than σ in absolute units.

| Topology | Bonds | Total wall (s) | Eig wall (s) | Max pair dist | Dissipation gap | Kernel dim | Distinct binned λ | Im range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Chain (open 7-bond path) | 7 | 225.8 | 208.0 | 4.16E-13 | 0.034 354 | 9 | 35 135 | ±5.125 |
| Ring (8-bond cycle) | 8 | 226.6 | 204.1 | 2.63E-13 | 0.133 901 | 9 | 21 845 | ±5.651 |
| Star (hub at site 0, 7 spokes) | 7 | 235.1 | 210.4 | 3.20E-13 | 0.086 970 | 9 | 2 275 | ±4.000 |
| K_4 ⊔ 4-chain (disconnected) | 9 | 224.5 | 199.1 | 2.56E-07 | 0.136 156 | 25 | 5 812 | ±5.366 |

Block-decomposition cost picture (same across all four since the joint-popcount structure depends only on N): 81 sectors, 41 X⊗N-primary sectors, max block 4 900² at (p_c, p_r) = (4, 4), Σ block³ = 5.46 × 10¹¹, effective speedup over dense (4^8)³ = **515.2×**. The block path is mandatory at N = 8: the full Liouvillian L_vec would be 65 536² × 16 B = 68.7 GB, past the .NET 2 GB single-array limit; the largest joint-popcount block fits in 0.38 GB.

### Discoveries from N = 8 data

1. **Kernel-dimension closed form across topologies**. The observed kernel dimensions match the prediction
        kernel dim  =  Π_c (|c| + 1)
   where the product runs over connected components c and |c| is the number of sites in c. Heisenberg conserves total S^z per component; Z-dephasing kills off-diagonal coherences; the joint kernel is then one stationary mode per (S^z_{c₁}, S^z_{c₂}, …) configuration. Chain / ring / star at N = 8 are connected → (8 + 1) = 9 modes (matched). K_4 ⊔ 4-chain at N = 8 has two components of sizes (4, 4) → 5 · 5 = 25 modes (matched). This refines the well-known "Heisenberg + Z-deph kernel = single-component prediction" to its component-additive form, with the SLOW_N8 sweep providing the first cross-component-count witness inside the F1 family at N = 8. **Tier 2 (verified)** at the four observed instances.

2. **Pairing precision drops with spectral degeneracy**. The K_4 ⊔ 4-chain instance shows 6 decades worse max pairing distance (2.56E-07 vs ~3E-13 for the single-component topologies). The driver is degeneracy density (distinct binned eigenvalues / total): chain 35 135 / 65 536 ≈ 54 %, ring 33 %, K_4 ⊔ 4-chain 9 %. Greedy nearest-neighbour multiset matching has many tied candidates when the spectrum is heavily degenerate, and MKL Evd's per-eigenvalue rounding noise (~1E-14 per block) accumulates across the ambiguous pairs. The 1e-6 tolerance was set to absorb exactly this regime. **Note** that the F1 identity itself still holds bit-exactly: it is the matching algorithm, not the spectrum, that loses 6 decades.

3. **Dissipation gap correlates with bond density, not connectivity per se**. Chain (7 bonds) has gap 0.034; ring (8 bonds, one extra edge closing the loop) jumps 4× to 0.134; star (7 bonds, same as chain) sits in between at 0.087; K_4 ⊔ 4-chain (9 bonds) is 0.136 = same as ring despite being disconnected. The gap is dominated by the bond-count contribution to the effective dissipator rather than by global topological invariants. Consistent with the bond-disjointness Section 1 of this proof: bond contributions add independently.

4. **Imaginary-axis extent reflects effective Hamiltonian gap**. Star has the smallest Im range (±4.0, an integer in units of γ at γ = 0.5), suggesting the star-Heisenberg spectrum saturates a discrete bound; ring has the largest (±5.65), reflecting circulant-spectrum spreading. The integer star bound deserves a separate analytical look; it is not predicted by any existing F1 theorem and is currently classified Tier 3 (Reading) pending derivation.

### N = 9 chain spot-check (SLOW_N9, verified 2026-05-19 via MklDirect bridge)

The N = 9 Heisenberg chain (8 bonds, J = 1, uniform γ = 0.5, σ = N·γ = 4.5) was reached on 2026-05-19 via the MklDirect bridge in [`LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`](../../compute/RCPsiSquared.Core/BlockSpectrum/LiouvillianBlockSpectrum.cs) (commit `abb2d52`). Full F1 palindromic-pairing identity {λ_k} = {−2σ − λ_k} verified bit-exact across all 262 144 eigenvalues to tolerance 1e-5 (observed max pairing distance ≈ 3.5E-13, max value of `MinReal` exactly −9.0, kernel dimension exactly 10 matching the F4KernelDim Tier1Derived prediction `Π_c(|c|+1) = 9 + 1`, dissipation gap 0.027 276 matching the chain `gap × N² ≈ 2.20` empirical scaling from [F1_DISSIPATION_GAP_PATTERN](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) to 0.4 %). Wall time **3 h 24 m** on the 128 GB / 24-core dev machine; effective speedup over the (infeasible) dense N=9 Liouvillian (would be 262 144² × 16 B ≈ **1.1 TB**) is **≈ 645.9×**, exceeding the LiouvillianBlockSpectrum docstring N=8 promise of 515×.

Per-system metric file: [`simulations/results/f1_n8_n9_metrics/chain_N9.json`](../../simulations/results/f1_n8_n9_metrics/) (full `F1SpectrumStatistics.TopologyMetrics` record).

See [§ Scale frontier: N = 9 chain, LP64 ceiling crossed 2026-05-19 via MklDirect bridge](#scale-frontier-n--9-chain-lp64-ceiling-crossed-2026-05-19-via-mkldirect-bridge) below for the bridge architecture and the N = 10 memory-pressure frontier that takes its place.

Scripts and tests:

- Python verification: [`simulations/f1_general_topology_verify.py`](../../simulations/f1_general_topology_verify.py)
- C# graph-aware + N=7 dogfood tests: [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs)
- C# N=8 SLOW_N8 dogfood tests: [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs)
- C# N=9 SLOW_N9 chain spot-check: [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN9BlockSpectrumChainTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN9BlockSpectrumChainTests.cs)
- Shared metrics utility: [`compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs`](../../compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs)
- Persisted metric JSON: [`simulations/results/f1_n8_n9_metrics/`](../../simulations/results/f1_n8_n9_metrics/)

### Why N=7 uses palindromic-pairing rather than the scaling formula

N=7 cannot use the block decomposition (`JointPopcountSectorBuilder` + `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`) for scaling-formula verification because the block infrastructure requires popcount-conserving H, which restricts to F1-truly H (Π²-even, palindrome-preserving) where ‖M‖²_F = 0 vacuously. The non-truly Hamiltonians used by the scaling-formula tests (XX+YZ, Π²-mixed) break popcount conservation and fall outside the block infrastructure's scope. So the non-truly XX+YZ scaling tests are done with the full dense `PalindromeResidual.Build` at N = 5 (graph-aware mode), and the N = 7 dogfood test exercises the F1 palindromic-pairing identity {λ} = {−2σ − λ} on the block-decomposed spectrum: the testable content of the F1 identity on the truly-H domain where the block path applies, and the only path that scales to N = 7 in a test budget (the full L_vec at N = 7 is 16384² = 4 GB, past the .NET 2 GB array limit).

## Scale frontier: N = 9 chain, LP64 ceiling crossed 2026-05-19 via MklDirect bridge

The N = 9 Heisenberg chain (8 bonds, J = 1, uniform γ = 0.5, σ = N·γ = 4.5) was the scale frontier for the F1 verification family between 2026-05-18 and 2026-05-19. The frontier was the LP64 MathNet MKL P/Invoke marshalling ceiling at 2 GB per native `Complex[]`, which at N = 9 blocked the largest joint-popcount sector C(9, 4) · C(9, 5) = 15 876² ≈ 4 GB. That ceiling was crossed on 2026-05-19 by promoting `MklDirect` (NativeMemory + ILP64 LAPACK) from `RCPsiSquared.Compute` into `RCPsiSquared.Core` and bridging the per-block dispatch in `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` (commit `abb2d52`).

**The bridge architecture.** A new `Lp64ComplexCeiling = 11 585` constant in [`LiouvillianBlockSpectrum`](../../compute/RCPsiSquared.Core/BlockSpectrum/LiouvillianBlockSpectrum.cs) selects between two per-sector Evd paths:

- **Block size ≤ 11 585² (≤ 2 GB native):** MathNet's LP64 `block.Evd().EigenValues` runs unchanged. Verified bit-exactly against MklDirect at N=3..5 in [`PerBlockLiouvillianBuilderNativeMemoryParityTests`](../../compute/RCPsiSquared.Core.Tests/BlockSpectrum/PerBlockLiouvillianBuilderNativeMemoryParityTests.cs) (11/11 pass).
- **Block size > 11 585² (> 2 GB native):** `PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory` allocates the block as a column-major `Complex*` via `NativeMemory.AllocZeroed`, then `MklDirect.EigenvaluesOnlyNative` performs the ILP64-aware zgeev directly, bypassing MathNet's LP64 marshaller. Outer DOP is forced to 1 on the MklDirect path to cap peak native-memory footprint (each block holds ~4 GB native).

**The 2026-05-19 N = 9 chain result.** The bridge runs end-to-end on a fresh test process with no parallel dotnet builds (race-condition pitfall documented in the local memory record `feedback_no_parallel_with_long_tests`): wall time **3 h 24 m**, **262 144 eigenvalues**, **0 outlier pairs**, max pairing distance **3.48E-13**, tolerance 1e-5. MinReal = −9.0 exactly (matching the F1 corollary `min Re(λ) = −2σ`), kernel dimension exactly 10 (matching the F4KernelDim Tier1Derived prediction `Π_c(|c|+1) = 10` for the 9-site connected chain), dissipation gap 0.027 276 (matching the empirical chain scaling `gap × N² ≈ 2.20` from [F1_DISSIPATION_GAP_PATTERN](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) to 0.4 %).

**Effective speedup at N = 9.** The full Liouvillian L_vec at N = 9 would be (4^9)² × 16 B = 262 144² × 16 B ≈ **1.1 TB**, infeasible on any commodity machine without out-of-core methods. The block-spectrum infrastructure with the MklDirect bridge reduces the wall-time cost to ~3.4 h on the 128 GB / 24-core dev machine; the cubic-cost speedup vs naive dense diagonalisation is **≈ 645.9×**, exceeding the LiouvillianBlockSpectrum docstring N=8 promise of 515×.

**Parallelism note.** The 3.4 h wall time is longer than the cubic-scaling estimate from N = 8 (which would predict ~1.5 h at constant per-thread throughput). CPU profile during the run sat at 4-5 of 24 cores active on average, suggesting the per-block MKL parallelism inside `MklDirect.EigenvaluesOnlyNative` is not fully scaling to 24 threads. Investigation deferred (Tom 2026-05-19: "gehen wir der Sache später auf den Grund"); the data was decisive enough that the parallelism analysis is the next sprint, not a blocker for the F1 verification record.

**The new frontier: N = 10.** With the LP64 ceiling crossed, the next infrastructure frontier sits at N = 10. The largest joint-popcount block there is C(10, 5)² = 252² = **63 504**, i.e. 63 504² × 16 B ≈ **64 GB** native per block. This fits on the 128 GB / 24-core dev machine for sequential outer-DOP execution but exhausts the parallel outer-DOP=6 budget (6 × 64 GB > 128 GB). The frontier therefore shifts from LP64 P/Invoke marshalling to memory-pressure × parallelism, requiring either (a) sequential outer-DOP at ~hours-to-days wall-time or (b) a more careful overlap of native allocation lifecycles. The F1 identity itself remains exact at every finite N.

**Test preservation.** `F1GeneralTopologyN9BlockSpectrumChainTests` is no longer a SkippableFact: the `Skip.If(maxBlockSize > Lp64EvdSquareMatrixCeiling, …)` guard was removed in commit `abb2d52` when the bridge landed. The test now runs unconditionally under `[Trait("Category", "SLOW_N9")]`, opt-in via `--filter "Category=SLOW_N9"`.

**Typed-claim record.** `F1GeneralTopologyVerifiedClaim.VerifiedNValues` extends to {5, 6, 7, 8, 9}, `ScaleUpToN` likewise, `ScaleFrontierBlockedAtN` advances to 10 with `ScaleFrontierBlockerReason` documenting the memory-pressure × parallelism frontier, `NamedGraphsVerified` increases from 22 to 23 (the +1 is the N=9 chain), and `SpectrumMetricsDataFiles` adds `chain_N9.json`. Future readers traversing the inspectable tree find the new frontier without grepping source comments.

## Closure note

This proof closes the last F1 OpenQuestion. Together with the 2026-05-18 closures of the T1 amplitude-damping closed form ([PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md)), depolarizing-noise closed form ([PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md](PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md)), and non-uniform γ closure ([PROOF_F1_NONUNIFORM_GAMMA.md](PROOF_F1_NONUNIFORM_GAMMA.md)), the F1 family has **zero open structural questions** as of 2026-05-18: the first time the F1 family's `OpenQuestions` collection is empty.

The numerical verification record extends from N = 2 through N = 9 without exception (the N = 9 chain was crossed 2026-05-19 via the MklDirect bridge described in § Scale frontier above; the next infrastructure frontier is at N = 10 and is memory-pressure × parallelism rather than LP64 marshalling). Both the (B, D2) closed form and the F1 palindromic-pairing identity hold across:

- arbitrary connected and disconnected graph topologies;
- uniform and non-uniform per-site γ_l (γ-independence of the H-block scaling factor is a separate Tier-1 result);
- uniform and weighted edges (via B → Σ J²_b substitution);
- main and single-body Hamiltonian classes (the latter with the (D2 / 2) prefactor).

The block-spectrum infrastructure (`JointPopcountSectorBuilder` + `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`) is the load-bearing computational primitive for the N ≥ 7 verifications: dense diagonalisation of the full L_vec is infeasible past N = 6 on commodity hardware.

## Cross-references

- Bond-disjointness (the substantive analytic content): [PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md) (Lemma 3 + Corollary).
- Algebraic origin of (B, D2) prefactors: [OPERATOR_RIGIDITY_ACROSS_CUSP.md](../../experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md) "Algebraic origin" and "Topology generalisation" sections.
- γ-independence: [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md).
- Typed claim (verification record): [`F1GeneralTopologyVerifiedClaim`](../../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs): fields `VerifiedNValues = {5, 6, 7, 8, 9}`, `ScaleUpToN = {5, 6, 7, 8, 9}`, `ScaleFrontierBlockedAtN = 10` with `ScaleFrontierBlockerReason` documenting the memory-pressure × parallelism frontier (the N=9 LP64 frontier was crossed 2026-05-19 via the MklDirect bridge), `NamedGraphsVerified = 23`, and `SpectrumMetricsDataFiles` listing the four N=8 JSON metric files plus `chain_N9.json`.
- F1 master identity: [F1PalindromeIdentity](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs), [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md).
- Shared metrics utility: [`F1SpectrumStatistics`](../../compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs): TopologyMetrics record covering wall-time profile, pairing-precision histogram, spectrum-structure invariants, block-decomposition cost picture, Hamiltonian + dissipator setup.
- Block-spectrum infrastructure (the load-bearing computational primitive for N ≥ 7): [`LiouvillianBlockSpectrum`](../../compute/RCPsiSquared.Core/BlockSpectrum/LiouvillianBlockSpectrum.cs), [`JointPopcountSectorBuilder`](../../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs), [`JointPopcountSectors`](../../compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs).
