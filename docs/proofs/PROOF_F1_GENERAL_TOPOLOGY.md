# PROOF: F1 H-block residual closed form extends to arbitrary graph topology

**Status:** Tier 1 (synthesis). Closes the F1 OpenQuestion "general topology beyond chain/ring/star/K_N" by combining the existing analytic result PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary (universality of the bond-disjointness argument across any graph) with disconnected-component and weighted-edge extensions. Verification record: Python at N = 5, 6 across named graphs, random connected Erdős-Rényi, disconnected, weighted, single-body; C# at N = 5 (graph-aware mode of `PalindromeResidualScalingClaim`) and N = 7 (F1 palindromic pairing via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` on multiple graph topologies including disconnected).
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

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

For the main class each ‖M_e‖²_F equals the per-bond constant c_H · 4^(N−2) (Lemma 1 of PROOF_CROSS_TERM_FORMULA, multiplied by the Liouville-space extension factor 4^(N−2)). The sum over E(G) gives the bond-count B(G) prefactor. For the single-body class, the per-site contributions get a degree-weighted assembly (each Iσ + σI bond term contributes to the site at both endpoints), giving the D2(G) / 2 prefactor (see OPERATOR_RIGIDITY_ACROSS_CUSP.md "Algebraic origin" for the derivation of why D2 enters the single-body case). Both prefactors are linear in graph invariants — no higher-order topological dependence (cycle count, triangle count, hub presence) ever appears in the residual norm closed form.

This establishes universality of (B, D2) parameterisation across any **connected** graph.

### Section 2. Disconnected components

For G = G_1 ⊔ G_2 a disjoint union of two connected components, the bond set E(G) = E(G_1) ⊔ E(G_2) and the bond-additive identity from Section 1 gives

    ‖M(N, G)‖²_F  =  Σ_{e ∈ E(G_1)} ‖M_e‖²_F  +  Σ_{e ∈ E(G_2)} ‖M_e‖²_F
                  =  c_H · (B(G_1) + B(G_2)) · 4^(N − 2)
                  =  c_H · B(G) · 4^(N − 2),

using B(G) = B(G_1) + B(G_2) for disjoint unions. The degree-squared sum D2(G) is similarly component-additive: D2(G) = Σ_i deg_G(i)² = D2(G_1) + D2(G_2) since the degree of each site only counts edges to its own component.

This generalises to any finite number of components by induction. **The dimensional factor 4^(N − 2) uses the global N (all qubits), not per-component sizes** — the Liouville-space extension factor is a property of the ambient Hilbert space, not the connected component.

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

Scripts and tests:

- Python verification: [`simulations/_f1_general_topology_verify.py`](../../simulations/_f1_general_topology_verify.py)
- C# graph-aware + N=7 dogfood tests: [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs)

## Closure note

This proof closes the last F1 OpenQuestion. Together with the 2026-05-18 closures of the T1 amplitude-damping closed form ([PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md](PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md)), depolarizing-noise closed form ([PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md](PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md)), and non-uniform γ closure ([PROOF_F1_NONUNIFORM_GAMMA.md](PROOF_F1_NONUNIFORM_GAMMA.md)), the F1 family has **zero open structural questions** as of 2026-05-18 — the first time the F1 family's `OpenQuestions` collection is empty.

## Cross-references

- Bond-disjointness (the substantive analytic content): [PROOF_CROSS_TERM_FORMULA](PROOF_CROSS_TERM_FORMULA.md) (Lemma 3 + Corollary).
- Algebraic origin of (B, D2) prefactors: [OPERATOR_RIGIDITY_ACROSS_CUSP.md](../../experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md) "Algebraic origin" and "Topology generalisation" sections.
- γ-independence: [PROOF_F1_NONUNIFORM_GAMMA](PROOF_F1_NONUNIFORM_GAMMA.md).
- Typed claim (verification record): [`F1GeneralTopologyVerifiedClaim`](../../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs).
- F1 master identity: [F1PalindromeIdentity](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs), [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md).
