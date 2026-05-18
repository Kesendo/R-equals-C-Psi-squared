# PROOF: F4 kernel-dimension closed form across connected components

**Status:** Tier 1 candidate. The structural argument (popcount conservation + tensor-sum factorisation of the Heisenberg Liouvillian L_H across connected components) is sound but no fully formal write-up exists; promotion to Tier 1 derived requires the analytic step in § "Open analytic step" below. Four bit-exact empirical anchors at N = 8 (integer kernel dimension; no rounding tolerance), captured during the F1 SLOW_N8 sweep on 2026-05-18.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Statement

For a graph G with connected components {G_1, G_2, ..., G_k} where component G_c contains |c| qubits (so Σ_c |c| = N), the kernel dimension of the Heisenberg Liouvillian L_H = −i[H_Heisenberg, ·] at Σγ = 0 factorises as

    dim ker L_H(G)  =  Π_c (|c| + 1).

The product runs over connected components c of G. For a single connected component of size N (the chain / ring / star / K_N / arbitrary connected case) the formula collapses to N + 1, which matches the F4 popcount-sector count Σ_w 1 for w ∈ {0, 1, ..., N}.

## Empirical anchors (bit-exact at N = 8)

Captured by the F1 SLOW_N8 sweep on 2026-05-18 (commit 89f725e, four `simulations/results/f1_n8_n9_metrics/<topology>_N8.json` files; observed kernel dimension surfaced in the `KernelDimension` field of each JSON).

| Topology | Components | Component sizes | Predicted Π(|c|+1) | Observed |
|---|---|---|---|---|
| chain N=8       | 1 | {8}    | 9    | 9 |
| ring N=8        | 1 | {8}    | 9    | 9 |
| star N=8        | 1 | {8}    | 9    | 9 |
| K_4 + disjoint 4-chain N=8 | 2 | {4, 4} | 5·5 = 25 | 25 |

The first three rows verify the connected-component specialisation (kernel dim = N + 1); the fourth row is the new content the data revealed (multiplicative factorisation across disjoint components, observed 25-fold kernel where the per-component prediction is 5).

## Proof sketch

### Section 1. Connected case (popcount labelling)

The Heisenberg bond Hamiltonian on bond (i, j) is

    H_(i,j) = J · (X_i X_j + Y_i Y_j + Z_i Z_j) = 2J · (S⃗_i · S⃗_j)  −  (J/2) · I.

All three terms preserve the joint Z-popcount w(σ) := #{i : σ_i ↑} on any computational-basis state, hence the global popcount

    Ŵ := Σ_l (I − Z_l) / 2

commutes with every Heisenberg bond Hamiltonian, [H, Ŵ] = 0 for any subgraph. The Liouvillian L_H = −i[H, ·] therefore acts within each pair of popcount sectors (w_left, w_right) of the operator algebra. For the kernel of L_H restricted to the diagonal-popcount sector (w_left = w_right = w), there is at least one mode per allowed popcount value w ∈ {0, 1, ..., N} (the projector onto popcount w, which commutes with H by construction).

For a single connected component of size N the popcount labels exhaust the kernel of L_H: this is the F4-trivial sector content (J = N/2 fully-symmetric multiplet contributes 1 mode per popcount, see [F4StationaryModeCountPi2Inheritance.cs](../../compute/RCPsiSquared.Core/Symmetry/F4StationaryModeCountPi2Inheritance.cs) for the full Clebsch-Gordan decomposition where the popcount-only modes are the m(J = N/2) = 1 contribution and the remaining (2J+1)² − N − 1 dimensions of F4's Stat(N) live in non-popcount-projector kernel modes). The popcount labels alone give N + 1 distinct kernel modes per connected component.

### Section 2. Tensor-sum factorisation across components

For G = G_1 ⊔ G_2 with sites partitioned into the two component supports, every Heisenberg bond Hamiltonian H_e is supported on a single component (e ∈ E(G_1) XOR e ∈ E(G_2)). Hence the full Hamiltonian decomposes as a tensor sum

    H(G) = H(G_1) ⊗ I_{G_2}  +  I_{G_1} ⊗ H(G_2),

and the Liouvillian inherits the tensor-sum structure on the operator algebra A(G) = A(G_1) ⊗ A(G_2):

    L_H(G) = L_H(G_1) ⊗ I_{A(G_2)}  +  I_{A(G_1)} ⊗ L_H(G_2).

The kernel of a tensor-sum operator factorises as the tensor product of per-summand kernels (standard result; see e.g. any textbook on operator equations of Lyapunov / Sylvester type):

    ker L_H(G) = ker L_H(G_1)  ⊗  ker L_H(G_2).

Combining with Section 1, dim ker L_H(G_c) ≥ |c| + 1 for each component, and the product of per-component popcount-labelled kernel modes gives

    dim ker L_H(G) ≥ Π_c (|c| + 1).

The four N = 8 anchors verify the equality side: the popcount-labelled kernel modes saturate the L_H kernel for the connected pieces tested (chain {8}, ring {8}, star {8}, K_4 {4}, 4-chain {4}). Induction on the number of components extends the formula to any finite k.

### Section 3. Open analytic step

What is rigorously established by Sections 1-2 is the lower bound

    dim ker L_H(G) ≥ Π_c (|c| + 1),

via the popcount-projector basis on each component plus tensor-product factorisation of kernels under the tensor-sum decomposition of L_H. The matching upper bound (popcount-projector basis exhausts ker L_H(G_c) for any connected component c) is observed bit-exactly at every (N, G) tested in this and prior sweeps (chain / ring / star / K_4 + disjoint at N=5..8 in [F1GeneralTopologyVerifiedClaim.cs](../../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs)), but a fully formal proof of equality requires either (a) the SU(2) Clebsch-Gordan decomposition restricted to the (w_left, w_right) = (w, w) operator sector showing that the J = N/2 multiplet's m = 1 popcount-projector basis is the unique kernel basis, or (b) an explicit count of the joint-J = 0 cohomology of the Heisenberg L_H acting on a single connected component of arbitrary topology. Either path closes Tier 1 candidate → Tier 1 derived.

The Clebsch-Gordan route is the natural one given F4's existing structure (Σ_J m(J)·(2J+1)² already partitions the full stationary count and the J = N/2 contribution m = 1 isolates the symmetric kernel basis cleanly). The tensor-sum kernel factorisation in Section 2 is already standard and rigorous; the open piece is purely the connected-case upper bound.

## Cross-references

- Parent: [F4 stationary mode count](../ANALYTICAL_FORMULAS.md#f4-stationary-mode-count-tier-1-clebsch-gordan-decomposition) (Clebsch-Gordan closed form for the connected case, Tier 1 derived in [F4StationaryModeCountPi2Inheritance.cs](../../compute/RCPsiSquared.Core/Symmetry/F4StationaryModeCountPi2Inheritance.cs))
- Sister Tier-2 verification record: [F1GeneralTopologyVerifiedClaim](../../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs) (the SLOW_N8 sweep that produced the four anchor JSON files; the dim-ker numbers are recorded there per topology)
- Typed claim: [F4KernelDimensionByComponentsClaim](../../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 candidate; `Predict(componentSizes)` returns Π(|c|+1))
- Data anchors: `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json` (`KernelDimension` field of each)
- Related: [PROOF_F1_GENERAL_TOPOLOGY](PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) parameterisation of the F1 residual norm under the same disconnected / weighted graph extensions used by this proof's Section 2)
