# PROOF: F4 kernel-dimension closed form across connected components

**Status:** Tier 1 derived. Connected-case upper bound `dim ker ≤ N+1` closed by [`DEGENERACY_PALINDROME.md`](../../experiments/DEGENERACY_PALINDROME.md) Result 2 (magnetization conservation: the only operators satisfying both `[H, Q] = 0` and `D(Q) = 0` are the N+1 diagonal functions of total S_z). Multi-component factorisation `Π_c(|c|+1)` follows from per-component closure plus tensor-sum factorisation of L_H across disjoint components. Four bit-exact empirical anchors at N=8 (integer kernel dimension; no rounding tolerance) and 12 more at N=3..6 across chain/ring/star via Python-anchor scripts captured 2026-05-19.
**Date:** 2026-05-18 (Tier 1 candidate); promoted Tier 1 derived 2026-05-19.
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The dephased Heisenberg Liouvillian (Heisenberg coupling under uniform Z-dephasing, at any Σγ > 0) has a steady-state kernel: the operator subspace left untouched by the dynamics. For a connected chain of N qubits the kernel dimension is N + 1, matching the number of total-spin-z eigenvalues. The dephasing is essential here: at Σγ = 0 the bare commutator kernel is the full commutant of H, which is strictly larger (10 at N = 2, 24 at N = 3 for the chain); the Z-dephasing restricts the kernel to the diagonal functions of total S_z, leaving N + 1. The question this proof answers is what happens when the graph has multiple connected components.

The answer is a clean factorization. If the graph splits into components of sizes |c_1|, |c_2|, ..., |c_k|, the kernel dimension factorizes as the product of (|c_i| + 1) across components. A single component of size N reproduces the chain answer N + 1; two components of sizes 4 and 4 give (4 + 1) · (4 + 1) = 25; the multiplicative structure is the same as if each component lived in its own independent Heisenberg dynamics, which it effectively does because L_H tensor-sums across disjoint components.

The per-component statement (each connected component contributes a factor |c| + 1) is the deeper result. It comes from a magnetization-conservation argument: the operators left untouched by the dephased dynamics on a connected component (those that both commute with H and are annihilated by the Z-dephasing dissipator) are exactly the diagonal functions of the total z-spin in that component, and there are |c| + 1 such functions (one per total-z eigenvalue from −|c|/2 to +|c|/2 in unit steps). Disjoint components have independent total-z conservations, so their kernels multiply.

The diagnostic upshot is that the F1 kernel dimension is an immediate graph topology readout. Measuring the operator subspace fixed by the Heisenberg Hamiltonian tells you the component structure of the underlying connectivity, and the factor pattern in the kernel dimension reveals whether the system is one piece or several. F4 turns the kernel into an observable for graph reconstruction, which is occasionally a useful inverse problem.

## Statement

For a graph G with connected components {G_1, G_2, ..., G_k} where component G_c contains |c| qubits (so Σ_c |c| = N), the kernel dimension of the dephased Heisenberg Liouvillian L_H (here the full Liouvillian −i[H_Heisenberg, ·] + L_D under uniform Z-dephasing at any Σγ > 0, whose kernel is the joint kernel of the commutator and the dephasing dissipator, i.e. the steady-state space) factorises as

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

The matching upper bound `dim ker L_H(G_c) ≤ |c| + 1` for any single connected component follows from [DEGENERACY_PALINDROME.md](../../experiments/DEGENERACY_PALINDROME.md) Result 2, which proves via magnetization conservation that the only operators commuting with H_{G_c} and annihilated by the dephasing dissipator D are the identity plus the |c| popcount-sector projectors. The lower and upper bounds match, giving `dim ker L_H(G_c) = |c| + 1` for any connected component.

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

### Section 3. Upper-bound closure (resolved 2026-05-18)

What Sections 1-2 establish constructively is the lower bound

    dim ker L_H(G) ≥ Π_c (|c| + 1),

via the popcount-projector basis on each component plus tensor-product factorisation of kernels under the tensor-sum decomposition of L_H. The matching upper bound `dim ker L_H(G_c) ≤ |c| + 1` for any single connected component is closed by [DEGENERACY_PALINDROME.md](../../experiments/DEGENERACY_PALINDROME.md) Result 2 (magnetization conservation), which proves that the only operators satisfying both `[H, Q] = 0` (Hamiltonian conservation) and `D(Q) = 0` (annihilated by the Z-dephasing dissipator) are the N+1 diagonal functions of total S_z: the identity operator and the N projectors onto magnetization sectors with total S_z = m for m ∈ {−|c|/2, ..., +|c|/2}. Equivalently in our basis these are exactly the |c|+1 popcount projectors P_n constructed in Section 1. There is no further independent kernel element to exhaust, so the lower bound saturates and `dim ker L_H(G_c) = |c| + 1` exactly. Combined with Section 2's tensor-sum kernel factorisation across disjoint components, this gives the full equality `dim ker L_H(G) = Π_c (|c| + 1)`, completing the Tier 1 derived chain.

The connected-case upper bound was already proved (April 2026) inside DEGENERACY_PALINDROME's Result 2 derivation; the contribution of this proof file is the tensor-sum factorisation across multiple components (the new content surfaced by the four 2026-05-18 N=8 anchors with `K_4 + disjoint 4-chain`).

### Alternative proof routes

The Clebsch-Gordan route is the natural one given F4's existing structure (Σ_J m(J)·(2J+1)² already partitions the full stationary count and the J = N/2 contribution m = 1 isolates the symmetric kernel basis cleanly). The SU(2) Clebsch-Gordan decomposition restricted to the (w_left, w_right) = (w, w) operator sector would also show that the J = N/2 multiplet's m = 1 popcount-projector basis is the unique kernel basis. A second alternative is an explicit count of the joint-J = 0 cohomology of the Heisenberg L_H acting on a single connected component of arbitrary topology. Both routes converge on the same N+1 connected-case count and the magnetization-conservation route used above is the most economical; the alternatives are listed here for cross-checking richness, not because the Tier 1 derived chain depends on them.

## Cross-references

- Parent: [F4 stationary mode count](../ANALYTICAL_FORMULAS.md#f4-stationary-mode-count-tier-1-clebsch-gordan-decomposition) (Clebsch-Gordan closed form for the connected case, Tier 1 derived in [F4StationaryModeCountPi2Inheritance.cs](../../compute/RCPsiSquared.Core/Symmetry/F4StationaryModeCountPi2Inheritance.cs))
- Connected-case upper-bound closure: [experiments/DEGENERACY_PALINDROME.md](../../experiments/DEGENERACY_PALINDROME.md) Result 2 (magnetization conservation: `d_real(0) = N+1`, proven for connected Heisenberg components under Z-dephasing; this is the analytic anchor that saturates the Section 1 lower bound to equality)
- Sister Tier-2 verification record: [F1GeneralTopologyVerifiedClaim](../../compute/RCPsiSquared.Core/F1/F1GeneralTopologyVerifiedClaim.cs) (the SLOW_N8 sweep that produced the four anchor JSON files; the dim-ker numbers are recorded there per topology)
- Typed claim: [F4KernelDimensionByComponentsClaim](../../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (Tier 1 derived 2026-05-19; `Predict(componentSizes)` returns Π(|c|+1))
- Data anchors: `simulations/results/f1_n8_n9_metrics/{chain,ring,star,k4_plus_disjoint_4chain}_N8.json` (`KernelDimension` field of each)
- Corroborating per-weight ker breakdown: [docs/proofs/PROOF_WEIGHT1_DEGENERACY.md](PROOF_WEIGHT1_DEGENERACY.md) § Appendix 2026-05-17 (per-weight ker(w) decomposition across topologies; the w=0 row pins the kernel-projector count to N+1 for every connected graph tested at N=3..5, corroborating the boundary upper-bound used here)
- Related: [PROOF_F1_GENERAL_TOPOLOGY](PROOF_F1_GENERAL_TOPOLOGY.md) (the (B, D2) parameterisation of the F1 residual norm under the same disconnected / weighted graph extensions used by this proof's Section 2)

## The kernel, drawn live

![The live Liouvillian spectrum of an N = 5 chain at Q = 1.5 (Symphony export). The N + 1 = 6 frozen modes (Re λ = 0, the dephased kernel of this proof) sit on the right edge, highlighted; every other mode fades.](../../simulations/results/symphony_reel/without_t_axis_spectrum.png)

This proof's count, seen at a glance: the 6 modes on the Re = 0 edge are exactly N + 1 = 6, the dephased kernel for the connected N = 5 chain. The same figure is the shared anchor of the [F1 palindrome](MIRROR_SYMMETRY_PROOF.md) (the mirror about −σ) and the [absorption rungs](PROOF_ABSORPTION_THEOREM.md) Re λ = −2γ·n_XY. Exported by `inspect --root symphony --N 5 --J 0.075 --gamma 0.05 --export`, drawn by `simulations/reel_and_projector.py`.
