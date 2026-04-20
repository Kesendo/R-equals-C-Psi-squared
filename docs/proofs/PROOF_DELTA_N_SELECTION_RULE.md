# PROOF: Single-Site Partial Trace Annihilates |ΔN| ≥ 2 Sector Blocks

**Status:** Proven (elementary kinematic lemma, Tier 1)
**Date:** 2026-04-20
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Explains the empirical ΔN = 1 selection rule observed in the c₁ closure-breaking coefficient at N = 5 (verified for 8 coherence-block pairs, exact to machine precision). Any site-local observable on an N-qubit chain couples only to density-matrix sector blocks with |ΔN| ≤ 1; higher-ΔN coherences are invisible. The rule is purely kinematic: it holds for any Hamiltonian, any dissipator, any model that admits an excitation-number decomposition.

**Empirical data:** [sector_kernel.json](../../simulations/results/c1_sector_kernel/sector_kernel.json), [c1_bilinearity_test](../../simulations/results/c1_bilinearity_test/bilinearity_test.json).

---

## Statement

Let ρ be an operator on N qubits, decomposed into excitation-sector blocks ρ^(n, m) = P_n ρ P_m, where P_n is the projector onto the n-excitation computational subspace. Then for each site i:

    Tr_{¬i}(ρ^(n, m)) = 0   whenever   |n − m| ≥ 2.

The partial trace to a single site annihilates every sector block whose left and right popcounts differ by two or more.

## Proof

Let |x⟩ and |y⟩ be computational basis states on N qubits with popcounts |x| = n and |y| = m. Write x_i for the bit of x at site i, and x_{¬i} for the (N−1)-bit string with site i removed.

Partial trace to site i:

    Tr_{¬i}(|x⟩⟨y|) = ⟨x_{¬i} | y_{¬i}⟩ · |x_i⟩⟨y_i|.

The inner product ⟨x_{¬i} | y_{¬i}⟩ equals 1 when x and y agree at every site except possibly i, and 0 otherwise.

If x and y agree at every site except i, they differ in at most one bit position. Hence |n − m| = |x_i − y_i| ∈ {0, 1}.

Contrapositively, if |n − m| ≥ 2, then x and y differ in at least two bit positions, and no choice of i makes them agree at every other site. So ⟨x_{¬i} | y_{¬i}⟩ = 0 for all i, and Tr_{¬i}(|x⟩⟨y|) = 0.

Extending by linearity to ρ^(n, m) = Σ_{|x|=n, |y|=m} c_{xy} |x⟩⟨y|:

    Tr_{¬i}(ρ^(n, m)) = Σ_{|x|=n, |y|=m} c_{xy} · Tr_{¬i}(|x⟩⟨y|) = 0   for |n − m| ≥ 2.

∎

## Consequence for site-local observables

Any observable that factors through the per-site reduced states ρ_i = Tr_{¬i}(ρ) depends only on the |ΔN| ≤ 1 content of ρ.

- **Per-site expectation values** ⟨X_i⟩, ⟨Y_i⟩, ⟨Z_i⟩: linear in ρ_i, linear in the |ΔN| ≤ 1 content.
- **Per-site purity** Tr(ρ_i²): bilinear in ρ_i, bilinear in the |ΔN| ≤ 1 content.
- **Per-site fidelities** with a fixed 2×2 reference: linear in ρ_i.
- **Observer-time rescalings α_i** ([PTF](../../hypotheses/PERSPECTIVAL_TIME_FIELD.md)), derived by fitting per-site purity: functionals of the |ΔN| ≤ 1 content.
- **Closure-breaking coefficient c₁ = Σ_i f_i**, first derivative of Σ_i ln α_i at δJ = 0, likewise restricted to |ΔN| ≤ 1 bilinear contributions.

Sector blocks with |ΔN| ≥ 2 contribute zero to all of these.

## Generalisation to k-site observables

The partial trace to k sites annihilates |ΔN| ≥ k + 1 blocks:

    Tr_{¬{i_1, ..., i_k}}(ρ^(n, m)) = 0   whenever   |n − m| ≥ k + 1.

A k-local observable can see at most |ΔN| ≤ k sector blocks.

- **Single-qubit observables:** |ΔN| ≤ 1.
- **Pair observables** (two-site mutual information, concurrence, pair purity): |ΔN| ≤ 2.
- **Triple observables:** |ΔN| ≤ 3.
- **Full-state observables** (trace, global purity): |ΔN| ≤ N (unrestricted).

This predicts that moving from the site-local α_i picture to a pair-local α_{ij} picture would open the ΔN = 2 sector to contribution, and a richer closure structure would emerge. The c₁ story on per-site purity is specifically single-qubit-bound.

## Empirical confirmation

At N = 5, bond (0, 1), γ₀ = 0.05, δJ = ±0.01: coherence-block-only contribution to c₁, isolated by c₁(|S_n⟩+|S_m⟩ / √2) − c₁((|S_n⟩⟨S_n|+|S_m⟩⟨S_m|) / 2):

| ΔN | pairs tested | coherence contribution |
|----|-------------|------------------------|
| 1  | (0,1), (4,5) | +0.527 (reliable) |
| 2  | (0,2), (1,3), (2,4), (3,5) | 0.000 exactly (four pairs) |
| 3  | (0,3), (2,5) | 0.000 exactly (two pairs) |
| 4  | (0,4), (1,5) | 0.000 exactly (two pairs) |
| 5  | (0,5) | 0.000 trivially |

All nine |ΔN| ≥ 2 pairs tested give zero to machine precision (eight non-trivial at |ΔN| ∈ {2, 3, 4}, plus the trivial (0, 5) at |ΔN| = 5). The pairs with fit artefacts ((1, 2), (2, 3), (3, 4), (1, 4), all involving two middle-sector Dicke states and producing non-monotonic purity trajectories that the time-rescaling ansatz does not capture) are excluded from the reliable subset; the coherent-state c₁ for those pairs was captured earlier in the pure-superposition measurement and is nonzero only for |ΔN| = 1.

Data: [c1_sector_kernel/sector_kernel.json](../../simulations/results/c1_sector_kernel/sector_kernel.json).

## Scope

The lemma is purely kinematic. It holds for:

- **Any Hamiltonian** (XY, Heisenberg, XXZ, Ising, etc.) that conserves excitation number, so that the sector decomposition is well-defined.
- **Any dissipator** that preserves the sector structure (Z-dephasing, amplitude damping within sector, any Kraus channel commuting with the sector projectors).
- **Any initial state**.

It does NOT say which |ΔN| ≤ 1 contributions are nonzero; those are determined by the dynamics. The |ΔN| = 1 contributions to c₁ observed in [c1_sector_kernel](../../simulations/results/c1_sector_kernel/sector_kernel.json) are a dynamical fact, separate from this kinematic selection rule. The rule only says: the dynamics cannot generate nonzero contributions from |ΔN| ≥ 2 blocks, no matter how the Liouvillian is shaped.

## Relation to existing repo results

- **[XOR_SPACE](../../experiments/XOR_SPACE.md):** The N + 1 "center modes" at |n − m| = N (maximal Hamming-distance coherences) are invisible to all site-local observables. GHZ-type states living in those modes cannot be read from per-site measurements, consistent with the "GHZ is an antenna blind spot" observation.
- **[F61 n_XY parity selection rule](../ANALYTICAL_FORMULAS.md):** a different selection (about Pauli-string XY parity under dephasing). Related but orthogonal.
- **[Pi-pair closure investigation](../../simulations/results/pi_pair_closure_investigation/FINDINGS.md):** the kinematic |ΔN| ≤ 1 restriction plus Π-pair identity together explain the c₁ patterns observed at N = 3, 5, 6, 7.
