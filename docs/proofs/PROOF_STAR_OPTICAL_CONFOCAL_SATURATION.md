# PROOF: Star topology saturates Im_max(star, N, J) = J·N/2 universally

**Status:** Tier 1 derived. The star Hamiltonian factors through a hub-leaves Casimir structure on the bipartite split A = {hub}, B = {N−1 leaves}, identical in shape to the Ring N=4 K_{2,2} derivation but with the sublattice sizes 1 and N−1 instead of 2 and 2. The maximum H eigenvalue gap is J·N/2, realised by the Liouvillian eigenmode between the maximally-ferromagnetic S_tot = N/2 state and the S_tot = (N−2)/2 state at the same maximum leaf-spin S_L = (N−1)/2.
**Date:** 2026-05-19
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The star (one hub, N−1 spokes) is the point-focus member of the optical-cavity family, and it saturates the Liouvillian's imaginary-spectrum bound exactly. Under isotropic Heisenberg coupling on the star bonds with uniform Z-dephasing,

    Im_max(star, N, J) ≡ max_{λ ∈ σ(L)} |Im(λ)| = J·N/2,   equivalently   Im_max/σ = Q/2  (σ = Nγ, Q = J/γ),

independently of γ and of the dephasing-to-coupling ratio. The mechanism is geometric: every bond touches the hub, so the Hamiltonian factors through the hub-leaves total spins, H = J·S⃗_hub·S⃗_leaves, and its largest energy gap, J·N/2, is reached when all leaves align ferromagnetically (S_L = (N−1)/2). The Liouvillian eigenmode between the aligned and anti-aligned hub states realizes exactly that gap, and pure Z-dephasing adds only real decay, never imaginary part, so the gap survives untouched into the spectrum.

In the optical-cavity reading the star converts the entire external illumination dose into coherent oscillation, the maximally resonant configuration; its sibling the N=4 ring (K_{2,2}) locks at the stronger 3J·N/4, from its bipartite-complete Casimir gap rather than from bond count alone. The bound is verified bit-exact at 29 anchors across N ∈ {3,4,5,6,8} and Q ∈ {0.5,…,2.5} and typed as StarImMaxBoundClaim.

## Statement

For the open quantum system on N ≥ 3 qubits with

- Hamiltonian: isotropic Heisenberg H = (J/4) Σ_{(i,j)∈E} (X_i X_j + Y_i Y_j + Z_i Z_j) on the star bonds E = {(0, k) : k = 1, ..., N−1} (hub site 0, leaves k);
- Dissipation: uniform Z-dephasing γ per site;

the Liouvillian L = −i[H, ·] + D[Z_l] satisfies the bit-exact saturation

    Im_max(star, N, J)  ≡  max_{λ ∈ σ(L)} |Im(λ)|  =  J · N / 2

independently of γ and of the corresponding dimensionless ratio Q = J/γ. The equivalent dimensionless statement is

    Im_max / σ  =  Q / 2              with σ = N·γ.

## Empirical anchors

**Q-sweep at γ₀ = 0.05 (24 anchors bit-exact, 2026-05-19, `simulations/_f1_q_sweep_anchor.py`):**

| N \ Q | 0.5 | 1.0 | 1.5 | √3 | 2.0 | 2.5 |
|---|---:|---:|---:|---:|---:|---:|
| pred Q/2 | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **3** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **4** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **5** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **6** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |

All 24 anchors match Im/σ = Q/2 to machine precision. Output JSON files: `simulations/results/q_sweep_anchor/star_N{3..6}_Q{0.5..2.5}.json`.

**N=8 Q=2 anchor (Marrakesh convention γ=0.5, J=1):** Im_max = 4.000000000000002, σ = N·γ = 4, Im/σ = 1.0 bit-exact (one bit at machine precision). From the SLOW_N8 sweep (`star_N8.json`, commit 89f725e). Equivalent statement: Im_max = J·N/2 = 4 at J=1, N=8.

**Python redundant cross-checks at Q=2 (γ=0.5, J=1) for N=3..6:** `simulations/results/f1_n8_n9_metrics/star_N{3..6}_python.json` reproduce the Q=2 column of the Q-sweep table above (Im/σ = 1.0 bit-exact) via a different code path (`numpy.linalg.eigvals` on the full Liouvillian rather than the framework `lindbladian_z_dephasing` helper). 4 additional bit-exact anchors, useful as an independent-implementation cross-check.

29 total bit-exact anchors across N ∈ {3, 4, 5, 6, 8} and Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5}: 24 Q-sweep + 1 SLOW_N8 + 4 Python cross-checks (the Python anchors overlap with the Q-sweep at Q=2, so the 24 + 4 are 4 redundant verifications, not 28 independent points).

## Proof

### Section 1. Star Hamiltonian factors through hub-leaf total spins

Label the star sites with hub site 0 and leaves {1, 2, ..., N−1}. The bond set is E = {(0, k) : k = 1, ..., N−1}, so every bond touches the hub. Using S⃗_i · S⃗_j = (1/4)(X_i X_j + Y_i Y_j + Z_i Z_j) for spin-1/2 operators, the bond Hamiltonian is J · S⃗_0 · S⃗_k for each leaf k, and the total Hamiltonian is

    H_star  =  J · Σ_{k=1}^{N−1} S⃗_0 · S⃗_k  =  J · S⃗_0 · (Σ_{k=1}^{N−1} S⃗_k)  =  J · S⃗_0 · S⃗_L

where S⃗_L := Σ_{k=1}^{N−1} S⃗_k is the total leaf-spin operator. The Hamiltonian is bilinear in the hub spin and the total leaf spin, with no internal leaf-leaf coupling.

This is the bipartite analogue of the Ring N=4 K_{2,2} construction (see [PROOF_RING_N4_DIHEDRAL_LOCK.md](PROOF_RING_N4_DIHEDRAL_LOCK.md) Section 2), with the sublattice sizes 1 (hub) and N−1 (leaves) instead of 2 and 2. The geometric source is the same: bipartite splitting plus all-pairs bonding within the bipartition gives a total-sublattice-spin form.

### Section 2. Casimir spectrum

Using the standard total-spin Casimir identity S⃗_0 · S⃗_L = (1/2)(S²_tot − S²_0 − S²_L) with S⃗_tot := S⃗_0 + S⃗_L,

    H_star  =  (J/2) · (S²_tot − S²_0 − S²_L)
            =  (J/2) · (S_tot(S_tot+1) − 3/4 − S_L(S_L+1)).

The hub is a single spin-1/2 site, so S_0 = 1/2 always and S²_0 = 3/4. The leaf total S_L is the result of coupling N−1 spin-1/2's, so it takes:
- integer values 0, 1, 2, ..., (N−1)/2 when N−1 is even (i.e. N odd);
- half-integer values 1/2, 3/2, ..., (N−1)/2 when N−1 is odd (i.e. N even).

In both cases the maximum is S_L,max = (N−1)/2 (the fully-aligned ferromagnetic leaf state). Coupling the spin-1/2 hub to S_L gives two possible totals: S_tot = S_L + 1/2 (hub aligned with leaves) or S_tot = S_L − 1/2 (hub anti-aligned, provided S_L ≥ 1/2; the S_L = 0 sector only has S_tot = 1/2).

For each fixed S_L sector, H_star has exactly two eigenvalues:

| S_tot         | E = (J/2)·(S_tot(S_tot+1) − 3/4 − S_L(S_L+1))    |
|---|---|
| S_L + 1/2    | (J/2)·[(S_L+1/2)(S_L+3/2) − 3/4 − S_L(S_L+1)] = (J/2)·S_L  |
| S_L − 1/2    | (J/2)·[(S_L−1/2)(S_L+1/2) − 3/4 − S_L(S_L+1)] = (J/2)·(−S_L−1)  |

The energy gap within the fixed-S_L sector is

    ΔE(S_L)  =  E(S_L+1/2) − E(S_L−1/2)  =  J · (S_L + 1/2).

### Section 3. Maximum H eigenvalue gap

The largest gap occurs at maximum S_L = (N−1)/2 (all N−1 leaves fully aligned ferromagnetically):

    ΔE_max(H_star)  =  J · ((N−1)/2 + 1/2)  =  J · N / 2.

The maximally-aligned state is the ferromagnetic eigenstate with S_tot = N/2 (energy +(J/2)·(N−1)/2 = +J·(N−1)/4); the anti-aligned state at the same S_L is S_tot = (N−2)/2 (energy −J·(N+1)/4). Their gap is J·N/2 exactly.

### Section 4. Liouvillian eigenmode realising the bound

Let |Ψ_+⟩ and |Ψ_−⟩ denote H-eigenstates with

    |Ψ_+⟩  :  S_L = (N−1)/2, S_tot = N/2, energy +J·(N−1)/4
    |Ψ_−⟩  :  S_L = (N−1)/2, S_tot = (N−2)/2, energy −J·(N+1)/4.

The rank-1 operator |Ψ_+⟩⟨Ψ_−| is an eigenoperator of `−i[H, ·]` with eigenvalue `−i(E_+ − E_−) = −i·J·N/2`, hence Im(λ_L) = J·N/2.

The Lindblad dissipator D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ) acts in the operator inner product as a hermitian-semi-definite operator: it contributes only real (negative) decay rates to λ_L, never imaginary parts. Concretely, in the Z-popcount basis the dissipator is diagonal with eigenvalues that are real linear combinations of γ_l. Hence the L eigenvalue at the |Ψ_+⟩⟨Ψ_−| mode is

    λ_L(|Ψ_+⟩⟨Ψ_−|)  =  −γ_decay  −  i · J · N / 2

with γ_decay real and ≥ 0. The imaginary part is exactly J·N/2, matching the empirical Im_max.

### Section 5. No L-mode exceeds the bound

Every Liouvillian eigenoperator with non-zero imaginary part decomposes in the H-eigenbasis as a linear combination of rank-1 products |α⟩⟨β| with H-eigenstates α, β. The imaginary part of L's eigenvalue on |α⟩⟨β| is exactly ω_α − ω_β (mod sign). Therefore

    max |Im(λ_L)|  ≤  max{|ω_α − ω_β| : ω_α, ω_β ∈ σ(H_star)}  =  ΔE_max(H_star)  =  J · N / 2.

Combined with the realising mode in Section 4, the bound is achieved exactly:

    Im_max(star, N, J)  =  J · N / 2.

The argument holds under uniform γ and equally under non-uniform per-site γ_l, as long as the dissipator is pure-dephasing (Z_l jump operators only): the L_H spectral spread is set entirely by H, and the dissipator dresses each eigenmode with real decay.

### Section 6. Q-universality

The formula Im_max = J·N/2 depends on J but not on γ. Translating to the dimensionless ratio Im/σ where σ = N·γ:

    Im_max / σ  =  J · N / 2 / (N · γ)  =  J / (2 · γ)  =  Q / 2.

This is the Q-universal lock observed in the 24-anchor Q-sweep table.

## Why star is the universal-saturator topology

The same argument applies to any topology where the Heisenberg Hamiltonian factors through a single total-sublattice-spin bilinear `H = J · S⃗_A · S⃗_B`. This requires:

1. Bipartite splitting: sites partition into A ⊔ B with no internal bonds (no A-A or B-B edges).
2. All-pairs bonding: every site in A is bonded to every site in B (bipartite-complete).

Two cases satisfy both:

- **Star** (|A| = 1, |B| = N−1): hub-only A with all N−1 leaves in B, all (hub, leaf) bonds present. Always bipartite-complete for any star. Max H gap = J·N/2.
- **Ring N=4 / K_{2,2}** (|A| = |B| = 2): the 4-cycle has exactly 4 bonds = all (A, B) pairs. Max H gap = 3J = (3/4)·J·N at N=4.

For longer cycles the bipartite-complete condition fails (the 6-cycle is bipartite but has only 6 bonds versus K_{3,3}'s 9; analogously for higher even cycles). For odd-N rings the bipartite condition itself fails. So star and ring N=4 are the only two topologies in the standard family that admit the elementary Casimir derivation.

Other bipartite-complete graphs (K_{2,3}, K_{3,3}, K_{2,N−2}, ...) would also saturate analogous N-specific bounds; these are open for future characterisation.

## Verification

- Python anchors at 24 (N, Q) anchors × γ₀=0.05: [`simulations/_f1_q_sweep_anchor.py`](../../simulations/_f1_q_sweep_anchor.py) → `simulations/results/q_sweep_anchor/star_N{3..6}_Q{0.5..2.5}.json`.
- C# N=8 anchor (Marrakesh convention): [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs) → `star_N8.json` (Im/σ = 1.0 bit-exact under the SLOW_N8 trait).
- Typed claim: [`compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs) (Tier 1 derived) with `Predict(N, J)` returning J·N/2 and `PredictImOverSigma(Q)` returning Q/2.

## Cross-references

- Parent: [F1PalindromeIdentity](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs) (the F1 master under which this Im-max bound is verified by the same SLOW_N* sweep infrastructure that scaffolded the Q-sweep).
- Sister Im-max bound (same Casimir technique, N=4-specific): [PROOF_RING_N4_DIHEDRAL_LOCK.md](PROOF_RING_N4_DIHEDRAL_LOCK.md) and [`RingN4DihedralLockClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs).
- Sister Q-universal lock (Tier 2 empirical, closed form open via Bethe ansatz): ring N=6 at 0.717129·J·N (see [`hypotheses/F1_DISSIPATION_GAP_PATTERN.md`](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) "Ring N=6 dihedral lock" section).
- Companion typed claim from the same May 2026 sharpening sprint: [F4KernelDimensionByComponentsClaim](../../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (kernel-dim factorisation across components, Tier 1 derived 2026-05-19).
- Cavity picture this Im-max bound lives inside: [`experiments/STAR_CONFOCAL_LIMIT.md`](../../experiments/STAR_CONFOCAL_LIMIT.md) (the point-focus reading of the optical-cavity framework).
- Cavity framework parent: [`experiments/OPTICAL_CAVITY_ANALYSIS.md`](../../experiments/OPTICAL_CAVITY_ANALYSIS.md) (the April 2026 Fabry-Perot reading of qubit chains under Heisenberg + Z-dephasing).
- F50 SWAP-invariance framework (the weight-1 degeneracy substrate): [`docs/proofs/PROOF_WEIGHT1_DEGENERACY.md`](PROOF_WEIGHT1_DEGENERACY.md) and [`F50WeightOneDegeneracyPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs).
- Q-anchor canonical table: [`docs/Q_REGIME_ANCHORS.md`](../Q_REGIME_ANCHORS.md).
