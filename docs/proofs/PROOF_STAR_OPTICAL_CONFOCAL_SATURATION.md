# PROOF: Im_max(star, N, J) = J·N/2 at uniform γ

**Status:** Tier 1 derived. The star Hamiltonian factors through a hub-leaves Casimir structure on the bipartite split A = {hub}, B = {N−1 leaves}, identical in shape to the Ring N=4 K_{2,2} derivation but with the sublattice sizes 1 and N−1 instead of 2 and 2. The maximum H eigenvalue gap is J·N/2, and it is realised in the Liouvillian spectrum by the coherences between a fully polarised computational state and the E_min states at each Hamming rung.
**Date:** 2026-05-19 (Sections 1 to 3 and the anchors), 2026-07-31 (Sections 4 to 6, the scope fences and the minimiser search)
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The star (one hub, N−1 spokes) has a closed-form Hamiltonian spectral spread, and among connected graphs on N ≤ 6 sites it is the unique smallest. Under isotropic Heisenberg coupling on the star bonds with uniform Z-dephasing,

    Im_max(star, N, J) ≡ max_{λ ∈ σ(L)} |Im(λ)| = J·N/2,   equivalently   Im_max/σ = Q/2  (σ = Nγ, Q = J/γ),

independently of γ and of the dephasing-to-coupling ratio. The mechanism is geometric: every bond touches the hub, so the Hamiltonian factors through the hub-leaves total spins, H = J·S_hub·S_leaves, and its largest energy gap, J·N/2, is reached when all leaves align ferromagnetically (S_L = (N−1)/2).

Three statements have to be kept apart. (i) ΔE_max(H_star) = J·N/2 is star-specific and exact at every N. (ii) The saturation `max|Im λ_L| = ΔE_max(H)` is **not** star-specific: it holds for every isotropic Heisenberg graph, connected or not, because the realising modes need only that a fully polarised state is an extreme of H, which Section 5b derives in general. (iii) What is star-specific besides the closed form is that J·N/2 is the **minimum** over connected graphs, searched exhaustively at N ≤ 6 with the star as the unique minimiser.

The bound is verified at 25 distinct anchors across N ∈ {3,4,5,6,8} and Q ∈ {0.5,…,2.5}, over 29 runs (four are Python re-runs of points the Q-sweep already covers), agreeing to a worst relative deviation of 1.98e-14, and typed as StarImMaxBoundClaim. The sibling N=4 ring (K_{2,2}) locks at the larger 3J·N/4, from its bipartite-complete Casimir gap rather than from bond count alone.

## Statement

For the open quantum system on N ≥ 3 qubits with

- Hamiltonian: isotropic Heisenberg H = (J/4) Σ_{(i,j)∈E} (X_i X_j + Y_i Y_j + Z_i Z_j) on the star bonds E = {(0, k) : k = 1, ..., N−1} (hub site 0, leaves k);
- Dissipation: uniform Z-dephasing γ per site;

the Liouvillian L = −i[H, ·] + D[Z_l] satisfies

    Im_max(star, N, J)  ≡  max_{λ ∈ σ(L)} |Im(λ)|  =  J · N / 2

independently of γ and of the corresponding dimensionless ratio Q = J/γ. The equivalent dimensionless statement is

    Im_max / σ  =  Q / 2              with σ = N·γ.

## Empirical anchors

**Q-sweep at γ₀ = 0.05 (24 anchors, 2026-05-19, `simulations/f1_q_sweep_anchor.py`):**

| N \ Q | 0.5 | 1.0 | 1.5 | √3 | 2.0 | 2.5 |
|---|---:|---:|---:|---:|---:|---:|
| pred Q/2 | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **3** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **4** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **5** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **6** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |

All 24 anchors match Im/σ = Q/2 to machine precision. Output JSON files: `simulations/results/q_sweep_anchor/star_N{3..6}_Q{0.5..2.5}.json`. Five of the six Q values are canonical anchors from [`docs/Q_REGIME_ANCHORS.md`](../Q_REGIME_ANCHORS.md); Q = 0.5 is a below-onset probe that is not in that table.

**N=8 Q=2 anchor (Marrakesh convention γ=0.5, J=1):** Im_max = 4.000000000000002, σ = N·γ = 4, so Im/σ = 1.0 to machine precision. From the SLOW_N8 sweep (`star_N8.json`, commit 89f725e). Equivalent statement: Im_max = J·N/2 = 4 at J=1, N=8.

**Python re-runs at Q=2 (γ=0.5, J=1) for N=3..6:** `simulations/results/f1_n8_n9_metrics/star_N{3..6}_python.json` reproduce the Q=2 column of the Q-sweep table above. These are re-runs of the same code, not an independent implementation: `f1_q_sweep_anchor.py` imports its `run` from `f1_topology_heisenberg_small_n_anchor.py`, which builds the Liouvillian with the framework helper `framework.lindblad.lindbladian_z_dephasing`. Both paths therefore call the same function, so these four runs confirm reproducibility, not cross-implementation agreement.

25 distinct (N, Q) anchors across N ∈ {3, 4, 5, 6, 8} and Q ∈ {0.5, 1.0, 1.5, √3, 2.0, 2.5}, over 29 runs: 24 Q-sweep + 1 SLOW_N8 + 4 Python re-runs of Q=2 points the sweep already covers.

**Precision.** No stored anchor equals J·N/2 exactly in floating point, so these agreements are machine-precision, not bit-exact. The worst relative deviation across all 29 stored runs is **1.98e-14** (largest at N=6, where the Liouvillian is 4096 × 4096).

## Proof

### Section 1. Star Hamiltonian factors through hub-leaf total spins

Label the star sites with hub site 0 and leaves {1, 2, ..., N−1}. The bond set is E = {(0, k) : k = 1, ..., N−1}, so every bond touches the hub. Vector operators carry no arrow here: S_i = (S^x_i, S^y_i, S^z_i) and the leaf total S_L are three-component operators and the dot is their scalar product, while the same letters read as quantum numbers wherever they carry no dot and no square, as in S_L = (N−1)/2 or S_tot(S_tot+1) from Section 2 on. Using S_i · S_j = (1/4)(X_i X_j + Y_i Y_j + Z_i Z_j) for spin-1/2 operators, the bond Hamiltonian is J · S_0 · S_k for each leaf k, and the total Hamiltonian is

    H_star  =  J · Σ_{k=1}^{N−1} S_0 · S_k  =  J · S_0 · (Σ_{k=1}^{N−1} S_k)  =  J · S_0 · S_L

where S_L := Σ_{k=1}^{N−1} S_k is the total leaf-spin operator. The Hamiltonian is bilinear in the hub spin and the total leaf spin, with no internal leaf-leaf coupling.

This is the bipartite analogue of the Ring N=4 K_{2,2} construction (see [PROOF_RING_N4_DIHEDRAL_LOCK.md](PROOF_RING_N4_DIHEDRAL_LOCK.md) Section 2), with the sublattice sizes 1 (hub) and N−1 (leaves) instead of 2 and 2. The geometric source is the same: bipartite splitting plus all-pairs bonding within the bipartition gives a total-sublattice-spin form.

### Section 2. Casimir spectrum

Using the standard total-spin Casimir identity S_0 · S_L = (1/2)(S²_tot − S²_0 − S²_L) with S_tot := S_0 + S_L,

    H_star  =  (J/2) · (S²_tot − S²_0 − S²_L)
            =  (J/2) · (S_tot(S_tot+1) − 3/4 − S_L(S_L+1)).

The hub is a single spin-1/2 site, so S_0 = 1/2 always and S²_0 = 3/4. The leaf total S_L is the result of coupling N−1 spin-1/2's, so it takes:
- integer values 0, 1, 2, ..., (N−1)/2 when N−1 is even (i.e. N odd);
- half-integer values 1/2, 3/2, ..., (N−1)/2 when N−1 is odd (i.e. N even).

In both cases the maximum is S_L,max = (N−1)/2 (the fully-aligned ferromagnetic leaf state). Coupling the spin-1/2 hub to S_L gives two possible totals: S_tot = S_L + 1/2 (hub aligned with leaves) or S_tot = S_L − 1/2 (hub anti-aligned), the second requiring S_L ≥ 1/2.

For each fixed S_L ≥ 1/2 sector, H_star has exactly two eigenvalues:

| S_tot         | E = (J/2)·(S_tot(S_tot+1) − 3/4 − S_L(S_L+1))    |
|---|---|
| S_L + 1/2    | (J/2)·[(S_L+1/2)(S_L+3/2) − 3/4 − S_L(S_L+1)] = (J/2)·S_L  |
| S_L − 1/2    | (J/2)·[(S_L−1/2)(S_L+1/2) − 3/4 − S_L(S_L+1)] = (J/2)·(−S_L−1)  |

The energy gap within the fixed-S_L sector is

    ΔE(S_L)  =  E(S_L+1/2) − E(S_L−1/2)  =  J · (S_L + 1/2).

The S_L = 0 sector, which exists whenever N is odd, is the exception: it admits only S_tot = 1/2 and so carries a single level, E = 0, with no gap defined. It is not the maximum and does not enter Section 3. (The measured star spectra confirm this: at N = 3, 5, 7 the value E = −J/2 that a two-levels-in-every-sector reading would predict is absent.)

### Section 3. Maximum H eigenvalue gap

The largest gap occurs at maximum S_L = (N−1)/2 (all N−1 leaves fully aligned ferromagnetically):

    ΔE_max(H_star)  =  J · ((N−1)/2 + 1/2)  =  J · N / 2.

The maximally-aligned state is the ferromagnetic eigenstate with S_tot = N/2 (energy +(J/2)·(N−1)/2 = +J·(N−1)/4); the anti-aligned state at the same S_L is S_tot = (N−2)/2 (energy −J·(N+1)/4). Their gap is J·N/2 exactly.

### Section 4. No L-mode exceeds the bound

Write L = K + D with K = −i[H, ·] and D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ), both as superoperators on the Hilbert-Schmidt space of operators.

Because each jump operator Z_l is Hermitian, D is **self-adjoint** with respect to the Hilbert-Schmidt inner product, and it is negative semidefinite. K is anti-Hermitian and normal, with spectrum {−i(ω_α − ω_β) : ω_α, ω_β ∈ σ(H)}.

Let v be an eigenvector of L with eigenvalue λ. Then λ = ⟨v, Lv⟩/‖v‖², and since ⟨v, Dv⟩ is real,

    Im(λ)  =  ⟨v, Kv⟩ / (i·‖v‖²)

which lies in the convex hull of spec(K)/i = {ω_β − ω_α}, i.e. in [−ΔE_max, +ΔE_max]. Hence

    max |Im(λ_L)|  ≤  max{|ω_α − ω_β| : ω_α, ω_β ∈ σ(H_star)}  =  ΔE_max(H_star)  =  J · N / 2.

This is a field-of-values argument and it holds for any γ_l ≥ 0, uniform or not.

**The hypothesis that carries it is the Hermiticity of the jump operators**, not the informal statement that the dissipator "only adds real decay". The informal version is false. Counterexample, verified in the gate: one qubit, H = 0 (so ΔE_max = 0), single jump operator c = I + iY. The generator is trace preserving, unital, has all Re λ ≤ 0 and no Hamiltonian part whatsoever, yet its spectrum is {0, 0, −2+2i, −2−2i}, so max|Im λ| = 2 > 0 = ΔE_max. What that dissipator lacks is self-adjointness: ‖D − D†‖ = 5.6569.

### Section 5. Liouvillian eigenmodes realising the bound

Equality in Section 4 forces v into the extremal eigenspace of K **and** forces v to be an eigenvector of D. That is a joint condition, and it is not automatic: the extremal H-coherence block is not D-invariant. A generic rank-1 product |Ψ_+⟩⟨Ψ_−| between the two extremal multiplets of Section 3 is therefore **not** an L-eigenoperator; the gate measures worst residuals of 0.94, 1.15 and 1.20 at N = 3, 4, 5 over all such products.

The realisers are a specific sub-family, and what selects them is the Z-eigenbasis, not the Casimir structure. Let

    |ferro⟩  :  a fully polarised computational state, |0…0⟩ or |1…1⟩
    |β_k⟩    :  the E_min eigenstate in the Hamming-weight-k sector, k = 1, …, N−1.

Two facts make this well posed. First, |ferro⟩ is the H **maximum**: each bond term J·S_i·S_j has largest eigenvalue +J/4 on the triplet, and |ferro⟩ is triplet on every bond simultaneously, so it attains all of them at once, at energy J·B/4 = J(N−1)/4 for the star's B = N−1 bonds. (With the convention J > 0 the fully polarised state is the top of the spectrum, not the ground state.) Second, |β_k⟩ exists at E_min for **every** k in 1..N−1: the E_min multiplet has S_tot = (N−2)/2 at S_L = (N−1)/2, and its S_z values run over −(N−2)/2 … +(N−2)/2, which relative to |0…0⟩ is exactly the Hamming rungs 1 … N−1. This is a star-specific fact and it is what fixes the mode count below; on the chain, ring and complete graph the global minimum sits only in the middle rung, so only that rung realises.

|ferro⟩ is a Z-product state, so Z_l|ferro⟩ = ±|ferro⟩ for every l, and the weight-k sector is an eigenspace of Σ_l Z_l with eigenvalue N − 2k. Under **uniform** γ the dissipator therefore acts on |β_k⟩⟨ferro| as the scalar

    D[ |β_k⟩⟨ferro| ]  =  γ·( (N − 2k) − N )·|β_k⟩⟨ferro|  =  −2γk · |β_k⟩⟨ferro|

and since |β_k⟩ sits at E_min while |ferro⟩ sits at E_max,

    L |β_k⟩⟨ferro|  =  ( −2γk  +  i · J·N/2 ) · |β_k⟩⟨ferro|

exactly. Combined with Section 4, the bound is achieved:

    Im_max(star, N, J)  =  J · N / 2.

Counting: two polarised extremes × (N−1) rungs × two signs of Im gives **4(N−1)** eigenvalues at |Im| = J·N/2 **counted with multiplicity**, measured as 8, 12, 16, 20 at N = 3, 4, 5, 6. The distinct values number 2(N−1), each of multiplicity 2.

**Uniform γ is required here.** Under site-dependent γ_l the dissipator acts on the weight-1 block as diag(−2γ_l) rather than as a scalar, no longer commutes with the one-magnon Hamiltonian, and the saturation fails. Measured (star, J = 1):

| N | γ profile | max\|Im\| | J·N/2 |
|---|---|---:|---:|
| 3 | (2.0, 0.5, 0.5) | 1.032100016 | 1.5 |
| 4 | (2.0, 0.5, 0.5, 0.5) | 1.592689809 | 2.0 |
| 4 | (0.1, 0.2, 0.3, 0.4) | 1.981222484 | 2.0 |
| 5 | (0.11, 0.29, 0.47, 0.65, 0.83) | 2.440907393 | 2.5 |

Always strictly below, never equal. The Section 4 bound survives non-uniform γ; the Section 5 saturation does not.

### Section 6. What the saturation is, and is not, a property of

The equality `max|Im(λ_L)| = ΔE_max(H)` is **not** special to the star, and Sections 4 and 5 already prove it in general. Neither used the star geometry. Restating the two steps without it:

- **Bound.** Section 4 used only that the jump operators are Hermitian. Any graph, any γ_l ≥ 0.
- **Attainment.** Section 5 needed |ferro⟩ to be an extreme of H, that it be a Z-product state, and that the opposite extreme lie in a single Hamming rung. For isotropic Heisenberg on **any** graph: |ferro⟩ attains the maximum +J/4 of every bond term simultaneously, so it is the global maximum at J·B/4; H commutes with Σ_l Z_l, so an E_min eigenvector may always be chosen inside one rung k; and then |β⟩⟨ferro| has definite Hamming weight k and the dissipator acts on it as the scalar −2γk.

So `max|Im λ_L| = ΔE_max(H)` for isotropic Heisenberg on any graph under uniform dephasing. Note the difference from the star case: in general only the rung containing the global minimum realises, whereas for the star every rung 1..N−1 does. Measured for chain, ring and complete at N = 3..6, an asymmetric six-vertex graph (|Aut| = 1) and a disconnected K₃ ⊔ K₃ at N = 6, agreeing to the last printed digit in every case.

What is star-specific is the **value**: ΔE_max(H_star) = J·N/2, and by exhaustive search the star is the unique minimiser of ΔE_max at N = 4, 5, 6 (38, 728 and 26704 connected labelled graphs searched; minimum exactly J·N/2; minimisers exactly the N stars). That minimality is searched, not proved, and is open past N = 6.

The hypothesis really in play is the polarised extreme, not "dephasing plus any H". Under the same uniform dephasing the XY chain does not saturate (max|Im| = 1.4398736 against spread 2.2360680 at N=4) and neither does the transverse-field Ising chain (3.7344319 against 4.1883993 at N=4). Whether the polarised extreme is also *necessary* is open.

### Section 7. The dimensionless form

The formula Im_max = J·N/2 depends on J but not on γ. Translating to the dimensionless ratio Im/σ where σ = N·γ:

    Im_max / σ  =  J · N / 2 / (N · γ)  =  J / (2 · γ)  =  Q / 2.

This is the Q-universal lock observed in the 24-anchor Q-sweep table.

## Which topologies admit the elementary Casimir derivation

This section is about which graphs give a *closed form* for ΔE_max by the argument of Sections 1 to 3. It is not about which graphs saturate: by Section 6, all of them do.

The Casimir derivation applies to any topology where the Heisenberg Hamiltonian factors through a single total-sublattice-spin bilinear `H = J · S_A · S_B`. This requires:

1. Bipartite splitting: sites partition into A ⊔ B with no internal bonds (no A-A or B-B edges).
2. All-pairs bonding: every site in A is bonded to every site in B (bipartite-complete).

Two cases satisfy both:

- **Star** (|A| = 1, |B| = N−1): hub-only A with all N−1 leaves in B, all (hub, leaf) bonds present. Always bipartite-complete for any star. Max H gap = J·N/2.
- **Ring N=4 / K_{2,2}** (|A| = |B| = 2): the 4-cycle has exactly 4 bonds = all (A, B) pairs. Max H gap = 3J = (3/4)·J·N at N=4.

For longer cycles the bipartite-complete condition fails (the 6-cycle is bipartite but has only 6 bonds versus K_{3,3}'s 9; analogously for higher even cycles). For odd-N rings the bipartite condition itself fails. So star and ring N=4 are the only two topologies in the standard family that admit the elementary Casimir derivation.

Other bipartite-complete graphs (K_{2,3}, K_{3,3}, K_{2,N−2}, ...) admit analogous N-specific closed forms; these are open for future characterisation.

Among all connected graphs, the star's closed form is the **minimum**. Exhaustive search over every connected graph on N vertices (38 at N=4, 728 at N=5, 26704 at N=6) gives a minimum ΔE_max of exactly J·N/2 in each case, attained by the N stars and by nothing else. Whether the star remains the unique minimiser at every N is open.

## Verification

- The gate: [`simulations/star_saturation_gate.py`](../../simulations/star_saturation_gate.py), 129 checks covering the Casimir closed form (N=3..7), the star law across (J, γ), the universality of the saturation across topologies, the realising modes and their count, all three scope fences (non-uniform γ, non-ferromagnetic H, non-Hermitian jumps), the minimiser search, the additive composition of disconnected components, and the stored-anchor precision.
- Python anchors at 24 (N, Q) anchors × γ₀=0.05: [`simulations/f1_q_sweep_anchor.py`](../../simulations/f1_q_sweep_anchor.py) → `simulations/results/q_sweep_anchor/star_N{3..6}_Q{0.5..2.5}.json`.
- C# N=8 anchor (Marrakesh convention): [`compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs`](../../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs) → `star_N8.json`. That test asserts the F1 palindromic-pairing identity; `MaxImag` is one of the statistics it logs, not something it asserts, so it anchors this claim as recorded data rather than as a test.
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
