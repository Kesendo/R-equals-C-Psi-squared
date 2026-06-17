# PROOF: Ring N=4 dihedral lock, Im_max(ring, N=4, J) = (3/4)·J·N = 3J

**Status:** Tier 1 derived. The 4-cycle is the bipartite-complete graph K_{2,2}; its isotropic Heisenberg Hamiltonian factors through total sublattice spins via SU(2) Casimir, yielding a 6-eigenvalue spectrum {−2J, −J, 0, 0, 0, +J} with max gap 3J. The Liouvillian eigenmode realising this gap is the standard transition between the maximally-antialigned (S_A=1, S_B=1, S_tot=0) ground and the maximally-aligned (S_tot=2) ferromagnet, and the Z-dephasing dissipator only adds real decay so the gap survives into Im(λ).
**Date:** 2026-05-19
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Distinct from:** [`PROOF_RING_GAP_DOMINANCE.md`](PROOF_RING_GAP_DOMINANCE.md), the *XY* ring gap-dominance result (max|Im| = 2J = J·ρ, the dihedral lock). Same words "ring N=4 dihedral lock", different Hamiltonian and different result; this proof is the *isotropic-Heisenberg* ring (max|Im| = 3J via the K_{2,2} Casimir gap).

## Abstract

The 4-cycle is the one ring size that coincides with the bipartite-complete graph K_{2,2}: two sublattices of two qubits, every inter-sublattice bond present, none within. On it the isotropic Heisenberg Hamiltonian factors through total sublattice spins, H = J·S⃗_A·S⃗_B, and the SU(2) Casimir spectrum {−2J, −J, 0, 0, 0, +J} has maximal gap 3J. Under uniform Z-dephasing the Liouvillian's largest imaginary eigenvalue is pinned exactly to that gap,

    Im_max(ring, N=4, J) = (3/4)·J·N = 3J,   equivalently   Im/σ = 3Q/4,

independently of γ and Q = J/γ. The dissipator is pure dephasing, adding only real decay, so no eigenmode can exceed the Hamiltonian's spread: the bound is saturated, not approached.

This is a finite-N exact, Q-universal lock with an exactly rational coefficient, the dihedral-symmetry sibling of the star's point-focus bound Im_max = J·N/2 (3/4 vs 1/2, traced to the bipartite-complete K_{2,2} Casimir gap 3J versus the star's hub-spoke gap, not to bond count alone). What does not carry to larger even rings is the rational closed form: they keep a Q-universal lock but at transcendental constants (ring N=6 at 0.7171·J·N, descending toward ln 2), because the bipartite-complete structure is special to N=4. Typed as RingN4DihedralLockClaim.

## Statement

For the open quantum system on N=4 qubits with

- Hamiltonian: isotropic Heisenberg H = (J/4) Σ_{(i,j)∈E} (X_i X_j + Y_i Y_j + Z_i Z_j) on the 4-cycle bonds E = {(0,1), (1,2), (2,3), (3,0)};
- Dissipation: uniform Z-dephasing γ per site;

the Liouvillian L = −i[H, ·] + D[Z_l] satisfies the bit-exact saturation

    Im_max(ring, N=4, J)  ≡  max_{λ ∈ σ(L)} |Im(λ)|  =  (3/4) · J · N  =  3 · J         (for N=4)

independently of γ and of the corresponding dimensionless ratio Q = J/γ. The equivalent dimensionless statement is

    Im_max / σ  =  3Q/4              with σ = N·γ.

## Empirical anchors (bit-exact at 6 Q-values × γ₀=0.05)

Q-sweep on 2026-05-19 (`simulations/_f1_q_sweep_anchor.py`, output under `simulations/results/q_sweep_anchor/ring_N4_Q*.json`):

| Q | predicted 3Q/4 | observed Im/σ | rel. error |
|---:|---:|---:|---:|
| 0.5    | 0.375000 | 0.375000 | 1.5e-16 |
| 1.0    | 0.750000 | 0.750000 | 1.0e-15 |
| 1.5    | 1.125000 | 1.125000 | 7.9e-16 |
| √3     | 1.299038 | 1.299038 | 5.1e-16 |
| 2.0    | 1.500000 | 1.500000 | 3.3e-15 |
| 2.5    | 1.875000 | 1.875000 | 5.1e-15 |

All six anchors hit the prediction to within machine precision (relative error < 1e-14). The lock is Q-universal: the absolute Im_max value scales with J, but the dimensionless `Im/σ = 3Q/4` ratio is exact at every Q.

## Proof

### Section 1. The 4-cycle is the bipartite-complete graph K_{2,2}

The ring on N=4 sites has bond set E = {(0,1), (1,2), (2,3), (3,0)} (the 4-cycle). Partition the sites into sublattices A = {0, 2} and B = {1, 3} (two-colouring of the bipartite 4-cycle). The bond set then reads E = {(0,1), (0,3), (2,1), (2,3)}: every pair (a, b) with a ∈ A and b ∈ B is a bond, and no within-sublattice bonds exist. This is exactly the bipartite-complete graph K_{2,2} on 2+2 sites.

The 4-cycle C_4 and the bipartite-complete K_{2,2} are isomorphic as graphs. This is a coincidence specific to N=4: for N=6 the 6-cycle C_6 is bipartite (sublattices {0,2,4} and {1,3,5}) but has only 6 bonds versus K_{3,3}'s 9, so the bipartite-complete decomposition does not apply.

### Section 2. Total-sublattice-spin factorisation

Let S⃗_i = (S^x_i, S^y_i, S^z_i) be the spin-1/2 operator on site i. Using S⃗_i · S⃗_j = (1/4)(X_i X_j + Y_i Y_j + Z_i Z_j), the bond Hamiltonian is J · S⃗_i · S⃗_j and the total Hamiltonian on K_{2,2} is

    H  =  J · Σ_{a ∈ A, b ∈ B} S⃗_a · S⃗_b
       =  J · (S⃗_0 + S⃗_2) · (S⃗_1 + S⃗_3)
       =  J · S⃗_A · S⃗_B

with total sublattice spins S⃗_A := S⃗_0 + S⃗_2 and S⃗_B := S⃗_1 + S⃗_3. The Hamiltonian is bilinear in the sublattice total spins, with no within-sublattice term.

### Section 3. Casimir spectrum

Using the standard total-spin Casimir identity S⃗_A · S⃗_B = (1/2)(S²_tot − S²_A − S²_B) with S⃗_tot := S⃗_A + S⃗_B,

    H  =  (J/2) · (S²_tot − S²_A − S²_B).

Each sublattice contains two spin-1/2 sites, so S_A, S_B ∈ {0, 1}. Each (S_A, S_B) pair has its own internal Clebsch-Gordan multiplicity m_inner from coupling 1/2 ⊗ 1/2 → S_A (and 1/2 ⊗ 1/2 → S_B): m_inner = 1 for the singlet S = 0, m_inner = 1 for the triplet S = 1 (there is exactly one way to make each from two spin-1/2). Then coupling S⃗_A and S⃗_B gives S_tot ∈ {|S_A − S_B|, ..., S_A + S_B}, with the (2S_tot + 1) M_tot states inside each S_tot multiplet. The full eigenvalue list and the dimensions are:

| (S_A, S_B) | S_tot | E = (J/2)·(S_tot(S_tot+1) − S_A(S_A+1) − S_B(S_B+1)) | dimension = m_inner(A) · m_inner(B) · (2S_tot+1) |
|---|---:|---:|---:|
| (0, 0) | 0 | (J/2)·(0 − 0 − 0) = **0**     | 1 · 1 · 1 = 1 |
| (0, 1) | 1 | (J/2)·(2 − 0 − 2) = **0**     | 1 · 1 · 3 = 3 |
| (1, 0) | 1 | (J/2)·(2 − 2 − 0) = **0**     | 1 · 1 · 3 = 3 |
| (1, 1) | 0 | (J/2)·(0 − 2 − 2) = **−2J**   | 1 · 1 · 1 = 1 (outer-CG singlet of two inner triplets) |
| (1, 1) | 1 | (J/2)·(2 − 2 − 2) = **−J**    | 1 · 1 · 3 = 3 (outer-CG triplet) |
| (1, 1) | 2 | (J/2)·(6 − 2 − 2) = **+J**    | 1 · 1 · 5 = 5 (outer-CG quintuplet) |

The inner-CG multiplicities m_inner(A) and m_inner(B) are both 1 because two spin-1/2's couple into each of S = 0, S = 1 in exactly one way. The dimensions tracked are then just (2S_tot + 1) per row. Total Hilbert space dimension: 1 + 3 + 3 + 1 + 3 + 5 = 16 = 2^4 ✓.

The eigenvalue multiset of H is therefore

    σ(H)  =  { −2J, −J·(triplet, mult 3), 0·(seven-fold), +J·(quintuplet, mult 5) }.

Maximum eigenvalue is +J (the ferromagnetic S_tot = 2 multiplet); minimum is −2J (the perfect singlet of two anti-aligned triplet dimers). Maximum H eigenvalue gap is therefore

    ΔE_max(H_K22)  =  E_max − E_min  =  J − (−2J)  =  3J  =  (3/4) · J · N         for N=4.

### Section 4. Liouvillian eigenmode realising the bound

The Lindblad Liouvillian L = −i[H, ·] + D where D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ) is the pure-dephasing dissipator. For any pair of H-eigenstates |α⟩, |β⟩ with eigenvalues ω_α, ω_β, the rank-1 operator |α⟩⟨β| is an eigenoperator of `−i[H, ·]` with eigenvalue `−i(ω_α − ω_β)` (so Im(λ_L) = −(ω_α − ω_β)). The dissipator D acts diagonally in the joint-popcount basis and contributes only real decay (it commutes with the popcount projection and is hermitian-semi-definite in the operator inner product), hence the dephasing only adds Re(λ_L) corrections, never Im(λ_L).

In particular, the rank-1 operator |Ψ_+⟩⟨Ψ_−| with

    |Ψ_+⟩  ∈  the S_tot = 2 ferromagnetic multiplet (E = +J),
    |Ψ_−⟩  ∈  the (S_A=1, S_B=1, S_tot=0) singlet (E = −2J),

is a Liouvillian eigenoperator with eigenvalue λ = γ_decay − i · (J − (−2J)) = γ_decay − 3i·J. The imaginary part is exactly 3J, matching the empirical Im_max = (3/4)·J·N.

### Section 5. No mode exceeds the bound

Every Liouvillian eigenoperator with non-zero imaginary part is a linear combination of products |α⟩⟨β| with H-eigenstates α, β and Im(λ_L) entries in the set {−(ω_α − ω_β) : ω_α, ω_β ∈ σ(H)}. The maximum |Im(λ_L)| over the L-spectrum is therefore bounded by

    max |Im(λ_L)|  ≤  max{|ω_α − ω_β| : ω_α, ω_β ∈ σ(H)}  =  ΔE_max(H).

For K_{2,2} this is 3J. Combined with the realising mode in Section 4, the bound is achieved exactly:

    Im_max(ring, N=4, J)  =  ΔE_max(H_K22)  =  3J  =  (3/4) · J · N.

The same bound argument applies under non-uniform γ_l per site, as long as the dissipator is pure-dephasing (Z_l jump operators only). The L_H spectral spread is set entirely by H, and the dissipator dresses it with real decay.

### Section 6. Q-universality

The formula Im_max = (3/4)·J·N depends on J but not on γ. Translating into the dimensionless ratio Im/σ where σ = N·γ:

    Im_max / σ  =  (3/4) · J · N / (N · γ)  =  (3/4) · (J/γ)  =  (3/4) · Q.

This is the Q-universal lock observed in the Q-sweep table (Section "Empirical anchors").

## Why this is N=4-specific

The bipartite-complete structure C_4 = K_{2,2} relies on the 4-cycle having exactly 4 bonds (one per (A, B) pair). For longer cycles:

- 6-cycle has 6 bonds, K_{3,3} has 9: a 6-cycle is bipartite but NOT bipartite-complete.
- 8-cycle has 8 bonds, K_{4,4} has 16: same story.

For odd N (3-cycle, 5-cycle, ...) the ring is not even bipartite. So the K_{2,2} = C_4 coincidence is unique to N=4. Ring N=6, ring N=8 etc. show Q-universal locks too (empirically, [`hypotheses/F1_DISSIPATION_GAP_PATTERN.md`](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md): ring N=6 = 0.717129·J·N at 6 Q-anchors), but the per-N closed form requires Bethe ansatz on the cyclic dispersion rather than the simple Casimir argument that closes N=4.

## The N → ∞ limit: c_∞ = ln 2 (resolved 2026-06-04)

Although the per-N value needs Bethe ansatz, the LIMIT is closed. Section 5 reduces the lock to the Hamiltonian alone, Im_max(L) = ΔE_max(H) = E_max − E_min (the dephasing adds only real decay), so the dimensionless constant is

    c_N  ≡  Im_max / (J·N)  =  (E_max − E_min) / (J·N)  =  1/4 − E₀(N)/(J·N),

with E_max = J·N/4 (the ferromagnet, exact) and E₀(N) the antiferromagnetic Heisenberg-ring ground state. The per-bond ground energy of the spin-½ Heisenberg ring has the Bethe/Hulthén thermodynamic limit E₀/(J·N) → 1/4 − ln 2, hence

    c_∞  =  1/4 − (1/4 − ln 2)  =  ln 2  =  0.693147…     (NOT 1/√2 = 0.707107).

Computing E₀(N) directly ([`simulations/ring_dihedral_lock_limit.py`](../../simulations/ring_dihedral_lock_limit.py), sparse ground state, N = 4..16) confirms it:

| N | c_N = 1/4 − E₀/N | c_N − ln 2 |
|---:|---:|---:|
| 4 | 0.75000 | +0.05685 |
| 6 | 0.71713 | +0.02398 |
| 8 | 0.70639 | +0.01324 |
| 10 | 0.70154 | +0.00840 |
| 12 | 0.69895 | +0.00580 |
| 14 | 0.69740 | +0.00425 |
| 16 | 0.69639 | +0.00325 |

The N=6 value 0.71713 reproduces the F1_DISSIPATION_GAP empirical 0.717129, validating the Im_max = ΔE_max(H) reduction. The sequence decreases monotonically to ln 2 with a clean 1/N² finite-size approach (c_N − ln 2 quarters as N doubles: 0.01324 at N=8 → 0.00325 at N=16). It crosses 1/√2 at N=8; 1/√2 is only a value the sequence passes through, not the limit (the same red-herring lesson as the birth-canal s* = 0.709). So the Q-universal ring dihedral lock, left open for general N, has the exact limit c_∞ = ln 2.

## Relationship to the star saturation

[STAR_CONFOCAL_LIMIT.md](../../experiments/STAR_CONFOCAL_LIMIT.md) shows the analogous result for the star topology:

    Im_max(star, N, J)  =  J·N/2

via a parallel hub-spoke Casimir construction H = J · S⃗_0 · S⃗_L with S⃗_L = Σ leaf spins. The Casimir spectrum has max gap J·N/2 (the maximally-ferromagnetic-leaves S_L = (N-1)/2 sector flipping the hub). The ring N=4 result is structurally the same kind of object: a topology where the SU(2)-invariant Heisenberg Hamiltonian factors through a sum-of-Casimirs form, giving a closed Im-max spectral bound.

The ratio between the two is `(3/4)·J·N / (J·N/2) = 3/2`: ring N=4 carries 50% more imaginary spread than star N=4. The geometric reason is the larger bond count (4 in K_{2,2} vs N−1 = 3 in the N=4 star) combined with the bipartite-complete structure that maximises the inter-sublattice Casimir gap.

## Verification

- Python anchor at 6 Q-values × γ₀=0.05: [`simulations/_f1_q_sweep_anchor.py`](../../simulations/_f1_q_sweep_anchor.py) → `simulations/results/q_sweep_anchor/ring_N4_Q*.json`.
- Typed claim: [`compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs) (Tier 1 derived) with `Predict(J)` returning `(3/4) · J · N = 3J` at N=4.

## Cross-references

- Parent: [F1PalindromeIdentity](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs) (the F1 master under which this Im-max bound is verified by the same SLOW_N* sweep infrastructure).
- Sister Im-max bound (open for formal Tier 1 derivation in the same style): [STAR_CONFOCAL_LIMIT.md](../../experiments/STAR_CONFOCAL_LIMIT.md).
- Sister Q-universal lock (Tier 2 empirical, closed form open): ring N=6 at 0.717129·J·N (see [F1_DISSIPATION_GAP_PATTERN.md](../../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) "Ring N=6 dihedral lock" section).
- Companion typed claim from the same May 2026 sharpening sprint: [F4KernelDimensionByComponentsClaim](../../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (kernel-dim factorisation across components, Tier 1 derived 2026-05-19).
- Q-anchor canonical table: [`docs/Q_REGIME_ANCHORS.md`](../Q_REGIME_ANCHORS.md).
