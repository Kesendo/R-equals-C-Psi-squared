# PROOF: Y-phase Π²-odd lens combinatorial theorem

**Status:** Tier 1 derived. Combinatorial proof + bit-exact verification at N=2, 3, 4, 5 across all M.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Abstract

The Y-phase lens theorem characterizes which product states can probe the Π²-odd (memory) sector of the operator algebra. Take an N-qubit tensor product of X- and Y-eigenstates (each site one of |±⟩ or |±i⟩) and let M be the number of Y-basis sites. The Pauli strings in the support of ρ = |ψ⟩⟨ψ| split under Π²_Z parity by a clean combinatorial rule: an X-only state (M=0) is entirely Π²-even, while any state with M ≥ 1 splits exactly in half, 2^(N−1) even and 2^(N−1) odd, independent of all sign choices.

The consequence is practical. X-only states are blind to the odd sector that F80/F81 dynamics populate (and are F88b-lens-blind); a single Y-basis site is enough to turn a product state into a perfect 50/50 probe of Π²-parity, which is why Y-basis and mixed X-Y product states are the canonical test probes for the polarity cube. It is structurally parallel to F102 (Y-parity as an independent polarity axis), and it supplies the canonical probe states for testing F107's truly-sector y-parity-zero purity and the F108 Π²-even palindrome family.

## Statement

Consider an N-qubit tensor-product state

    |ψ⟩ = ⊗_{i=1}^{N} |basis_i⟩

where each |basis_i⟩ is an eigenstate of σ_X or σ_Y (i.e., one of |±⟩ or |±i⟩). Let

    M := number of sites in the Y-basis,    so M ∈ {0, 1, ..., N}.

The Pauli strings in supp(ρ = |ψ⟩⟨ψ|) split under the Π²_Z parity as:

    M = 0 (X-only):   2^N Π²-even, 0 Π²-odd        [Π²-classical state class]
    M ≥ 1:            2^(N−1) Π²-even, 2^(N−1) Π²-odd

The split is independent of the sign choices ε_i ∈ {±1} per site; sign flips only change per-string amplitudes, not supp membership.

## Conventions

- **σ_X, σ_Y, σ_Z** are the standard 2×2 Pauli matrices; basis-state eigenvectors: |+⟩, |−⟩ for σ_X (eigenvalues +1, −1); |+i⟩, |−i⟩ for σ_Y.
- **Π** is the F1 palindrome conjugation operator (see [PROOF F1 / MIRROR_SYMMETRY](MIRROR_SYMMETRY_PROOF.md) and [F1 in ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)). **Π²_Z** is its square under Z-dephasing, diagonal on every Pauli string σ_α with eigenvalue (−1)^(Σ_i bit_b(α_i)), where the per-letter parity is `bit_b(I) = bit_b(X) = 0` and `bit_b(Y) = bit_b(Z) = 1`. See [F88a in ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md#f88a) for the two-axis Klein decomposition this parity belongs to.
- **supp(ρ)** denotes the set of Pauli strings α with non-zero expectation `Tr(ρ σ_α)` ≠ 0; i.e., the Pauli-basis support of the density matrix. For an N-qubit state ρ, supp(ρ) is a subset of the 4^N Pauli strings.
- **L** is the Lindbladian for Heisenberg + uniform Z-dephasing on the chain; its kernel (the static sector) is the span of popcount-sector projectors `{P_n : n = 0, ..., N}`. Pauli strings outside the kernel form the dynamic ("memory") sector.

## Empirical anchor (motivation, observed 2026-05-03)

At N = 3 with M ≥ 1, the F88b-Lens reads **Π²-odd-fraction-within-memory = 4/7 ≈ 0.5714 exactly**, the same value for Y-only |+i, +i, +i⟩ and for any X-Y mix with at least one Y-site. The theorem below explains this exactness: ρ has 8 supp Pauli strings, the theorem splits them as 4 Π²-even + 4 Π²-odd; the identity component III is one of the 4 Π²-even strings and lies entirely in static (kernel of L), leaving 7 strings in the memory sector that split 3 Π²-even + 4 Π²-odd, giving Π²-odd/memory = 4/7.

## Proof

**Per-site Pauli expansion of the four basis-state density matrices:**

    |+⟩⟨+|   = (I + X)/2,    |−⟩⟨−|   = (I − X)/2
    |+i⟩⟨+i| = (I + Y)/2,    |−i⟩⟨−i| = (I − Y)/2

So for a tensor-product state on N sites with X-sites and Y-sites, the density matrix expands as

    ρ = ∏_{i=1}^{N} (I + ε_i · B_i) / 2     where B_i ∈ {X, Y}, ε_i ∈ {±1}.

Multiplying out the product gives 2^N Pauli-string terms, each determined by a binary choice per site (whether site i contributes the I-factor or the ε_i·B_i factor). The Pauli string in supp(ρ) is therefore characterised by a subset

    S ⊆ {1, ..., N}    "active sites" where B_i appears

and the bijection (subset S) ↔ (supp string) is exact: every supp string is reached by exactly one S, and every S yields a distinct string (since the B_i ∈ {X, Y} per site differ from I).

**Π²_Z parity of the string indexed by S:**

    Π²_Z eigenvalue = (−1)^(number of Y-sites in S)

because X contributes bit_b = 0 and Y contributes bit_b = 1.

**Case M = 0 (no Y-sites):** the bit_b sum is identically 0 for every S; all 2^N strings are Π²-even. The state lies entirely in the Π²-classical class.

**Case M ≥ 1:** within the M Y-sites, the count of Y in S ranges 0 to M. The number of strings with Y-count k at Y-positions is C(M, k). The standard binomial identity (valid for M ≥ 1, follows from (1+1)^M ± (1−1)^M)

    Σ_{k odd}  C(M, k) = 2^(M−1)
    Σ_{k even} C(M, k) = 2^(M−1)

splits the Y-active patterns half-odd and half-even. The X-sites contribute 2^(N−M) independent I/X choices, all with bit_b = 0, so they multiply both halves equally without affecting parity. Total counts:

    Π²-even strings = 2^(N−M) · 2^(M−1) = 2^(N−1)
    Π²-odd strings  = 2^(N−M) · 2^(M−1) = 2^(N−1)

The sign factors ε_i affect only the per-string amplitudes (sign of `Tr(ρ σ_α)`), not which subsets S yield a non-zero amplitude; hence supp(ρ) is sign-independent. ∎

## Verification

[`simulations/y_phase_pi2_odd_verify.py`](../../simulations/y_phase_pi2_odd_verify.py) enumerates all 4^N Pauli strings, computes ⟨ψ|σ_α|ψ⟩ for each, and counts the supp(ρ) by Π²_Z parity. Tested at N ∈ {2, 3, 4, 5} across all M ∈ {0..N}, plus a sign-independence check at N=3, M=2 across all 2^3 sign patterns. All 18 (N, M) cases plus the 8 sign-pattern cases match the predicted split exactly.

## Connection to F88b state-level lens

This theorem characterises a clean class of test states for the F88b-Lens ([`MemoryAxisRho.cs`](../../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs)): any X-Y mixed product state with M ≥ 1 surfaces Π²-odd content; X-only states (M = 0) are Π²-blind and cannot probe F80 cluster dynamics or F81 operator-shift dynamics. For canonical Π²-odd-driving probes, use Y-basis tensor products (|+i⟩ per site).

The **pair-state companion** is the F88b popcount-coherence Krawtchouk closed form ([`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md) §"F88b: popcount-coherence Π²-odd / memory closed form"): pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 have a continuous Π²-odd-fraction parametrised by popcount-pair and Hamming distance; product states have the discrete 2-anchor structure proven here (0 Π²-odd at M=0, 2^(N−1) Π²-odd at M≥1).

## Remark: self-similar half-split structure

The 2^(N−1) + 2^(N−1) split on the X-Y product-state sub-algebra reproduces F88a's 4^N/2 + 4^N/2 split on the full Pauli operator space, on a smaller sub-algebra; both invoke the same Π²_Z parity operating on different domains. Every Π²-stable sub-algebra inherits the "half + half" structure from the parent F88a Klein decomposition.

> *The 1/2 is the structural operation the framework applies to itself, not a value.* (Tom 2026-05-03)

## Cross-references

- [`compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs`](../../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs) — the state-level diagnostic that surfaces Π²-odd content
- [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md) F88b — popcount-coherence Krawtchouk closed form (the pair-state companion)
- [`simulations/y_phase_pi2_odd_verify.py`](../../simulations/y_phase_pi2_odd_verify.py) — bit-exact verification at N=2..5
- Memory `project_y_phase_pi2_odd_lens` — original empirical observation at N=3 and cockpit application notes (user memory, outside repo)
