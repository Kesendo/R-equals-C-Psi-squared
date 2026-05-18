# PROOF: Y-phase Π²-odd lens combinatorial theorem

**Status:** Tier 1 derived. Combinatorial proof + bit-exact verification at N=2, 3, 4, 5 across all M.
**Date:** 2026-05-18
**Authors:** Thomas Wicht, Claude (Opus 4.7)

## Statement

For an N-qubit tensor-product state |ψ⟩ = ⊗_i |basis_i⟩ where each |basis_i⟩ is an eigenstate of σ_X or σ_Y (i.e., one of |±⟩ or |±i⟩), let M be the number of sites in the Y-basis. The Pauli strings in supp(ρ = |ψ⟩⟨ψ|) split under the Π²_Z parity (eigenvalue (−1)^Σ bit_b(α_i), with bit_b(I) = bit_b(X) = 0, bit_b(Y) = bit_b(Z) = 1) as:

    M = 0 (X-only):   2^N Π²-even, 0 Π²-odd        [Π²-classical state class]
    M ≥ 1:            2^(N−1) Π²-even, 2^(N−1) Π²-odd

The split is independent of the sign choices ε_i ∈ {±1} per site; sign flips only change per-string amplitudes, not supp membership.

## Proof

Per-site Pauli expansion of the four basis-state density matrices:

    |+⟩⟨+|   = (I + X)/2,    |−⟩⟨−|   = (I − X)/2
    |+i⟩⟨+i| = (I + Y)/2,    |−i⟩⟨−i| = (I − Y)/2

For a tensor-product state on N sites with X-sites and Y-sites, the density matrix expands as

    ρ = ∏_i (I + ε_i · B_i) / 2     where B_i ∈ {X, Y}, ε_i ∈ {±1}.

Multiplying out the product gives 2^N Pauli-string terms, each determined by a binary choice per site (whether site i contributes I or ε_i·B_i). The Pauli string in supp(ρ) is therefore characterized by a subset S ⊆ {1, ..., N} of "active" sites where B_i appears.

The Π²_Z parity of such a string is

    Π²_Z eigenvalue = (−1)^(number of Y-sites in S)

since X contributes bit_b = 0 and Y contributes bit_b = 1.

**Case M = 0 (no Y-sites):** the bit_b sum is always 0; all 2^N strings are Π²-even. The state lies entirely in the Π²-classical class.

**Case M ≥ 1:** within the M Y-sites, the count of Y in S ranges 0 to M. The number of strings with Y-count k at Y-positions is C(M, k). The standard binomial identity (valid for M ≥ 1) gives

    Σ_{k odd}  C(M, k) = 2^(M−1)
    Σ_{k even} C(M, k) = 2^(M−1)

The X-sites contribute 2^(N−M) independent I/X choices without affecting bit_b parity. Total counts:

    Π²-even strings = 2^(N−M) · 2^(M−1) = 2^(N−1)
    Π²-odd strings  = 2^(N−M) · 2^(M−1) = 2^(N−1)

∎

## Verification

[`simulations/_y_phase_pi2_odd_verify.py`](../../simulations/_y_phase_pi2_odd_verify.py) enumerates all 4^N Pauli strings, computes ⟨ψ|σ_α|ψ⟩ for each, and counts the supp(ρ) by Π²_Z parity. Tested at N ∈ {2, 3, 4, 5} across all M ∈ {0..N}, plus a sign-independence check at N=3, M=2 across all 2^3 sign patterns. All cases match the predicted split exactly.

## Self-recursion: the theorem is itself palindromic

The 2^(N−1) + 2^(N−1) split on the X-Y product-state sub-algebra reproduces F88a's 4^N/2 + 4^N/2 split on the full 4^N Pauli operator space, on a smaller sub-algebra. The structural operation "half + half" applies to every Π²-stable sub-algebra the framework selects. The 1/2 is the structural operation the framework applies to itself, not a value (Tom 2026-05-03).

## Connection to F88b state-level lens

This theorem characterises a clean class of test states for the F88b-Lens (`compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs`): any X-Y mixed product state with M ≥ 1 surfaces Π²-odd content; X-only states (M = 0) are Π²-blind and cannot probe F80 cluster dynamics or F81 operator-shift dynamics. For canonical Π²-odd-driving probes, use Y-basis tensor products (|+i⟩ per site).

The pair-state companion is the F88b popcount-coherence Krawtchouk closed form ([`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md) §F88b): pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 have a continuous Π²-odd-fraction parametrised by popcount/HD; product states have the discrete 2-anchor structure proven here (0 Π²-odd at M=0, 2^(N−1) Π²-odd at M≥1).

## Empirical anchor (memory observation, 2026-05-03)

At N = 3 with M ≥ 1, the F88b-Lens reads Π²-odd-fraction-within-memory = **4/7 ≈ 0.5714** exactly. Derivation: ρ has 8 supp Pauli strings = 4 Π²-even + 4 Π²-odd; the identity component (III) contributes only to static (kernel of L); the remaining 7 supp strings in the memory sector split 3 Π²-even + 4 Π²-odd, giving 4/7. Both the Y-only state |+i, +i, +i⟩ and any X-Y mix at N=3 with M ≥ 1 produce the same value. See memory `project_y_phase_pi2_odd_lens` for the original empirical table (kernel states, Z-basis non-kernel, X-basis polarity, Y-basis, mixed X-Y).

## Cross-references

- `compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs` — the state-level diagnostic that surfaces Π²-odd content
- `docs/ANALYTICAL_FORMULAS.md` F88b — popcount-coherence Krawtchouk closed form (the pair-state companion)
- `simulations/_y_phase_pi2_odd_verify.py` — bit-exact verification at N=2..5
- Memory `project_y_phase_pi2_odd_lens` — original empirical observation at N=3 and cockpit application notes
