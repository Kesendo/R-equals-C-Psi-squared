# PROOF: Chromaticity of Single-Step Coherence Blocks

**Status:** Proven (elementary combinatorial theorem, Tier 1)
**Date:** 2026-04-22
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Formalizes the chromaticity formula c(n, N) = min(n, N−1−n) + 1 observed in the (n, n+1) popcount coherence blocks of U(1)-conserved systems under uniform Z-dephasing. The formula labels the number of distinct pure dephasing rates present in each block at J = 0, and remains the structural sector label once H is turned on. Primary empirical ground: [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 3, verified at J = 0 for all tested (N, n) with N = 3..8.

**F-entry:** [F74 in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md).

---

## Statement

Let the N-qubit system have a Hamiltonian H with [H, N_total] = 0 (total excitation number conserved) and uniform Z-dephasing at rate γ₀ on all sites. Decompose the operator space into popcount blocks via the projectors P_n onto the n-excitation computational subspace. The (n, n+1) coherence block is spanned by basis elements |x⟩⟨y| with |x| = n, |y| = n+1.

At J = 0, this block contains exactly

    c(n, N) = min(n, N−1−n) + 1

distinct pure dephasing rates, given by

    α(HD) = 2γ₀ · HD,    HD ∈ {1, 3, 5, ..., 2c−1},

where HD is the Hamming distance between the basis states |x⟩ and |y⟩.

## Proof

The proof has two lemmas and a final combination.

### Lemma 1 (combinatorial). Distinct Hamming distances in an (n, n+1) block.

**Claim.** For basis states |x⟩, |y⟩ ∈ {0, 1}^N with |x| = n and |y| = n+1, the Hamming distance HD = |{i : x_i ≠ y_i}| takes values in the set

    {1, 3, 5, ..., 2·min(n, N−1−n) + 1},

yielding exactly c(n, N) = min(n, N−1−n) + 1 distinct values.

**Proof.** Partition the N sites by the (x_i, y_i) pattern into four types:

- **match:** both 1. Count = m.
- **type-a:** x_i = 1, y_i = 0. Count = a.
- **type-b:** x_i = 0, y_i = 1. Count = b.
- **zero:** both 0. Count = N − m − a − b.

The popcounts give |x| = m + a = n and |y| = m + b = n + 1, so b = a + 1.

Hamming distance HD = a + b = 2a + 1. Every HD value is odd.

Constraints on a:

- a ≥ 0 (trivially).
- b = a + 1 ≤ N − n, since type-b sites must come from the N − n sites with x_i = 0. So a ≤ N − n − 1.
- m = n − a ≥ 0, so a ≤ n.

Therefore a ∈ {0, 1, ..., min(n, N − n − 1)}, giving HD ∈ {1, 3, ..., 2·min(n, N − n − 1) + 1}.

The count of distinct HD values is min(n, N − n − 1) + 1 = min(n, N − 1 − n) + 1. ∎

### Lemma 2 (algebraic). Pauli representation of |x⟩⟨y| has ⟨n_XY⟩ = HD.

**Claim.** For any basis pair (|x⟩, |y⟩) with Hamming distance HD, every Pauli string with nonzero coefficient in the decomposition of |x⟩⟨y| has exactly HD factors from {X, Y} and N − HD factors from {I, Z}. Equivalently, ⟨n_XY⟩(|x⟩⟨y|) = HD.

**Proof.** The operator |x⟩⟨y| tensor-factorizes across sites:

    |x⟩⟨y| = ⊗_{i=1}^N |x_i⟩⟨y_i|.

Each single-qubit factor has a closed Pauli expansion:

    |0⟩⟨0| = (I + Z)/2,     |1⟩⟨1| = (I − Z)/2,
    |0⟩⟨1| = (X + iY)/2,    |1⟩⟨0| = (X − iY)/2.

Sites with x_i = y_i (both 0 or both 1) contribute a factor (I ± Z)/2, supported on {I, Z}. These are the N − HD "match" sites.

Sites with x_i ≠ y_i contribute (X ± iY)/2, supported on {X, Y}. These are the HD "differ" sites.

Expanding the tensor product over all sites: every Pauli string in the expansion of |x⟩⟨y| has X or Y on exactly the HD differ-sites, and I or Z on exactly the N − HD match-sites. Every nonzero term has n_XY = HD.

For an operator v with all Pauli coefficients living at a single n_XY value, the Absorption Theorem expectation

    ⟨n_XY⟩(v) = Σ_P |c_P|² · n_XY(P) / ||v||²

collapses to the single value HD. ∎

### Combining the lemmas.

At J = 0, the Liouvillian reduces to the dephasing dissipator L = L_D. In the Pauli basis, L_D acts diagonally:

    L_D(P) = −2γ₀ · n_XY(P) · P

for any Pauli string P, by the Absorption Theorem (AT) applied to single Pauli strings.

Let (|x⟩, |y⟩) be a basis pair in the (n, n+1) block with Hamming distance HD. By Lemma 2, every Pauli string in the decomposition of |x⟩⟨y| has n_XY = HD, so all terms share the eigenvalue −2γ₀·HD. By linearity, |x⟩⟨y| is a right eigenvector of L_D with eigenvalue −2γ₀·HD, and the decay rate is

    α(|x⟩⟨y|) = 2γ₀ · HD.

By Lemma 1, HD ranges over {1, 3, ..., 2c−1} within the block, where c = min(n, N−1−n) + 1. Each value is attained (pick a with 0 ≤ a ≤ c − 1 and set HD = 2a + 1). So the block contains exactly c distinct pure dephasing rates 2γ₀ · {1, 3, ..., 2c−1}. ∎

---

## Consequences

**C1. Monochromatic endpoint blocks (c = 1).** At n = 0 and n = N − 1, min(n, N−1−n) = 0 so c = 1. Only HD = 1 is accessible, and the single rate is 2γ₀. Every basis pair in these blocks decays at the same rate. [F73](../ANALYTICAL_FORMULAS.md) (spatial-sum coherence closure for vac-SE probes) is the c = 1, n = 0 case.

**C2. Maximal chromaticity at the center.**

- **Odd N:** Unique center block at n = (N − 1)/2, with c_max = (N + 1)/2.
- **Even N:** Two adjacent center blocks at n = N/2 − 1 and n = N/2, both with c_max = N/2.

At small N: c values by (N, n) are

| N | c values (n = 0 to N−1) |
|---|-------------------------|
| 3 | (1, 2, 1) |
| 4 | (1, 2, 2, 1) |
| 5 | (1, 2, 3, 2, 1) |
| 6 | (1, 2, 3, 3, 2, 1) |
| 7 | (1, 2, 3, 4, 3, 2, 1) |
| 8 | (1, 2, 3, 4, 4, 3, 2, 1) |
| 9 | (1, 2, 3, 4, 5, 4, 3, 2, 1) |

The pattern generalises: c_max grows as ⌈(N+1)/2⌉, reached at the center(s).

**C3. Chromaticity labels H-mixing bands.** When J > 0, the Hamiltonian H couples the c pure-rate channels within each block (H preserves popcount but changes n_XY by ±2 via XY-site hopping). Dressed eigenmodes appear at fractional ⟨n_XY⟩ values between the integer rungs {1, 3, ..., 2c−1}. The dressed-mode weight W(Q) measured in [Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) quantifies this H-mixing and shows that the observable peak location is c-specific and N-invariant (Q_peak(c) = 1.5, 1.6, 1.8 for c = 2, 3, 4). The c = 1 blocks have no inner channels to mix: W ≡ 0 for all Q (consistent with F73's J-independence).

## Scope

The theorem is purely combinatorial plus one application of the Absorption Theorem. The minimal assumptions are:

- **Hamiltonian:** any H with [H, N_total] = 0 (XY, XXZ, Heisenberg plus optional Z-field, any sector-preserving interaction). Only needed so that the (n, n+1) sector decomposition is dynamically invariant; at J = 0 the theorem holds purely kinematically.
- **Dissipator:** uniform Z-dephasing γ₀ on every site, so that L_D acts diagonally in the Pauli basis with eigenvalue −2γ₀ · n_XY.
- **Graph topology:** irrelevant. The counting in Lemma 1 and the Pauli decomposition in Lemma 2 use only the site label, not any connectivity structure. Holds for chain, ring, star, complete graph, arbitrary graph.
- **N:** any N ≥ 1. The formula is trivially satisfied for N = 1 (only c = 1 possible), becomes structurally informative for N ≥ 3.

## Breaks for

- **Non-uniform γ_i (site-dependent dephasing).** The dissipator eigenvalue on a basis pair with differ-site set D becomes −2 · Σ_{i ∈ D} γ_i rather than −2γ₀ · HD. Pairs with the same HD but different differ-site patterns then decay at different rates, breaking the clean c-label structure. The block still decomposes into HD-classes, but each HD-class further splits by differ-set, and the number of distinct rates can exceed c(n, N).
- **Non-U(1) Hamiltonians** (transverse fields X_i, Y_i, odd-popcount interactions). The (n, n+1) sector decomposition is no longer H-invariant, so the H-mixing picture in C3 breaks down. At J = 0 the theorem still holds (Lemma 1 and Lemma 2 are U(1)-independent), but the dynamical interpretation as a stable sector label requires U(1).
- **Dissipators with off-diagonal Pauli action** (amplitude damping, X-dephasing, depolarizing). L_D no longer acts with eigenvalue −2γ₀ · n_XY on single Pauli strings, so the basis pair |x⟩⟨y| is no longer an eigenvector of L_D at rate 2γ₀ · HD. Chromaticity remains a useful Pauli-counting label but loses its dynamical rate interpretation.
- **Coherences across popcount gaps larger than 1** (|ΔN| ≥ 2 blocks). Lemma 1 counts HD in the (n, n+1) block specifically. For (n, n+k) blocks the analogous counting gives HD ∈ {k, k+2, ..., min(2n+k, 2N−2n−k)}, with min(n, N−n−k) + 1 distinct values. This generalises the formula but the site-local invisibility result from [F70](../ANALYTICAL_FORMULAS.md) applies: |ΔN| ≥ 2 blocks contribute zero to any single-site observable.

## Empirical verification

**Block structure c-values for N = 3..8:** exact match to the combinatorial prediction, verified by enumerating popcount sectors and counting distinct HD values in each (n, n+1) block ([Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 3).

**J = 0 spectrum verification:** for each tested (N, n) block at J = 0, the block-restricted Liouvillian L_D has eigenvalues −2γ₀ · {1, 3, ..., 2c−1} with multiplicities consistent with the HD-class sizes. No extra rates appear.

**H-mixing prediction verified at J > 0:** Q_SCALE_THREE_BANDS Result 2 shows that dressed modes at fractional ⟨n_XY⟩ appear continuously between the integer rungs as J increases from 0, saturating toward c-specific peaks at Q_peak(c) ∈ {1.5, 1.6, 1.8} for c ∈ {2, 3, 4}. The c = 1 blocks remain at W ≡ 0 for all Q.

Data: `simulations/results/q_scale_n_scaling/` (block-L computations for N = 4..8).
Script: [`q_scale_n_scaling.py`](../../simulations/q_scale_n_scaling.py).

## Relation to existing repo results

- **[Absorption Theorem](PROOF_ABSORPTION_THEOREM.md):** central prerequisite. The identity α = 2γ₀ · ⟨n_XY⟩ is the single non-combinatorial input in the proof.
- **[F70 ΔN selection rule](PROOF_DELTA_N_SELECTION_RULE.md):** site-local observables couple only to |ΔN| ≤ 1 blocks. The (n, n+1) coherence blocks studied here are exactly the |ΔN| = 1 case that survives to single-site observables. F70 gives the dynamical relevance of chromaticity: only these blocks enter per-site purity, per-site expectations, or any c₁-style closure-breaking coefficient.
- **[F73 spatial-sum closure](../ANALYTICAL_FORMULAS.md#f73):** the c = 1 monochromatic case. F73's universal decay (1/2) · exp(−4γ₀ t) for the vac-SE coherent probe in the (0, 1) block is the single rate 2γ₀ · 1 at HD = 1, times a factor of 2 from the probe normalization.
- **[F61 n_XY parity selection rule](PROOF_PARITY_SELECTION_RULE.md):** partitions the Liouvillian spectrum by n_XY mod 2. In the (n, n+1) block, all HD are odd, so all basis pairs have odd n_XY parity. F61 confirms this is a stable sector under dynamics.
- **[Q_SCALE_THREE_BANDS](../../experiments/Q_SCALE_THREE_BANDS.md) Result 3:** primary empirical verification. Chromaticity first named and formalised there; this proof gives the algebraic backing.
