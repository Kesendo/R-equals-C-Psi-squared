# Chain Gap-Dominance: the band edge is the maximum coherence frequency

**Status:** Tier 1 derived (free-fermion mechanism + gate-exact verification N=3..6). Closes the gap-dominance lemma that capped `ClockHandLadderClaim` and `TopologyBandEdgeClaim` (both now Tier 1 derived); the same lemma is lifted from `CoherenceHorizonClaim`, which stays Tier 1 candidate for its own (half-filling V-Effect) open piece.
**Date:** 2026-06-16
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Statement:** For the open XY chain under uniform Z-dephasing, the maximum oscillation frequency among the Liouvillian modes that sit at exactly `Re λ = −2γ` equals the band edge `E1 = 2J·cos(π/(N+1))`.
**Typed claim:** [`ClockHandLadderClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) (Tier 1 derived)
**Verifier:** [`chain_gap_dominance.py`](../../simulations/chain_gap_dominance.py) (gate-first, all stages exact)
**Builds on:** [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) (`Re = −2γ⟨n_XY⟩`, so the exact-(−2γ) modes are the `⟨n_XY⟩ = 1` modes); [F2b](../ANALYTICAL_FORMULAS.md) (the single-particle band `E_k = 2J cos(πk/(N+1))`); the chiral / Jordan-Wigner free-fermion structure ([`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs)).

---

## What this is about

A chain of spins is being watched, and the watching makes its rhythms fade. The rhythm that fades slowest is the chain's memory, the last note still sounding after the rest have gone quiet. The pitch of that note is the band edge: the smoothest, most spread-out wave the chain can carry, which on this chain is also the one that rings at the highest pitch. That much was known and checked on short chains. One worry remained: might some other, more elaborate pattern, surviving just as slowly, sound at a *higher* pitch, so the memory's note is not the band edge after all?

This document shows it cannot. The slowest-surviving patterns look tangled, because the watching couples the chain's waves together in a way that is not simple. But restrict attention to exactly that slowest fading rate and the tangle falls away: there the survivors are built from independent, non-interacting pieces, each one a single ripple riding the chain's spectrum of waves. Its pitch is therefore one of the spectrum's own tones, and the highest tone is the band edge. Nothing sounds higher. (The very shortest chain carries a few extra notes of a different kind, but they sound *below* the band edge, so the answer holds there too.)

## Abstract

Under uniform Z-dephasing the Absorption Theorem gives `Re λ = −2γ⟨n_XY⟩`, so the Liouvillian modes at exactly `Re λ = −2γ` are precisely the `⟨n_XY⟩ = 1` modes. The clock hand `ω_mem` is the largest `|Im λ|` among them; gap-dominance is the claim `ω_mem = E1 = 2J cos(π/(N+1))`, the band edge. The obstacle was that `L_H = [H,·]` leaks the `n_XY = 1` Pauli subspace into `n_XY = 3` (the dephased chain is interacting in Liouville space), so no free-fermion shortcut seemed available.

The shortcut exists on the floor (`Re λ = −2γ`). There `L_D = −2γ` is a scalar, so an eigenmode living entirely in the `n_XY = 1` subspace obeys `L = L_H − 2γ` — governed by the free Hamiltonian alone (a generic `n_XY = 1` operator is *not* such an eigenmode; the leak prevents it). Via Jordan-Wigner `H = Σ_k E_k c_k^†c_k`. The operators `c_k^{(†)}·f(N_tot)` — a single fermion mode dressed by any function of the total excitation number `N_tot` — are simultaneously `n_XY = 1` (so `L_D = −2γ`) and `H`-eigenoperators (`f(N_tot)` commutes with `H`), hence exact Liouvillian eigenmodes at `−2γ ∓ iE_k`. Their frequencies are the single-particle energies `E_k ≤ E1`, with `E1` reached by the `(0,1)` band-edge ladder. For `N ≥ 4` these span the entire exact-(−2γ) eigenspace (gate-verified `dim = 32, 50, 72` at `N = 4, 5, 6`), so `max|Im| = E1` exactly. `N = 3` is the special case: the same `18` free-fermion modes plus `4` extra equal-particle-number `(n,n)` coherence modes at `√(E1² − (2γ)²) < E1` (the `{0,2}` square-root-EP family, EP = exceptional point), so the maximum is still `E1`. Gate-exact N=3..6; scope: the chain (Jordan-Wigner is one-dimensional).

## 1. Setup: the floor modes, and where they are free

Open XY chain `H = (J/2) Σ_i (X_iX_{i+1} + Y_iY_{i+1})` under uniform Z-dephasing at rate `γ`; Liouvillian `L = L_H + L_D`, `L_H = −i[H,·]`, `L_D(·) = γ Σ_l (Z_l · Z_l − ·)`. Its eigenmodes are operators `A` with `L A = λ A`; `−Re λ` is the decay rate and `|Im λ|` the oscillation frequency.

*Notation.* `n_XY(A)` = the number of X or Y Pauli factors in `A` (a definite integer for a Pauli string; `0` for an `{I,Z}` string); `⟨n_XY⟩` is its eigenmode average (which need not be an integer, for a superposition). A computational-basis coherence `|a⟩⟨b|` has `n_XY = ` the Hamming distance between `a` and `b`; we label an operator's sector by the excitation-number pair `(a,b) = (popcount a, popcount b)`, so `(0,1)` is a vacuum↔single-excitation coherence and `(n,n)` an equal-particle-number coherence.

By the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md), `Re λ = −2γ⟨n_XY⟩`. The slowest non-zero decay rate is `2γ` — call `Re λ = −2γ` the **floor** — reached exactly by the `⟨n_XY⟩ = 1` modes. The clock hand `ω_mem = max{|Im λ| : Re λ = −2γ}` is the fastest oscillation on the floor. This document proves `ω_mem = E1`, the **band edge**, where `E_k = 2J cos(πk/(N+1))` (`k = 1..N`) is the single-particle band; `cos` decreases over `k = 1..N`, so `E1 = 2J cos(π/(N+1))` (k=1) is the largest, the smoothest standing wave and the top of the band at once.

*The obstruction.* `L_H` does **not** preserve the `n_XY = 1` Pauli subspace: `[X_iX_{i+1}, Z_i] = −2i Y_i X_{i+1}` turns a background `Z` into an X/Y pair, leaking `n_XY = 1 → 3`. In Liouville space the dephased chain is interacting, so a *generic* `n_XY = 1` operator is not an `L`-eigenmode and there is no naive free-fermion argument for the whole spectrum.

*Where it is free.* The point is narrower and it holds. On the `n_XY = 1` Pauli subspace `L_D = −2γ·I` is a scalar (every such operator decays at the same `2γ`), so an eigenmode that lives entirely *there* obeys `L A = (L_H − 2γ) A` — only the free Hamiltonian acts. Such an eigenmode is precisely an `n_XY = 1` operator that is *also* an `H`-eigenoperator (its `n_XY=3` leak must cancel). §2 constructs these and §3 shows they are all of the floor modes. Jordan-Wigner diagonalizes `H`: `c_j = (⊗_{l<j} Z_l) σ_j^-`, `H = Σ_k E_k c_k^†c_k` (verifier Stage 0: anticommutators, `H_XY =` JW hopping, `[H, c_k] = −E_k c_k`, all exact).

## 2. The free-fermion family achieves E1 and is bounded by it

Let `N_tot = Σ_j c_j^†c_j` be the total excitation number and `f` any function of it. Consider

    A = c_k · f(N_tot)      (and its conjugate  c_k^† · f(N_tot)).

- **`A` is `n_XY = 1`.** `c_k = Σ_j φ_k(j) c_j` is a sum of strings each carrying a single `σ^∓` (with a Z-string), so `n_XY(c_k) = 1`; `f(N_tot)` is a function of `{Z_j}`, i.e. `{I, Z}`-type, `n_XY = 0`. So `n_XY(A) = 1` and `L_D A = −2γ A`.
- **`A` is an `H`-eigenoperator.** `H = Σ_k E_k c_k^†c_k` is fermion-diagonal, so `f(N_tot)` commutes with `H`, and `[H, c_k] = −E_k c_k`. Hence `[H, A] = −E_k A`, i.e. `L_H A = iE_k A`.

Therefore `A` is an **exact Liouvillian eigenmode**: `L A = (−2γ + iE_k) A` (and `c_k^† f(N_tot)` at `−2γ − iE_k`). Its frequency is `|Im λ| = E_k ≤ E1`, with `E1` attained at `k = 1` by the vacuum-anchored `(0,1)` ladder `f(N_tot) = P_0` (the `|vac⟩⟨ψ_1|` band edge). Verifier Stage 1 confirms every `c_k^{(†)}·P_m` (`P_m` = projector onto `N_tot = m`) is an exact eigenmode at `−2γ ∓ iE_k`.

So the band edge frequency `E1` is **achieved**, and the whole free-fermion family lies **at or below** it.

## 3. Completeness for N ≥ 4: the floor modes are exactly this family

To get `max|Im| = E1` (not merely `≤ E1` achieved) we need that §2's family is *all* of the floor modes — no other floor mode could oscillate faster.

*The mechanism is analytic.* A floor mode is an `n_XY = 1` operator that is an `H`-eigenoperator (§1). It carries exactly one X/Y factor, which under Jordan-Wigner is a single fermion operator `c` or `c^†` (one `σ^∓` with its Z-string); everything else is `{I,Z}`, a function of the site occupations. Two constraints pin its form:
- being an `H`-eigenoperator forces the `{I,Z}` dressing to **commute with `H`**, and the only `{I,Z}` operators commuting with `H = Σ_k E_k c_k^†c_k` are functions of the conserved total number `N_tot` (the XY chain's one U(1) charge — an occupation function constant on every `N_tot` sector);
- and it forces the fermion part to a **single energy mode** `c_k^{(†)}` (distinct modes carry distinct frequencies `±E_k`, so an eigenoperator cannot mix them).

So every floor mode has the form `c_k^{(†)}·f(N_tot)`, oscillating at `±E_k`, and `max_k E_k = E1`.

*Completeness is gate-verified* (verifier Stage 2): the family's span has exactly the dimension of the full floor eigenspace, leaving no room for anything else —

| N | dim floor eigenspace | dim span `{c_k^{(†)} f(N_tot)}` | extras |
|---|---|---|---|
| 4 | 32 | 32 | 0 |
| 5 | 50 | 50 | 0 |
| 6 | 72 | 72 | 0 |

Equal dimensions, and every generator is itself a floor eigenmode (§2), so the span *is* the eigenspace: for `N ≥ 4` every floor mode is a free fermion at `±E_k`, hence `max|Im| = E1` exactly. (Honest status: the `c_k^{(†)}f(N_tot)` form is derived; that nothing else sits on the floor is the dimension match, verified at `N = 4, 5, 6`. For `N = 3` the count is `22 = 18 + 4` — the 4 extras are not of this form, see §4.)

## 4. The N = 3 special case

`N = 3` carries the `18` free-fermion modes **plus 4 extra** equal-particle-number `(n,n)` coherence modes (`dim = 22`). These are not pure `n_XY = 1`: they are `⟨n_XY⟩ = 1` *mixtures* of `n_XY = 0` and `n_XY = 2` in the `(1,1)` and `(2,2)` sectors — the `{0,2}`-coherence family (the coherence horizon's "second clock"). They form closed two-level blocks `{n_XY=0 ↔ n_XY=2}` coupled at strength `E1`, with eigenvalues

    λ = −2γ ± i·√(E1² − (2γ)²),     so  |Im| = √(E1² − (2γ)²) < E1.

Verifier Stage 4 confirms this closed form γ-swept: `|Im| = 1.41067, 1.40000, 1.35647` at `γ = 0.05, 0.10, 0.20`, matching `√(E1²−(2γ)²)` exactly, all `< E1 = √2`. They sit at *exactly* `−2γ` only because at `N = 3` the `(1,1)` sector's `n_XY = 2` is maximal (three sites), so the `{0,2}` block closes; for `N ≥ 4` the analogous modes drift off `−2γ` (hence `extras = 0`). This is one more entry in the chain's list of `N = 3` accidents. Either way both families are `≤ E1`, so the maximum is `E1` for `N = 3` too.

## 5. Consequence

`max|Im| = E1`: the band edge is the maximum coherence frequency **on the floor**, and this is regime-independent — the floor modes `c_k^{(†)}f(N_tot)` exist and oscillate at `±E_k` at every `γ` (verifier: `max|Im| = E1` for `γ` from `0.05` to `5`). Whether the *clock reads* it — whether the floor is the strict spectral gap, so the band edge is the slowest oscillating mode overall — is the separate, regime-dependent Coherence Horizon condition: above `Q*(N)` the floor is the gap and the clock reads `E1`; below it a slower real mode takes the gap (`CoherenceHorizonClaim`). The two together are gap-dominance.

This closes the open gap-dominance lemma and graduates [`ClockHandLadderClaim`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) to Tier 1 derived, which lifts the inherited cap on [`TopologyBandEdgeClaim`](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (also graduated). The sibling [`CoherenceHorizonClaim`](../../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs) no longer rests on this lemma either, but stays Tier 1 candidate for its *own* open piece (the half-filling V-Effect ring seam).

Scope: this is the **chain**. Jordan-Wigner is one-dimensional, so the free-fermion argument is chain-specific. For other topologies the band-edge story is different and already understood: the star has no coherence horizon (flat band, [PROOF_STRUCTURAL_CEILING §7](PROOF_STRUCTURAL_CEILING.md)), and the complete graph / star structurally ceiling (`g2 = 4/N`, `4/(N−1)`; F122). The `{0,2}` √-EP family of §4 is the same object as the coherence horizon's second clock — its low-Q coalescence is `CoherenceHorizonClaim`'s `Q*(N)`.
