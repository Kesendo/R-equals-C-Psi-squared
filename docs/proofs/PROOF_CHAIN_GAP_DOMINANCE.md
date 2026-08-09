# Chain Gap-Dominance: the band edge is the maximum coherence frequency

**Status:** Tier 1 derived FOR N ≥ 3 (free-fermion mechanism + gate-exact verification N=3..6). At N = 2 the statement is FALSE for Q > 2/√3 (§4.2, found 2026-08-02 while hardening the §4 verifier); the N ≥ 3 result is unaffected, and so are the two claims that rest on it, but not because they avoid N = 2: `ClockHandLadderClaim` explicitly carries the N = 2 exceptional point and already states the pulled hand `2√(J²−γ²)` there. What this document did not cover, and now does, is that at N = 2 that hand is the floor MAXIMUM above `Q = 2/√3`. Closes the gap-dominance lemma that capped `ClockHandLadderClaim` and `TopologyBandEdgeClaim` (both now Tier 1 derived); the same lemma is lifted from `CoherenceHorizonClaim`, whose own open piece (the ring 2-excitation `(2,2)/(N−2,N−2)` doublet V-Effect seam) was later resolved by the reviewed [Ring Handover Slope proof](PROOF_RING_HANDOVER_SLOPE.md); all three claims are now Tier 1 derived.
**Date:** 2026-06-16
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Statement:** For the open XY chain under uniform Z-dephasing, at `N ≥ 3`, the maximum oscillation frequency among the Liouvillian modes that sit at exactly `Re λ = −2γ` equals the band edge `E1 = 2J·cos(π/(N+1))`. At `N = 2` it holds only for `Q = J/γ ≤ 2/√3`; above that the floor's maximum is `2√(J²−γ²) > E1` and the equality fails (§4.2). Bond uniformity is NOT required for the maximality (Jordan-Wigner is blind to the hopping profile, so §§2-3 carry over and `E1` becomes the largest single-particle energy of the weighted hopping matrix, gated at `N = 4, 5` under random bond weights by `mixed02_block_threshold.py` G7); it is required only for §4.1's `N = 3` extras to sit exactly on the floor (§4.3). It is vacuous at `N = 2`, which has a single bond, so the counterexample of §4.2 cannot be detuned away and is unconditional in bond space.
**Typed claim:** [`ClockHandLadderClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) (Tier 1 derived)
**Verifier:** [`chain_gap_dominance.py`](../../simulations/chain_gap_dominance.py) (gate-first, all stages exact) + [`chain_gap_dominance_skeptic.py`](../../simulations/chain_gap_dominance_skeptic.py) (the independent-convention adversarial companion: a from-scratch column-stack re-derivation, a basis-independent purity gate, and a γ≪J/γ~J/γ≫J regime sweep: two vec conventions agree, so the result is not a stacking artifact) + [`mixed02_block_threshold.py`](../../simulations/mixed02_block_threshold.py) (§4: where the {0,2} extras are, the model-dependent rung, the N = 2 crossover at Q = 2/√3, and the bond-detuning fence; six gates, each two-sided)
**Builds on:** [the Absorption Theorem](PROOF_ABSORPTION_THEOREM.md) (`Re = −2γ⟨n_XY⟩`, so the exact-(−2γ) modes are the `⟨n_XY⟩ = 1` modes); [F2b](../ANALYTICAL_FORMULAS.md) (the single-particle band `E_k = 2J cos(πk/(N+1))`); the chiral / Jordan-Wigner free-fermion structure ([`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs)).

---

## What this is about

A chain of spins stands in light, and the light makes its rhythms fade. The rhythm that fades slowest is the chain's memory, the last note still sounding after the rest have gone quiet. The pitch of that note is the band edge: the smoothest, most spread-out wave the chain can carry, which on this chain is also the one that rings at the highest pitch. That much was known and checked on short chains. One worry remained: might some other, more elaborate pattern, surviving just as slowly, sound at a *higher* pitch, so the memory's note is not the band edge after all?

This document shows it cannot. The patterns that could survive that long look tangled, because the chain's own motion mixes a pattern that pays one unit of the light's price with patterns that pay three. But restrict attention to exactly that slowest fading rate and the tangle falls away: there the survivors are built from independent, non-interacting pieces, each one a single ripple riding the chain's spectrum of waves. Its pitch is therefore one of the spectrum's own tones, and the highest tone is the band edge. Nothing sounds higher. (Short chains carry a few extra notes of a different kind. At three spins they sound *below* the band edge, so the answer holds there too. At two spins they sound *above* it whenever the coupling is strong against the light, and there the answer is no: the memory's note is one of the extras, not the band edge. §4 is where that is worked out.)

## Abstract

Under uniform Z-dephasing the Absorption Theorem gives `Re λ = −2γ⟨n_XY⟩`, so the Liouvillian modes at exactly `Re λ = −2γ` are precisely the `⟨n_XY⟩ = 1` modes. The clock hand `ω_mem` is the largest `|Im λ|` among them; gap-dominance is the claim `ω_mem = E1 = 2J cos(π/(N+1))`, the band edge. The obstacle was that `L_H = −i[H,·]` leaks the `n_XY = 1` Pauli subspace into `n_XY = 3` (the dephased chain is interacting in Liouville space), so no free-fermion shortcut seemed available.

The shortcut exists on the floor (`Re λ = −2γ`). There `L_D = −2γ` is a scalar, so an eigenmode living entirely in the `n_XY = 1` subspace obeys `L = L_H − 2γ`, governed by the free Hamiltonian alone (a generic `n_XY = 1` operator is *not* such an eigenmode; the leak prevents it). Via Jordan-Wigner `H = Σ_k E_k c_k^†c_k`. The operators `c_k^{(†)}·f(N_tot)`, a single fermion mode dressed by any function of the total excitation number `N_tot`, are simultaneously `n_XY = 1` (so `L_D = −2γ`) and `H`-eigenoperators (`f(N_tot)` commutes with `H`), hence exact Liouvillian eigenmodes at `−2γ ∓ iE_k`. Their frequencies are the single-particle energies `E_k ≤ E1`, with `E1` reached by the `(0,1)` band-edge ladder. For `N ≥ 4` these span the entire exact-(−2γ) eigenspace (gate-verified `dim = 32, 50, 72` at `N = 4, 5, 6`), so `max|Im| = E1` exactly. `N = 3` carries the same `18` free-fermion modes plus `4` extra equal-particle-number `(n,n)` coherence modes at `√(E1² − (2γ)²) < E1` (the `{0,2}` square-root-EP family, EP = exceptional point), so the maximum is still `E1`. `N = 2` carries the same kind of extras and they are NOT below the band edge: they sit at `2√(J²−γ²)`, which exceeds `E1 = J` for `Q > 2/√3`, so the equality fails there (§4.2). Gate-exact N=3..6; the low-N extras and the `N = 2` failure are gated separately by `mixed02_block_threshold.py`. Scope: the chain (Jordan-Wigner is one-dimensional) and `N ≥ 3`, unqualified in `γ` and in the bond profile; uniform bonds are needed only for §4's extras, not for the maximality (§4.3).

## 1. Setup: the floor modes, and where they are free

Open XY chain `H = (J/2) Σ_i (X_iX_{i+1} + Y_iY_{i+1})` under uniform Z-dephasing at rate `γ`; Liouvillian `L = L_H + L_D`, `L_H = −i[H,·]`, `L_D(·) = γ Σ_l (Z_l · Z_l − ·)`. Its eigenmodes are operators `A` with `L A = λ A`; `−Re λ` is the decay rate and `|Im λ|` the oscillation frequency.

*Notation.* `n_XY(A)` = the number of X or Y Pauli factors in `A` (a definite integer for a Pauli string; `0` for an `{I,Z}` string); `⟨n_XY⟩` is its eigenmode average (which need not be an integer, for a superposition). A computational-basis coherence `|a⟩⟨b|` has `n_XY = ` the Hamming distance between `a` and `b`; we label an operator's sector by the excitation-number pair `(a,b) = (popcount a, popcount b)`, so `(0,1)` is a vacuum↔single-excitation coherence and `(n,n)` an equal-particle-number coherence.

By the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md), `Re λ = −2γ⟨n_XY⟩`. In this document's regime (coupling ratio above the N-dependent threshold Q*_gap of [that proof](PROOF_ABSORPTION_THEOREM.md) §4.3; note the unit conversion, since that proof tabulates Q*_gap in Pauli-J while this document's H is `(J/2)Σ(XX+YY)`, so its own thresholds are twice the tabulated XY values: 1.414, 1.879, 2.372 at N=3,4,5) the slowest non-zero decay rate is `2γ`; call `Re λ = −2γ` the **floor**, reached exactly by the `⟨n_XY⟩ = 1` modes. The clock hand `ω_mem = max{|Im λ| : Re λ = −2γ}` is the fastest oscillation on the floor. This document proves `ω_mem = E1`, the **band edge**, for `N ≥ 3` (§4.2 has the `N = 2` exception and its `Q` threshold), where `E_k = 2J cos(πk/(N+1))` (`k = 1..N`) is the single-particle band; `cos` decreases over `k = 1..N`, so `E1 = 2J cos(π/(N+1))` (k=1) is the largest, the smoothest standing wave and the top of the band at once.

*The obstruction.* `L_H` does **not** preserve the `n_XY = 1` Pauli subspace: `[X_iX_{i+1}, Z_i] = −2i Y_i X_{i+1}` turns a background `Z` into an X/Y pair, leaking `n_XY = 1 → 3`. In Liouville space the dephased chain is interacting, so a *generic* `n_XY = 1` operator is not an `L`-eigenmode and there is no naive free-fermion argument for the whole spectrum.

*Where it is free.* The point is narrower and it holds. On the `n_XY = 1` Pauli subspace `L_D = −2γ·I` is a scalar (every such operator decays at the same `2γ`), so an eigenmode that lives entirely *there* obeys `L A = (L_H − 2γ) A`: only the free Hamiltonian acts. Such an eigenmode is precisely an `n_XY = 1` operator that is *also* an `H`-eigenoperator (its `n_XY=3` leak must cancel). §2 constructs these and §3 shows they are all of the floor modes. Jordan-Wigner diagonalizes `H`: `c_j = (⊗_{l<j} Z_l) σ_j^-`, `H = Σ_k E_k c_k^†c_k` (verifier Stage 0: anticommutators, `H_XY =` JW hopping, `[H, c_k] = −E_k c_k`, all exact).

## 2. The free-fermion family achieves E1 and is bounded by it

Let `N_tot = Σ_j c_j^†c_j` be the total excitation number and `f` any function of it. Consider

    A = c_k · f(N_tot)      (and its conjugate  c_k^† · f(N_tot)).

- **`A` is `n_XY = 1`.** `c_k = Σ_j φ_k(j) c_j` is a sum of strings each carrying a single `σ^∓` (with a Z-string), so `n_XY(c_k) = 1`; `f(N_tot)` is a function of `{Z_j}`, i.e. `{I, Z}`-type, `n_XY = 0`. So `n_XY(A) = 1` and `L_D A = −2γ A`.
- **`A` is an `H`-eigenoperator.** `H = Σ_k E_k c_k^†c_k` is fermion-diagonal, so `f(N_tot)` commutes with `H`, and `[H, c_k] = −E_k c_k`. Hence `[H, A] = −E_k A`, i.e. `L_H A = iE_k A`.

Therefore `A` is an **exact Liouvillian eigenmode**: `L A = (−2γ + iE_k) A` (and `c_k^† f(N_tot)` at `−2γ − iE_k`). Its frequency is `|Im λ| = E_k ≤ E1`, with `E1` attained at `k = 1` by the vacuum-anchored `(0,1)` ladder `f(N_tot) = P_0` (the `|vac⟩⟨ψ_1|` band edge). Verifier Stage 1 confirms every `c_k^{(†)}·P_m` (`P_m` = projector onto `N_tot = m`) is an exact eigenmode at `−2γ ∓ iE_k`.

So the band edge frequency `E1` is **achieved**, and the whole free-fermion family lies **at or below** it.

## 3. Completeness for N ≥ 4: the floor modes are exactly this family

To get `max|Im| = E1` (not merely `≤ E1` achieved) we need that §2's family is *all* of the floor modes: no other floor mode could oscillate faster.

*The mechanism is analytic.* A floor mode is an `n_XY = 1` operator that is an `H`-eigenoperator (§1). It carries exactly one X/Y factor, which under Jordan-Wigner is a single fermion operator `c` or `c^†` (one `σ^∓` with its Z-string); everything else is `{I,Z}`, a function of the site occupations. Two constraints pin its form:
- being an `H`-eigenoperator forces the `{I,Z}` dressing to **commute with `H`**, and the only `{I,Z}` operators commuting with `H = Σ_k E_k c_k^†c_k` are functions of the conserved total number `N_tot` (the XY chain's one U(1) charge, an occupation function constant on every `N_tot` sector);
- and it forces the fermion part to a **single energy mode** `c_k^{(†)}` (distinct modes carry distinct frequencies `±E_k`, so an eigenoperator cannot mix them).

So every floor mode has the form `c_k^{(†)}·f(N_tot)`, oscillating at `±E_k`, and `max_k E_k = E1`.

*Completeness is gate-verified* (verifier Stage 2): the family's span has exactly the dimension of the full floor eigenspace, leaving no room for anything else:

| N | dim floor eigenspace | dim span `{c_k^{(†)} f(N_tot)}` | extras |
|---|---|---|---|
| 4 | 32 | 32 | 0 |
| 5 | 50 | 50 | 0 |
| 6 | 72 | 72 | 0 |

Equal dimensions, and every generator is itself a floor eigenmode (§2), so the span *is* the eigenspace: for `N ≥ 4` every floor mode is a free fermion at `±E_k`, hence `max|Im| = E1` exactly. (Honest status: the `c_k^{(†)}f(N_tot)` form is derived; that nothing else sits on the floor is the dimension match, verified at `N = 4, 5, 6`. For `N = 3` the count is `22 = 18 + 4`; the 4 extras are not of this form, see §4.)

## 4. The low-N extras, and the one place the statement fails

`N = 3` carries the `18` free-fermion modes **plus 4 extra** equal-particle-number `(n,n)` coherence modes (`dim = 22`). These are not pure `n_XY = 1`: they are `⟨n_XY⟩ = 1` *mixtures* of `n_XY = 0` and `n_XY = 2` in the `(1,1)` and `(2,2)` sectors, the `{0,2}`-coherence family (the coherence horizon's "second clock"). They form closed two-level blocks `{n_XY=0 ↔ n_XY=2}` coupled at strength `E1`, with eigenvalues

    λ = −2γ ± i·√(E1² − (2γ)²),     so  |Im| = √(E1² − (2γ)²) < E1.

Verifier Stage 4 confirms this closed form γ-swept: `|Im| = 1.41067, 1.40000, 1.35647` at `γ = 0.05, 0.10, 0.20`, matching `√(E1²−(2γ)²)` exactly, all `< E1 = √2`. So at `N = 3` both families are `≤ E1` and the maximum is `E1` there too.

### 4.1 It is not an N = 3 accident: the rung moves with the model

The extras sit at *exactly* `−2γ` when their two-level `{n_XY = 0 ↔ n_XY = 2}` block closes, and it is worth being precise about what closes it, because the obvious answer is not the operative one.

The obvious answer is combinatorial: the block closes when every equal-particle-number coherence in play is at Hamming distance `≤ 2`, i.e. `2·min(n, N−n) ≤ 2`, which holds iff `N ≤ 3`. That is a true statement about the XY chain and it is what bounds the phenomenon here. But it contains no Hamiltonian, so on its own it can say nothing about *which model* shows the extras, and the models differ.

They differ by one rung. Adding the `ZZ` term in this document's own normalisation (Heisenberg, `H = (J/2) Σ (XX+YY+ZZ)`; the factor matters, since with `J Σ(·)` the hopping doubles and the frequency below reads 3.998750 instead) puts a diagonal on the single-excitation sector, and whether that diagonal is CONSTANT decides whether it can break the block. (Constant is the operative word, not non-degenerate: `(0, −J, 0)` repeats `0` and is degenerate, yet it splits the sector perfectly well. A constant diagonal is a multiple of the identity and commutes with everything; a merely degenerate one does not.)

| N | SE `ZZ` diagonal (this document's units) | effect | XY extras | Heisenberg extras |
|---|---|---|---|---|
| 2 | `(−J/2, −J/2)` | **constant** on the sector, so `ZZ` acts as a multiple of the identity there and cannot break anything | 2 | 2 |
| 3 | `(0, −J, 0)` | **not constant**, so `ZZ` splits the sector and the block opens | 4 | 0 |
| ≥ 4 | not constant (at `N = 4` it reads `(J/2, −J/2, −J/2, J/2)`: degenerate in pairs, but not flat) | block already open combinatorially | 0 | 0 |

The two count columns are read at high `Q` (γ = 0.05, J = 1). They are functions of `(N, Q)` and not of `N` alone: each block has its own `Q*` below which the pair coalesces off the floor (`Q*(2) = 1`, `Q*(3) = √2`), so at `γ = 0.8` the `N = 3` extras have already gone and at `γ = 1.5` all of them have. §4.2's regime table is the same statement for `N = 2`.

So the `ZZ` term costs exactly one rung: what the XY chain shows at `N = 3`, Heisenberg shows at `N = 2`, and the Heisenberg `N = 2` frequency reproduces the XY `N = 2` one to machine precision (`1.997498435544` against `1.997498435544` at `γ = 0.05`) for exactly the reason in the table. A low-N structural threshold whose rung depends on the model is not an accident and is not `N = 3`'s.

### 4.2 N = 2, where this document's statement is false

`N = 2` also carries `2` extras, and there they are not below the band edge. With `E1 = 2J cos(π/3) = J` and the block coupling `2J` rather than `E1`, the closed form reads

    |Im| = √((2J)² − (2γ)²) = 2√(J² − γ²),

The crossover itself is not new here: the repository has recorded since 2026-06-12 that the live `max|Im|` clock surfaces the pulled `N = 2` mode only above `Q = 2/√3` (the `clock_hand_ladder` arc, `ClockHandLadderClaim`, and the F2b corollary in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)). What is new is what it means for THIS document: that above that ratio the pulled mode is the floor's MAXIMUM, so gap-dominance fails there. The frequency exceeds `E1 = J` exactly when `2√(J²−γ²) > J`, i.e. when

    Q = J/γ > 2/√3 ≈ 1.1547.

Above that ratio the floor's maximum is the extras' frequency and **not** the band edge, measured directly on the full floor: at `Q = 20` the floor max is `1.997498` against `E1 = 1.000000`. So the headline statement of this document is false at `N = 2` in the high-Q regime, which is the regime the rest of the repository works in.

At or below `Q = 2/√3` it is true, and below `Q = 1` it is true for a second, different reason: the square root turns real, the pair coalesces and leaves the floor altogether (the floor dimension drops from `10` to `8`), which is the coherence horizon's `Q*(2) = 1` exceptional point seen from this side. Three regimes, one statement:

| regime | extras on the floor | floor max | statement |
|---|---|---|---|
| `Q > 2/√3` | yes, at `2√(J²−γ²)` | the extras | **false** |
| `1 < Q ≤ 2/√3` | yes, at `2√(J²−γ²)` | `E1` | true |
| `Q ≤ 1` | none (coalesced off the floor) | `E1` | true |

`N = 2` is therefore a second band-edge counterexample beside the ring's 4-cycle, in the `Q > 2/√3` regime. [`PROOF_RING_GAP_DOMINANCE`](PROOF_RING_GAP_DOMINANCE.md) calls the 4-cycle the unique one; that uniqueness holds among the graphs it sweeps, which start at `N = 4`, and not across `N`.

### 4.3 The extras need uniform bonds; the maximality does not

The `{0,2}` block closes only on the uniform chain. Detuning a single bond by `δ` takes every extra off the floor: at `N = 3, γ = 0.05` the count goes `4 → 0` at `δ = 10⁻²`, `10⁻³` and `10⁻⁴`. This is the same fine-tuning the ring sibling records for its own case, and it is a fence on §4 only, not on §§2-3: the free-fermion family survives non-uniform bonds, since Jordan-Wigner diagonalizes any bond profile.

It also has a moral about tolerances, and the first version of that moral was wrong in the way it was warning about. The first verifier admitted the floor at `|Re λ + 2γ| < 10⁻⁸`, which at `δ = 10⁻⁴` still counted all four detuned modes as floor modes: the gate passed on input where the property it tests is false. Tightening to `10⁻¹¹` does not fix that, it moves it one decade, because the drift off the floor is **quadratic in the detuning**, measured `O(δ^1.999)`:

| δ | 10⁻² | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ |
|---|---|---|---|---|---|
| `\|Re λ + 2γ\|` | 9.9·10⁻⁶ | 1.0·10⁻⁷ | 1.0·10⁻⁹ | 1.0·10⁻¹¹ | 1.0·10⁻¹³ |

So `10⁻¹¹` is a window at `δ ≈ 10⁻⁵`, not a machine-precision test, and any fixed tolerance is a window on `δ`. The quadratic law is the statement without one, and it is what the verifier now gates.

**Verifier:** [`mixed02_block_threshold.py`](../../simulations/mixed02_block_threshold.py) gates all of §4: where the extras are (γ-swept), their `{0: ½, 2: ½}` histogram, the closed form, the `Q = 2/√3` crossover with the facing flipping there and nowhere else, the bond-detuning sensitivity, and the `ZZ`-diagonal mechanism.

## 5. Consequence

`max|Im| = E1` for `N ≥ 3`: the band edge is the maximum coherence frequency **on the floor**, and at those `N` this is regime-independent: the floor modes `c_k^{(†)}f(N_tot)` exist and oscillate at `±E_k` at every `γ` (verifier: `max|Im| = E1` for `γ` from `0.05` to `5`). The `γ`-independence is what fails at `N = 2`, where the extras of §4.2 overtake `E1` above `Q = 2/√3`; there the floor maximum is `γ`-dependent. Whether the *clock reads* it (whether the floor is the strict spectral gap, so the band edge is the slowest oscillating mode overall) is the separate, regime-dependent Coherence Horizon condition: above `Q*(N)` the floor is the gap and the clock reads `E1`; below it a slower real mode takes the gap (`CoherenceHorizonClaim`). The two together are gap-dominance.

This closes the open gap-dominance lemma and graduates [`ClockHandLadderClaim`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) to Tier 1 derived, which lifts the inherited cap on [`TopologyBandEdgeClaim`](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (also graduated). The sibling [`CoherenceHorizonClaim`](../../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs) no longer rests on this lemma either; its *own* open piece (the ring 2-excitation `(2,2)/(N−2,N−2)` doublet V-Effect seam) was resolved by [the Ring Handover Slope proof](PROOF_RING_HANDOVER_SLOPE.md), and it is now Tier 1 derived as well.

Scope: this is the **chain**. Jordan-Wigner is one-dimensional, so the free-fermion argument is chain-specific. For other topologies the band-edge story is different and already understood: the star has no coherence horizon (flat band, [PROOF_STRUCTURAL_CEILING §7](PROOF_STRUCTURAL_CEILING.md)), and the complete graph / star structurally ceiling (`g2 = 4/N`, `4/(N−1)`; F122). The `{0,2}` √-EP family of §4 is the same object as the coherence horizon's second clock; its low-Q coalescence is `CoherenceHorizonClaim`'s `Q*(N)`. The cyclic sibling is [`PROOF_RING_GAP_DOMINANCE.md`](PROOF_RING_GAP_DOMINANCE.md): the ring max is `2J = J·ρ` as well, with N=4 the lone half-filling `(2,2)` √-EP exception that *exceeds* the band top, the twin of this chain's N=2 exception, which also exceeds its band top; the chain's N=3 exception is the other-facing one, sitting below `E1` and changing nothing. `max|Im| = J·ρ` on the floor is not a chain or ring property at all: it was measured on all 38 connected labelled graphs at N=4, on star, complete and asymmetric graphs at N=5 and 6, and under random bond weights, with the three labellings of the 4-cycle the only violations **in that sweep**. That sweep starts at `N = 4`. Across `N` the 4-cycle is not the only violation: the two-site chain `P₂` is a second one for `Q > 2/√3` (§4.2), and it is invisible to every gate listed above because none of them builds `N = 2`. What is topology-specific is the closed form for ρ, and which exceptions a given sweep can see.
