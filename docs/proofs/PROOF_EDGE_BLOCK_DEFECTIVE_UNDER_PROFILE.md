# PROOF: Under a Rate Profile the Edge Block Carries a Defective EP, and on the XY Chain It Takes One Knob

**Status:** Tier 1 derived for the structural half, which is an exact algebraic identity with no eigensolver in it. The existence half is certified by a sign change of a real function, an intermediate-value argument rather than a tolerance. Both are gated. **The fencing of the sites in §(i) is deliberately NOT done here**: the arc that owns them asks for it as its own pass with its own review rounds, and §(i) is the inventory that pass should start from.
**Date:** 2026-08-20
**Authors:** Thomas Wicht, Claude (Opus 5)
**Script:** [`simulations/edge_block_defective_ep_gate.py`](../../simulations/edge_block_defective_ep_gate.py) (147 gates, ~3 min)
**Builds on:**
- [`PROOF_CODIM1_BY_ADDITIVITY.md`](PROOF_CODIM1_BY_ADDITIVITY.md), whose **Edge lemma** is the statement refuted here and whose **window-edge lemma** and **rate window** survive and do the confining.
- **F152** in [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md), the vacuum-block generator in closed form. F152 is fenced to |Δ| = 1; the Edge lemma is about Δ = 0. Both books appear below, and §(a) says which is which.
- [`simulations/d10_block_closure_verify.py`](../../simulations/d10_block_closure_verify.py), the committed reader that takes the block off the full 4^N Liouvillian with the leak measured rather than assumed.
- [`D10_W1_DISPERSION.md`](derivations/D10_W1_DISPERSION.md), which owns the block and warns that "the w=1 sector" is the wrong name for it.

## What the repo already held, and what this adds

**Sweep record.** `compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs` returned the arc `site_resolved_vacuum_block` (line 3819). The **first item of its (b2) inventory** (line 4118) is exactly this question and says of it, verbatim, *"'no defective EP can live on an edge block' is no longer forbidden a priori ... Settle it before fencing, do not assume either way."* Its **second** item (line 4136) is a separate and larger inventory, the unfenced whole-sector `Re = −2γ` claims, which this document does not address and does not close. `docs/ANALYTICAL_FORMULAS.md` returned **F152**, which already owns the operator, already records the non-normality under a profile, already carries the N = 2 real-part coalescence threshold, and already states that *"the flat-γ argument 'the dissipator is scalar, so no Jordan block' does not transfer to the profile pencil"*; **F153**, the correctly fenced sibling, whose header carries "SCOPE: uniform γ REQUIRED"; **F125**, which repeats the Edge lemma verbatim inside an entry scoped to "arbitrary rates γ_j"; and **F2**, whose N−1 distinct frequencies are the uniform-J special case of §(g). `docs/proofs/` returned `PROOF_CODIM1_BY_ADDITIVITY.md`; `PROOF_R90_FROZEN_DIVISOR.md`, whose section 9 exhibits Jordan blocks under γ profiles but of a different operator, see §(h); and `PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md:145`, which already fences **this document's own premise** in one sentence, *"Under site-dependent γ_l the dissipator acts on the (0,1) coherence block as diag(−2γ_l) rather than as a scalar, no longer commutes with the one-magnon Hamiltonian"*. `hypotheses/` returned `INHERITED_RULES_AND_THE_OWN.md:88`, which **already states §(f)'s window** on this exact block, and `DIABOLIC_BY_INTEGRABILITY.md`, the repo's defective-versus-semisimple discriminator. `experiments/` returned measurements of this block under profiles (`CONCENTRATOR_OPTICS.md` Result 3, `ANALYTICAL_SPECTRUM.md`, `VEFFECT_CAVITY_MODES.md`) and, on its **defectiveness**, nothing in either direction, with one caution recorded in §(h); `THE_MIRROR_TRANSVERSAL_CERTIFICATE.md:59-66` returned a standing warning that reading §6's window under a profile is *"this note's extrapolation, not §6's own statement"*, which §(f) answers rather than ignores. `fw.Confirmations` and the C# `ConfirmationsRegistry` returned no entry bearing on edge-block normality, one that bears on the object itself (`block_cpsi_saturation_kingston_may2026`, a hardware measurement of the (0,1) block), and one EP-adjacent entry that explicitly disclaims this question. `docs/GLOSSARY.md` returned no entry for defective versus diabolic on a coherence block, and one alias trap named in §(h). `docs/CAUGHT_ERRORS.md` returned the 2026-06-21 Petermann retraction and its 2026-07-07 partial reversal, which is why no Petermann quantity appears below, and the 2026-08-20 entry, whose lesson about uniform fixtures is the blindness that let this question stand open.

**What is new, stated narrowly.** Not the operator, not the observation that the old argument fails, and not the N = 2 point: F152 owns the first two, and `BlockSpectrumWitnessTests.cs:517-535` together with `OpenArcsRegistry.cs:3968` already call that point *"a defective 2x2"* in a different J book. What is new is **Lemma A**, which appears nowhere in the repo; the **symbolic certificate** for the N = 2 reading the repo had asserted and corroborated but never proved; the **codimension** result of §(e), which turns the question from fine-tuning into one knob; and the **answer** the arc asked for.

## (a) Two operators, and which one the lemma is about

On the (0,1) coherence block, spanned by the N coherences `|0⟩⟨j|` with a single excitation at site j, the generator is an N × N matrix, and it takes two forms:

  **Δ = 1 (XXX):  M = −2i·𝓛_J − 2·diag(γ)**, the Laplacian form, which is F152's and is fenced there to |Δ| = 1.
  **Δ = 0 (XY):   M = +2i·A_J − 2·diag(γ)**, the adjacency form.

`PROOF_CODIM1_BY_ADDITIVITY.md:24` and F125 are scoped to *"the N-site XY chain (Δ=0, open ends, arbitrary bond profile J_b)"*. **The Edge lemma is about Δ = 0**, so every witness below is on the adjacency form. This is not pedantry: a first version of this work put its witnesses on the Laplacian form and was fully green while measuring an operator the refuted lemma is not about. On the adjacency form those same profiles have eigenvalue splits of 1.54, 1.67 and 0.107, which are not EPs at all.

Gate **G1** anchors both forms against the full 4^N Liouvillian through `d10_block_closure_verify.measure_block`, which measures the leak over the whole complement. The leak is exactly 0.0 at every size and in both books. The residual against the closed form obeys a measured noise law, `residual ≤ C·eps·‖block‖·N` with C flat across N = 3 to 7 and no trend; the linear factor N is the block dimension, and it was measured rather than assumed.

Write `M = A + iB` with `A = −2·diag(γ)` in either book. At uniform γ that term is scalar and commutes with everything, which is the whole content of the old argument.

## (b) Lemma A: on a path, M is non-derogatory

**Lemma A.** Let the graph be a path, vertices in their natural order, with every bond coupling `J_b ≠ 0`. Then `rank(M − λI) ≥ N − 1` for every λ ∈ ℂ, in either Δ book. Equivalently, the geometric multiplicity of every λ is at most one, hence exactly one at every eigenvalue.

*Proof.* On a path, M is tridiagonal, and its off-diagonal entries are `M[b, b+1] = M[b+1, b] = 2i·J_b` in both books, the Δ term being diagonal. Delete row 0 and column N−1 from `M − λI`. The remaining (N−1) × (N−1) matrix is triangular, with those off-diagonal entries on its diagonal, so its determinant is `∏_b (2i·J_b)`: free of λ, and nonzero whenever every coupling is nonzero. A matrix with a nonvanishing (N−1) × (N−1) minor has rank at least N−1. Rank is invariant under relabelling, so the natural order is a convenience and not a hypothesis. ∎

Gate **G2** carries this symbolically for both Δ and N = 2 to 6, comparing exactly, with no tolerance anywhere in it.

## (c) The corollary, and why no EP instrument appears here

**Corollary.** On a path edge block with all `J_b ≠ 0`, **any** repeated eigenvalue is defective, and an eigenvalue of algebraic multiplicity m heads a single Jordan block of size m.

Geometric multiplicity one with algebraic multiplicity m ≥ 2 leaves room for nothing else. There is no diabolic alternative to exclude, which is why neither this proof nor its gate invokes `EpCharacter`, a Petermann factor, a phase rigidity, or any other defective-versus-diabolic discriminator: Lemma A closes by an exact rank argument what those instruments would have to decide numerically. The repo's history is why that matters, the real-axis Petermann magnitudes having been retracted as grid artifacts on 2026-06-21. Gate **G3** measures the rank statement where it can bite, at the EPs themselves, and prints the minimum eigenvalue gap over random draws to show why a random λ cannot test it.

## (d) N = 2 in closed form

With `M = [[−2γ₀, 2iJ], [2iJ, −2γ₁]]` the discriminant of the characteristic polynomial is

  **disc = 4·(γ₁ − γ₀)² − 16·J²,**

so the eigenvalues coalesce exactly on `|γ₁ − γ₀| = 2J`. On the branch `γ₁ = γ₀ + 2J` the doubled eigenvalue is `tr M / 2`, the nilpotent `K = M − (tr M/2)·I` equals `[[2J, 2iJ], [2iJ, −2J]]`, which is nonzero, and `K² = 0`. A genuine 2 × 2 Jordan block at real positive rates and a real coupling. The two Δ books differ at N = 2 by the scalar `−2iJ·I`, so the condition is the same in both. Gate **G4** carries all of it symbolically.

The repo already read this point as a Jordan point: `BlockSpectrumWitnessTests.cs:517-535` calls it *"an exact coalescence (a defective 2x2)"* and sqrt-scales its tolerance accordingly, and `OpenArcsRegistry.cs:3968` says the same. Those are in the C# spin book, where `J_gate = 4·J_proof`. What §(d) adds is the symbolic certificate for a reading that was asserted and corroborated but never proved.

## (e) Existence, and the codimension that makes it cheap

**On the XY chain the discriminant is real.** A path is bipartite, so with `S = diag((−1)^k)` one has `S·conj(M)·S⁻¹ = M` **exactly** at Δ = 0, since `S·A_J·S⁻¹ = −A_J` while S leaves the rate diagonal alone. M is therefore similar to its own conjugate, its characteristic polynomial has real coefficients, and the discriminant is real. Gate **G5** asserts the similarity against literal zero, and asserts that the same similarity **fails** at Δ = 1, where the degree diagonal breaks it (`S·D_J·S⁻¹ = D_J`, not `−D_J`).

**So the EP set has codimension one.** `disc = 0` is one real equation in the rate space rather than two, and its solution set is a curve rather than a set of isolated points. Existence then needs nothing more than a **sign change**: a real continuous function taking opposite signs at the ends of an interval has a zero between them. No topological degree, no radius, no point count, and no way for a near-coalescence to masquerade as a zero. Gate **G6** asserts the two signs, that the imaginary part is at the float noise, and that the bracketed root is a coalescence.

Seven witnesses, all at the canonical **J = 0.075** (`docs/Q_REGIME_ANCHORS.md`, the Q = 1.5 anchor) with the canonical base rate γ₀ = 0.05, each reached by turning **one** rate:

| Witness | γ profile at the root | Re λ | rate window | max/min |
|---|---|---|---|---|
| N = 4 | [0.05, 0.220959, 0.07, 0.09] | −0.21236 | (−0.4419, −0.1000) | 4.42 |
| N = 5 (a) | [0.05, 0.169674, 0.343164, 0.07, 0.09] | −0.27044 | (−0.6863, −0.1000) | 6.86 |
| N = 5 (b) | [0.05, 0.195675, 0.343164, 0.07, 0.09] | −0.47648 | (−0.6863, −0.1000) | 6.86 |
| N = 5 (c) | [0.05, 0.493074, 0.343164, 0.07, 0.09] | −0.80437 | (−0.9861, −0.1000) | 9.86 |
| N = 6 (a) | [0.05, 0.230093, 0.191302, 0.07, 0.09, 0.11] | −0.18362 | (−0.4602, −0.1000) | 4.60 |
| N = 6 (b) | [0.05, 0.342699, 0.191302, 0.07, 0.09, 0.11] | −0.47842 | (−0.6854, −0.1000) | 6.85 |
| N = 6 (c) | [0.05, 0.846400, 0.191302, 0.07, 0.09, 0.11] | −0.27912 | (−1.6928, −0.1000) | 16.93 |

The gate stores **brackets**, not these coordinates: the certificate is the sign change across the bracket and the root is found from it, so the printed profiles are outputs rather than inputs.

The block also obeys `M(γ, J) = J·M(γ/J, 1)` exactly (gate **G11**), so the EP condition depends on the ratio γ/J alone, and the table is a statement about that ratio and not about a scale.

**What this would cost on a device, in the units that matter.** Codimension one means one knob and not two, which removes the fine-tuning objection. What remains is the **contrast**: the mildest witness needs a rate ratio of 4.4 across four sites at Q = 1.5. A factor of 4.4 in T₂ between neighbouring qubits is ordinary; a Q of order one is not. And the claim is emphatically **not** supported by pointing at the IBM Torino profile, whose 0.134 is a **sacrifice-zone** qubit, deliberately degraded and roughly thirty times the device median; reading it as "the kind of spread a real device delivers" was an error in an earlier draft of this section. Whether the regime is reachable is open, and this document does not claim that it is.

## (f) What falls, and what survives

**The Edge lemma's conclusion is false under a rate profile.** §(e) contradicts it at N = 4, 5 and 6, on the operator it is about. Its premise fails first, `A = −2γ·I` being scalar only at uniform γ, but a failed premise merely removes a guarantee; §(e) removes the conclusion.

**The window-edge lemma survives, and it does the confining.** It needs only that the Hermitian part be Hermitian, which holds for any profile, and it forbids a defective eigenvalue only **at an edge** of the rate window `[−2γ_max, −2γ_min]`. Bendixson together with that lemma give, as a theorem and with no measurement:

> A defective eigenvalue on a path edge block **must** sit strictly inside `[−2γ_max, −2γ_min]`.

"Must", not "anywhere": the witnesses are seven points, and nothing here shows the interior is exhausted. Gate **G7** corroborates the containment; it does not establish it. At uniform γ the interval collapses to a point, there is no interior, and the old conclusion returns as the zero-width case, which is how `PROOF_CODIM1_BY_ADDITIVITY.md:115` already describes the relation. The window statement is not new either: `hypotheses/INHERITED_RULES_AND_THE_OWN.md:88` states it for this block, Tier 3 labelled.

Of the two lemmas it is the weaker structural one that fails and the stronger derived one that holds. To be fair to the source: the Edge lemma carries its scope inside its own symbol (`A = −2γ·I`), and its proof runs at uniform γ = 1 throughout. What went wrong is that its **verbatim string was transported** into entries and claims scoped to arbitrary rates, and those are what §(i) lists.

## (g) Why uniform γ never showed this

Lemma A applies to the Hermitian Jacobi matrix `A_J` as well, so a path's hopping spectrum is **simple for any bond profile**. Hence at uniform γ the block's spectrum is simple: there is no repeated eigenvalue at all, and the Edge lemma's conclusion holds there for a reason that has nothing to do with normality. Gate **G9** measures 2700 random J profiles at uniform γ and never finds a degeneracy, and confirms that M is normal there, with a commutator that is exactly zero.

An earlier draft argued this from `2J(1 − cos(mπ/N))`, which is F2's closed form and holds only at **uniform J**, so it did not cover the case it was invoked for. The Jacobi argument is both correct and free, being Lemma A again.

This gives the statement its final shape:

> On a path edge block, the per-site RATE profile is what can create a degeneracy at all, and Lemma A then leaves that degeneracy no form but a Jordan block.

A bond profile cannot do it. A rate profile can, and once it has, there is no choice left.

## (h) Scope, and three neighbours this is not

**Lemma A is about paths.** Gate **G10** is the control that isolates that hypothesis rather than a convenient one: on a **star with a rate profile**, M is genuinely non-normal, `‖[M, M†]‖` of order five to seven, and the repeated eigenvalue is nonetheless **semisimple**, nullity 2. So non-normality does not force defectiveness; tridiagonality does. An earlier version used a star at uniform γ, where M is normal and semisimplicity is forced by normality alone, so tridiagonality never entered and the fixture could not have failed.

**Not the frozen divisor.** `PROOF_R90_FROZEN_DIVISOR.md` section 9 exhibits Jordan blocks at explicit (γ profile, J) pairs, including a size-3 block at N = 3. That is `M̃` on the **(1,1)** block on the R₉₀ locus, and it is tuned in **J** at a fixed profile, where this document tunes a **rate** at fixed J.

**Not a contradiction of the atmosphere null result.** `simulations/atmosphere_cluster_ep.py` scans a γ-profile shape toward a coalescence and finds it semisimple. It is the **same geometry**, an open path at the same canonical J and γ₀, so "different geometry" is not the distinction, and an earlier draft was wrong to say so. Two things separate them and the second is decisive. Its profile is **palindromic**, and `ATMOSPHERE_CLUSTER.md` Finding 4 records that the spatial mirror commutes with L exactly then, which protects the crossing, while the witnesses here are non-palindromic. And its coalescence carries geometric multiplicity 32, which Lemma A forbids on any path edge block at any profile, so whatever it reports cannot be this block.

**Not the (1,1) Haken-Strobl block.** Four experiment files call a *"single-excitation block"* defective: `COHERENCE_HORIZON_EP_SENSOR_DEBATE.md`, `FOLD_AND_CUSP_TWO_SEAMS.md`, `THE_FLOW_BETWEEN_TWO_SINGULARITIES.md`, `THE_HUB_KILLS_THE_HORIZON.md`. That is the (1,1) block, of dimension N² and rate 4γ, not this one, of dimension N and rate 2γ. The alias is a genuine trap and the two must not be merged.

**Not a statement about the full spectrum.** The block is closed, so its eigenvalues are eigenvalues of the full generator; nothing here says this EP is the slowest mode, or that it is visible in any observable.

## (i) The inventory for the fencing pass

The fencing is **not** done in this document, and that is deliberate: the arc that owns these sites asks for *"its own pass with its own review rounds, not as a rider"*, and its own second inventory item names a further class this document has not touched. What follows is the starting inventory, and it is known to be incomplete.

The two sites where the Edge lemma is a **premise** rather than a restatement. §(f)'s window does not repair these, because the band cap needs no-Jordan and not a location:

| Site | What it does |
|---|---|
| `PROOF_CODIM1_BY_ADDITIVITY.md:101` | uses the lemma to cap the band chain |
| `PROOF_CODIM1_BY_ADDITIVITY.md:127` | the generalized form, "only the boundary blocks' normality forbids Jordan blocks outright" |

The restatements, each asserting unconditionally what §(e) contradicts:

| Site | Note |
|---|---|
| `PROOF_CODIM1_BY_ADDITIVITY.md:105` | the lemma itself, the source |
| `docs/ANALYTICAL_FORMULAS.md:5170` (F125) | verbatim, inside an entry scoped at `:5152` to "arbitrary rates γ_j", with the proof's units clause dropped |
| `SpectatorIntertwinerClaim.cs:49` and `:126` | in a claim that advertises "for any rates γ_j" in four places |
| `MultiSectorMonodromyVerdictClaim.cs` | five sites; the file contains no occurrence of γ at all, so the restriction is not merely unstated but unstatable in its vocabulary |
| `SectorBraidWitness.cs` | three sites; γ = 1 is hard-coded in its builder and never said in prose |
| `SpectatorIntertwinerGateTests.cs:572-598` | the only site that EXECUTES the implication, passing legitimately at γ = 1 |
| `StructuralCeilingClaim.cs:60` | the claim string dropped the word "uniform" its own docstring carries; this claim asserts a rate and a bound, which a profile breaks by a different route, so it is not a member of the same class |

The fence to write is not "uniform γ assumed" alone, which would leave a reader thinking the profile case is merely unproven. It is §(f)'s window: under a profile the conclusion is false, and what replaces it is a location statement. `BlockSpectrumWitness.cs:33-36` is the model wording already in the repo, with one caution: that file writes the generator in a different J book and as the conjugate partner, and F152 says outright *"never copy a sign or an imaginary value between the two"*.
