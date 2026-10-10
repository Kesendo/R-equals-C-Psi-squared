# The singlet under a shared field: correlation moves a pure spin pair's first loss by its spin-spin covariance, and a field shared in every axis never touches the singlet

**Date:** 2026-10-10
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Verifier:** [`singlet_shared_field.py`](../simulations/singlet_shared_field.py) (exact propagation against every closed form below; all green)
**Grew from:** a pointer an NMR expert gave in answer to our question about internal correlation, recorded on the banner of [Internal and External Observers](../docs/historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md): TROSY, where two relaxation mechanisms interfere because the same molecular motions drive both, and long-lived nuclear singlet states as a related idea.
**A word fence:** a *shared* field acts on both sites at once, one jump operator for the pair; the repository's other words for it are *collective* ([Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §2, [Noise Robustness](NOISE_ROBUSTNESS.md) Q1) and *correlated* ([Quantum Sonar](QUANTUM_SONAR.md)). η here is the correlation coefficient of the two sites' fields; it is not the η-pairing of the frozen-band proofs, not the Glossary's watcher factor η(r), and *correlation* here is not the retired tool's correlation bridge. *Covariance* is the state's spin-spin covariance Cov_a; the field's own covariance, the off-diagonal ηγ, appears in the price-pair page's book as c_ij = 2ηγ. *Dark* means annihilated by a jump.

## What this is about

Two spins that dephase are usually treated as if each sat in its own weather. Keep the weather at each spin exactly as it is and change only how alike the two spins' weather is: a correlation η, from −1 (opposite) through 0 (independent) to 1 (the same field at both). Each spin's own dephasing is the same at every η; the pair is what can tell.

For a pure pair at the first moment, how much correlation changes its loss is a single number: how its two spins co-vary, summed over the noise axes. Where that spin-spin covariance is negative, a positively correlated field starts out doing less harm than independent ones; where it is positive, more; where it is zero, as for a product state, the first loss does not depend on η; an anti-correlated field reverses the two. Later the bond and the decay move the covariance, and that first ordering need not last. The singlet is the extreme: its covariance is as negative as it gets, and a field shared in every direction cannot touch it at all, because the singlet has total spin zero and a shared field acts on the pair through its total spin. The three triplets are the opposite extreme and, in the field shared in every direction, cross the quarter boundary sooner the more it is shared (in the concurrence book).

On the quarter boundary this reads simply. The singlet's CΨ starts at ⅓ (in either book), above ¼, and falls through ¼ after a dose set only by the part of the field it can feel: the crossing time stretches as 1/(1 − η); under z alone the dose is the concurrence constant of [F14](../docs/ANALYTICAL_FORMULAS.md#f14-k-invariance-tier-2-lindblad-scaling). At η = 1 it never crosses.

## Abstract

Two qubits, Heisenberg bond H = J(XX + YY + ZZ). For every axis a ∈ S (S = {z} or {x, y, z}) each site feels σ_a at rate γ and the two sites' fields have correlation η ∈ [−1, 1]: per axis the Kossakowski matrix is γ[[1, η], [η, 1]], realised by a shared jump √(|η|γ)(σ_a⊗1 ± 1⊗σ_a) and local jumps √((1 − |η|)γ)σ_a. A pure state loses fidelity initially at γ Σ_{a∈S} [Var A_a + Var B_a + 2η Cov_a], A_a = σ_a⊗1, B_a = 1⊗σ_a, so the sign of Σ_a Cov_a decides whether correlation lowers or raises that initial loss. On the Bell states Var = 1 and Cov_a = s_a, their parity under σ_a⊗σ_a, giving 2γ Σ_a(1 + s_aη): 6γ(1 − η) for the singlet under x, y, z and (6 + 2η)γ for each triplet. The shared sum jump annihilates exactly the Bell states with s_a = −1, so only the singlet is dark to all three. The singlet's trajectory is the uncorrelated one at the dose (1 − η)γt; in the concurrence book (concurrence × l₁/3) its CΨ crosses ¼ at t× = K/((1 − η)γ), K = ln(4/3)/8 under z (F14's concurrence constant) and ln(6/(1 + √19))/8 under x, y, z. At η = 1 under x, y, z the generator's kernel is span{1, |s⟩⟨s|}. Exact propagation meets every closed-form trajectory to about 10⁻¹⁵ and every crossing time to 10⁻¹⁴.

## What the repo already holds

Swept 2026-10-10 by a survey agent and two reviewers, stores by name.

- `docs/ANALYTICAL_FORMULAS.md`: [F14](../docs/ANALYTICAL_FORMULAS.md#f14-k-invariance-tier-2-lindblad-scaling), K = γ·t_cross constant within a fixed book, shown for Bell+ under local Z where H cannot touch the state and claiming no universality across states or channels, with the concurrence constant ln(4/3)/8, which is this page's crossing at η = 0 under z (the word dose is this page's); the Absorption Theorem's registry paragraph names collective dephasing. No entry treats a correlation coefficient between sites.
- `docs/proofs/`: [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) §2, the correlated-jump algebra and the collective law (a coherence pays 2γ(Δpopcount)² under one jump √γΣ_kZ_k); [PROOF_MONOTONICITY_CPSI](../docs/proofs/PROOF_MONOTONICITY_CPSI.md) Test B, Bell+ dark to the anti-correlated Z₁ − Z₂; [PROOF_ASYMPTOTIC_SECTOR_PROJECTION](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md), the fixed-point algebra as decoherence-free subalgebra and commutant of H, the frame of the η = 1 kernel.
- `experiments/`: [Noise Robustness](NOISE_ROBUSTNESS.md) Q1, [Decoherence Relativity](DECOHERENCE_RELATIVITY.md) and [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md), Bell+ under the additive collective jump, whose |00⟩⟨11| coherence decays at 8γ (here Φ+ at η = 1 under z, whose fidelity loss starts at 4γ, half the coherence rate); [Mathematical Findings](MATHEMATICAL_FINDINGS.md) §9, the split γ_A against γ_B at a fixed sum, the other axis of the same 2×2 Kossakowski matrix, decided by the Z⊗Z parity sectors (its Bell+ baseline crossing 0.719 at γ_total = 0.1 is this page's η = 0.5 z value, the same dose); [When Psi Matters](WHEN_PSI_MATTERS.md), the Werner-line values of CΨ that this page's x, y, z trajectory runs along; [Quantum Sonar](QUANTUM_SONAR.md), a correlated-bath knob on another observable, normalised so that a single site's own rate falls to (1 − η/2)γ, unlike here; [Operator Feedback](OPERATOR_FEEDBACK.md), Bell states decoherence-free under σ_x⊗σ_x as eigenstates of the jump; [PRICE_PAIR_HARDWARE_PREDICTION](PRICE_PAIR_HARDWARE_PREDICTION.md), the field-covariance form of the same Kossakowski structure, correlated dephasing adding 2·Cov_ij to a decay rate and the aligned pair coherence paying γ_i + γ_j + 2c_ij, and on hardware a 5.1σ covariance that did not recur, beside an anti-correlated reading re-attributed to coherent ZZ; [ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION](ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md), the bond jump D[Z_a + Z_b] with an excess of 4d on the pair it acts on and none on another pair, converting c_ij = 2d. These two hold the z axis of this page's Bell rule in coherence form.
- `docs/carbon/`: [Painter Alternation NMR Bridge](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md), the repository's other TROSY page, which fences TROSY as needing a molecular Hamiltonian, bath and readout. `docs/quantum/`: [Schrödinger's Cat Translated](../docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md), superdecoherence under common-mode dephasing, the N-qubit face of the collective law.
- `hypotheses/`: [The Other Side](../hypotheses/THE_OTHER_SIDE.md) §21 tests a ZZ collective dephasing across a bridge. `docs/CAUGHT_ERRORS.md`: the entry on the retired tool's collective σ_z ratio. `docs/GLOSSARY.md`: the 0.036/γ crossing in the concurrence book, this page's η = 0 z value, the rule that the CΨ book be named at the point of use, and the table of spent senses of *dark*, which lacks this page's sense, already used by Test B. `review/OPEN_QUESTIONS_INDEX.md`: OQ-146, the collective z jump, resolved for Bell+, other states not classified.
- `fw.Confirmations` / ConfirmationsRegistry: the price-pair entry's common-mode decoding. The Claim graph: `AbsorptionTheoremClaim` names collective dephasing, with no correlation coefficient. [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md) points to the literature name, decoherence-free subspaces. The OpenArcs registry, the Diagnostics witnesses, `compute/MirrorWorld/`: nothing on a correlation between sites or on singlet protection.
- Not in the repository before this page: the covariance law for a pure state's initial fidelity loss, the Bell-state rule across the x and y axes and the three-axis sum (its z axis is in the two hardware pages above in coherence form), the singlet's crossing time as a function of η, and the triplets' earlier crossing.

## The model

Per axis a, each site's own dephasing is fixed at γ and the correlation η sets the off-diagonal of the Kossakowski matrix γ[[1, η], [η, 1]]. For η ≥ 0 the jumps are √(ηγ)(σ_a⊗1 + 1⊗σ_a) and √((1 − η)γ)σ_a on each site; for η < 0 the shared jump carries the difference σ_a⊗1 − 1⊗σ_a with weight |η|. Traced over the partner, the cross terms cancel and every site dephases at γ for every η and every state (checked on random mixed states). The shared sum σ_a⊗1 + 1⊗σ_a is twice the total spin along a.

## The covariance law

For a pure state |v⟩ and Hermitian jumps the initial fidelity loss is Σ_J (⟨J²⟩ − ⟨J⟩²). With the Kossakowski matrix above this is

  γ Σ_{a∈S} [Var A_a + Var B_a + 2η Cov_a],   Cov_a = ⟨A_aB_a⟩ − ⟨A_a⟩⟨B_a⟩.

Correlation enters only through the summed covariance: at t = 0 and for η > 0 a negative one loses less, a positive one more (an anti-correlated field reverses the two), and a product state's initial loss does not depend on η. This is a statement about the first moment; the bond and the decay move the covariance, and at finite times a state with negative initial covariance can end up with a lower fidelity under correlation than without it (checked). Among random pure states both signs occur. A state need not be dark to be spared: cos 0.3|01⟩ − sin 0.3|10⟩ is gripped by the shared field and still loses less, 0.464 → 0.319 → 0.174 at η = 0, 0.5, 1 (γ = 0.1).

## The Bell states

A Bell state has variance 1 for every σ_a on either site and covariance s_a, its parity under σ_a⊗σ_a, so its loss is

  2γ Σ_{a∈S} (1 + s_aη):

each axis costs 2γ(1 − η) where the parity is odd and 2γ(1 + η) where it is even. The sum jump annihilates a Bell state where s_a = −1, the difference jump where s_a = +1, since (σ_a⊗1 + s'·1⊗σ_a)|v⟩ = σ_a⊗1(1 + s's_a)|v⟩. The singlet is odd in all three axes; each triplet is odd in exactly one (Ψ+ in z, Φ+ in y, Φ− in x). Under a field correlated alike in all three axes the singlet loses 6γ(1 − η) and every triplet (6 + 2η)γ, for every η. Under z alone Ψ+ is the singlet's twin, spared alike. With one η for all axes, only the singlet is dark to every shared jump (η = 1 under x, y, z); under z alone Ψ+ is dark at η = 1 and Φ± at η = −1, the darkness of Bell+ that Test B names.

## The singlet's trajectory and its crossing

Write the Kossakowski matrix as (1 − η)·1 + η·[[1, 1], [1, 1]]; the dissipator is then (1 − η) times the uncorrelated one plus η times the dissipator of the shared sum jump at rate γ, for every sign of η. The singlet's trajectory stays on states built from eigenprojectors of the bond: under x, y, z the Werner line, span{|s⟩⟨s|, 1}; under z the span of |s⟩⟨s| and |Ψ+⟩⟨Ψ+|. The uncorrelated dissipator keeps these spans (the Werner weight decays at 8γ, the z coherence at 4γ), and the bond commutes with every such state, so it acts as zero on the trajectory (F14's Hamiltonian-dead condition), and the shared sum jump annihilates |s⟩ (and, under z, |Ψ+⟩) and leaves the identity fixed, so it acts as zero there as well. So the whole trajectory is the uncorrelated one at the dose (1 − η)γt, ρ_η(t) = ρ₀((1 − η)t):

**z only.** ρ(t) = ½(|01⟩⟨01| + |10⟩⟨10|) − ½f(|01⟩⟨10| + |10⟩⟨01|), f = e^(−4(1−η)γt). Concurrence f and l₁ = f, so in the concurrence book CΨ = f²/3.

**x, y, z.** ρ(t) = p|s⟩⟨s| + (1 − p)·1/4, p = e^(−8(1−η)γt). Concurrence max(0, (3p − 1)/2) and l₁ = p, so CΨ = p(3p − 1)/6.

CΨ = ¼ is met at

  t× = K / ((1 − η)γ),   K = ln(4/3)/8 (z),   K = ln(6/(1 + √19))/8 (x, y, z),

the first K being F14's concurrence constant for Bell+ under local z; here the same dose is read at the only rate the singlet can feel. At γ = 0.1: 0.360 and 0.141 at η = 0, 0.719 and 0.283 at η = 0.5, 3.596 and 1.413 at η = 0.9. The 1/(1 − η) holds in either CΨ book, since the whole state rescales in time; the constants are the concurrence book's. Under x, y, z the triplets go the other way: in the concurrence book a triplet crosses at 0.141, 0.129 and 0.119 at η = 0, 0.5 and 1.

At η = 1 the singlet is a fixed point of the generator. Under x, y, z the kernel is span{1, |s⟩⟨s|}, so every Werner state is stationary and the singlet is the only stationary pure state: the decoherence-free subspace of collective decoherence, in the literature's name.

## What this does and does not say

It says that a correlated field sets a pair's first loss through the pair's spin-spin covariance, axis by axis, and that the singlet, whose covariance is −1 in every axis, is the state a field shared in every direction leaves alone. Whether the singlet crosses ¼ is decided by η = 1 alone; when it crosses is set by the uncorrelated rate (1 − η)γ.

It does not model TROSY, which is a cross-correlation between two different relaxation mechanisms of one spin pair (dipole-dipole and chemical shift anisotropy) and selects one multiplet line; the shared jump here is an analogue of the same kind, a Kossakowski off-diagonal, between two sites. It does not model long-lived nuclear singlet states; it shares with them only that the singlet is blind to whatever acts on both spins alike. It does not model a radical pair, whose relaxation runs through many hyperfine couplings and g-anisotropy. And it does not say how large η is anywhere. The expert found a TROSY-like cancellation hard to imagine in cryptochrome, whose many hyperfine couplings are, in a dynamic protein, probably little correlated; in this model's terms that would mean an η near zero and a singlet that gains almost nothing, but that is our mapping, not a measurement.

## Reproduction

`python simulations/singlet_shared_field.py` (a few seconds): the per-site rate on random mixed states, the covariance law on random pure states (both signs, and a finite-time reversal), the Bell-state rule over η ∈ {−1, −0.5, 0, 0.5, 0.9, 1}, the annihilation table, every closed-form trajectory and the dose rescaling at t = 0.5, 2, 7, every crossing root-found on exact propagation, the η = 1 kernel, and the triplet crossings.
