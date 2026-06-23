# The Asymmetry Is the Chirality of the Un-recycled Drain

**Status:** Tier 3 — **Direction 1 (the headline reframe) REFUTED 2026-06-20** by gate-first probe `simulations/f112_nojump_cancellation_gate.py` (controls at 1e-15). The spine reframe ("chirality of the un-recycled drain") is WRONG. See "REFUTED + what was actually found" below; the directions are kept as the failed hypothesis, for the record. This is the generative wing working as designed: the move proposed, the gate disposed.
**Date:** 2026-06-20
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Origin:** the **generative** pass of the `reviewing-before-it-lands` review workflow, run on the F112 scope-correction (see `docs/CAUGHT_ERRORS.md` 2026-06-20). The defensive wing fixed the overclaim; this is what the sharpening *opened*.

## REFUTED 2026-06-20 (gate-first) — and what was actually found

The headline reframe is dead. Probe `simulations/f112_nojump_cancellation_gate.py`, gate-first, controls at 1e-15:
- The un-recycled drain has **no chirality**: `DRAIN-ONLY asym = 0`, and the full no-jump generator of any physical spin chain (bond H + σ⁻ drain) is **balanced** (asym = 0 at N=2,3) — identical to the full Lindbladian. Nothing for the jump to recycle (G2 fired).
- Direction 2 ("Δ = f(B), blind to A") is the **opposite** of true: the no-jump asymmetry is a linear functional of **H**, zero for chain H, nonzero for random Hermitian H.
- Direction 3 (cooling/heating sign) is **vacuous**: η = 0 for the physical chain.

**What is actually true (verified to 1e-15, with a positive control `chain + 0.37·Z₁ = 0.37·contrib(Z₁)`):** the no-jump asymmetry is the **cross-term between the commutator and the drain** — an exact **linear functional of H whose only single-Pauli carriers are the single-site Z_l strings** (a local detuning on the same axis as the drain's c†c = n = (I−Z)/2), closed form `contrib(Z_l) = −16·γ·4^{N−2}`. All NN bond terms {XX,YY,ZZ} are in its kernel, so **every standard topology (chain, ring, star, Heisenberg) gives asym = 0**. The spine's "132/270 for random H" is a **random-H artifact** (random H carries Z_l content); reality of H is irrelevant.

**Consequence for the F112 scope-correction** (`docs/CAUGHT_ERRORS.md`, commit `a4d7eb7`): its "balance FAILS for the physical generator (132/270)" illustration is **misleading** — true for random H, but a physical uniform-damping bond-Hamiltonian chain (the named PT/gain-loss spin systems) has a **balanced** no-jump generator. The core (commutator ≠ physical generator) stands; the illustration needs sharpening. The full boundary for general/structured gain-loss (random B) is open — not yet grounded.

---

## The spine (grounded facts) — the (refuted) original hypothesis, kept for the record

| Object | polarity asymmetry ‖M₊ᵢ‖²−‖M₋ᵢ‖² |
|---|---|
| commutator −i[H,·], any matrix H | 0 — F112 proven, universal N |
| bit_b-homogeneous dephasing | 0 — F112 proven |
| full amplitude-damping Lindbladian (drain + jump c⊗c*) | 0 — measured (`POLARITY_COORDINATES.md`; re-confirm in Dir. 1) |
| physical no-jump generator −i(Hρ−ρH†) = −i[A,ρ]+{B,ρ} (drain, NO jump) | ≠0 — 2026-06-20 (mean\|·\| ≈ 132/270 at N=2/3, fixed-norm random ensemble) |

## The reframe (plausible, to confirm)

The asymmetry is the **chirality of the un-recycled drain** — the gain-loss that the recycling jump term normally feeds back. It is the conditional / no-click / post-selected / PT generator's signature. Structurally it is the **circular Stokes component (V) of M**: Π is the order-4 (90°) operator whose ±i eigenspaces are the two circular senses; F83's anti-fraction is one Stokes axis of M, this asymmetry is the orthogonal (circular) axis. The project built two Stokes parameters of M without yet naming the Poincaré sphere.

## Directions (gate-first; each can fail)

1. **No-jump cancellation (strongest, cheapest).** *Hyp:* the jump term exactly cancels the drain's asymmetry — full Lindbladian = 0, no-jump piece ≠ 0, and they sum to 0; this resurrects the T1-asymmetry that `POLARITY_COORDINATES.md` recorded "refuted" (refuted only for the *full* Lindbladian). *Probe:* N=2,3, c=σ⁻, build full / no-jump / jump-only L, feed `polarity_coordinates_from_L`; assert full=0, no-jump≠0, (no-jump + jump + cross)=0 to machine precision. Touches F84, F82, F113.
2. **Closed form Δ = f(B), blind to A (most likely typed result).** *Hyp:* Δ depends only on the anti-Hermitian part B, is independent of A, quadratic in B, with a kernel diagonal in the Π-basis (only B's bit_b-odd content carries it) — same family as F82/F84. *Probe:* A-blindness (fix B, sweep A → Δ const); quadratic (B→λB → Δ∝λ²); kernel vanishes on bit_b-even strings; fit Δ(N) vs c·4^{N−1}. Must first show it is distinct from F113.
3. **Circular dichroism / complete M's Poincaré sphere (borrowing crown).** *Hyp:* η = (‖M₊ᵢ‖²−‖M₋ᵢ‖²)/(‖M₊ᵢ‖²+‖M₋ᵢ‖²) ∈ [−1,1] is the scale-free chirality; cooling η>0, heating η<0, detailed balance η=0 (racemic), tracking F84's net γ_↓−γ_↑. *Probe:* η of the no-jump generator at cooling/heating/balance, N=2,3; sign(cool)=−sign(heat), η_balance=0, |η|≤1. Makes [[reference_diopter_as_polarity_bridge]] load-bearing — "polarity" was the borrowed polarization term all along.
4. **Petermann discriminator across the EP (highest connectivity).** *Hyp (to refute the naive identity):* Δ (gain-loss injected, finite) ≠ K=1/r² (spectral consequence, diverges at EP) — complementary, not equal. *Probe:* sweep the F86 c=2 doublet through Q_EP=1.5; assert Δ/η finite-smooth while K→∞. Joins the F112-polarity cluster to the F86-EP cluster ([[project_everything_is_one_object]]).
5. **PT "dose not phase."** *Hyp:* for a PT family, η grows smoothly (∝γ²) through the PT-breaking threshold (spectrum-blind), so η and the spectrum are complementary PT diagnostics. *Probe:* PT (+iγ/−iγ) family sweep across threshold; η analytic while K diverges at it.

**Honest-empty:** topology-dependence (likely null if Dir. 2's A-blindness holds; only structured/boundary loss could differ); an RMT closed form of the raw 132/270 (ensemble variance, not a law — the physics is in *structured* B, not random H).

**Strongest single move: Direction 1** — cheap (N=2,3, existing `polarity_coordinates_from_L`), fails three ways, and if it passes it converts today's scope-correction into a positive statement: *the asymmetry is the chirality of the un-recycled drain.*

## Anchors

`docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md` · `docs/CAUGHT_ERRORS.md` (2026-06-20) · `reflections/POLARITY_COORDINATES.md` · F83 (anti-fraction) · F84 (`PROOF_F84_AMPLITUDE_DAMPING.md`) · F113.
