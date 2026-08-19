# The Drain Has No Chirality of Its Own

**Status:** Tier 3, the record of a refuted direction. The reframe this document was written to propose ("the F112 polarity asymmetry is the chirality of the un-recycled drain") was killed the same day by its own gate-first probe, `simulations/f112_nojump_cancellation_gate.py`. Directions 1 and 2 are refuted, Direction 3 turned out vacuous rather than tested, and **Directions 4 and 5 were never gated and are still open**. What the probe found instead was already ours: it is F113, registered three and a half weeks earlier.
**Date:** 2026-06-20; corrected and renamed 2026-08-19.
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Origin:** the **generative** pass of the `reviewing-before-it-lands` review workflow, run on the F112 scope-correction (see `docs/CAUGHT_ERRORS.md` 2026-06-20). The defensive wing fixed the overclaim; this is what the sharpening opened.
**Naming note (2026-08-19):** filed as `ASYMMETRY_IS_THE_UNRECYCLED_DRAIN.md`, whose title asserted the claim the document's own gate had already killed. Renamed to what the episode actually established, with the scope in the title: the drain **by itself** carries no chirality. The generator it sits in certainly can.

## The objects, since the rest of the page uses them

- **the asymmetry** := ‖M₊ᵢ‖² − ‖M₋ᵢ‖², where M is the polarity decomposition of a Liouvillian against Π, the order-4 palindrome operator, and M₊ᵢ / M₋ᵢ are its ±i eigenspace parts. Definitions: [`reflections/POLARITY_COORDINATES.md`](../reflections/POLARITY_COORDINATES.md) and the primitive `simulations/framework/diagnostics/polarity_coordinates.py`. Δ below is the same quantity, in the directions' original wording.
- **η** := asymmetry / (‖M₊ᵢ‖² + ‖M₋ᵢ‖²), the scale-free version, in [−1, 1].
- **the drain** := the no-jump part of an amplitude-damping channel, −½{c†c, ρ}; the **jump** is the recycling term c ⊗ c*. Full Lindbladian = commutator + drain + jump; **no-jump generator** (also called the physical generator of PT / gain-loss / post-selection dynamics) = −i(Hρ − ρH†) = commutator + drain. Writing H = A + iB with A and B Hermitian, that generator is −i[A, ρ] + {B, ρ}: A is the ordinary Hamiltonian and iB is the gain/loss part. Every asymmetry quoted below is of the no-jump generator unless the sentence says otherwise.
- **bit_b-homogeneous** (used in the table) is F112's hypothesis on the collapse operators: every Pauli string in c_k has the same #Y + #Z parity.
- Gate labels G1, G2, G3, B1, B2a, B2b below are the FIRST probe's own (`f112_nojump_cancellation_gate.py`, thirty gates, and there a FIRED gate is the finding, so a failure is the result). The second probe cited here, `f112_gain_loss_carrier_check.py`, reuses G0..G9 for unrelated gates that must all pass. The two label sets have nothing to do with each other.

## The refutation

The probe, gate-first, N = 2 and 3, H = open XY chain at J = 1, c = σ⁻ per site at γ = 0.1. Six gates fired, at both N: **G2**, **B1** and **B2a**.

- **The drain has no chirality of its own** (G2): the drain term alone has asymmetry exactly 0, and so does the full no-jump generator of that chain, exactly as the full Lindbladian does. There was nothing for the recycling jump to feed back. Re-checked 2026-08-19 and the claim came back wider than it was made: the drain alone is exactly 0.0 not only for mixed per-site σ⁻/σ⁺ profiles at N = 2, 3, 4 but for random per-site 2×2, full random and non-normal collapse operators too (gate G9 of `simulations/f112_gain_loss_carrier_check.py`). The shape is the reason, and it is the transpose pairing rather than Hermiticity: −½(K ⊗ I + I ⊗ Kᵀ), the same operator lifted to both sides with a transpose, has zero asymmetry for every K whatsoever, and the drain is of that form because K = c†c. Break the pairing, K ⊗ I + I ⊗ conj(K) with non-Hermitian K, and the asymmetry is macroscopic; that broken arm is also the checker's positive control, since a gate that only ever reports zeros has not shown it can see anything (gate G9b). The chirality was never the drain's to have.
- **Direction 2's guess is backwards** (B1): the no-jump asymmetry is not blind to the Hermitian part A. It is a linear functional of **H**, and B1 fired because the value ranged over 8.406 at N = 2 and 9.549 at N = 3 across its six Hamiltonians, while the B-only prediction was 0.
- **Direction 3 was not tested, it was vacuous** (B2a): the cooling/heating sign comparison ran on a configuration where η = 0 on both sides.

## What was actually found, and whose it already was

The asymmetry is the cross-term between the commutator and the drain: a linear functional of H whose only single-Pauli carriers are the single-site Z_l strings, the axis the drain's own c†c = n = (I−Z)/2 sits on, and it is same-site, a Z on site l pairing only with damping on site l. The bond bilinears (XX, YY, ZZ, XY) and the X- and Y-drives are all in its kernel, so a bond-only H gives asymmetry 0 on any topology.

**That is F113**, Tier 1 derived on 2026-05-26, closed form asymmetry = (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l) for H ⊃ Σ_l (ω_l/2)·Z_l, with its own proof and registry entry. The `contrib(Z_l) = −16·γ·4^{N−2}` this document reported as its own finding is that formula at one site with a unit Pauli Z (ω = 2), uniform loss γ, no pump: −1.6 at N = 2 and −6.4 at N = 3 for γ = 0.1, which is what the probe prints. Direction 2 below says, in its own words, *"Must first show it is distinct from F113."* The check was never run; it is not distinct. Because it was never run, the generalisation walked out of here into two committed proofs, `docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md` and `docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md`, where it was repaired on 2026-08-19 (`docs/CAUGHT_ERRORS.md`).

One thing here was genuinely new, and it now lives on the F113 entry rather than in this document: the recycling term c ⊗ c* is asymmetry-inert, so F113's closed form reads the no-jump generator exactly as it reads the full Lindbladian. That is a measurement at N = 2, 3, 4 over the σ⁻/σ⁺ family, not a derivation, and F113's own theorem is stated for Hermitian H.

**What F113 does not cover**, and this document originally read as settled: an anti-Hermitian part B outside the per-site n = (I−Z)/2 family. The boundary is mapped inside that family, not beyond it.

**On the sign.** `contrib(Z_l)` is negative in the pipeline's pairing and positive in the other one; F113's entry pins this. The magnitude and the zero-versus-nonzero verdict are convention-free, the direction is not, and this document originally quoted the signed number without the pin.

## The spine, corrected

Each row holds only for the family named in it. Rows 1 and 2 are F112's own scope; rows 3 to 6 are the amplitude-damping family (per-site σ⁻ and σ⁺), which is OUTSIDE F112's bit_b-homogeneous hypothesis, σ⁻ = (X + iY)/2 being bit_b-mixed. That is why those rows can be nonzero at all. The original table read rows 3 and 4 as one contrast, and that contrast was an artifact of measuring the two rows on different Hamiltonians.

| Object | family | polarity asymmetry ‖M₊ᵢ‖²−‖M₋ᵢ‖² | status |
|---|---|---|---|
| commutator −i[H,·] | any matrix H | 0 | F112, proven, universal N |
| Lindbladian with bit_b-homogeneous c | any matrix H, in the commutator | 0 | F112, proven, universal N |
| drain term alone, −½{c†c, ·} | ANY collapse operator | 0 | measured, N = 2, 3, 4; the transpose-paired shape is the reason |
| full amplitude-damping Lindbladian | bond-only H | 0 | F113 at ω_l = 0 |
| full amplitude-damping Lindbladian | H with single-site Z content | F113's closed form, generically nonzero | F113, derived |
| no-jump generator −i(Hρ−ρH†) = −i[A,ρ]+{B,ρ} | either | the same value as the full Lindbladian | measured, N = 2, 3, 4 |

The same no-jump generator under X- or Y-dephasing instead of amplitude damping is a different question, and it is open: see the Klein-V₄ proof's scope note.

The "mean ≈ 132/270 at N = 2/3 over a fixed-norm random ensemble" that stood in the old row 4 is not reproducible: no specification of that ensemble exists anywhere in the repo, and no probe computes it. It should not be quoted as a measurement, and a mean over random H is not a boundary in any case, since a random H generically carries Z_l content.

## What the probe's passes do and do not show

Six of its thirty gate instances (G1, G3, B2b, each at N = 2 and 3; G_sigma too) compare zero against zero on the bond-only chain: full asymmetry = 0, |full| ≪ |no-jump| with both sides 0, and detailed balance η ≈ 0 where η is 0 regardless. They pass, and the document originally read them as confirmations "at 1e-15". They are consistent with the finding, not evidence for it. What carries the carrier statement is the mechanism block instead: M2 (the functional is linear in H, homogeneous and additive), M3 (it is not identically zero, some strings carry it), and Mscan (a random H reconstructed from the single-Pauli scan, residual 2.89e-15 at N = 2 and 1.33e-15 at N = 3). Since 2026-08-19 there is also `simulations/f112_gain_loss_carrier_check.py`, which gates the closed form against a direct build, gets the no-carrier zeros exactly in all 144 swept configurations, and reads the cancelling ones against a measured noise floor instead of pretending they are exact everywhere.

## The directions, as they were proposed

1. **No-jump cancellation.** *Hyp:* the jump term exactly cancels the drain's asymmetry, so the full Lindbladian is 0 while the no-jump piece is not. REFUTED: on the chain both are 0, and off it they are equal, the jump being inert.
2. **Closed form Δ = f(B), blind to A.** *Hyp:* Δ depends only on the anti-Hermitian part, quadratic in B, with a Π-diagonal kernel. REFUTED, and backwards: it is linear in H and reads A through the single-site-Z moment. Its own precondition, "must first show it is distinct from F113", is where the episode should have stopped.
3. **Circular dichroism, completing M's Poincaré sphere.** *Hyp:* η is a scale-free chirality, positive for cooling, negative for heating, zero at detailed balance. NOT TESTED, and vacuous as posed: on the bond-only chain η = 0 on both sides of the comparison. F113's γ_pump − γ_T1 is what carries the cooling/heating sign, and it needs a Z-carrier to show at all, so the direction is re-runnable on a configuration that has one.
4. **Petermann discriminator across the EP.** *Hyp:* Δ (gain-loss injected, finite) is complementary to K = 1/r² (the Petermann factor, r the eigenvector overlap, divergent at an exceptional point), not equal to it. NEVER GATED. **Caveat added 2026-08-19:** the F86a instrument this direction reaches for was corrected in two halves (`docs/CAUGHT_ERRORS.md`, entries 2026-06-21 and the 2026-07-07 sequel; current verdict in `docs/proofs/PROOF_F86A_EP_MECHANISM.md` §The real-axis EP). The Petermann-K magnitudes of the original real-axis sweep were grid artifacts and are dropped; a real-axis defective seed does exist at every odd N, by F89's nullity count. So the direction stands, but it must be built on the surviving half.
5. **PT "dose not phase".** *Hyp:* for a PT family, η grows smoothly through the PT-breaking threshold, making η and the spectrum complementary PT diagnostics. NEVER GATED. F113 sharpens what to ask: since gain and loss enter with opposite signs, a PT configuration is balanced or broken depending on the field profile it meets, so the question is whether η moves at all as the threshold is crossed at fixed carrier content.

**Honest-empty:** an RMT closed form of the raw "132/270" is not a thing to chase; that number has no ensemble behind it.

## Anchors

`docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md` · `docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md` · `docs/CAUGHT_ERRORS.md` (2026-06-20, and the 2026-08-19 sequel) · `reflections/POLARITY_COORDINATES.md` · F83 (anti-fraction, `docs/ANALYTICAL_FORMULAS.md`) · F84 (`docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md`) · F113 (`docs/ANALYTICAL_FORMULAS.md`, `docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md`) · probes `simulations/f112_nojump_cancellation_gate.py` and `simulations/f112_gain_loss_carrier_check.py`.
