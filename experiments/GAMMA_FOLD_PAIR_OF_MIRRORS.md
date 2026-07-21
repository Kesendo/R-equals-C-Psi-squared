# The gamma fold: the pair of mirrors on the gamma axis

*2026-07-21. The home-side move of the F134/F139 arc: after [PROOF_F139_SEAM_IDENTITY](../docs/proofs/PROOF_F139_SEAM_IDENTITY.md) closed the character-side mechanism (the wall as a Chebyshev divisor), the arc's plan asked whether the pair of mirrors exists at home, on the gamma axis of the project's own Lindblad family. It does, and it is a one-line corollary of committed facts; this note states it, places it against its siblings, and records the adoption into MirrorWorld.*

## The statement

Take the usual family: an XXZ chain under local Z-dephasing with per-site rates γ_l, σ = Σ_l γ_l. Two involutions act on the watching parameter:

- **s, the gain turn**: γ_l ↦ −γ_l (reflection through the unwatched zero; negative rates amplify),
- **s₀, the anti-watch turn**: agreement watched instead of disagreement (the turned rule −2γ(N−k) of [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md), site-resolved: −2·Σ_{l agrees} γ_l).

They chain through one exact identity:

  **L_anti(γ⃗) = L(−γ⃗) − 2σ·Id**,

and the trajectory wears the shift as a scalar veil:

  **ρ_anti(t) = e^(−2σt) · ρ_gain(t)**  (gain = the same seed evolved at −γ⃗).

*Derivation (one line, from committed pieces):* the Hamiltonian leg never sees γ, and on a coherence cell |i⟩⟨j| the turned rate is plain arithmetic: −2·Σ_{agree} γ_l = +2·Σ_{differ} γ_l − 2σ = (rate at −γ⃗) − 2σ. Both inputs are owned: the site-resolved rate −2·Σ_l γ_l·(bit l of i⊕j) is the Absorption Theorem's cost identity, and the turned rule is the Lattice adoption. ∎

On the rate functions r the two turns read s: r ↦ −r and s₀: r ↦ −r − 2σ, so composing gain after the turn (s∘s₀) is the **translation r ↦ r + 2σ**: two mirrors make a translation, ⟨s, s₀⟩ is the infinite dihedral group (σ ≠ 0) with the full price 2σ as its step, and the fixed locus of s₀ is r = −σ, the palindrome center. This is exactly the two-mirror shape [F134](../docs/proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md) carries on the character side (s: μ₁ ↦ −μ₁, s₀: μ₁ ↦ 22−μ₁, translation 22), now sitting on the home γ axis.

## What it is not (the placement)

- **Not the F1 fold.** The palindrome Π L Π⁻¹ = −L − 2σ ([MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)) sends the whole operator to −L − 2σ, in particular flipping the Hamiltonian leg (Π L_H Π⁻¹ = −L_H); the gamma fold keeps H untouched and flips only γ. The two are siblings sharing the price constant 2σ, not the same map.
- **Not R₉₀.** The parameter-side anti-palindromic reshuffle of F91/F92/F93 ([PROOF_F91_GAMMA_NINETY_DEGREES](../docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md) on this axis) is γ_l ↦ 2·γ_avg − γ_{N−1−l}, which preserves σ. The gain turn negates σ. Different involutions on the same axis.
- **The gain world was already sighted, informally.** [THE_OTHER_SIDE](../hypotheses/THE_OTHER_SIDE.md) names the laser regime (Σγ < 0) as the physical realization of the other side. The identity above is its exact algebraic form: the anti-watched world IS the gain world in the price veil. One caution on the label: negated pure dephasing is time-reversed dephasing, not a CP channel like laser gain (that one inverts populations through a jump operator); the identity is plain algebra and does not depend on the physical reading.

Composing with the committed one-sided X^N bridge (the Lattice reading L(t)[i,j] = e(t)[~i,j]) gives the third face: **X^N ρ(t) = e^(−2σt) · [gain evolution of X^N ρ(0)]**, the watched world read through the complement is the gain world wearing the price.

## Verification (the adopted pins)

Adopted into MirrorWorld as `GammaFold` (`compute/MirrorWorld/GammaFold.cs`, run mode `gammafold N`; suite 240 via `GammaFoldTests`):

- the generator identity per cell at machine zero (2.2e−16) for a non-uniform site profile, ZZ on or off;
- the dihedral closure exactly (involution, translation step 2σ);
- the veil law over twin RK4 at 2.6e−8 (dt = 0.02, the RK4 truncation of the scalar shift), with the amplification witness gain/anti novelty = e^(+2σt) (novelty = the summed off-diagonal coherence weight; the law is analytically exact, the twin-RK4 witness matches it to 1e−4; the trace is blind to γ in every world, the growth lives in the coherences);
- the discriminator: the veil against the unflipped (+γ) world misses at 0.057, so the gain flip is load-bearing;
- the X^N cross-dock at 4.5e−8.

## Open: the divisor question this sharpens

F139's lesson was that the wall is a **divisor, not a symmetry**: the reflection law fell out of dividing by the vanishing polynomial of an angle lattice, after every symmetry reading had been closed off. The gamma fold gives the home side the same cast of characters: σ plays the shifted level, the dihedral plays ⟨s, s₀⟩, and the home angle lattice already exists (the single-excitation spectrum 2J·cos(rπ/(N+1)), whose vanishing polynomial is the Chebyshev S_N; the [Niven rationality root](../docs/carbon/F99_NIVEN_COMPLETENESS.md) claim owns its arithmetic). What is missing is the home polynomial that the wall factor divides: the candidates assembled so far are the characteristic polynomial of the effective Liouvillian on the Q = J/γ₀ axis (where a naive similarity factorization already failed, the same false start F134's centroid theorem closed on the character side) and the F89 discriminant factorization disc = C·w^v·A₁·A₂² (a committed home-side "factor divides out exactly"). That hunt is the arc's next round, not this note's claim.
