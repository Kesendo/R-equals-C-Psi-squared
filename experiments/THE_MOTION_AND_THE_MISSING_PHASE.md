# The Motion and the Missing Phase: one lit seat leaves a motion untouched and fades a phase that no small page can see

**Status:** Exploration note. The encoding identities, the F70 invisibility, the negativity and the veil identity B_γ(t) = e^{−2γt}·A_{−γ}(t) of §8.1 are exact algebra assembled from F157, F70, F158 and the gamma fold, at every Δ. The short-time defect coefficients (Δ = 0) are derived in §8; the exact-arithmetic checks §9 reports at N = 5, 7, 9, 11, 13 have no committed artifact and are re-derivable from its recurrence only. The N = 7 long-time run (Δ = 0, t ≤ 6000, two grids) is floating point with its diagnostics beside the data. No hardware run, no F number; C# live-object adoption open, specified in §9.
**Date:** 2026-09-05
**Authors:** Thomas Wicht; Codex (OpenAI): §1 to §12 and the producer; Claude (Anthropic, Fable 5.1): the entry layers and the sweep record
**Producer:** [`missing_phase_long_time.py`](../simulations/missing_phase_long_time.py), with [five tests](../simulations/tests/test_missing_phase_long_time.py) on the output maps and the reduced generators
**Data:** [`simulations/results/missing_phase_long_time/`](../simulations/results/missing_phase_long_time/) (two grids × two defect signs, CSV + JSON)
**Grew from:** F157 and F158 read together. The blind seat says what one lit site cannot touch; the two-end count says the untouched has a partner that pays the whole price. The question was what a single state looks like that carries both at once. Tom's clock hint, that the two blocks of §8.1 are adjoint partners, L_B = −2γ·Id − L_A†, joined the motion to the fading phase and became §11.
**A word fence:** two letters and five words here carry more than one sense. *σ* is the internal state's density matrix, σ_0 and σ(t), throughout; the repository's σ = Σ_l γ_l is written out as Σ_l γ_l on this page, and equals γ. *q* is the copy label in the Abstract and §1 to §6 and the bright weight 1 − p in §10.1. Five words are spent elsewhere in the repository. *Dark* in this note means γ = 0, a chain with no light on it; it is not the glossary's dark, an exact zero in a record, and not the dark state of [The Seat That Cuts](THE_SEAT_THAT_CUTS.md), in its quantum-optics sense a state the dephasing cannot touch, although in that last sense the blind motion D of §1 is dark too. *Mirror* in "breaks the mirror" means the chain's left-right reflection, not F1's palindrome mirror and not the one-sided relation of §8.1 and §11, L_B = −2γ·Id − L_A†, which the detuned chain keeps. *Light* is the arriving γ throughout; in the Absorption Theorem's own words, which §11 borrows, it is also the {X, Y} letter content a mode exposes, and w_c is that content's share at the centre, the proof's light_l(v) at l = c read on the one-excitation block. *Clock* in §11 means a mode's two hands, its decay and its turning rate, the sense of [`Clock.cs`](../compute/MirrorWorld/Clock.cs) and the clock-hand ladder, not the cracked bell's clock, which is a time a crack sets. *Block* is mostly a plain matrix block here, the A and B of the 2N × 2N reduced form in §8.1 and the centre/outer partition of h in §10.1; in the joint-popcount sense of Liouville space (bra popcount first), A is the (1,1) and the (N−1,N−1) block and B the (1,N−1) and the (N−1,1), the pairing F1's corollary Π² = X^{⊗N} names.

## What this is about

Put a chain of spins in the dark and let light fall on one seat only, the middle one. That is all dephasing at a single site is, in this repository's words: the chain stands in [light](../docs/quantum/DEPHASING_TRANSLATED.md) at one place, and what arrives there tells two possibilities apart only where they differ at that seat, and passes through where they agree. Now start one excitation moving along the chain in the pattern whose left and right halves are always in anti-phase. Such a motion has a node at the middle. The seat under the light is never occupied, so the light finds nothing there to tell apart, and the motion runs as if the whole chain were dark: a mode with a node under the damping pad, the picture [the Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)'s preface sets aside for uniform rates and its Theorem 2 brings back exactly for site-dependent ones. On an odd chain the middle seat is [the blindest seat there is](THE_SEAT_THAT_CUTS.md), because the two halves it leaves behind agree about everything.

Then take the same motion in a second copy of the chain, with every spin flipped. Flipping every spin at once is a symmetry of the hopping, so the flipped copy moves exactly as the first; but where the first keeps the middle seat empty, the flipped one keeps it full. Put the chain into a single state that holds both copies, the way [the cat](../docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md) holds both its branches, except that the cat's two branches stand still and these two keep moving. The motion is shared. What is new is a phase between the two copies, and that phase is exactly what the light at the middle CAN tell apart, since at every instant one version shows it an empty seat and the other a full one. So the phase between the copies fades at the plain rate one lit seat charges, while the motion inside each copy goes on untouched. One state, two fates.

The question of this note is where that fading phase is. Look at any one spin, any two, any four of seven: the lit chain and the dark chain give identical answers, although entanglement between the middle and the rest is really being lost. The phase lives in a joint pattern of the middle with all but one mirror pair of the outer spins, N−2 of the N at once, and no [page](../docs/quantum/SPOOKY_ACTION_TRANSLATED.md), the repository's word for what a chosen handful of sites shows on its own ([`Marginal.cs`](../compute/MirrorWorld/Marginal.cs)), carries a trace of it below N−2 spins. That is the missing phase of the title: missing from every small page, not from the state. Nor is it lost for good. A reversible re-encoding, the middle spin flipping every outer spin, brings the whole phase onto the middle spin alone, where one measurement reads it and the motion on the outer spins is left undisturbed; reading every spin at once in the flipped basis gives the same number and destroys the motion in doing so. A correlation order and a number of measurement settings are two different things.

The second half of this note breaks the mirror on purpose. One end bond is detuned, the two halves no longer agree about everything, and the light can now reach the motion: the middle seat is occupied a little, and the light begins to write into the motion itself. We follow how soon each way of looking notices, and the answer counts hops: the more hops from an end to the middle, so the longer the chain, the later any reading feels the light's touch, and the decoded outer spins feel it no later than any two-spin page does. Then we follow seven spins for a long time and look at what the light leaves standing: one blind direction that even the detuned chain keeps when the middle seat is odd-indexed, as it is at seven spins, forced by a parity rather than by the mirror; a plateau the decoded motion settles onto; and a pair of clocks, the slowest fading mode and its partner in the fading phase, turning at the same rate and splitting one fixed sum of decays, 2γ, between them. So the note's three questions, what survives, what is visible, and what remains after reading, get three answers about one state. What survives is the motion while the mirror holds; once it is broken, one parity-forced standing weight, and nothing that moves. What is visible on any small page is nothing of the phase. What remains after reading depends on how one reads: the decoder leaves the motion standing, the all-spin readout destroys it. (γ as light is this repository's Tier-4 canvas, [GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md), and nothing below leans on it.)

## Abstract

On an odd open XY/XXZ chain, N = 2m+1 ≥ 5 with reflection-symmetric bonds and Z-dephasing at rate γ on the centre site c only, the reflection-odd single-excitation space D (dimension m, zero centre amplitude: F157's blind seat, forced by the chain's reflection) and its global-flip image FD, F = X^{⊗N}, carry the same internal Hamiltonian h and opposite centre parity. The isometry V: |q⟩ ⊗ |j⟩_in ↦ F^q|d_j⟩ intertwines exactly, HV = V(I ⊗ h), Z_cV = V(Z ⊗ I), FV = V(X ⊗ I). For the encoded product |+⟩⟨+|_q ⊗ σ_0 the trajectory factorises, ρ(t) = V[ρ_q(t) ⊗ σ(t)]V† with ρ_q = ½[[1, η], [η, 1]], η = e^{−2γt}, and σ(t) unitary under h: the purity falls to Tr σ_0² · (1+η²)/2 and the centre/rest negativity is η/2, while by F70 every reduced state on fewer than N−2 sites is independent of η. The first sufficient correlation order N−2 is attained by the weight-(N−2) observable O with V†OV = X ⊗ I, so ⟨O⟩ = ⟨F⟩ = η at Var O = m − η² against Var F = 1 − η². The decoder U = ∏_{l≠c} CNOT(c→l) satisfies U†X_cU = F and maps the family to ρ_q on the physical centre tensor the internal state on the outer sites. In the 2N-dimensional block form ρ = ½W[[A, B], [B, A]]W† (A the one-excitation block, B the cross-copy block, §8.1; in Liouville space A is the (1,1) and the (N−1,N−1) joint-popcount block and B the (1,N−1) and the (N−1,1), bra popcount first) the cross-copy block is the Lattice's one-sided reading (X^N ρ or ρ X^N under the turned rule, [`Lattice.cs`](../compute/MirrorWorld/Lattice.cs) and [the opening law](LATTICE_OPENING_LAW.md)), here a genuine block of a physical state rather than a relabelling of the normal trajectory, and it obeys B_γ(t) = e^{−2γt} A_{−γ}(t), [the gamma fold](GAMMA_FOLD_PAIR_OF_MIRRORS.md)'s veil identity with its Σ_l γ_l equal to γ, equivalently L_B = −2γ·Id − L_A†, the one-sided relation of [the two-end proof](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md); every A mode −α + iω has a B partner −(2γ−α) + iω.

At Δ = 0 and general odd N with one end bond detuned to 1+ε the encoding leaks, p_out = ε²t² + …, and the light reaches the motion. Exact short-time laws at fixed ε and γ: the decoded outer trace distance d_out opens as [√2·m·2^m/(m+1)!]·|ε|γ·t^{m+1}; the best two-site page d_2 opens at t^{m+1} for even m and one order later, t^{m+2}, for odd m, because its leading centre/end coefficient is imaginary and cancels between the copies; the dephasing-induced leakage difference opens as (−1)^m·[2^{2m+1}/(2m+1)!]·ε²γ·t^{2m+1}, so at N = 7, 11, … dephasing first reduces the leakage population. §9 records exact Gaussian-rational checks at N = 5, 7, 9, 11, 13, re-derivable from its recurrence. At N = 7, J = 1, γ = 0.3, ε = ±0.1, through t = 6000 on two grids: the leakage difference crosses sign once, negative to positive (t ≈ 1.25 and 1.51; the producer's search starts at t = 0.5, past the onset's rounding flips), ‖B‖_F is 9.0 × 10⁻¹⁷ and 6.0 × 10⁻¹⁷ by t = 100 while ‖A − A_∞‖_F is still 0.68 and 0.69, one blind vector v with hv = 0 survives the defect (F157's parity-forced kind, forced by the odd size of the two zero-diagonal blocks the seat leaves behind, present because the centre index m = 3 is odd, so at N ≡ 3 mod 4; blind weight p = 441/926 and 361/686), A_∞ = p·vv† + (1−p)·(I − vv†)/(N−1), and d_out plateaus at [5(1−p)/6 + √((5(1−p)/6)² + 4p(1−p))]/2 = 0.76327 and 0.73432 (sampled 0.76324 and 0.73417 at t = 6000). The slowest non-stationary A modes have α = 2γ·w_c with w_c ≈ 0.27 % and 0.22 % (Theorem 2 of the Absorption Theorem, its vector form with site-dependent rates; w_c is its light_l(v) at l = c on the one-excitation block), lifetimes 608 and 742, and B partners of lifetime 1.67.

Not claimed: no hardware run, no C# live object, no F number minted; the sampled ordering d_out > d_2 is a finite observation, not a theorem; the short-time coefficients are derivation checks at five N, not an all-N proof; N = 3 is excluded because the two copies are no longer orthogonal on the outer sites.

## What the repo already held

The repository search followed the primitives through the stores. [The formula registry](../docs/ANALYTICAL_FORMULAS.md) supplies F157's blind space, F70's local selection rule, F158's commuting/anticommuting operator conditions, and F1's corollary that the global flip X^{⊗N} is Π², the square of the palindrome conjugation, pairing the joint-popcount sectors (p_c, p_r) ↔ (N−p_c, N−p_r): it is why the two copies of A, (1,1) and (N−1,N−1), and the two copies of B, (1,N−1) and (N−1,1), coincide in §8.1, while the A ↔ B relation is the one-sided fold below.

In `docs/proofs/`, [the blind-seat proof](../docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md) supplies the invariant Krylov complement, [the local-selection proof](../docs/proofs/PROOF_DELTA_N_SELECTION_RULE.md) supplies the partial-trace argument, and [the two-end proof](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) supplies the one-sided U; [the structural-ceiling proof](../docs/proofs/PROOF_STRUCTURAL_CEILING.md) §2 already uses particle-hole conjugation to carry the (1,1) sector onto (N−1, N−1).

In `experiments/`, [The Blind Site](THE_BLIND_SITE.md) already holds moving blind states, controls with the light away from the centre, and a mirror-broken arm (bond (1,2) detuned to J = 1.6 at N = 11, Heisenberg, the centre filling to a few percent) that its own page keeps as an arm whose effect cannot be attributed to the broken mirror alone, the trusted control being the reflection-keeping double detune; the end-bond detuning here is that arm's XY neighbour with exact coefficients; [The Blind Seat on the Road](THE_BLIND_SEAT_ON_THE_ROAD.md) separates forced from accidental blindness under knobs that keep the chain reflection, and [The Seat That Cuts](THE_SEAT_THAT_CUTS.md) §7 records that the XY book's parity-forced kind survives a signed-profile sweep whole and leaves the Heisenberg half open, whether a mirror-forced and a met blindness part as the palindrome, the chain's reflection symmetry, detunes; §10.1 is that question's XY neighbour, a bond knob that breaks the reflection, under which the centre seat's blind space drops from m to one at every detuning with r² ≠ 1, r = 1+ε, the one that stays, the centre being odd-indexed at N = 7, being the parity-forced kind; [the Lattice opening law](LATTICE_OPENING_LAW.md) holds the one-sided reading X^N ρ under the turned rule for any site-resolved profile (its gate T0; the C# `Lattice` runs the uniform case), which is the rule the B block of §8.1 obeys, there as a reading of the normal trajectory and here as a block of a physical state; [Gamma Fold](GAMMA_FOLD_PAIR_OF_MIRRORS.md) holds the veil identity ρ_anti(t) = e^{−2Σ_l γ_l t}·ρ_gain(t) and its cross-dock with the Lattice reading, X^N ρ(t) = e^{−2Σ_l γ_l t}·[gain run of X^N ρ(0)], written there with σ for Σ_l γ_l, the face the B identity of §8.1 wears at Σ_l γ_l = γ (typed as MirrorWorld `GammaFold`); [the price-pair flight](PRICE_PAIR_HARDWARE_PREDICTION.md) already uses one all-X setting for several subset correlators; [Survivor Flip and Reflection-Odd](SURVIVOR_FLIP_AND_REFLECTION_ODD.md) reads the same two symmetries, X^{⊗N} and the chain reflection, as parities of a different object, the slowest non-stationary density mode under uniform dephasing.

The [glossary](../docs/GLOSSARY.md) separates the dephasing rate from this repository's reading of γ as light, a Tier-4 canvas laid on the same algebra, and [Caught Errors](../docs/CAUGHT_ERRORS.md) records, in the first of its 2026-08-30 entries, that a guard which seeds a fixed point has tested the integrator and not the property.

The [open-arcs ledger](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs) returned, among the arcs that bear on this note: `compressed_density_laws`, named with the siblings below; `gamma_is_the_sender_not_the_watching`, which names the second reading X^N ρ, the B block's nearest neighbour there, as one of three outside objects with no eye in them; `site_resolved_vacuum_block`, the (0,1) block under a γ profile, the nearest sibling of a block under centre-only light; `handshake_decoder` (retired), where the other decoder named below was built; `the_gate_that_does_not_gate`, whose annotation already holds this note's preparation (the reflection-odd single-excitation state at N = 11 is moved by the centre seat at γ = 50 by 6.2 × 10⁻¹⁶, by its neighbour by 0.40) and whose remaining items, the level form of the mediator bridge and the walk-time corollary, this note does not answer; and `the_forced_and_the_met`, which records that F157's count splits by seat, a reflection-fixed seat keeping it at every u and Δ, the knobs there keeping the reflection where the knob here breaks it.

The [hardware registry](../simulations/framework/confirmations.py) returned two flights reading an X string on every qubit, the price-pair flight (all seven subsets of a three-qubit line, its flight page reading them from one all-X setting) and the F120 moment tower on Kingston (an XXX string, |⟨XXX⟩| ≤ 0.02 throughout), and no flight of a flipped copy or a decoder.

Core and Diagnostics contain [the blind-seat claim](../compute/RCPsiSquared.Core/Symmetry/SeatCutBlindnessClaim.cs), [its witness](../compute/RCPsiSquared.Diagnostics/Foundation/SeatCutBlindnessWitness.cs), [the two-end claim](../compute/RCPsiSquared.Core/Symmetry/PalindromeTwoEndCountClaim.cs) and [its witness](../compute/RCPsiSquared.Diagnostics/Foundation/PalindromeTwoEndCountWitness.cs), and [`F89CrossFoldSimilarityClaim`](../compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs) types the antiunitary form of the same one-sided fold; [`XGlobalChargeConjugationPairing`](../compute/RCPsiSquared.Core/SymmetryFamily/XGlobalChargeConjugationPairing.cs) types F1's corollary.

Nearer siblings of the block pair of §8.1: MirrorWorld `Mirror` runs the (1,1)/(1,N−1) fold pair with the price veil at uniform γ past the wall, the A/B pair at price 2Nγ where here the price is 2γ because only the centre is lit; the arc `compressed_density_laws` records the undressed one-sided fold F D F = −D − 2Σ_l γ_l, exact for any Δ and profile; [Palindromic Partner Mode](PALINDROMIC_PARTNER_MODE.md) with F67 and F68 holds, at endpoint dephasing, the bonding mode's palindromic partner and a Bell-pair-like two-branch encoding whose off-diagonal decays at the partner rate, §11's pairing already used as an encoding at another seat; F64 is the single-site rate −Re λ = 2γ_B·|v(B)|², the Absorption Theorem applied to the single-excitation sector, which is what §11 reads; and [Slow-Mode R-Parity](SLOW_MODE_R_PARITY.md) named the R-odd half of the slow-mode landscape as what F86's channel-uniform probe leaves out.

What none of them assemble is the coherent two-copy superposition with a moving internal state as a prepared state (the cat pair of the opening law is its motionless analogue, internal space of dimension one, under uniform light), its negativity η/2, and the CNOT decoder with U†X_cU = X^{⊗N}; the repository's other decoder, [`DefectDecoder.cs`](../compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs), reads defect locations. No hardware run of this construction exists. The ingredients already have homes; this note puts them together.

## 1. The physical system and its two copies

Let N = 2m + 1 ≥ 5, with sites numbered 0 through N−1 and centre c = m. Use Pauli matrices, ℏ = 1, and

    H = Σ_{l=0}^{N−2} J_l (X_l X_{l+1} + Y_l Y_{l+1} + Δ Z_l Z_{l+1}),
    J_l = J_{N−2−l},       J_l and Δ real,
    L(ρ) = −i[H,ρ] + γ(Z_c ρ Z_c − ρ),       γ > 0.

There are no longitudinal fields and no other dissipation channels. The hopping matrix element is 2J_l in this convention. Zero bonds are allowed algebraically, although they may remove the motion of interest.

Write |j⟩ for the computational state with its only excitation at site j. The reflection-odd single-excitation space D has the orthonormal basis

    |d_j⟩ = (|j⟩ − |N−1−j⟩)/√2,       j = 0,…,m−1.

H preserves excitation number and commutes with spatial reflection. Thus D is H-invariant. Every vector in D has zero centre amplitude: Z_c acts as +I there. This is F157's blind seat at the centre.

Let F = X_0 X_1 … X_{N−1}, the global spin flip, the repository's Π², the square of the palindrome conjugation (F1's corollary). Since [H,F] = 0, the space FD has the same internal Hamiltonian. Its states have one hole rather than one excitation, and Z_c acts as −I there. D and FD are orthogonal.

Define an isometry V from a two-level copy label q and an m-dimensional internal coordinate space into the physical chain:

    V(|0⟩_q ⊗ |j⟩_in) = |d_j⟩,
    V(|1⟩_q ⊗ |j⟩_in) = F|d_j⟩.

The label q is an encoded degree of freedom, not yet a separate physical qubit. With h_{jk} = ⟨d_j|H|d_k⟩,

    H V = V(I_2 ⊗ h),
    Z_c V = V(Z ⊗ I_m),
    F V = V(X ⊗ I_m).

These are invariant-subspace identities, not just equal compressed matrix elements. Both H and Z_c keep the state inside the image of V.

For a uniform five-site chain the internal Hamiltonian is

    h = J [[2Δ, 2],
           [2,  0]].

At Δ = 0, |d_0⟩ evolves into cos(2Jt)|d_0⟩ − i sin(2Jt)|d_1⟩. The centre stays empty throughout. The ZZ term changes this internal motion without exposing it to centre dephasing.

## 2. The motion and the fading phase separate exactly

Prepare the encoded product |+⟩⟨+|_q ⊗ σ_0, where |+⟩ = (|0⟩+|1⟩)/√2 and σ_0 is any internal density matrix. Physically this means a coherent superposition of the two copies with the same internal state.

The generator identities give

    ρ(t) = V[ρ_q(t) ⊗ σ(t)]V†,
    ρ_q(t) = ½ [[1, η], [η, 1]],       η = exp(−2γt),
    σ(t) = exp(−iht) σ_0 exp(+iht).

Consequently,

    Tr ρ(t)² = Tr σ_0² · (1 + η²)/2.

For pure σ_0 the global purity falls from 1 towards 1/2 while the internal state remains pure and continues its unitary motion. Protection here belongs to that internal degree of freedom, not to the entire physical state.

The F158 connection is at the operator level: F is [the two-end proof](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)'s one-sided U, with [F, H] = 0 and F Z_c F⁻¹ = −Z_c, so one-sided multiplication by F gives ℒ_F L ℒ_F⁻¹ = −L† − 2Σ_l γ_l, here −L† − 2γ. If |u⟩ and |v⟩ are energy eigenvectors in D, with ω = E_u−E_v, then

    L(|u⟩⟨v|)   = −iω |u⟩⟨v|,
    L(F|u⟩⟨v|)  = (−2γ−iω) F|u⟩⟨v|.

Applying F on one side after taking the adjoint gives the point-reflected eigenvalue −2γ+iω. The damped partners connect the one-excitation and one-hole sectors. They are operator directions, not density matrices by themselves. The physical superposition above makes such cross-sector coherence part of a valid state.

## 3. What every small reduced state misses

The difference between two values of η consists entirely of operators connecting excitation numbers 1 and N−1. Their difference in excitation number is N−2. By F70,

    Tr_{outside S}[ρ(η,t) − ρ(η′,t)] = 0       if |S| < N−2.

Here the internal state σ(t) is the same in both arms. These reduced states may change with time; they simply do not depend on η, hence on γ in this experiment.

| N | Internal dimension m | Every smaller marginal misses the phase | First sufficient correlation order |
|---|---:|---|---:|
| 5 | 2 | one and two sites | 3 |
| 7 | 3 | up to four sites | 5 |
| 9 | 4 | up to six sites | 7 |
| 11 | 5 | up to eight sites | 9 |

The last column is attained, not merely a lower bound; the observable in §5 supplies it. This statement concerns reduced states and their functions. It does not exclude collecting high-order correlations by measuring individual sites, or converting the phase into a local signal through an intervention.

## 4. Entanglement is actually lost

For pure internal input, write its physical state with centre removed as |a(t)⟩. Let |b(t)⟩ = X^{⊗(N−1)}|a(t)⟩ on the outer sites. Their outer excitation numbers are 1 and N−2, so they are orthogonal for N ≥ 5.

Across the centre/rest partition, the state is

    ρ(t) = ½[|0a⟩⟨0a| + |1b⟩⟨1b|
             + η|0a⟩⟨1b| + η|1b⟩⟨0a|].

The partial transpose on the centre has nonzero eigenvalues

    1/2, 1/2, η/2, −η/2.

Thus the negativity, defined as (‖ρ^{T_c}‖₁−1)/2, is η/2. Its decay is a loss of entanglement across this partition even though all the marginals specified in §3 agree between the γ > 0 and γ = 0 arms.

For mixed σ_0, diagonalize it with eigenvalues p_a. Orthogonal internal eigenvectors give orthogonal outer blocks, each carrying p_a times these four eigenvalues. The negativity is still η/2, while the purity carries the factor Tr σ_0² in §2.

N = 3 is excluded deliberately: the outer excitation numbers then coincide, the orthogonality argument fails, and the internal space has only one dimension.

## 5. A correlation order is not a number of measurement settings

Define products with identities at all unlisted sites:

    O = −Σ_{j=0}^{m−1} ∏_{l ∉ {j,N−1−j}} X_l.

Each summand has weight N−2. Its action on |d_j⟩ is −F|d_j⟩; on another basis pair it leaves the encoded space. Consequently,

    V† O V = X ⊗ I_m,       ⟨O⟩ = η.

At N = 5 this is O = −(X_0 X_2 X_4 + X_1 X_2 X_3).

One all-X measurement setting supplies every summand from the same outcome strings. The global product F uses those same strings and also gives ⟨F⟩ = η. Discarding their joint structure and retaining only low-order summaries removes the very information being sought.

The two estimators have different fluctuations. Two distinct summands of O flip four outer positions; they take a one-excitation or one-hole code vector outside the encoded space. Therefore

    V† O² V = m I_{2m},
    Var(O) = m − η²,       Var(F) = 1 − η².

These are ideal single-shot variances, without readout errors. Equal compressed observables give equal means, not identical measurement processes. In particular O need not preserve the encoded space. Measuring every physical X also destroys the unknown internal state; ignoring the individual outcomes afterward does not undo that measurement.

The global-F exponential alone does not establish the protected subsystem: [H,F] = 0 and Z_c F Z_c = −F give the same decay law for any initial state under this dissipator. Internal motion and its γ-independence are separate parts of the question.

## 6. Bring the phase onto the centre

At a chosen readout time apply the ideal unitary

    U = ∏_{l ≠ c} CNOT(c → l).

The centre controls a flip of every outer qubit. These gates commute and U² = I. Directly on the basis,

    U|d_j⟩  = |0⟩_c ⊗ |a_j⟩,
    UF|d_j⟩ = |1⟩_c ⊗ |a_j⟩,

where |a_j⟩ is |d_j⟩ with its empty centre removed. Hence the product family becomes ρ_q(t) on the physical centre tensor the internal state encoded on the outer sites. Measuring only X_c now reads η without measuring the internal state.

The operator identity U† X_c U = F identifies this as reversible parity extraction. It gives the same readout mean as multiplying all the X outcomes, with a different remaining state.

This is an ideal operation at the readout time, not a gate schedule under the still-running Hamiltonian. Gate duration, available connections and errors are not included. For initially correlated copy/internal inputs, measuring the copy can condition the internal state; the undisturbed-product statement uses the preparation in §2.

## 7. Break protection while the phase stays locally hidden

An asymmetric bond breaks reflection protection but still preserves excitation number. Centre dephasing does too. Therefore a difference supported only between the one-excitation and one-hole sectors stays in those operator sectors, even after this perturbation. F70 still makes it invisible to every marginal with fewer than N−2 sites. This comparison is between coherent and incoherent copy preparations evolved under the same generator. It does not say that the γ > 0 and γ = 0 trajectories still have identical small marginals: the light can now reach the motion within each sector.

The next question is how that internal damage becomes visible while the intercopy coherence remains locally hidden.

Start with N = 5, J = 1, Δ = 0, σ_0 = |0⟩_in⟨0| and the copy in |+⟩. Change only the left end bond from 1 to 1+ε; retain centre-only dephasing. Compare each perturbed trajectory with its γ = 0 counterpart at the same ε and the same time.

Write ρ_{ε,γ}(t) for this trajectory, P = VV† for the original symmetric code projector, and δρ = ρ_{ε,γ}(t) − ρ_{ε,0}(t). Use three readings:

    p_out(ε,γ,t) = 1 − Tr[P ρ_{ε,γ}(t)],
    d_2(ε,γ,t) = max_{|S|=2} ½ ‖Tr_{outside S}(δρ)‖₁,
    d_out(ε,γ,t) = ½ ‖Tr_c(U δρ U†)‖₁.

Report p_out for both γ > 0 and γ = 0; coherent leakage alone can make it nonzero. The other two readings are trace distances to the matched Hamiltonian-only trajectory. The decoded outer comparison includes the full outer density matrices, so it includes output outside the intended internal subspace without postselection. Its reference is the perturbed γ = 0 output, not the symmetric h evolution of §2.

These readings do not share a common onset. At N = 5, for fixed ε and γ as t → 0⁺,

    p_out(ε,γ,t) = ε²t² − ε²(ε²+2)t⁴/3 + (4/15)ε²γt⁵ + O(t⁶),
    d_2(ε,γ,t)   = (2/3)|ε|γt³ + O(t⁴),
    d_out(ε,γ,t) = (4√2/3)|ε|γt³ + O(t⁴).

In particular, the leakage difference between the γ > 0 and γ = 0 arms begins with (4/15)ε²γt⁵, although the two state distances begin at t³. The t² leakage itself is already present without dephasing. This distinguishes a population reading from a change in coherence; it is not a finite waiting time before an effect exists.

Write |s_j⟩ = (|j⟩+|4−j⟩)/√2. The first steps are visible directly:

    H|d_0⟩ = (2+ε)|d_1⟩ + ε|s_1⟩,
    H|s_1⟩ = (2+ε)|s_0⟩ + ε|d_0⟩ + 2√2|2⟩.

The reflection-even amplitude starts at order εt. The centre amplitude starts at −√2 εt². Its coherence with an initially occupied end therefore exists at order εt², before its own population, which begins at order ε²t⁴. Applying dephasing to that coherence supplies the t³ response.

## 8. Distance and the real part of the readout

The method has existing homes: [F96's Dyson argument](../docs/proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md) separates the outcomes of one measurement by their leading time order, [the frozen-divisor proof §8](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) grades cells by their walking distance to the anti-diagonal, and [the chiral trajectory proof](../docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) supplies the bipartite sublattice structure; the (−2i)^m alternation below is its hop-parity face, this note's own reading. Their observables and normalizations differ; the coefficients below follow from this preparation and partial trace, not from transferring their numbers.

Keep the setup of §7, now at N = 2m+1 ≥ 5, centre c = m, and initial internal state |d_0⟩. Only J_0 = 1+ε changes; every other bond remains 1. Here m is the distance from an initially occupied end to the centre, not the distance from the nearer endpoint of the changed bond. The following are short-time expansions at fixed finite ε and γ > 0, with ε ≠ 0. They are not finite-time fits or statements about arbitrary preparations.

    d_out = [√2 m 2^m/(m+1)!] |ε|γ t^(m+1) + O(t^(m+2)),

    d_2 = [2^m/(m+1)!] |ε|γ t^(m+1) + O(t^(m+2))                 (m even),
    d_2 = [2^(m+1)/(m+1)!] |ε|γ max(1,|1+ε|) t^(m+2)
          + O(t^(m+3))                                           (m odd),

    δp_out := p_out(ε,γ,t) − p_out(ε,0,t)
            = (−1)^m [2^(2m+1)/(2m+1)!] ε²γ t^(2m+1)
              + O(t^(2m+2)).

Both leakage arms begin with ε²t². At ε = 0 the original encoding is invariant and all three differences vanish exactly; at γ = 0 the matched differences vanish by definition. To restore a positive common J, replace t by Jt and γ by γ/J in these J = 1 expressions.

| N | m | First order of d_out | First order of d_2 | First order and sign of δp_out |
|---|---:|---:|---:|---|
| 5 | 2 | t³ | t³ | +t⁵ |
| 7 | 3 | t⁴ | t⁵ | −t⁷ |
| 9 | 4 | t⁵ | t⁵ | +t⁹ |
| 11 | 5 | t⁶ | t⁷ | −t¹¹ |
| 13 | 6 | t⁷ | t⁷ | +t¹³ |

Thus at N = 7, 11, … dephasing initially reduces the leakage population relative to the same asymmetric Hamiltonian without dephasing. This does not make the decoded states identical: d_out is positive at its earlier order. A smaller leakage population is not a certificate of an unchanged internal state.

### 8.1 An exact reduced generator, including the coherence between copies

Let W map |0⟩⊗|j⟩ to the physical one-excitation state |j⟩ and |1⟩⊗|j⟩ to F|j⟩, for all j = 0,…,N−1. Unlike V, W includes the centre and the reflection-even states. Excitation conservation makes this 2N-dimensional space invariant even with the defect. Nothing in this subsection uses Δ = 0: a ZZ term adds a diagonal to h, and [H, F] = 0 carries it into both blocks alike, so the block form and the two identities below hold at every Δ; only the shortest-path coefficients of §8.2 to §8.4 use the zero diagonal.

Let h be the N × N real matrix with zero diagonal and off-diagonal entries

    h_{0,1} = h_{1,0} = 2(1+ε),
    h_{j,j+1} = h_{j+1,j} = 2       for j = 1,…,N−2.

Let z = I−2|c⟩⟨c|. The physical density matrix has the form

    ρ = ½ W [[A, B], [B, A]] W†,
    A(0) = B(0) = |d_0⟩⟨d_0|,
    Ȧ = −i[h,A] + γ(zAz−A),
    Ḃ = −i[h,B] − γ(zBz+B).

Both A and B are Hermitian; A is a normalized density matrix and B is an operator block, not a density matrix. On a cell (i,j) the dissipative entry of the B equation is −γ(z_i z_j + 1): −2γ where i and j agree at the centre and 0 where they differ, the Lattice's turned rule read site by site. In particular B_γ(t) = exp(−2γt) A_{−γ}(t). The negative-rate expression is an algebraic solution, not a physical dephasing experiment. This is the one-sided flip relation already used in §2.

This reduction does not discard the intercopy coherence: B carries it. Retaining only the ordinary single-excitation density matrix A would not reproduce the decoded readout.

### 8.2 Why the two-site readout loses one order at odd m

The shortest paths from the two ends to c have m hops. Their amplitudes cancel at ε = 0; with the defect their difference is

    ⟨c|exp(−iht)|d_0⟩ = ε(−2i)^m t^m/(√2 m!) + O(t^(m+2)).

For δA = A_γ−A_0, the first dephasing-sensitive entries are

    δA_{c,0} = −γε(−2i)^m t^(m+1)/(m+1)! + O(t^(m+2)),
    δA_{c,N−1} = −δA_{c,0} + O(t^(m+2)),

and their Hermitian conjugates. A dephasing action must follow the m hops needed to reach the centre row or column. At that order there are no other sensitive entries. One extra dephasing action cannot turn an imaginary coefficient into a real one.

For any two distinct sites a,b, put s = δA_{aa}+δA_{bb} and v = Re δA_{ab}. Directly tracing both copies gives, in the pair's computational basis,

    δρ_{ab} = [[−s/2, 0,   0,   0],
               [0,    s/2, v,   0],
               [0,    v,   s/2, 0],
               [0,    0,   0,  −s/2]],
    ½‖δρ_{ab}‖₁ = |s|/2 + max(|s|/2,|v|).

The B blocks vanish under this two-site trace by F70. The imaginary part of A_{ab} cancels between the two copies; this is a different cancellation from F70. At even m the displayed c/end coefficient is real and yields the first branch of d_2, attained by pairs {0,c} and {c,N−1}.

At odd m it is imaginary. At the next order the largest real coefficients are

    δA_{c,1} = −γε(1+ε)(−2i)^m(2i) t^(m+2)/(m+1)!
               + O(t^(m+3)),

and the c,N−2 counterpart whose magnitude replaces |1+ε| by 1. Moving the centre index one step instead gives coefficients smaller by a factor m+2 than the corresponding end-index move. Diagonal changes need at least 2m hops and one dephasing action, so cannot compete at this order. Taking the largest pair gives the odd-m branch. This argument uses zero diagonal XY hopping; adding a ZZ term or a longitudinal field requires a separate calculation.

### 8.3 The decoder keeps the missing component

After U and tracing the centre, the occupied outer basis consists of the N−1 one-excitation states |j⟩_out (j ≠ c) and the all-ones state |f⟩_out. The exact output is

    ρ_out = Σ_{i,j≠c} A_{ij}|i⟩_out⟨j|
            + A_{cc}|f⟩_out⟨f|
            + Σ_{j≠c}(B_{cj}|f⟩_out⟨j| + B_{jc}|j⟩_out⟨f|).

The identity B_γ = exp(−2γt)A_{−γ}, expanded on the first centre/end entry, gives

    δB_{c,0} = −mγε(−2i)^m t^(m+1)/(m+1)! + O(t^(m+2)),
    δB_{c,N−1} = −δB_{c,0} + O(t^(m+2)).

At this order δρ_out is a Hermitian coupling between |f⟩_out and the difference of the two end states. Its nonzero eigenvalues are ±√2|δB_{c,0}| to leading order, giving d_out for either parity of m. No real-part projection occurs here.

### 8.4 Why the leakage difference changes sign

In the one-excitation space, let R|j⟩ = |N−1−j⟩ and Q = (I+R)/2, the complement of the original reflection-odd space. Then p_out = Tr(QA). Define U_t = exp(−iht), ψ_t = U_t d_0 and C = |c⟩⟨c|. Since z = I−2C, the term with one dephasing insertion is

    δp_out = γ ∫₀ᵗ [−4 Re(ψ_c(s)⟨ψ_t|Q U_{t−s}|c⟩)
                      +4|ψ_c(s)|²⟨c|U_{t−s}† Q U_{t−s}|c⟩] ds
              + terms with at least two dephasing insertions.

For r = t−s, the required shortest-path terms are

    ψ_c(s) = ε(−2i)^m s^m/(√2 m!) + higher orders,
    ⟨ψ_t|Q U_r|c⟩ = ε(2i)^m [s^m−(−r)^m]/(√2 m!) + higher orders,
    ⟨c|U_r† Q U_r|c⟩ = 1 + O(r²).

For the middle identity, before the averaging over the two parts of Q, the identity part contributes ε(2i)^m s^m/(√2 m!). In its reflection part, expand exp(iht)R exp(−ihr). Among terms with a+b=m hops, the end-bond defect lies on the first segment when a ≥ 1 and on the reflected segment when a = 0. The endpoint difference changes sign in that last case. The binomial sum is therefore s^m−2(−r)^m; averaging the two parts of Q gives the displayed expression.

The s^(2m) contributions in the integrand cancel. What remains at the first possible order is

    (−1)^m ε²γ 2^(2m+1)/(m!)² · s^m(t−s)^m.

Using ∫₀ᵗ s^m(t−s)^m ds = (m!)² t^(2m+1)/(2m+1)! gives δp_out in §8. Another dephasing insertion adds a time power without supplying a hop, so cannot contribute at this order. The sign is a property of this fixed projector and preparation, not a general rule that dephasing increases or decreases leakage.

## 9. Reproduce and resume

The exact coefficient calculation needs only the N × N matrices in §8.1. It does not need the full 4^N Liouvillian. For rational ε and γ every Taylor coefficient lies in the Gaussian rationals. Write A(t) = Σ A_k t^k and B(t) = Σ B_k t^k, with A_0 = B_0 = d_0d_0†. Iterate

    A_{k+1} = [−i(hA_k−A_kh) + γ(zA_kz−A_k)]/(k+1),
    B_{k+1} = [−i(hB_k−B_kh) − γ(zB_kz+B_k)]/(k+1).

Run a third sequence with γ = 0 for the matched reference. Apply the linear projector/partial-trace/decoder maps to the coefficient differences BEFORE taking a norm. Locate the first nonzero coefficient matrix exactly, then take its trace norm; do not take a norm of the full matrix series coefficient by coefficient and call that a finite-time distance.

Exact-arithmetic checks were run in the writing session at N = 5, 7, 9, 11, 13 with ε = 1/3, γ = 2/5, and at N = 7 also with ε = −1/3, −1, −3 (a severed and a sign-reversed end bond among them) plus ε = 0 and γ = 0 controls; the scripts were not kept, so the checks are re-derivable from the recurrence and output maps above and not re-locatable in the repository. They were finite checks of the derivation, not its all-N proof.

Sections 10–12 carry the finite-time N = 7 comparison at ε = ±1/10, γ = 3/10. The short-time recurrence remains the entry point for exact coefficients; the finite-time producer evaluates matrix exponentials through their spectral representation and checks that route independently.

C# live-object adoption is not yet implemented. Its starting specification is §8.1 plus the recurrence and output maps in this section. The ordinary single-excitation propagator alone is insufficient: the B block and the decoded output must be retained. No OpenArcs item or hardware confirmation is closed by this note.

## 10. Seven spins over long times

The [experiment producer](../simulations/missing_phase_long_time.py) evolves A_γ, B_γ and the matched A_0 at J = 1, γ = 0.3, ε = ±0.1. The runs cover 0 ≤ t ≤ 6000, first with spacing 0.1, then 0.025. The fine runs have 240001 points each. Their parameter records, selected readings, numerical diagnostics and crossing brackets are in [plus_fine.json](../simulations/results/missing_phase_long_time/plus_fine.json) and [minus_fine.json](../simulations/results/missing_phase_long_time/minus_fine.json); the [plus](../simulations/results/missing_phase_long_time/plus_coarse.json) and [minus](../simulations/results/missing_phase_long_time/minus_coarse.json) coarse records hold the same readings at spacing 0.1.

| t | d_2, ε = +0.1 | d_out, ε = +0.1 | d_2, ε = −0.1 | d_out, ε = −0.1 |
|---:|---:|---:|---:|---:|
| 1 | 0.00279837 | 0.00585356 | 0.00308074 | 0.00781037 |
| 100 | 0.11391515 | 0.17920605 | 0.09746687 | 0.15256187 |
| 1000 | 0.41394906 | 0.66976453 | 0.37503932 | 0.60345248 |
| 6000 | 0.62778098 | 0.76324206 | 0.31520079 | 0.73417439 |

On both grids, every sampled point from t = 0.5 through 6000 has d_out > d_2. This is a finite sampled ordering, not a theorem of pointwise ordering: the two maps retain different centre coherences, so a general data-processing argument does not compare them at finite time.

The early reduction of leakage does not last. Both grids locate one negative-to-positive crossing of δp_out after t = 0.5: near t = 1.250075 for ε = +0.1 and 1.513442 for ε = −0.1. These times are numerical roots refined inside the sign-change brackets, not closed-form constants. No further sign change is resolved through 6000 on either grid.

By t = 100, ‖B‖_F is about 9.0 × 10⁻¹⁷ or 6.0 × 10⁻¹⁷, respectively. Yet ‖A−A_∞‖_F is still about 0.68 or 0.69. The late output distance therefore reads changed motion relative to the γ = 0 reference, not a remaining intercopy phase. A floating zero for B at very late times is underflow, not finite-time disappearance.

### 10.1 The stationary state and the plateau

The single lit centre is an odd-indexed site of this odd XY chain: c = (N−1)/2 is odd exactly when N ≡ 3 mod 4. F157's parity-forced blindness survives the asymmetric bond. Put r = 1+ε. Its remaining zero-energy vector is

    v = (1,0,−r,0,r,0,−r)ᵀ / √(1+3r²),
    hv = 0,       z v = v,
    p = |⟨v|d_0⟩|² = (1+r)²/[2(1+3r²)].

This weight is conserved: hv = 0 and zv = v give d⟨v|A|v⟩/dt = 0, so p is the blind weight at every t.

For the two defects here, the left and right three-site principal characteristic polynomials are x[x²−4(r²+1)] and x(x²−8). Their gcd is exactly x because r² ≠ 1. Thus v spans the whole blind space, not merely one vector inside it. The same argument applies in a punctured neighbourhood of ε = 0, where all bonds remain nonzero and r² ≠ 1. The blind-seat commutant result fixes the kernel of the reduced generator at 1 + blind = 2, the zero-free chain supplying the simple spectrum on the Krylov complement that the identity needs, spanned by vv† and the scalar on its complement; which mixture the preparation relaxes to is argued below, and the stationary density is

    A_∞ = p vv† + (1−p)(I−vv†)/6.

There are no additional undamped oscillations in A. An imaginary-axis eigenoperator must commute with z. At nonzero frequency its centre entry vanishes, and its outer block annihilates the centre coupling vector and every vector in that vector's outer-Hamiltonian Krylov span. What remains is supported on the one-dimensional blind space, whose Hamiltonian is zero, contradicting nonzero frequency. The six-dimensional complement therefore relaxes to its scalar identity part.

B has no undamped mode either. Zero dissipative cost in B requires support only between c and the outer sites. Writing the Hamiltonian in centre/outer blocks as [[0,u†],[u,h_o]], the eigenoperator equations force its off-diagonal vectors proportional to u and u†. The remaining equations require h_o u = 0; here h_o u ≠ 0. Consequently B → 0.

At γ = 0 the decoder maps A_0 isometrically to the occupied outer basis. The same isometry maps A_∞ because its centre off-diagonal entries vanish. Since [h,A_∞] = 0, the limiting distance to the still-moving pure A_0(t) is constant. With q = 1−p,

    lim_{t→∞} d_out(t) = [5q/6 + √((5q/6)²+4pq)]/2.

To obtain this expression, split A_0's pure vector into its blind component and its normalized bright component. On their span A_∞−A_0 has diagonal entries 0 and −5q/6, with off-diagonal magnitude √(pq). The other five bright directions each have eigenvalue q/6. Their trace norm gives the displayed result.

| ε | Conserved blind weight p | Limiting d_out | d_out at t = 6000 |
|---:|---:|---:|---:|
| +0.1 | 441/926 ≈ 0.47624190 | 0.76326556 | 0.76324206 |
| −0.1 | 361/686 ≈ 0.52623907 | 0.73431605 | 0.73417439 |

The two-site readout can retain oscillations against the unitary reference. Its asymptotic curve is bounded above by this plateau: its page is a quantum channel applied to A, so trace-distance contraction applies to A_∞−A_0(t). That asymptotic argument does not establish the finite-time sampled ordering above.

## 11. The clock connection: same rotation, complementary fading

The clock search returned two applicable anchors: [the site-resolved Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), including its live [slow-light reading](../compute/RCPsiSquared.Diagnostics/Ptf/SlowLightDistribution.cs), and [the gamma-fold relation](GAMMA_FOLD_PAIR_OF_MIRRORS.md). The uniform-dephasing [clock-hand ladder](../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) supplies the two-hand vocabulary, a coherence hand and a Takt hand, but its 2γ gap and protected band-edge formula do not transfer to this centre-only profile. Its carrier convention also has hopping J, whereas this note has hopping 2J.

Taking the Hilbert-Schmidt adjoint of the block generator in §8.1 gives an exact identity:

    L_B = −2γ Id − L_A†,
    λ_B = −2γ − conjugate(λ_A).

The Hamiltonian commutator is skew-adjoint, and conjugation by z is self-adjoint. Hence a mode with λ_A = −α+iω has a B partner with

    λ_B = −(2γ−α)+iω.

The angular frequency is the same. The radial decay rates divide the fixed sum 2γ = 0.6. These are spectral partners; an arbitrary prepared trajectory is a sum of modes, not one pair of exponentials.

The slow nonstationary A modes give the following clock readings, computed from the 49 × 49 generators:

| ε | α_A | Frequency magnitude | 1/α_A | α of its B partner | Partner lifetime |
|---:|---:|---:|---:|---:|---:|
| +0.1 | 0.00164538 | 2.90132577 | 607.7622 | 0.59835462 | 1.67125 |
| −0.1 | 0.00134810 | 2.76074395 | 741.7834 | 0.59865190 | 1.67042 |

The oscillation periods are about 2.17 and 2.28. These modes turn quickly and fade slowly; their partners turn at the same rate and fade quickly. The listed B partner is not B's slowest mode. Nor is a spectral lifetime automatically the fitted lifetime of d_2 or d_out: preparation and readout determine which modes contribute.

The Absorption Theorem identifies the physical source of the small α. For an A eigenoperator M,

    α = 2γ · w_c,
    w_c = Σ_{i,j: exactly one of i,j equals c} |M_{ij}|² / ‖M‖_F².

The slow modes have w_c ≈ 0.00274230 and 0.00224684 for the two defects. Thus only about 0.27% or 0.22% of their squared operator weight lies on entries the centre dephasing charges. These are operator weights, not probabilities that an excitation occupies the centre. The producer computes that share from an eigenvector separately from its decay rate and records their residual.

## 12. Executable resumption

Run from the repository root:

```powershell
python -m pytest simulations/tests/test_missing_phase_long_time.py -q
python simulations/missing_phase_long_time.py --epsilon .1 --dt .1 --output simulations/results/missing_phase_long_time/plus_coarse.csv
python simulations/missing_phase_long_time.py --epsilon -.1 --dt .1 --output simulations/results/missing_phase_long_time/minus_coarse.csv
python simulations/missing_phase_long_time.py --epsilon .1 --dt .025 --output simulations/results/missing_phase_long_time/plus_fine.csv
python simulations/missing_phase_long_time.py --epsilon -.1 --dt .025 --output simulations/results/missing_phase_long_time/minus_fine.csv
```

The defaults are γ = 0.3 and t_max = 6000. Each run writes a full CSV and its JSON summary beside it. This is a one-shot Python experiment producer, not the still-open C# live-object adoption. The framework supplies the Hamiltonian commutator; the two reduced dissipators are the diagonal rules of §8.1, applied directly in γ without a square-root round trip: passing √γ·z as a jump operator and squaring would put A and B at rates differing by rounding, and the fifth test asserts the adjoint identity L_B + L_A† = −2γ·Id exactly.

The [tests](../simulations/tests/test_missing_phase_long_time.py) distinguish real from imaginary two-site coherence, distinguish A from B at the decoder's centre row, retain noncentre entries, and compare both reduced generators with a full doubled-space Lindblad construction. Numerical propagation is checked separately: the JSON reports eigenvector conditioning, direct-matrix-exponential comparisons through t = 6000, a direct ODE comparison at three requested accuracies through t = 10, trace/Hermiticity residuals and the minimum eigenvalue of the full physical density matrix. No positivity clipping or state renormalization is applied. These floating diagnostics are not exact-arithmetic certificates.

The next question is how the long relaxation scale depends on the end defect as ε approaches zero. The limits already differ: at ε = 0, d_out(t) = 0 for all t; taking t → ∞ first and then ε → 0 through nonzero defects gives p → 1/2 and d_out → 3/4 by §10.1. Determine the slowing rate from the generator, rather than treating this singular order of limits as a finite-time jump. The current executable's long-time stationary formula is scoped to ε = ±0.1; extend its domain and controls before using it for that sweep.
