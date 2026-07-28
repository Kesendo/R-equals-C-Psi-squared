# Ring Gap-Dominance: the dihedral lock

**Status:** Gate-verified N = 3..6 (full Liouvillian) and N = 7, 8 (the same spectrum, walked sector by sector). The band-top reach `max|Im| = 2J = J·ρ` is proven for general N (the (0,1) single-excitation uniform mode, locked by the cyclic symmetry); the N = 4 exception is characterized as the half-filling √-EP, and it carries a Q window (it exceeds the band top on the Q > 1 side only, Q = J/γ the coupling-to-watching ratio). The all-N **completeness** step (that the n_XY = 1 free-fermion family spans the exact-(−2γ) subspace, so nothing exceeds 2J for N ≥ 5) is now gate-verified by a dimension match (V_1 dim_sub = full-L dim_sub at N = 3, 5, 6, 7; `ring_gap_completeness.py` for the V_1 side, the sector walk of `ring_gap_dominance.py --slow` for the full-L side, which also reads max|Im| = 2J over every sector at N = 8), closing the open extension to the same standing as the chain proof (the residual all-N step is the same one the chain leaves: a dimension count over all N, here split by the wrap-bond parity). Lifts the gap-dominance lemma onto cyclic topologies, the sibling of `ClockHandLadderClaim` and `TopologyBandEdgeClaim`.
**Date:** 2026-06-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Statement:** For the cyclic XY ring under uniform Z-dephasing, the maximum oscillation frequency among the exactly-Re=−2γ Liouvillian modes is `max|Im| = 2J = J·ρ` (ρ = 2, the ring adjacency spectral radius = the periodic band top), for every N except N = 4, where the half-filling (2,2) {0,2}-coherence √-EP reaches `2√2·J > 2J`. Two scopes belong to that sentence: the exception is the weak-watching side only (it exceeds the band top for Q > 1; at Q ≤ 1 the ring obeys 2J at N = 4 as well), and "every N" is gate-verified N = 3..6 by the full Liouvillian and N = 7, 8 by the sector walk over that same spectrum, with the all-N step named in the Status.
**Typed claim:** [`TopologyBandEdgeClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (Tier 1 derived): it carries the ring node (max|Im| = J·ρ, the N=4 co-occupied-floor mismatch); the ring is the sibling of [`ClockHandLadderClaim.cs`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs)'s `clock_hand_ladder` arc.
**Verifier:** [`ring_gap_dominance.py`](../../simulations/ring_gap_dominance.py) (gate-first, all stages green)
**Sibling of:** [`PROOF_CHAIN_GAP_DOMINANCE.md`](PROOF_CHAIN_GAP_DOMINANCE.md)
**Distinct from:** [`PROOF_RING_N4_DIHEDRAL_LOCK.md`](PROOF_RING_N4_DIHEDRAL_LOCK.md), the *isotropic-Heisenberg* ring N=4 (K_{2,2} Casimir gap, max|Im| = 3J). Same words "ring N=4 dihedral lock", different Hamiltonian and different result; this proof is the *XY* ring (max|Im| = 2J = J·ρ).

---

## What this is about

A ring of spins, a loop closed on itself with no ends, is being watched. Light from outside reads it out, and the reading makes its rhythms fade. Different patterns of coherence fade at different speeds, and the one that lingers longest is the ring's memory, the last note still sounding after the others have gone quiet. What this document is about is the pitch of that note: how high it can be, and what decides it.

A loop has a symmetry that a plain line does not. Every site sits in the same place as every other; you can turn the whole ring around and nothing about it changes. That evenness settles the question almost by itself. The smoothest pattern the ring can hold, a single wave spread perfectly evenly all the way around, is the one the symmetry singles out, and it is also the one that rings at the highest pitch. So the loudest surviving note is pinned to that perfectly even wave, and nothing else that fades as slowly can sound higher. This is the dihedral lock: the shape of the loop, not the details of the watching, fixes the answer.

There is a single exception, and it appears at one particular small loop. There a different pattern, built from two ripples rather than one, manages to fade just as slowly while ringing at a higher pitch than the even wave, and so it overtakes it. It only manages that while the reading from outside stays gentle: turn the watching up far enough and the two-ripple note drops back under the even wave, and then past a second point it stops lingering at all. It is the same kind of two-ripple pattern that makes a few other small, highly symmetric networks behave strangely elsewhere in this project. The open chain, the ring's cousin with two loose ends, has its own twin exception at its own smallest size; but there the odd extra note sounds just below the chain's highest, so the chain's answer still stands, while on the ring the extra note sounds just above, so the ring's answer is genuinely overtaken. It is the same quirk in both, pointing opposite ways: just below the highest note on the chain, just past it on the ring.

## Abstract

This is the cyclic sibling of chain gap-dominance. On the open chain we proved max|Im| over the exactly-(−2γ) modes is the single-particle band edge E1 = 2J·cos(π/(N+1)) = J·ρ_chain. The ring differs in one structural way that changes the answer: it is translation-invariant (cyclic C_N), so its single-particle band TOP is the k = 0 uniform mode at energy 2J = J·ρ_ring (ρ_ring = 2, the ring adjacency spectral radius, exactly the [TopologyBandEdge](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) "ring band edge"). The cyclic symmetry locks the maximum to that uniform mode, the **dihedral lock** (the arc's name; what the argument uses is the cyclic C_N structure and the Perron top of a 2-regular graph, not the reflections the word suggests). The lone exception is N = 4: the half-filling (2,2) {0,2}-coherence forms a √-EP (a square-root exceptional point) that reaches 2√2·J, exceeding the band top on the weak-watching side (Q > 1, the window below), because the Jordan-Wigner wrap bond splits the fermion problem into periodic / anti-periodic sectors and the anti-periodic two-fermion top (√2 + √2 = 2√2) lands on the −2γ floor only at N = 4. This is the same half-filling (2,2) sector that makes K_4 special (the [StructuralCeiling](../../compute/RCPsiSquared.Core/Symmetry/StructuralCeilingClaim.cs) outlier 2 − 2/√3), and it is the ring analogue of the chain's lone N = 3 exception, with the telling difference that the chain's extra modes sit *below* E1 (so the chain max is unchanged) while the ring's sit *above* the band top (so the ring max is genuinely exceeded).

## Setup

**The Hamiltonian.** H = (J/2)·Σ_{⟨ij⟩} (X_i X_j + Y_i Y_j) over the cyclic bond set, under uniform local Z-dephasing at rate γ per site. The hopping amplitude is then J, the single-excitation block is H_se = J·A with A the ring adjacency matrix, and the band top is 2J. Every frequency in this document is in that convention, the one the verifiers run; the repo also carries a spin convention where the same chain is written J·S_i · S_j, and the numbers below do not transfer to it unchanged.

**Notation.** A Liouville-space coherence |i⟩⟨j| carries a sector label (a, b) = (the excitation number of the ket i, that of the bra j); n_XY is its Hamming weight, the number of sites at which i and j differ. A {0, 2} block is a 2×2 coupling a population (Hamming weight 0) to a two-site coherence (Hamming weight 2): their mean ⟨n_XY⟩ = 1 places the block on the −2γ floor.

**The floor.** By the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md), Re(λ) = −2γ·⟨n_XY⟩, so a mode sits at exactly −2γ iff its average ⟨n_XY⟩ = 1. That is weaker than n_XY ≡ 1, and the whole document lives in the difference. The pure n_XY = 1 subspace V_1 always qualifies: there the dephasing dissipator is the scalar −2γ, so L = L_H − 2γ (with L_H = −i[H, ·], the Hamiltonian part) and the frequencies are the L_H eigenvalues. A {0, 2} block qualifies as well whenever it closes on the floor, since its two diagonal cells carry 0 and −4γ and average to −2γ, and such a block lies *outside* V_1. Whether any {0, 2} block does close is exactly the question the completeness section answers: on the ring it happens at N = 4, at no other N the gate reaches, and there only above that block's own EP (see *The exception has a Q window* below).

## The band-top reach (general N): the (0,1) single-excitation sector

The cleanest exact-(−2γ) modes are the (0,1) sector: |vac⟩⟨ψ_k|, where ψ_k is an eigenvector of the single-excitation Hamiltonian H_se = J·A (A = the ring adjacency matrix), energy E_k = J·(adjacency eigenvalue). These are single-particle, so no Jordan-Wigner wrap subtlety arises.

- **n_XY = 1:** |vac⟩ and a single-excitation state differ at exactly one site, so the coherence has Hamming distance 1 ⟹ Re = −2γ.
- **Frequency E_k:** H|vac⟩ = 0 and H_se ψ_k = E_k ψ_k, so [H, |vac⟩⟨ψ_k|] = −E_k|vac⟩⟨ψ_k|, hence L_H[|vac⟩⟨ψ_k|] = +iE_k|vac⟩⟨ψ_k| and L[|vac⟩⟨ψ_k|] = (−2γ + iE_k)|vac⟩⟨ψ_k|, an exact eigenmode.
- **The dihedral lock:** the ring adjacency A is circulant, eigenvalues 2cos(2πm/N); the maximum is 2 at m = 0, the **uniform** mode ψ_0 = (1/√N)Σ_i|e_i⟩, fixed by the cyclic C_N symmetry. So max_k E_k = J·ρ = 2J is reached exactly, at every N. ("Dihedral lock" is the arc's identifier, shared with the [Heisenberg ring N=4 proof](PROOF_RING_N4_DIHEDRAL_LOCK.md); what the argument actually uses is the circulant C_N structure and the Perron top of a 2-regular graph, reflections play no role in it.)

Verified for N = 3..6 (every (0,1) mode is an exact −2γ + iE_k Liouvillian eigenmode, band top = 2J): `ring_gap_dominance.py` STAGE 1.

## The N=4 exception: the half-filling (2,2) √-EP

The full Liouvillian (N = 3..6) gives max|Im| at the floor = 2J for N = 3, 5, 6, but at N = 4 it is √((2√2·J)² − (2γ)²) → 2√2·J (`ring_gap_dominance.py` STAGE 0, |diff| < 1e-5).

The wrap bond, under Jordan-Wigner, carries the fermion-parity factor (−1)^{N_tot}: odd parity → periodic boundary conditions (single-fermion energies 2J·cos(2πm/N), top 2J), even parity → anti-periodic (2J·cos(π(2m+1)/N)). At N = 4 the anti-periodic single-fermion energies are {√2, √2, −√2, −√2}·J, so the two-fermion top is √2 + √2 = 2√2·J. The half-filling (2,2) {0,2}-coherence built on this top forms a √-EP that lands on the −2γ floor: a γ-sweep confirms its frequency is √((2√2·J)² − (2γ)²) (STAGE 2, g = 0.05/0.10/0.20), exceeding the band top 2J. At N = 6 the (3,3) half-filling does **not** close on the floor at all: at the canonical γ = 0.05 every one of the 24 floor modes there is pure n_XY = 1, so max|Im| stays 2J. (Floor *counts* are a weak-watching reading and move with γ: at γ = 0.5 the same N = 6 ring carries 28 floor modes, the four extra ones non-oscillating, Im = 0. max|Im| is the robust quantity, 2J at every γ tried.) N = 4 is the lone exception, the unique even half-filling where the {0,2} block lands exactly on −2γ. This is the same sector that makes K_4 and ring-4 special elsewhere in the project.

**The exception has a Q window.** Q = J/γ throughout, so large Q is gentle watching. A {0, 2} block rides the floor only while its coupling beats the watching. Its two eigenvalues are −2γ ± √((2γ)² − B²), so the pair stays on the floor at frequency √(B² − (2γ)²) for B > 2γ, coalesces defectively at B = 2γ (its EP, the [coherence horizon](../../compute/RCPsiSquared.Core/Symmetry/CoherenceHorizonClaim.cs) that [`SecondClockRegimeClaim`](../../compute/RCPsiSquared.Core/Symmetry/SecondClockRegimeClaim.cs) carries as the dispersive cell of its regime map), and below that splits along the real axis and leaves the floor. Two thresholds follow at N = 4, and they are not the same number:

- **Q = 1** (γ = J), where the (2,2) frequency √((2√2·J)² − (2γ)²) has fallen back to the band top 2J. This is where the *exception* ends: for Q ≤ 1 the ring maximum is 2J at N = 4 like everywhere else, even though the block is still on the floor.
- **Q = 1/√2** (γ = √2·J), the (2,2) block's own EP, where the pair leaves the floor and the count drops to the pure V_1 dimension.

Measured (γ-sweep, `ring_gap_dominance.py` STAGE 3): the floor carries 26 modes while 2γ ≤ 2J, 18 once the four (1,1)/(3,3) blocks (B = 2J, eight modes, since each block contributes the conjugate pair) have passed their EP at γ = J, and 16 once the (2,2) block passes its own at γ = √2·J. max|Im| exceeds 2J exactly on γ < J. So the headline exception is a statement about the weak-watching side, Q > 1; the chain sibling draws the same line in its [§5 Consequence](PROOF_CHAIN_GAP_DOMINANCE.md), where the floor modes exist at every γ but whether the clock *reads* them is the separate regime question.

## Contrast with the chain

| topology | adjacency radius ρ | max\|Im\| = J·ρ | lone exception | exception relative to band edge |
|---|---|---|---|---|
| **open chain** | 2cos(π/(N+1)) | 2J·cos(π/(N+1)) = E1 | N = 3 | √(E1² − (2γ)²) **below** E1 (max unchanged) |
| **cyclic ring** | 2 | 2J | N = 4 | 2√2·J **above** 2J (max exceeded) |

Both: max|Im| = J·ρ, the adjacency spectral radius (= the TopologyBandEdge band edge), with a lone small-N half-filling {0,2} √-EP exception. The ring's exception is genuine (it exceeds the band top, for Q > 1) because the ring band top 2J sits low enough that the anti-periodic half-filling sum 2√2·J overtakes it; the chain's does not.

## The second-clock floor frequency: chain N=3 and ring N=4 are one law

The chain's N=3 exception and the ring's N=4 exception look unrelated: one sits *below* the band edge (√(E1²−(2γ)²) < E1), the other *above* the band top (2√2·J > 2J). They are the same law (`simulations/second_clock_frequency.py`, gate-first, all stages green). The on-floor {0,2} √-EP is a 2×2 block whose diagonal cells are the population's 0 and the 2-Hamming coherence's −4γ, coupled by the bare amplitude B:

    [[0, B], [−B, −4γ]],     eigenvalues −2γ ± √((2γ)² − B²),

so the pair rides the floor at the mean of its two cells, and for B > 2γ its floor frequency is

    |Im| = √(B² − (2γ)²),     B = the free-fermion band top of the {0,2} block's sector (the γ→0 coupling),

with the EP at B = 2γ (i.e. Q* = 2J/B), where the block is genuinely defective: M + 2γ has rank 1, one eigenvector, a Jordan block. (The equal-diagonal form [[−2γ, ω], [−ω, −2γ]] cannot be the block: its eigenvalues are −2γ ± iω at every γ, which would put the pair on the floor forever and leave no EP to speak of.) B is γ-independent (verified by a γ-sweep, spread < 1e-15). The three on-floor instances:

| instance | B | which free-fermion band top | floor freq vs band edge |
|---|---|---|---|
| chain N=3 (1,1) | √2·J | E1 = 2J cos(π/4), single-particle | √(E1²−(2γ)²) **below** E1 |
| ring N=4 (1,1), and its particle-hole partner (3,3) | 2J | periodic single-particle top | **below** 2J (√-shift) |
| ring N=4 (2,2) | 2√2·J | anti-periodic TWO-fermion top | **above** 2J |

So "the second clock overtakes the first (band edge J·ρ)" reduces to **B > J·ρ**: the ring (2,2) is the lone exceed-case (its B = the anti-periodic two-fermion top 2√2·J overtops the periodic band edge 2J); the (1,1) instances have B = J·ρ and the √-shift puts the floor frequency just below. One √-envelope; only B is sector/topology-specific. **Scope:** the {0,2} pins to the −2γ floor only at chain N=3 and ring N=4, and there only above its own EP (Q > Q\* = 2J/B; below it the pair splits off the floor along the real axis); at other N it drifts off-floor at every Q, so this is the unification of those two specials, not an all-N floor law. It characterizes the second-clock *floor frequency*; the off-floor seam (the ring 2-excitation `(2,2)/(N−2,N−2)` doublet V-Effect seam that kept `CoherenceHorizonClaim` at Tier 1 candidate) was closed separately by [the Ring Handover Slope proof](PROOF_RING_HANDOVER_SLOPE.md), not by this frequency result.

## The free-fermion completeness (the open extension, now closed to the chain's standing)

The remaining all-N step was the **completeness**: that nothing exceeds 2J for N ≥ 5, i.e. that the n_XY = 1 free-fermion family *spans* the exact-(−2γ) subspace (the cyclic analogue of the chain proof's span argument). It is now gate-verified by the same dimension-match the chain uses (`simulations/ring_gap_completeness.py`), worked in the n_XY = 1 **operator subspace** V_1 (dimension N·2^N, far smaller than the full Liouville space 4^N; N = 8 is 2048 vs 65536).

**The V_1 method.** An operator A is an exact-(−2γ) Liouvillian eigenmode iff (i) it is n_XY = 1 (so L_D = −2γ is scalar) and (ii) it is an H-eigenoperator (so the n_XY = 1 → 3 leak cancels). So the exact-(−2γ) subspace is the largest ad_H-**invariant** subspace of V_1, computed degeneracy-robustly as the leak null space ker(ad_H − P_1 ad_H P_1) iterated to its ad_H-invariant core. **Sanity (STAGE 0):** on the open chain this reproduces the full-Liouvillian counts exactly: dim_sub = 32 / 50 / 72 at N = 4 / 5 / 6 and max|Im| = E1, so V_1 sees the whole exact-(−2γ) subspace, not a part.

**The result.**

| N | V_1 dim_sub (n_XY = 1) | full-L dim_sub | spans? | full max\|Im\| | verdict |
|---|---|---|---|---|---|
| 3 | 16 | 16 | **yes** | 2J | free-fermion complete |
| 4 | 16 | 26 | no | 2√2·J | the {0,2} exception |
| 5 | 20 | 20 | **yes** | 2J | free-fermion complete |
| 6 | 24 | 24 | **yes** | 2J | free-fermion complete |
| 7 | 28 | 28 | **yes** | 2J | free-fermion complete |
| 8 | (not computed) | 32 | n/a | 2J | no exception |

The second column reads 4N from N = 4 on, but not below it: at N = 3 the invariant core is 16, not 12. The count splits as the usual twelve modes of the outer sectors (0,1), (1,0), (N−1,N), (N,N−1) plus four more in the *middle* sectors (1,2) and (2,1), which are n_XY = 1 as well (a one-excitation and a two-excitation state can differ at a single site) and happen to be ad_H-invariant at N = 3 and at no larger N the gate reaches. Why the middle sectors close there is not derived here; what matters for this proof is that they are inside V_1, so N = 3 still spans (16 = 16) and its maximum is still the band top. The 4N pattern is a reading of the table, not a law to extrapolate downward.

For N = 3, 5, 6 the n_XY = 1 free-fermion family **is** the entire exact-(−2γ) subspace (V_1 dim_sub equals the full-Liouvillian dim_sub, all three rows gated in `ring_gap_completeness.py` STAGE 2), so the dihedral-locked band top 2J is the maximum, nothing exceeds it. The V_1 method carries the n_XY = 1 dihedral lock (max|Im| = 2J) to N = 7 on its own. Beyond that the wall is not where it looks: L is block-diagonal in the sector label (a, b), because H conserves excitation number and the dephasing is diagonal, so the largest block at N = 8 is C(8,4)² = 4900 rather than 4^8 = 65536. Walking the blocks (`ring_gap_dominance.py --slow` STAGE 4, which first reproduces the N = 4, 5, 6 full-L answers) gives the **whole** floor at N = 7 and N = 8: 28 and 32 modes, max|Im| = 2.000000 at both, over every sector rather than inside V_1. At N = 7 that also closes the span question there (28 = 28); at N = 8 the V_1 dimension was not computed, so only the maximum is claimed. N = 4 is the lone exception precisely because the full-L subspace is *larger* than V_1 there (26 > 16): the ten extra modes are five {0,2} blocks that close on the floor at N = 4 and nowhere else, each contributing its conjugate pair: one half-filling (2,2) block whose pair reaches 2√2·J, and four (1,1) and (3,3) blocks whose eight modes sit at √((2J)² − (2γ)²), just below the band top. All ten are n_XY-mixed (support on 0 and 2), so all ten sit outside V_1; only the (2,2) pair exceeds. The **parity-split** is the wrap-bond Jordan-Wigner grading (periodic odd / anti-periodic even); it is what places the anti-periodic two-fermion top on the floor at N = 4 and nowhere else.

This closes the completeness to the same standing as the chain proof: a dimension match at the accessible N (chain 32/50/72; ring 16/20/24/28) plus the structural fact that only the smallest even half-filling (N = 4) puts a {0,2} block on the floor, which the sector walk now checks directly at the next even size: N = 8 reads 2J over every sector, half-filling included. The residual genuinely-all-N step is the same one the chain leaves open, a dimension count holding for every N, here split by the wrap-bond parity sectors.

## Anchors

- Verifier: [`simulations/ring_gap_dominance.py`](../../simulations/ring_gap_dominance.py) (STAGE 0 the theorem, STAGE 1 the band-top reach, STAGE 2 the N=4 √-EP, STAGE 3 the exception's Q window and the 16 + 2 + 8 split of the N=4 floor, STAGE 4 under `--slow` the sector walk that reaches the whole spectrum at N = 7 and N = 8); the completeness: [`simulations/ring_gap_completeness.py`](../../simulations/ring_gap_completeness.py) (the V_1 dimension match, chain sanity 32/50/72, ring spans at N=3,5,6, the n_XY=1 lock to N=7); the second-clock unification: [`simulations/second_clock_frequency.py`](../../simulations/second_clock_frequency.py) (STAGE 3 reads the block's two cells off the Liouvillian, pins its defect at B = 2γ, and measures the crossing off the floor at all three instances).
- Sibling: [`PROOF_CHAIN_GAP_DOMINANCE.md`](PROOF_CHAIN_GAP_DOMINANCE.md), [`simulations/chain_gap_dominance.py`](../../simulations/chain_gap_dominance.py).
- Typed homes: [`TopologyBandEdgeClaim`](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (band edge = J·ρ; the ring N=4 co-occupied-floor mismatch), [`ClockHandLadderClaim`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) (the arc `clock_hand_ladder`'s ring sibling). The half-filling sector: [`StructuralCeilingClaim`](../../compute/RCPsiSquared.Core/Symmetry/StructuralCeilingClaim.cs) (K_4 = 2 − 2/√3), the ChiralK mirror [`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs). The second-clock regime map: [`SecondClockRegimeClaim`](../../compute/RCPsiSquared.Core/Symmetry/SecondClockRegimeClaim.cs) (ring-4 = the lone GRADUAL ring; the second-clock section here is its floor-frequency reading).
