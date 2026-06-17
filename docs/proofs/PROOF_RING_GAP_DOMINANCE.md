# Proof of Ring Gap-Dominance: the dihedral lock

**Statement.** For the cyclic XY ring (N sites, bonds (i, i+1 mod N)) under uniform Z-dephasing, the maximum oscillation frequency among the exactly-Re=−2γ Liouvillian modes is

    max|Im| = 2J = J·ρ    (ρ = 2 = the ring adjacency spectral radius = the periodic band top)

for every N **except N = 4**, where the half-filling (2,2) {0,2}-coherence square-root EP mode reaches

    max|Im| = √((2√2·J)² − (2γ)²)  →  2√2·J > 2J.

**Status.** Gate-verified N = 3..6 (full Liouvillian) plus a general-N band-top-reach argument (the single-excitation sector); the N = 4 exception is characterized as the half-filling √-EP. The all-N completeness step (that no exact-(−2γ) mode exceeds 2J for N ≥ 5) is gate-verified through N = 6 and argued by the mechanism, not yet proven by a full ring free-fermion span (the open extension). Verifier: [`simulations/ring_gap_dominance.py`](../../simulations/ring_gap_dominance.py) (gate-first, all stages green).

**Date.** 2026-06-17. Sibling of [`PROOF_CHAIN_GAP_DOMINANCE.md`](PROOF_CHAIN_GAP_DOMINANCE.md).

---

## Abstract

This is the cyclic sibling of chain gap-dominance. On the open chain we proved max|Im| over the exactly-(−2γ) modes is the single-particle band edge E1 = 2J·cos(π/(N+1)) = J·ρ_chain. The ring differs in one structural way that changes the answer: it is translation-invariant (cyclic C_N), so its single-particle band TOP is the k = 0 uniform mode at energy 2J = J·ρ_ring (ρ_ring = 2, the ring adjacency spectral radius — exactly the [TopologyBandEdge](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) "ring band edge"). The cyclic symmetry locks the maximum to that uniform mode — the **dihedral lock**. The lone exception is N = 4: the half-filling (2,2) {0,2}-coherence forms a √-EP that reaches 2√2·J, exceeding the band top, because the Jordan-Wigner wrap bond splits the fermion problem into periodic / anti-periodic sectors and the anti-periodic two-fermion top (√2 + √2 = 2√2) lands on the −2γ floor only at N = 4. This is the same half-filling (2,2) sector that makes K_4 special (the [StructuralCeiling](../../compute/RCPsiSquared.Core/Symmetry/StructuralCeilingClaim.cs) outlier 2 − 2/√3), and it is the ring analogue of the chain's lone N = 3 exception — with the telling difference that the chain's extra modes sit *below* E1 (so the chain max is unchanged) while the ring's sit *above* the band top (so the ring max is genuinely exceeded).

## Setup

The exactly-Re=−2γ subspace is the n_XY = 1 coherence subspace: by the Absorption Theorem Re(λ) = −2γ·⟨n_XY⟩, so a mode sits at exactly −2γ iff ⟨n_XY⟩ = 1. On that subspace the dephasing dissipator is the scalar −2γ, so L = L_H − 2γ and the frequencies are the L_H eigenvalues there.

## The band-top reach (general N): the (0,1) single-excitation sector

The cleanest exact-(−2γ) modes are the (0,1) sector: |vac⟩⟨ψ_k|, where ψ_k is an eigenvector of the single-excitation Hamiltonian H_se = J·A (A = the ring adjacency matrix), energy E_k = J·(adjacency eigenvalue). These are single-particle, so no Jordan-Wigner wrap subtlety arises.

- **n_XY = 1:** |vac⟩ and a single-excitation state differ at exactly one site, so the coherence has Hamming distance 1 ⟹ Re = −2γ.
- **Frequency E_k:** H|vac⟩ = 0 and H_se ψ_k = E_k ψ_k, so [H, |vac⟩⟨ψ_k|] = −E_k|vac⟩⟨ψ_k|, hence L_H[|vac⟩⟨ψ_k|] = +iE_k|vac⟩⟨ψ_k| and L[|vac⟩⟨ψ_k|] = (−2γ + iE_k)|vac⟩⟨ψ_k|, an exact eigenmode.
- **The dihedral lock:** the ring adjacency A is circulant, eigenvalues 2cos(2πm/N); the maximum is 2 at m = 0, the **uniform** mode ψ_0 = (1/√N)Σ_i|e_i⟩, fixed by the cyclic C_N symmetry. So max_k E_k = J·ρ = 2J is reached exactly, at every N.

Verified for N = 3..6 (every (0,1) mode is an exact −2γ + iE_k Liouvillian eigenmode, band top = 2J): `ring_gap_dominance.py` STAGE 1.

## The N=4 exception: the half-filling (2,2) √-EP

The full Liouvillian (N = 3..6) gives max|Im| at the floor = 2J for N = 3, 5, 6, but at N = 4 it is √((2√2·J)² − (2γ)²) → 2√2·J (`ring_gap_dominance.py` STAGE 0, |diff| < 1e-5).

The wrap bond, under Jordan-Wigner, carries the fermion-parity factor (−1)^{N_tot}: odd parity → periodic boundary conditions (single-fermion energies 2J·cos(2πm/N), top 2J), even parity → anti-periodic (2J·cos(π(2m+1)/N)). At N = 4 the anti-periodic single-fermion energies are {√2, √2, −√2, −√2}·J, so the two-fermion top is √2 + √2 = 2√2·J. The half-filling (2,2) {0,2}-coherence built on this top forms a √-EP that lands on the −2γ floor: a γ-sweep confirms its frequency is √((2√2·J)² − (2γ)²) (STAGE 2, g = 0.05/0.10/0.20), exceeding the band top 2J. At N = 6 the (3,3) half-filling does **not** close on −2γ with excess (max|Im| stays 2J), so N = 4 is the lone exception — the unique even half-filling where the {0,2} block lands exactly on −2γ. This is the same sector that makes K_4 and ring-4 special elsewhere in the project.

## Contrast with the chain

| | adjacency radius ρ | max|Im| = J·ρ | lone exception | exception relative to band edge |
|---|---|---|---|---|
| **open chain** | 2cos(π/(N+1)) | 2J·cos(π/(N+1)) = E1 | N = 3 | √(E1² − (2γ)²) **below** E1 (max unchanged) |
| **cyclic ring** | 2 | 2J | N = 4 | 2√2·J **above** 2J (max exceeded) |

Both: max|Im| = J·ρ, the adjacency spectral radius (= the TopologyBandEdge band edge), with a lone small-N half-filling {0,2} √-EP exception. The ring's exception is genuine (it exceeds the band top) because the ring band top 2J sits low enough that the anti-periodic half-filling sum 2√2·J overtakes it; the chain's does not.

## The second-clock floor frequency: chain N=3 and ring N=4 are one law

The chain's N=3 exception and the ring's N=4 exception look unrelated — one sits *below* the band edge (√(E1²−(2γ)²) < E1), the other *above* the band top (2√2·J > 2J). They are the same law (`simulations/second_clock_frequency.py`, gate-first, all stages green). The on-floor {0,2}-coherence √-EP is a 2×2 block [[−2γ, ω],[−ω, −2γ]] (a population coupled to a 2-Hamming coherence), so its floor frequency is

    |Im| = √(B² − (2γ)²),     B = the free-fermion band top of the {0,2} block's sector (the γ→0 coupling),

with the EP at B = 2γ (i.e. Q* = 2J/B). B is γ-independent (verified by a γ-sweep, spread < 1e-15). The three on-floor instances:

| instance | B | = free-fermion band top | floor freq vs band edge |
|---|---|---|---|
| chain N=3 (1,1) | √2·J | E1 = 2J cos(π/4), single-particle | √(E1²−(2γ)²) **below** E1 |
| ring N=4 (1,1) | 2J | periodic single-particle top | **below** 2J (√-shift) |
| ring N=4 (2,2) | 2√2·J | anti-periodic TWO-fermion top | **above** 2J |

So "the second clock overtakes the first (band edge J·ρ)" reduces to **B > J·ρ**: the ring (2,2) is the lone exceed-case (its B = the anti-periodic two-fermion top 2√2·J overtops the periodic band edge 2J); the (1,1) instances have B = J·ρ and the √-shift puts the floor frequency just below. One √-envelope; only B is sector/topology-specific. **Scope:** the {0,2} pins to the −2γ floor (for all Q) only at chain N=3 and ring N=4 — at other N it drifts off-floor — so this is the unification of those two specials, not an all-N floor law. It characterizes the second-clock *floor frequency*; it does NOT close the CoherenceHorizon Tier1Candidate cap (the half-filling V-Effect seam / the dynamic C=0.5 survivor question), which is about the off-floor dynamics, not this frequency.

## Scope and the open extension

Proven for general N: 2J = J·ρ is reached (the (0,1) uniform mode). Gate-verified N = 3..6: it is the maximum, with N = 4 the exception. The remaining step for an all-N theorem is the ring free-fermion completeness (that the parity-resolved free-fermion family spans the exact-(−2γ) subspace for N ≥ 5, so nothing exceeds 2J) — the cyclic analogue of the chain proof's span argument, complicated by the periodic/anti-periodic parity split. The structural picture (N = 4 the unique even half-filling exception) and the N ≤ 6 verification make the all-N statement the strongly-expected extension.

## Anchors

- Verifier: [`simulations/ring_gap_dominance.py`](../../simulations/ring_gap_dominance.py) (STAGE 0 the theorem, STAGE 1 the band-top reach, STAGE 2 the N=4 √-EP).
- Sibling: [`PROOF_CHAIN_GAP_DOMINANCE.md`](PROOF_CHAIN_GAP_DOMINANCE.md), [`simulations/chain_gap_dominance.py`](../../simulations/chain_gap_dominance.py).
- Typed homes: [`TopologyBandEdgeClaim`](../../compute/RCPsiSquared.Core/Symmetry/TopologyBandEdgeClaim.cs) (band edge = J·ρ; the ring N=4 co-occupied-floor mismatch), [`ClockHandLadderClaim`](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs) (the arc `clock_hand_ladder`'s ring sibling). The half-filling sector: [`StructuralCeilingClaim`](../../compute/RCPsiSquared.Core/Symmetry/StructuralCeilingClaim.cs) (K_4 = 2 − 2/√3), the ChiralK mirror [`ChiralKClaim`](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs).
