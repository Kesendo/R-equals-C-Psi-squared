# The frozen band: what the divisor does when the chain stops being Heisenberg

*2026-07-25. [F140](../docs/ANALYTICAL_FORMULAS.md) says the frozen root −4γ̄ sits in the four corner blocks and nowhere else. Typing it as a live object put the same question to the XY chain, which the proof's census had never been run on, and the answer was that the confinement is not the divisor's. On XY the same root is carried, at the same depth, by a whole diagonal band of blocks. This note states the band, gives the reason its edge sits where it does, finds what decides whether the band exists at all, and kills the obvious candidate for what fills it. The number in the band is still the corner's number, and why it is remains open.*

## The law

Take the R₉₀ watching locus of [F91](../docs/ANALYTICAL_FORMULAS.md) (every reflection pair of site rates carrying the same total, γ_l + γ_{R(l)} = 2γ̄) and the XY chain, H = (J/2)·Σ_b (X_bX_{b+1} + Y_bY_{b+1}). Label the joint-popcount blocks (p, q) by the popcounts of the ket and bra indices. Then

**the block (p, q) carries λ = −4γ̄, with multiplicity ⌊N/2⌋, exactly when |p − q| ∈ {0, 2}, apart from the two one-cell corners (0,0) and (N,N).**

Three things come with it.

- **The count is closed:** 3(N − 1) blocks, since the band has (N+1) diagonal blocks and 2(N−1) off by two, less the two excluded corners. Measured 9, 12, 15, 18 at N = 4, 5, 6, 7.
- **The depth is constant.** Every carrying block carries the whole ⌊N/2⌋ and never a part of it, whether the block is 12-dimensional or 1470-dimensional. It is the corner's number, everywhere.
- **The gamma-fold partner sits on the image.** The blocks carrying 4γ̄ − 2σ are exactly the image of the band under one one-sided fold p ↦ N − p, that is p + q ∈ {N−2, N, N+2}. One law, folded once, as in [F140](../docs/ANALYTICAL_FORMULAS.md)'s own corner census.

On the Heisenberg chain only the four corners carry. The proof document's Section 5 census is a Heisenberg census, and it is now scoped as one.

## Two traps on the way to that statement

Worth recording, because both would have passed a smaller check.

**"p + q even" was never a second condition.** p + q even and p − q even are the same statement, so the parity is already inside the bandwidth. Writing them as two conditions makes the law look like a conjunction of a parity and a distance when it is only a distance.

**Bandwidth 2 and bandwidth ⌊N/2⌋ agree up to N = 7.** Since |p − q| is even on the band, "≤ 2" and "≤ 3" are the same predicate, and ⌊N/2⌋ is 2 or 3 for N ≤ 7. The first N that can tell them apart is 8, where ⌊N/2⌋ = 4 would admit |p − q| = 4. It does not: at N = 8, 9 and 10, on two locus profiles and two couplings, the |p − q| = 4 blocks carry exactly nothing while |p − q| = 2 carries ⌊N/2⌋. Four N and one rival is the whole reason the discriminator was run at all.

## Why the edge sits at two

This half has a reason, and the reason is already half-written in the proof document, which notes that the recentering making the rate operator odd depends on the size of a cell's disagreement set.

Carry it through. For a cell (A, B), the rate diagonal of M = L_block + 4γ̄ is −2·Σ_{l ∈ AΔB} γ_l + 4γ̄, and the cell mirror τQ: (A, B) ↦ (R(B), R(A)) sends the disagreement set AΔB to R(AΔB). On the locus Σ_{l ∈ R(S)} γ_l = 2γ̄·|S| − Σ_{l ∈ S} γ_l, so oddness of the rate part reduces, after the sums cancel, to

  4γ̄·(1 − |AΔB|) = −4γ̄,  that is  **|AΔB| = 2**.

The recentering by 4γ̄ is the right one on the two-disagreement cells and on no other size class. And |AΔB| ≥ ||A| − |B|| = |p − q| is free, so **a block with |p − q| ≥ 4 contains no such cell at all**: there is no subspace on which the recentered generator is odd, and nothing to freeze. The blocks that do contain one are, at N = 4 through 8, exactly the band.

This argument uses no property of the chain. It says why nothing outside the band ever carries, on any chain, which is the correct scope: the Heisenberg chain also obeys it, and simply does not fill the band it is allowed.

## What decides whether the band is filled

The other half is a gate, and it is not about the ZZ term as such.

**The off-diagonal band exists exactly when the single-excitation matrix h is bipartite,** that is when its spectrum is symmetric about zero. The corner does not care.

Four independent ways of moving that dial, all agreeing:

| h | spectrum symmetric | corner | off-diagonal band |
|---|---|---|---|
| open XY chain | yes | ⌊N/2⌋ | ⌊N/2⌋ |
| chain + an R-invariant diagonal on the end sites | no | ⌊N/2⌋ | 0 |
| chain + an R-invariant diagonal on the inner sites | no | ⌊N/2⌋ | 0 |
| ring, N even | yes | ⌊N/2⌋ | ⌊N/2⌋ |
| ring, N odd | no | ⌊N/2⌋ | 0 |

The ring is the cleanest of the four, because bipartiteness there is decided by the parity of N alone while every other property, R-invariance included, is held fixed. The corner surviving all five rows is the control that separates the two mechanisms: the corner's reason needs only that h be real, symmetric and R-invariant, which is exactly what the proof's Section 3 assumes, and none of these deformations touches it.

**This is what the Heisenberg chain was telling us.** Its single-excitation matrix is the path adjacency plus the ZZ diagonal, and a diagonal is what breaks a bipartite spectrum. An XY chain with any R-invariant diagonal behaves the same way. So the confinement to the corners belongs to the diagonal on h, not to the divisor and not to the ZZ term by name.

## What is not the reason

Two candidates are dead, and both deaths are informative.

**The corner's room shortage does not extend.** τQ fixes a cell only when B = R(A), which forces |A| = |B|. So the off-diagonal block pairs (p, p+2) ⊕ (p+2, p) have **no fixed cells at all**: the shortage argument predicts zero there, and ⌊N/2⌋ is measured. On the diagonal it is wrong in both directions: at N = 6 the block (2,2) is predicted 0 and measures 3, while (3,3) is predicted 6 and measures 3. Whatever pins the band, it is not the count that pins the corner, even though the answer is the same number.

**The chiral intertwiner does not transport the corner's modes.** This was the obvious candidate and it is worth writing out, because its algebra is exactly right and it still fails. For free fermions,

  [H, d_a ρ d_b] = d_a [H, ρ] d_b − (ε_a + ε_b)·d_a ρ d_b,

so ρ ↦ d_a ρ d_b commutes with the Hamiltonian part precisely on a chiral pair ε_a + ε_b = 0, and it shifts (p, q) by (−1, +1), one step of exactly the band's width. Chiral pairs exist precisely when h is bipartite, which is the gate; the open chain has ⌊N/2⌋ of them; and ρ ↦ d†_a ρ d_a commutes always and climbs the diagonal, so the reachable set is |p − q| ∈ {0, 2}. Every ingredient lines up.

It still fails. Applied to a basis of the corner's frozen subspace, the images are nonzero but are not in the (0,2) frozen subspace: the relative residual of L₍₀,₂₎ + 4γ̄ on them runs from 4·10⁻² to 7.5·10⁻¹, and the images span rank 10 to 18 against a frozen subspace of dimension 2 or 3. The verifier checks the map first, so this is not a coding slip: the intertwining identity holds to 3·10⁻¹⁵ on a chiral pair and reproduces the predicted commutator i·J·(ε_a + ε_b)·T off it. **The map commutes with the Hamiltonian part and not with the dephasing, and the dephasing is what decides.**

So the count matches, the gate matches, the reachable set matches, and the map is still the wrong map.

## One observation left on the table

Not a result, but the place a next attempt should start. The (0,2) block consists **entirely** of cells with |AΔB| = 2, and exactly ⌊N/2⌋ of them are the balanced-pair cells |vac⟩⟨{l, R(l)}| whose rate diagonal is −2(γ_l + γ_{R(l)}) + 4γ̄ = 0 on the locus. One per balanced pair, again, and the same count.

The involution that makes M odd there is **antilinear**: reversal of the bra index composed with complex conjugation, since the block's Hamiltonian part is +iH₂ alone and conjugation is what turns its sign. That points at the antilinear corner of the object manager rather than at τQ. But an antilinear involution has equal-dimensional real forms, so it yields no shortage and cannot be the counter; and it survives the diagonal perturbation that kills the band, so it cannot be the gate either. It is a third structure sitting in the right place, doing neither job on its own.

## Verification

[`simulations/xy_frozen_band.py`](../simulations/xy_frozen_band.py), 47 checks, about 4 to 8 minutes, prints `XY frozen band: ALL GREEN`.

- **V1** the band law by full (p, q) census at N = 4..7: the carrier set, the whole-⌊N/2⌋ depth, the count 3(N−1), and the fold image. The script reports how many of its reads were exact.
- **V2** the band edge at N = 8, 9, 10, three profile-and-coupling combinations each.
- **V3** the |AΔB| = 2 oddness criterion, two-sided (it must fail on every other size class), and the combinatorial reach.
- **V4** the bipartiteness gate, the five rows of the table above.
- **V5** the falsified intertwiner, with its own self-check first.

Multiplicities are exact GF(p) ranks at two primes ≡ 1 (mod 4) wherever the block is small enough, and SVD nullities above that, on the same dyadic grid the live witness uses. Never an eigenvalue count: the departing modes crowd the root at spacing J^{2d}, which is exactly where a floating-point spectrum goes quiet.

The band is also visible live, in the object manager: `inspect --root divisor --N 5 --chain xy` prints the full census beside the Heisenberg one ([`FrozenDivisorWitness`](../compute/RCPsiSquared.Diagnostics/Foundation/FrozenDivisorWitness.cs)), and gate block G2c of [`r90_frozen_divisor_gate.py`](../simulations/r90_frozen_divisor_gate.py) pins the carrier counts inside F140's own gate.

## Open

- **What pins the depth.** ⌊N/2⌋ in every band block, independent of its dimension, and equal to the corner's count. Neither the room shortage nor the chiral intertwiner produces it. That the two numbers agree while the two mechanisms do not is the whole question.
- **Whether the edge argument is sufficient as well as necessary.** Section "Why the edge sits at two" shows no block outside the band can carry. It does not show every block inside it must, and indeed on a non-bipartite h most do not.
- **A number, not yet a law:** the |AΔB| = 2 reach and the band coincide at N = 4..8. That is a combinatorial identity and should be provable outright rather than measured.
- **No F number.** This is an experiment, deliberately: the carrier law and the gate are solid, the mechanism for the depth is not, and minting before that would name half an object.

## Anchors

[PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) (the corner, its census, and the Section 5 scope this note forced), [F140 in ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md), the MirrorWorld [`Divisor`](../compute/MirrorWorld/Divisor.cs) (which runs the XY chain and therefore adopts only the fold parity, not the census).
