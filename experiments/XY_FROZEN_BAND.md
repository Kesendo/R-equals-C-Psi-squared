# The frozen band: what the divisor does when the chain stops being Heisenberg

*2026-07-25. [F140](../docs/ANALYTICAL_FORMULAS.md) said the frozen root −4γ̄ sits in the four corner blocks and nowhere else. Typing it as a live object put the same question to the XY chain, which the proof's census had never been run on, and the answer was that the confinement is not the divisor's. On XY the same root is carried, at the same depth, by a whole diagonal band of blocks. This note states the band, gives the reason its edge sits where it does, finds what decides whether the band exists at all, and buries the two candidates for what fills it. The number is the count of balanced pairs, which is visible directly once the coupling is switched off; what is still open is why every block collapses onto it from wildly different starting counts.*

## The law

Take the R₉₀ watching locus of [F91](../docs/ANALYTICAL_FORMULAS.md) (every reflection pair of site rates carrying the same total, γ_l + γ_{R(l)} = 2γ̄) and the XY chain, H = (J/2)·Σ_b (X_bX_{b+1} + Y_bY_{b+1}). Label the joint-popcount blocks (p, q) by the popcounts of the ket and bra indices. Then

**the blocks carrying λ = −4γ̄ are the band |p − q| ∈ {0, 2}, less the two one-cell blocks (0,0) and (N,N), each of them carrying multiplicity ⌊N/2⌋ at the couplings tested.**

That last qualifier is not decoration, and this note would be making the very mistake its sibling commit was written to correct if it dropped it: ⌊N/2⌋ is a generic depth, not an inventory. The corner block (1,1) is itself in the band, and at J = 0 it carries twice that, semisimply, which is F140's own corollary. Everything below is read at J = 3/4 and J = 3/2.

That is a measured law, and the two directions do not stand on the same footing. Nothing OUTSIDE the band carries, and the section below removes the only mechanism that could have let it; everything INSIDE the band does carry, and for that there is measurement and a gate, no argument. Scope of the measurement: a full census over every (p, q) at N = 4, 5, 6, 7, and band-edge probes at N = 8, 9, 10.

**Notation, because this note borrows from two directions.** Cells here are written (A, B) with A and B index SETS, so that the disagreement set AΔB has a size worth talking about; the proof document writes (a, b) for the single-site case and calls that set S. The **one-cell blocks** (0,0) and (N,N) excluded above are the vacuum and full popcount blocks, not the divisor's four corner blocks, which are always named as such. The fermion ladder operators below are d_a, d†_a, and the site distance the proof writes d_c appears here only inside the exponent J^{2d}. The transport map is written Ξ, not T, because the proof already spends T on the bordered cofactor matrix of its Section 8.

Three things come with it.

- **The count is closed:** 3(N − 1) blocks, since the band has (N+1) diagonal blocks and 2(N−1) off by two, less the two excluded one-cell blocks. Measured 9, 12, 15, 18 at N = 4, 5, 6, 7.
- **The depth is constant.** Every carrying block carries the whole ⌊N/2⌋ and never a part of it, across the range the census reaches: from the 6-dimensional (0,2) block at N = 4 to the 1225-dimensional (3,3) at N = 7. It is the corner's number, in every block the census reached.
- **The gamma-fold partner sits on the image.** The blocks carrying 4γ̄ − 2σ are exactly the image of the band under one one-sided fold p ↦ N − p. Written out that is p + q ∈ {N−2, N, N+2} **less the images of the two excluded one-cell blocks**, namely (N,0) and (0,N); without that exclusion the set would hold 3N − 1 blocks instead of the 3(N − 1) counted above, and the two extra ones do not carry. One law, folded once, as in [F140](../docs/ANALYTICAL_FORMULAS.md)'s own corner census.

On the Heisenberg chain only the four corners carry. The proof document's Section 5 census is a Heisenberg census, and it is now scoped as one.

## Two traps on the way to that statement

Worth recording, because both would have passed a smaller check.

**"p + q even" was never a second condition.** p + q even and p − q even are the same statement, so the parity is already inside the bandwidth. Writing them as two conditions makes the law look like a conjunction of a parity and a distance when it is only a distance.

**Bandwidth 2 and bandwidth ⌊N/2⌋ agree up to N = 7.** Since |p − q| is even on the band, "≤ 2" and "≤ 3" are the same predicate, and ⌊N/2⌋ is 2 or 3 for N ≤ 7. The first N that can tell them apart is 8, where ⌊N/2⌋ = 4 would admit |p − q| = 4. It does not: at N = 8, 9 and 10, in three profile-and-coupling combinations each (profile A at J = 3/4 and at J = 3/2, profile B at J = 3/4), the |p − q| = 4 blocks carry exactly nothing while |p − q| = 2 carries ⌊N/2⌋. Four N and one rival is the whole reason the discriminator was run at all.

## Why the edge sits where it does

This half has a reason. Read what it delivers before reading it: it shows that the ONE known freezing mechanism is unavailable outside the band, which is why the edge sits at two rather than anywhere else. It does not prove that −4γ̄ is absent out there for every coupling; that absence is measured. Mechanism-unavailability and a null are different statements, and only the first is argued here.

The ingredient is already half-written in the proof document, which notes that the recentering making the rate operator odd depends on the size of a cell's disagreement set.

Carry it through. For a cell (A, B), the rate diagonal of M = L_block + 4γ̄ is −2·Σ_{l ∈ AΔB} γ_l + 4γ̄, and the cell mirror τQ: (A, B) ↦ (R(B), R(A)) sends the disagreement set AΔB to R(AΔB). On the locus Σ_{l ∈ R(S)} γ_l = 2γ̄·|S| − Σ_{l ∈ S} γ_l, so oddness of the rate part reduces, after the sums cancel, to

  4γ̄·(1 − |AΔB|) = −4γ̄,  that is  **|AΔB| = 2**.

The recentering by 4γ̄ is the right one on the two-disagreement cells and on no other size class. And |AΔB| ≥ ||A| − |B|| = |p − q| is free, so **a block with |p − q| ≥ 4 contains no such cell at all**: there is no subspace on which the recentered generator is odd, and so no room shortage to freeze anything. The blocks that do contain one are, at N = 4 through 8, exactly the band.

This argument uses no property of the chain, and that is its scope: on ANY chain, a block outside the band is beyond the reach of this mechanism. The Heisenberg chain obeys it too, and simply does not fill the band it is allowed. What the argument leaves open is whether some other mechanism could put −4γ̄ outside the band; nothing found so far does, and the census says nothing does at N ≤ 10.

## What decides whether the band is filled

The other half is a gate, and it is not about the ZZ term as such.

**The off-diagonal band exists exactly when the spectrum of the single-excitation matrix h is symmetric about zero.** For a plain adjacency matrix that is bipartiteness of the hopping graph, and the word is used in that sense below; note that the deformed rows of the table carry a diagonal, so their h is not the adjacency matrix of any graph, and only the spectral reading applies to them. The corner does not care either way.

Four independent deformations, against the undeformed chain as the baseline, all agreeing:

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

**The corner's room shortage does not extend.** τQ fixes a cell only when B = R(A), which forces |A| = |B|. So the off-diagonal block pairs (p, p+2) ⊕ (p+2, p) have **no fixed cells at all**: the shortage argument predicts zero there, and ⌊N/2⌋ is measured. On the diagonal it is wrong in both directions at N = 6, where (2,2) is predicted 0 and measures 3 while (3,3) is predicted 6 and measures 3; at N = 5 it happens to match everywhere, which is worth saying, since a rule that is sometimes right is not refuted by one N. All of these counts are gated (V6), not hand-read. Whatever pins the band, it is not the count that pins the corner, even though the answer is the same number.

**The chiral intertwiner does not transport the corner's modes.** This was the obvious candidate and it is worth writing out, because its algebra is exactly right and it still fails. For free fermions,

  [H, Ξ(ρ)] = Ξ([H, ρ]) − (ε_a + ε_b)·Ξ(ρ),   Ξ(ρ) := d_a ρ d_b,

so Ξ commutes with the Hamiltonian part precisely on a chiral pair ε_a + ε_b = 0, and it shifts (p, q) by (−1, +1), one step of exactly the band's width. Chiral pairs exist precisely when h is bipartite, which is the gate; the open chain has ⌊N/2⌋ of them; and ρ ↦ d†_a ρ d_a commutes always and climbs the diagonal, so the reachable set is |p − q| ∈ {0, 2}. Every ingredient lines up.

It still fails. Applied to a basis of the corner's frozen subspace, the images are nonzero but are not in the (0,2) frozen subspace: the relative residual of L₍₀,₂₎ + 4γ̄ on them runs from 4·10⁻² to 7.5·10⁻¹, and the images span rank 10 to 18 against a frozen subspace of dimension 2 or 3. The verifier checks the map first, so this is not a coding slip: the intertwining identity holds to 3·10⁻¹⁵ on a chiral pair and reproduces the predicted commutator i·J·(ε_a + ε_b)·Ξ off it. **The map commutes with the Hamiltonian part and not with the dephasing, and the dephasing is what decides.**

So the count matches, the gate matches, the reachable set matches, and the map is still the wrong map.

## Where the number comes from, and what is still missing

The depth question has moved. It is no longer "why ⌊N/2⌋" but something sharper, and the sharpening is a result rather than a rewording. The (0,2) block consists **entirely** of cells with |AΔB| = 2, and exactly ⌊N/2⌋ of them are the balanced-pair cells |vac⟩⟨{l, R(l)}| whose rate diagonal is −2(γ_l + γ_{R(l)}) + 4γ̄ = 0 on the locus. One per balanced pair, again, and the same count.

The involution that makes M odd there is **antilinear**: reversal of the bra index composed with complex conjugation, since the block's Hamiltonian part is +iH₂ alone and conjugation is what turns its sign. That points at the antilinear corner of the object manager rather than at τQ. But an antilinear involution has equal-dimensional real forms, so it yields no shortage and cannot be the counter; and it survives the diagonal perturbation that kills the band, so it cannot be the gate either. It is a third structure sitting in the right place, doing neither job on its own. Hand it over as a fact, not as a lead.

**At J = 0 the count is manifest, and it is one per balanced pair.** Turn the coupling off and every band block is diagonal, so one can simply read which cells sit at the root: the entry is −2·Σ_{l ∈ AΔB} γ_l + 4γ̄, which vanishes exactly when the disagreement set is a balanced pair {l, R(l)}. Only |AΔB| = 2 cells qualify (a four-site set cannot sum to 2γ̄ on this locus), so the balanced pairs, and nothing else, put cells on the root.

**But their number is not ⌊N/2⌋, except in one block.** How many such cells a block holds depends on how much room there is beside the pair, and it varies enormously: at N = 7 the counts are 3, 6, 30, 15 and 60 for the blocks (0,2), (1,1), (2,2), (1,3) and (3,3). Every one of them collapses to exactly 3 the moment the coupling turns on. The (0,2) block is the exception where nothing departs, because a two-site bra beside an empty ket leaves no room: its J = 0 count is already ⌊N/2⌋.

**In that block the count is explained, but only in the limit.** The J = 0 kernel is exactly the span of the ⌊N/2⌋ balanced-pair cells, the frozen dimension stays ⌊N/2⌋ at every coupling, and as J shrinks the frozen subspace converges onto that span: the principal cosines reach 0.998, 0.988, 0.914 at N = 7, J = 1/256. So the pair cells index the modes **in the J → 0 limit**, which is where the number comes from. They are not what the modes look like at finite coupling: at J = 3/4 the same subspace has rotated almost entirely off them (cosines 0.059, 0.018, 0.002 at N = 7). Anyone re-deriving this should note that reading a null space basis out of a numerically zero matrix returns whatever basis the routine happens to emit; the pair alignment has to be measured against the pair-cell span, not read off a basis vector.

**What is missing is the collapse.** Everywhere else in the band, far more cells start on the root than end up frozen, and the survivors always number ⌊N/2⌋. Why a 60-cell start and a 3-cell start both land on 3 is the open question, and it is a better one than the count was: it asks about a departure rate, which is the same kind of object as the corner's own ladder J^{2d_c}. The natural next move is to measure the departure orders in a band block the way the proof measures them in the corner, and see whether ⌊N/2⌋ survivors is again a valuation statement.

Do NOT try to test the pair indexing by walking one pair off balance. It cannot be done: the pair deviations sum to zero, so no single pair can be unbalanced alone, and any perturbation moves γ̄ and with it the root, detuning every pair cell at once. The J = 0 reading above is the version of that experiment that is well posed.

## Verification

[`simulations/xy_frozen_band.py`](../simulations/xy_frozen_band.py), 67 checks, about 4 to 8 minutes, prints `XY frozen band: ALL GREEN`.

- **V1** the band law by full (p, q) census at N = 4..7: the carrier set, the whole-⌊N/2⌋ depth, the count 3(N−1), and the fold image. The script reports how many of its reads were exact.
- **V2** the band edge at N = 8, 9, 10, three profile-and-coupling combinations each.
- **V3** the |AΔB| = 2 oddness criterion, two-sided (it must fail on every other size class), and the combinatorial reach.
- **V4** the bipartiteness gate, the five rows of the table above.
- **V7** where the number comes from: the J = 0 census of cells on the root, block by block, and the continuation of the (0,2) pair cells into the J > 0 frozen subspace.
- **V6** the other falsified candidate as a count: the τQ fixed-cell census behind "the shortage does not extend", two-sided, since it must predict zero off the diagonal AND miss on it.
- **V5** the falsified intertwiner, with its own self-check first.

Multiplicities are exact GF(p) ranks at two primes ≡ 1 (mod 4) wherever the block is small enough, and SVD nullities above that, on the same dyadic grid the live witness uses. Never an eigenvalue count: the departing modes crowd the root at spacing J^{2d}, which is exactly where a floating-point spectrum goes quiet.

The band is also visible live, in the object manager: `inspect --root divisor --N 5 --chain xy` prints the full census beside the Heisenberg one ([`FrozenDivisorWitness`](../compute/RCPsiSquared.Diagnostics/Foundation/FrozenDivisorWitness.cs)), and gate block G2c of [`r90_frozen_divisor_gate.py`](../simulations/r90_frozen_divisor_gate.py) pins the carrier counts inside F140's own gate.

## Open

- **What pins the depth, now sharpened.** At J = 0 the number of cells on the root is one per balanced pair times whatever room the block leaves beside it, so it ranges from ⌊N/2⌋ to 60 at N = 7; at J > 0 every band block holds exactly ⌊N/2⌋. The open question is the COLLAPSE, not the count: why every starting number lands on the same survivor count. In the one block where nothing departs, (0,2), the modes are the pair cells continued and the number is understood.
- **Whether the edge argument is sufficient as well as necessary.** The section "Why the edge sits where it does" removes the known mechanism outside the band; it neither proves the null out there (that is measured) nor shows that every block inside must carry, and indeed on a non-bipartite h most do not.
- **A number, not yet a law:** the |AΔB| = 2 reach and the band coincide at N = 4..8. That is a combinatorial identity and should be provable outright rather than measured. Once it is, the whole edge argument stops being chain-specific evidence and becomes a lemma, at which point it belongs in [PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) beside the census it scopes, not in an experiment note.
- **No F number.** This is an experiment, deliberately. The gate is well controlled and the edge has an argument, but the carrier law's inside direction is measurement, and the mechanism for the depth is missing entirely. Minting now would name half an object.

## Anchors

[PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) (the corner, its census, and the Section 5 scope this note forced), [F140 in ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md), the MirrorWorld [`Divisor`](../compute/MirrorWorld/Divisor.cs) (which runs the XY chain and therefore adopts only the fold parity, not the census).
