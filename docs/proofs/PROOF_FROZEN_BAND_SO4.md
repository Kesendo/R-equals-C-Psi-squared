# What the frozen band contains: two ladders, an invariant disagreement, the seed rung in closed form, and a floor the other ladder sets

**Status:** Tier 1 derived, section by section, with what each imports named beside it.

- **Section 2**, both ladders, Σ-oddness as the exact condition for the SPIN ladder, and Lemma 2.5 deriving the off-diagonal floor from the corner's by injectivity. Corollary 2.4 records that Σ-oddness is sufficient for that floor and NOT necessary for the off-diagonal band to exist. The matching off-diagonal EQUALITY is the ceiling there and stays measured.
- **Section 3**, the band as ⌊N/2⌋ copies of one SO(4) irrep, derived as a FLOOR. The equality rests on the depth Section 7 bounds and Proposition 5.1 converts; j_spin ≤ 1 is imported from the measured band statement of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md).
- **Section 4**, F141, the invariant disagreement and the inheritance of a rung's spectrum from the rungs below. The criterion it yields is explicitly one-directional and lives on ker(ad_h), not on the full lowest-weight spaces.
- **Section 5**, the large-J reduction, every even order vanishing by transpose parity, and Proposition 5.3, which carries the bound from one coupling to all but finitely many by the pencil alone. Proposition 5.1 remains an upper bound at large J and claims nothing more. The finite exceptional set is not empty, and Section 5 says where it has been computed exactly.
- **Section 6**, F143, the seed rung's spectrum in closed form.
- **Section 7**, F144, the floor 𝒦 ≥ ℓ(N − ℓ)/(N + 1) on LW_ℓ ∩ V₀, from an exact identity whose whole deficit is the spin ladder. It supplies for all N at once the input the ceiling needs and that was previously computed one rung at a time, and with Proposition 5.1 it makes the depth exactly ⌊N/2⌋ in the LARGE-J regime at every N ≥ 6.

**Section 8 lists what is open**, and neither the all-N question nor the coupling question is on it any more. First is the exceptional set itself, whose size has no law and whose singlet character is proved only where the whole block has been enumerated; then the equality in the law above ℓ = 2; then the off-diagonal ceiling and the band statement, both unchanged.
**Date:** 2026-07-27
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the open XY chain under Z-dephasing, the band of blocks carrying the frozen root λ = −4γ̄ **contains** a single representation-theoretic object: ⌊N/2⌋ copies of the irreducible representation of Spin(4) = SU(2) × SU(2) with (η-spin, spin) = (N/2 − 1, 1), which descends to SO(4) only at even N, which places ⌊N/2⌋ in every one of the 3(N − 1) band blocks and is the proved floor. That the band is nothing more than those copies is the ceiling; its arithmetic input is proved for all N in Section 7, which with Proposition 5.1 settles the depth in the large-J regime and, through Proposition 5.3, at all but finitely many couplings; the exceptional couplings are real, and at one of them an extra mode joins the diagonal rungs, exactly one more and as a spin singlet wherever the whole block has been enumerated, which is why the equality is a generic statement and not a universal one. The two ladders that build it are the two commuting SU(2)s of the Hubbard model in disguise; the η one needs nothing of the Hamiltonian beyond number-conserving quadratic, and the spin one commutes exactly when the single-excitation matrix h is odd under the staggering Σ = diag((−1)^l). That condition is sufficient for the off-diagonal floor through the spin ladder and is not necessary for the off-diagonal band to be occupied, so it names a mechanism and does not replace the measured band-existence gate. The disagreement count is η-invariant, so it acts only on the multiplet-counting spaces, and the spectrum on rung p is the union of the lowest-weight spectra at rungs 0 ≤ ℓ ≤ min(p, N−p). At the seed rung the resulting operator is G = (1/M)(𝟏𝟏ᵀ + (I+R)/2) with M = N + 1 and 𝟏 the all-ones vector, whose spectrum is {0 with multiplicity ⌊N/2⌋, 1/M with multiplicity ⌈N/2⌉ − 1, 1 simple}: the frozen kernel is exactly the chiral-odd sector, and the gap above it is 1/(N+1). Above the seed the disagreement has a floor: on LW_ℓ ∩ V₀ it satisfies 𝒦 ≥ ℓ(N − ℓ)/(N + 1), by an exact identity in which every term subtracted from ℓ(ℓ+1) is a square and the leading one is the norm under the spin ladder, so the maximum of the double occupancy requires that BOTH ladders lower the state to nothing, and two further chiral conditions besides. For N ≥ 6 that floor exceeds 1, which is what the ceiling asks, and the single N where it lands exactly on 1 is N = 5.
**Verification:** [`simulations/eta_ceiling_reduction.py`](../../simulations/eta_ceiling_reduction.py) (must print "eta-ceiling reduction gate: ALL GREEN", 236 checks, about ten seconds; `--deep` adds the larger rungs and the wider sweeps, 286 checks in about four minutes, and carries F143 to N = 40) and, for Proposition 5.3 and the exceptional couplings, [`simulations/exceptional_couplings.py`](../../simulations/exceptional_couplings.py) (must print "exceptional couplings gate: ALL GREEN", 35 checks, about a minute; `--deep` adds N = 7 and the whole-block enumeration at N = 6)
**Depends on:** [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) (the frozen root, its corner, and the proved floor), [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) (the band |p − q| ∈ {0, 2}, the η ladder, and the per-N ceiling this document tries to make N-free)

---

## What this means

The frozen band was found one block at a time: a diagonal line of blocks all carrying the same undisturbable decay rate, plus two lines beside it, and a ladder that climbs the diagonal one excitation at a time. This document says the three lines are not three facts. They are one object seen from three angles, and the object is an old one: the chain under dephasing is the Hubbard model wearing a disguise, and a Hubbard model has two independent SU(2) symmetries rather than one. The disguise is worth naming precisely, because the two halves of a Hubbard model land in two different halves of a Liouvillian. The turning supplies only the hopping; what plays the part of the on-site repulsion is **the watching**, since on a diagonal block the dephasing is −4γ(ℓ − D̂) and so charges the double occupancy at rate 4γ. So this is not a Hubbard Hamiltonian with a dissipator attached. It is a generator whose kinetic term is the chain and whose interaction term is the measurement, which is why the two SU(2)s are exactly the ones a Hubbard model has and why the whole question below is about the double occupancy. The ladder already known is the first of them. The second is what steps sideways, off the diagonal onto the two neighbouring lines, and that is the whole reason the band is three lines wide and not one.

The second symmetry comes with a condition, and the condition is about the matrix rather than about its spectrum: the one-excitation matrix has to connect only sites of opposite parity. That is strictly stronger than having a spectrum symmetric about zero, which it implies and which does not imply it. What it buys is a mechanism: where it holds, the sideways ladder carries the corner's frozen modes onto the two neighbouring lines, and the floor there is the corner's floor for a reason rather than by measurement.

It does not, however, decide whether those lines are occupied. A star, one site joined to every other, fails the parity condition and fills them anyway. So the measured test the band note uses to predict occupancy is not superseded here, and the parity condition is one sufficient route to the floor rather than the gate itself. Getting this the wrong way round is easy, and this document did, until the star was asked.

There is a second thing the disguise gives, and it is the one that makes the remaining question small. Ask how much the two sides of a coherence disagree, count the sites where one has an excitation and the other does not, and that count does not notice the ladder at all: adding a matched pair to both sides adds no disagreement. So the whole dissipative half of the problem lives on the *count of multiplets* rather than on the multiplets themselves, and the question of how deep the freezing goes on a high rung reduces to the same question on the rungs below it. Everything the band does is then decided by what starts, and where.

At the bottom rung that decision can be written out completely. The operator that survives the reduction there is a two-line matrix built from the all-ones matrix and the mirror that swaps mode k with mode N+1−k, and reading its spectrum takes no more than noticing that the two commute. The frozen modes are exactly the combinations odd under that mirror, one per mirror pair, which is the ⌊N/2⌋ the whole arc has been counting, now in the mode picture rather than the site picture. And the next eigenvalue above them is 1/(N+1), which is not a chain constant but the length of the discrete sine transform the chain diagonalises under.

The rungs above the bottom one are what took the longest, and the answer came from the sideways ladder rather than from the climbing one. Ask how large the *agreement* between the two sides of a coherence can get, once the state is required both to start a multiplet and to be blind to the turning. Write that agreement out in the mode picture and it comes to one number, ℓ(ℓ+1), minus three quantities, each of which is a length squared and so cannot be negative. The largest of the three is how far the sideways ladder can still lower the state. So the agreement is greatest exactly where nothing can be lowered any further in either direction, and the disagreement, which is what the watching charges for, is correspondingly floored. That floor is above the frozen value at every chain of six sites or more, which is the whole ceiling in one sentence, where before it was one exact rank computation per rung. The one chain length where the floor lands exactly on the frozen value, rather than above it, is five.

All of that is an argument about how the chain behaves when the turning is fast. Carrying it to an ordinary coupling used to look like a formality, and the formality turned out to have contents. What carries it is not a limit at all but the shape of the problem: the block depends on the coupling linearly, and a linear family of matrices has one kernel dimension almost everywhere and a larger one on a finite set of exceptions. Since the floor holds at every coupling, a single coupling where the count is not too large already fixes the generic answer. The exceptions are then a question rather than a technicality, and asking it gives the surprise of this document: they are there. At isolated couplings, at ordinary values with no closed form, one more mode freezes. It is not one of the band's own; every one of them found so far carries no orientation in the sideways ladder, so it cannot leave the diagonal line it appears on, and the band the rest of this document describes is untouched by it. Whether that is true of every such mode at every chain length is open, and Section 8 says so. A chain at such a coupling holds one memory more than a chain at any coupling beside it.

## 1. Setting

Open chain of N sites, H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}), so the single-excitation matrix h is the N×N hopping matrix with h_{ab} = J on neighbours. Site-resolved Z-dephasing with rates γ_l, mean γ̄, and

  L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ).

Cells are |A⟩⟨B| with A and B site SETS; L preserves the joint-popcount block (p, q) = (|A|, |B|). The dephasing part is diagonal on cells with entry

  −2·Σ_{l ∈ AΔB} γ_l,

so the **disagreement set** AΔB is the only thing it sees. The frozen root is λ = −4γ̄, and the band is the set of blocks carrying it, |p − q| ∈ {0, 2} less the two blocks (0,0) and (N,N), as established in [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md).

**The doubled picture.** Read the ket index as one fermion species and the bra index as a second: a cell |A⟩⟨B| becomes a state with ↑-occupation A and ↓-occupation B. Right multiplication transposes, so H acting from the right becomes h^T on the ↓ species, and

  −i[H, ·]  becomes  −i·K,  K = Σ_{ab} h_{ab} c†_{a↑}c_{b↑} − Σ_{ab} h_{ba} c†_{a↓}c_{b↓}.

Write **𝒦** for the disagreement count as an operator: 𝒦 = Σ_l n_{l↑}(1 − n_{l↓}), which on a cell reads |A \ B|, and write 𝒦' for its mirror |B \ A|. At uniform rates the dissipator is −2γ(𝒦 + 𝒦'). **On a diagonal block, and only there**, the two agree, |A \ B| = |B \ A| = |AΔB|/2, and the dissipator collapses to −4γ·𝒦. On an off-diagonal band block it does not: a frozen cell there has |A \ B| = 2 and |B \ A| = 0, so the dissipator is −4γ while −4γ·𝒦 would be −8γ. Everything below that evaluates the spectrum of 𝒦 does so on a diagonal rung, where the collapse is valid.

**Symbols.** Collected here because the sections below use each other's, and a reader arriving at Section 5 should not have to hunt.

- **d_l, d†_l** are the Jordan-Wigner fermions of the chain, in which H is quadratic with coefficient matrix h; **c_{l↑}, c_{l↓}** are the two doubled species just introduced, ↑ carrying the ket index and ↓ the bra, with number operators **n_{l↑} := c†_{l↑}c_{l↑}** and likewise for ↓. The same two species in the MODE basis are written **f_{a↑}, f_{a↓}**, so that c_{l↑} = Σ_a v_a(l) f_{a↑} and likewise for ↓, the transform being real and orthogonal.
- **M := N + 1**, the length of the discrete sine transform the open chain diagonalises under. It is neither a popcount nor a prime. The modes of h are v_a(l) = √(2/M)·sin(πal/M) for a = 1, …, N, with energies ε_a = 2J·cos(πa/M); **P_a** is the spectral projector of h at ε_a, and **ā := M − a** is the **chiral partner** of mode a, the mode at −ε_a. (The mode index is written k where the chiral pairing is the point and a where it is a matrix index; they are the same index.)
- **α, β** denote ℓ-subsets of modes, so |α⟩⟨β| is a cell in the mode picture and **E_α := Σ_{a ∈ α} ε_a** is its Slater energy. The letters a, b, c, d are always SINGLE modes; α and β are always subsets, and E always carries a subset index where ε carries a single-mode one.
- **ad_h** is the map X ↦ [h, X], and **V₀ := ker(ad_h)**, spanned by the pairs (α, β) with E_α = E_β; **m_E** is the number of ℓ-subsets at energy E, so dim V₀ = Σ_E m_E².
- **LW_ℓ** always means the lowest-weight space of the **η ladder** Φ at rung ℓ, that is ker Ψ there, and never of the spin ladder S⁺; Section 2 introduces both, and Ψ := Φ†. Two indices name a rung and they agree where it matters: **ℓ** is the index of the multiplet that seeds it, **p** the popcount of the block it sits in.
- **[condition]** is the Iverson bracket, 1 when the condition holds and 0 otherwise, so [M | m] reads "M divides m".
- **Seed** here means the bottom rung of a ladder, (1,1) for the η one, together with the frozen modes P_a − P_ā that start there. The repo spends the same word on a different object, MirrorWorld's `Seed`, the within-block self-dual seed [F89](../ANALYTICAL_FORMULAS.md) holds as a count; neither statement depends on the other, and the word is kept because [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) already carries it in this sense.

## 2. The two ladders

Two maps on cells, both raising:

  **Φ**(ρ) := Σ_l d†_l ρ d_l = Σ_l c†_{l↑}c†_{l↓},   (p, q) ↦ (p+1, q+1)
  **S⁺**(ρ) := Σ_l (−1)^l d†_l ρ d†_l = Σ_l (−1)^l c†_{l↑}c_{l↓},   (p, q) ↦ (p+1, q−1)

Φ is Yang's η-pairing without the staggering, S⁺ is the spin raising with it. The swap of which one carries the stagger is the partial particle-hole transformation that turns this Liouvillian into a Hubbard Hamiltonian: the ↓ species hops with the opposite sign here, and that sign is exactly what moves the stagger from one SU(2) to the other.

**Lemma 2.1 (both are blind to the watching, at any profile).** Φ and S⁺ commute with the dissipator for every rate profile γ_l, and not only for the uniform one.

*Proof.* The dissipator is diagonal on cells with entry −2Σ_{l ∈ AΔB} γ_l, so it suffices that both maps leave AΔB unchanged. Φ acts on a cell only at a site l with l ∉ A and l ∉ B, so l ∉ AΔB before, and afterwards l lies in both A and B, so l ∉ AΔB again; no other site moves. S⁺ acts only at a site l with l ∉ A and l ∈ B, so l ∈ B \ A ⊆ AΔB before, and afterwards l ∈ A \ B ⊆ AΔB; again no other site moves. ∎

**Lemma 2.2 (Φ and the turning).** [K, Φ] = 0 for every single-excitation matrix h, symmetric or not.

*Proof.* [Σ_{ab} h_{ab} c†_{a↑}c_{b↑}, Σ_l c†_{l↑}c†_{l↓}] = Σ_{a,l} h_{al} c†_{a↑}c†_{l↓}, and the ↓ term contributes −Σ_{a,l} h_{la} c†_{l↑}c†_{a↓}, which after renaming (a ↔ l) is −Σ_{a,l} h_{al} c†_{a↑}c†_{l↓}. The two cancel. ∎

This is the ladder already carried by [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md); it is restated here because Lemma 2.3 is its sibling and reads best beside it. Note what the two lemmas together say: Φ commutes with the two halves of L **separately**, and for independent reasons.

**Lemma 2.3 (S⁺ and the turning).** Let Σ := diag((−1)^l). For real symmetric h,

  [K, S⁺] = Σ_{a,l} h_{al}·((−1)^a + (−1)^l)·c†_{a↑}c_{l↓},

which vanishes identically **iff h_{al} = 0 whenever a and l have the same parity**, that is, iff **Σ h Σ = −h**.

*Proof.* The ↑ term of K contributes Σ_{a,l} h_{al}(−1)^l c†_{a↑}c_{l↓}; the ↓ term contributes +Σ_{b,l} h_{bl}(−1)^l c†_{l↑}c_{b↓}, which after renaming is Σ_{a,l} h_{la}(−1)^a c†_{a↑}c_{l↓}. Adding and using h_{la} = h_{al} gives the displayed form. Each coefficient vanishes separately, and (−1)^a + (−1)^l is 0 for opposite parities and ±2 for equal ones. ∎

**Corollary 2.4 (what Σ-oddness does and does not decide).** Σ h Σ = −h implies that the spectrum of h is symmetric about zero, because Σ is an involution conjugating h to −h, and the converse fails: the **star**, one site joined to every other, is bipartite and so has a symmetric spectrum, but its two colours are centre-against-leaves rather than even-against-odd, so Σ h Σ ≠ −h and S⁺ does not commute there.

So Σ-oddness is **sufficient** for the off-diagonal floor through the spin ladder: where it holds, S⁺ carries the corner's frozen space onto the block (2,0) and the floor is inherited. It is **not necessary** for that block to be occupied. The star fills it regardless, carrying 2, 5 and 9 frozen modes at N = 5, 6 and 7 while failing Σ-oddness, and its (1,1) block carries far more than ⌊N/2⌋ besides, because its hopping matrix is highly degenerate.

Those three counts are for the **spin** sector, which is what L is built from, and a reader who reaches for free fermions will get different ones. The star's bonds are not adjacent in any site ordering, so the Jordan-Wigner strings do not cancel and the star is not a free-fermion model; pricing the same bonds by order instead of by occupation gives 4, 7 and 11 there, the Slater count 1 + C(N−2, 2). That is the arc's own pricing distinction turning up again, and it is worth stating because the free-fermion number is the one a quick recomputation produces. The measured band-existence test of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) is therefore not superseded by anything here, and Σ-oddness is not "the seed gate"; it is the condition for one mechanism.

This paragraph replaces an earlier one that had it backwards, and the way it went wrong is worth keeping. Gate check V1(b) verifies that Σ-oddness and S⁺-commutation agree, which is Lemma 2.3, and the star was added to it as the case separating Σ-oddness from the symmetric spectrum. That much is real. But S⁺-commutation is not band occupancy, and the check was then read as evidence for a claim about occupancy that it never touched. **A verification can be green, and be about the wrong thing.** V1(c) now asks the star the occupancy question directly, and the answer is the one above.

**Lemma 2.5 (S⁺ does not annihilate a seed, so the off-diagonal floor is derived and not measured).** Assume Σ h Σ = −h. Then S⁺ is injective on the frozen space of the corner block, and consequently each of (2,0) and (0,2) carries at least ⌊N/2⌋ frozen modes, at every N.

*Proof.* Three steps, each exact.

1. **S⁺ is the commutator with Σ.** On the corner block write X for the N×N matrix of the cell coefficients. S⁺ lands in the block (2,0), which is Λ²(C^N), and for **symmetric** X its image has Λ²-coordinates (S⁺X)_{ab} = [Σ, X]_{ab}. (Gate V2(b) pins this identity at residual exactly zero; it is the statement that acting with two creation operators antisymmetrises, and the staggering supplies the Σ.)
2. **The seeds are Σ-odd.** Σ conjugates h to −h, so it carries the ε_k eigenvector to an eigenvector of −ε_k, giving Σ P_k Σ = P_{k̄} for the chiral partner k̄. Hence Σ (P_k − P_{k̄}) Σ = −(P_k − P_{k̄}), and the same holds for every real combination of seeds, since the relation is linear.

   Worth spelling out, because it is the step a reader is most likely to doubt: these seeds are frozen at **every** coupling, not only in the large-J limit. Two facts do it, and neither mentions J. Each seed commutes with h, so the Hamiltonian part annihilates it; and each has zero diagonal, because chiral partners have identical squared amplitudes, so the rate part multiplies it by the single constant −4γ̄ that the two-disagreement cells carry. That the ⌊N/2⌋ of them are the WHOLE corner frozen space is a separate fact and comes from elsewhere: [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) Section 7 proves the corner multiplicity is exactly ⌊N/2⌋ at the uniform point for every J ≠ 0, and for all but finitely many J on any profile of the locus. Lemma 2.5 does not need that: it needs only that the seeds it names are symmetric and nonzero. The corner frozen space is accordingly a space of SYMMETRIC matrices at every coupling, which is what steps 1 and 3 need, and it is not the generic non-symmetric space one might expect from a non-Hermitian generator.
3. **An odd symmetric matrix has a nonzero commutator.** For Σ X Σ = −X and Σ² = I one has X Σ = −Σ X, so [Σ, X] = Σ X − X Σ = 2 Σ X, which vanishes only if X = 0 because Σ is invertible.

So S⁺X ≠ 0 for every nonzero X in the seed space, the image has dimension ⌊N/2⌋, and it consists of frozen vectors because S⁺ commutes with L. ∎

This is what upgrades the off-diagonal floor. Commutation alone gives only that S⁺ maps the corner's frozen space **into** (2,0); a seed annihilated by S⁺ would be a spin singlet and would contribute nothing there. Lemma 2.5 rules that out for all N, where gate check V2 only ever measured it at the N its census reaches. What stays measured is the matching **equality**, that (2,0) carries no more than ⌊N/2⌋, which is the ceiling and is not this section's.

## 3. The band as one SO(4) multiplet

**Lemma 3.1 (the seeds are η-lowest-weight).** Every frozen vector of the corner block (1,1) is traceless, hence annihilated by Ψ = Φ†.

*Proof.* Ψ maps (1,1) to the one-dimensional block (0,0) and is there proportional to the trace. The corner's frozen space is spanned by the differences P_k − P_{k̄} of spectral projectors of h at chiral partner pairs (Section 6 identifies them; each has trace 1 − 1 = 0). ∎

**Proposition 3.2.** Given the band statement |p − q| ∈ {0, 2}, the frozen space of the whole band is ⌊N/2⌋ copies of the single SO(4) irreducible representation with

  (j_η, j_spin) = (N/2 − 1, 1),  of dimension (2j_η + 1)(2j_spin + 1) = 3(N − 1),

(at odd N the η-spin is half-integer, so strictly this is a Spin(4) = SU(2) × SU(2) representation that does not factor through SO(4) = Spin(4)/Z₂; the name follows standard Hubbard usage and the dimension count, which is what the argument uses, is the same either way)

one dimension in each of the 3(N − 1) band blocks.

*Proof.* By Lemma 3.1 each of the ⌊N/2⌋ corner seeds is η-lowest-weight; its η-weight there is (2·1 − N)/2 = 1 − N/2, so its multiplet has j_η = N/2 − 1 and spans the diagonal rungs p = 1 up to p = N − 1, which is exactly the diagonal line of the band. For the spin: by Lemma 2.5 no seed is annihilated by S⁺, so none is a spin singlet and j_spin ≥ 1 (gate check V2 additionally measures that the image is all of the frozen space of (2,0), which is the matching ceiling and is not needed for the floor); and j_spin ≥ 2 would put frozen weight on a block with |p − q| = 4, which the band statement excludes, giving j_spin = 1. Counting dimensions, ⌊N/2⌋ copies of a 3(N − 1)-dimensional irrep fill each band block with ⌊N/2⌋, which is the measured depth. ∎

Two halves of this lean on measurement and both are named. First, j_spin ≤ 1 is the band statement, which [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) establishes by a census plus a mechanism-unavailability argument, not by a proof of the null. Second, and this is the sharper one: what the argument derives is the **floor**, that ⌊N/2⌋ copies are present. That each band block carries **only** ⌊N/2⌋, which is what makes "is" rather than "contains" the right verb in the proposition, is the ceiling; the proof imports it at the words "which is the measured depth". Read as a floor the proposition is derived; read as an equality it consumes the ceiling. What that import now costs has changed: Corollary 7.3 proves the arithmetic condition for every N and Proposition 5.3 carries it off the large-J limit, so the equality holds at all but finitely many couplings and no longer rests on a rank read per N. At an exceptional coupling it is the equality that fails, by one dimension and on the diagonal rungs alone, which is why Section 5 states the depth generically and this proposition should be read the same way.

## 4. F141: the disagreement is η-invariant, and what that buys

**Theorem 4.1 (F141).** [𝒦, Φ] = 0 and [𝒦, Ψ] = 0 exactly, at every rung. (𝒦 is the profile-independent integer count, so the quantifier over rate profiles has nothing to bite on here; it belongs to the stronger Lemma 2.1, that the DISSIPATOR commutes with both ladders at any profile.)

*Proof.* Lemma 2.1 for Φ; Ψ is its adjoint and the disagreement count is self-adjoint. In the fermion language the same statement reads 𝒦 = n_↑ − D̂ with D̂ = Σ_l n_{l↑}n_{l↓} the double occupancy, and [n_↑, Φ] = Φ and [D̂, Φ] = Φ cancel: adding a matched pair at an empty site raises the particle number by one and the double occupancy by one, so it raises their difference by nothing. ∎

**Corollary 4.2 (the spectrum on a rung is inherited from below).** Since Φ, Ψ and 𝒦 commute in the stated way, the space decomposes into sl(2) multiplets and 𝒦 acts only on the multiplicity spaces. Writing LW_ℓ for the lowest-weight space at rung ℓ,

  **spec(𝒦 on rung p) = ⋃_{0 ≤ ℓ ≤ min(p, N−p)} spec(𝒦 on LW_ℓ),  with multiplicities.**

Both limits are load-bearing. The upper one is min(p, N−p) and **not** p: a multiplet seeded at rung ℓ has η-spin N/2 − ℓ and so spans rungs ℓ through N − ℓ, which does not reach rung p once p > N − ℓ, so the union over ℓ ≤ p overcounts on the upper half of the ladder. The lower one is 0 and **not** 1: LW_0 is the vacuum cell, whose multiplet is the identity direction climbing every rung, and it contributes one mode at 𝒦 = 0. Dropping it leaves the union short by exactly one dimension at every rung. Gate check V3(b) pins both limits rung by rung, upper half included; it reads them on V₀, which is where the reduction lives, and the limits are the same on the full rung.

**What this does and does not reduce.** It does not yet say the ceiling is a statement about lowest weights, and two things stand in the way, both of which Section 5 removes rather than this section.

- **"Frozen" is not "𝒦 = 1" in general.** At finite coupling the frozen condition is (−i·ad_h + 4γ(1 − 𝒦))v = 0, and 𝒦 does not commute with ad_h, so a frozen mode is not usually a 𝒦-eigenvector at all. The two conditions coincide only on ker(ad_h), which is where the large-J reduction of Section 5 puts the question.
- **On the full lowest-weight spaces the criterion is vacuous.** 𝒦 is diagonal on cells with integer entries and preserves ker Ψ, so the eigenvalue 1 is present there with large multiplicity at essentially every N, including every N at which the ceiling holds. Read on the full spaces, "no LW_ℓ carries 𝒦 = 1" is simply false, and false in a way that would prove nothing if it were true.

So the usable form of the reduction is the one Section 5 supplies, on **V₀ := ker(ad_h)**, and it runs one way only:

  **no LW_ℓ ∩ V₀ with 2 ≤ ℓ ≤ min(p, N−p) carries the eigenvalue 1  ⟹  the ceiling holds at rung p.**

The converse fails, and N = 5 is the counterexample: there the reduced lowest-weight space at ℓ = 2 does carry 𝒦 = 1, the first-order certificate is not full rank, and the ceiling holds anyway, closed at third order. That is the same one-directional shape as Proposition 5.1, which bounds the multiplicity from above and never from below. The premise of that implication, that no such space carries the eigenvalue 1, is what Corollary 7.3 establishes for every N ≥ 6 and every ℓ ≥ 2.

Gate check V3(a) pins the commutators at machine zero on uniform and graded profiles; V3(b) pins the union identity.

## 5. The large-J reduction, and why the even orders are empty

Write the frozen question on a diagonal rung as a kernel. With h = J·h₀ and s = 1/J, and after dividing by J,

  **N(s) = A + s·B,  A = −i·ad_{h₀},  B = 4γ(1 − 𝒦),**

whose kernel for small s > 0 is the frozen space. A is anti-Hermitian and semisimple, so in the cell basis |α⟩⟨β| of the mode picture it is diagonal with entries −i(E_α − E_β), the projector **P_{V₀}** onto V₀ := ker A is orthogonal, and the reduced resolvent S is elementwise i/(E_α − E_β) off V₀. (Both index kinds appear in this section: α and β are ℓ-subsets carrying Slater energies E, while a, b, c, d below are single modes.)

**Proposition 5.1 (the reduction is an upper bound).** The number of eigenvalue branches of N(s) that vanish identically is at most the multiplicity of 0 in B₁ := P_{V₀} B P_{V₀} restricted to V₀.

*Proof.* Standard Rayleigh-Schrödinger with A semisimple: the branches emanating from 0 are s·λ_j(B₁) + O(s²), so a branch identically zero forces λ_j(B₁) = 0. ∎

That is an upper bound at large J and nothing more, which is all the proof gives and all the proposition claims. The rigorous per-N form of the same upper bound is the exact rank certificate of [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md): a full column rank over GF(q) forces a full rank over ℚ, which settles the non-existence outright and needs no limiting argument. Section 5 supplies the object on which that rank is read; Section 7 replaces the rank read, proving for all N and all ℓ ≥ 2 at once what the certificate confirms rung by rung.

What carries the step from one coupling to the others is not perturbation theory at all, and it is elementary once the pencil is written down.

**Proposition 5.3 (one coupling fixes the generic count).** Read the block at the frozen root as the pencil **M(J) = C + i·J·A₀**, with C the rate diagonal shifted by −4γ̄ and A₀ the Hamiltonian part at J = 1, both integer. Then:

1. dim ker M(J) is constant off a finite set of couplings and takes its minimum there. The entries of M are polynomials in J, so for r the rank over ℚ(i)(J), the rational functions of the coupling, the r×r minors are polynomials not all identically zero, the rank equals r off their common zero set, which is finite, and the rank can only drop on it.
2. Section 3's floor gives dim ker M(J) ≥ ⌊N/2⌋ at every J ≠ 0, so the generic count is at least ⌊N/2⌋.
3. The branches Proposition 5.1 counts are exactly the generic count, by the same finiteness, so with Corollary 7.3 the generic count is at most ⌊N/2⌋ for every N ≥ 6.

Hence **the depth is exactly ⌊N/2⌋ for all but finitely many couplings, at every N ≥ 6**, on a DIAGONAL block. The scope is worth saying twice: parts 1 and 2 hold on every band block, and only part 3 is diagonal, because Corollary 7.3 rests on the collapse of the dissipator to −4γ·𝒦 that Section 1 restricts to the diagonal. A single coupling with count ≤ ⌊N/2⌋ would do as well as a limit, and one is available at every N the exact GF(p) rank reaches, J = 1 included. On a side line the same rank supplies the input per N rather than for all N at once, so there the generic ceiling holds at each N where it has been read and the all-N statement stays the measured one of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md). ∎

**And "finitely many" is not a hedge.** The exceptional set is NOT empty. At N = 5, at the uniform watching point and J² = 3/2, the blocks (2,2) and (3,3) carry ⌊N/2⌋ + 1 frozen modes, as an exact rank over ℚ. At N = 5 the whole band has since been enumerated at every coupling at once, exactly: three exceptional couplings, on those two blocks and on no other, each a SIMPLE root of the deflated determinant and so raising the count by exactly one. At N = 6 the block (2,2) carries six, again by exact enumeration of the whole block and again simple. At N = 7 the number is 11, and there the exact polynomial is the one of the ℓ = 2 singlet subspace rather than of the whole block, so the count is exact while completeness and the multiplicity are not read. The extra mode is an η multiplet seeded at ℓ = 2 that is a spin SINGLET, so it occupies the diagonal rungs ℓ through N − ℓ and no side line, and neither the band statement nor the count on the two side lines moves. That last sentence is established over all couplings at N = 5, where the whole band was enumerated, and read at one coupling at N = 6; at larger N it is neither derived nor measured, and Section 8 lists it as open. [THE_EXCEPTIONAL_COUPLINGS](../../experiments/THE_EXCEPTIONAL_COUPLINGS.md) carries the certificates; gate [`simulations/exceptional_couplings.py`](../../simulations/exceptional_couplings.py) pins them. So no argument can promote this section's bound to every coupling, and the generic statement above is the sharp one.

That is what collapses the cost. The block has dimension C(N,p)², the reduced object has dimension dim V₀ = Σ_E m_E², the number of pairs of p-subsets of modes with equal Slater energy. At N = 12 on the middle rung that is 3584 against 853776.

**Proposition 5.2 (every even order vanishes).** Let T be the transpose involution (a, b) ↦ (b, a) on the operator basis. Then B is T-even and S is T-odd, so the n-th order reduced operator has T-parity (−1)^{n−1}. The frozen vectors are symmetric matrices, which is Lemma 2.5 step 2 for the corner seeds and is carried up the ladder by Φ, since Φ commutes with the transpose involution; hence they are of one T-parity, and a T-odd operator has zero matrix elements between two vectors of the same parity. Therefore **every even order vanishes identically on the surviving space**, and the next order after the first that can bite is the third.

*Proof.* In the mode basis the (a,b) row of B is R_{ab}[c,d] = (U^T diag(v_a) Dm diag(v_b) U)[c,d] with Dm the symmetric cell-level multiplier, so R_{ab}[d,c] = R_{ba}[c,d], which is T-evenness. S is elementwise i/(E_α − E_β) and so is manifestly T-odd. The n-th order operator is a product of n factors B and n − 1 factors S. ∎

Gate check V6 measures both halves: the second order is the zero matrix at every N tested, and where the first order is loose the third is what closes it.

## 6. F143: the seed rung in closed form

Let M := N + 1 and v_a(l) = √(2/M)·sin(πal/M) be the modes of h, with ε_a = 2J·cos(πa/M).

**Lemma 6.1 (the reduced operator at the seed).** On rung p = 1 the reduced operator P D̂ P is, in the mode basis, the matrix

  G_{ac} = Σ_l v_a(l)²·v_c(l)² = (W^T W)_{ac},  W_{lk} := v_k(l)²,

and since W is symmetric, G = W². G is positive semidefinite and never the zero matrix; the frozen vectors at the seed are the ones it annihilates, so what has to be found is ker G.

**Theorem 6.2 (F143).** Let R be the chiral involution a ↦ M − a on modes, and let **𝟏** be the all-ones vector in ℝ^N, so that 𝟏𝟏ᵀ is the all-ones matrix. (Written as an outer product on purpose: the letter this matrix usually carries is J, which in this document is the coupling.) Then

  **G = (1/M)·(𝟏𝟏ᵀ + (I + R)/2),**

and therefore, since 𝟏𝟏ᵀ and R commute,

  **spec(G) = { 0 with multiplicity ⌊N/2⌋,  1/M with multiplicity ⌈N/2⌉ − 1,  1 simple }.**

The kernel is exactly the R-odd sector, spanned by the chiral differences e_k − e_{M−k}, one per mirror pair. The gap above the frozen eigenvalue is **1/M = 1/(N+1)**, for N ≥ 3; at N = 2 the middle band is empty and the gap is 1.

*Proof.* Write W_{lk} = (2/M)sin²(πkl/M) = (1/M)(1 − cos(2πkl/M)); W is symmetric in k and l by inspection. Then

  G_{ac} = (1/M²)·Σ_{l=1}^{M−1} [1 − cos(2πal/M) − cos(2πcl/M) + cos(2πal/M)cos(2πcl/M)].

Use Σ_{l=1}^{M−1} cos(2πml/M) = M·[M | m] − 1. The first three terms give (M−1) + 1 + 1. The fourth is ½[(M·[M | a−c] − 1) + (M·[M | a+c] − 1)]. For a ≠ c and a + c ≠ M this is −1, giving G_{ac} = 1/M. For a + c = M it is M/2 − 1, giving 3/(2M). For a = c with 2a ≠ M it is likewise M/2 − 1, giving 3/(2M); for the self-paired mode a = M/2 both delta terms fire and the entry is 2/M. Those four cases are exactly the entries of (1/M)(𝟏𝟏ᵀ + (I+R)/2).

For the spectrum: R permutes coordinates and fixes 𝟏, so R·𝟏𝟏ᵀ = 𝟏𝟏ᵀ·R and the two are simultaneously diagonalisable. 𝟏𝟏ᵀ has eigenvalue N on 𝟏, which is R-even, and 0 on its orthogonal complement. (I+R)/2 is the orthogonal projector onto the R-even sector. Hence G is (N + 1)/M = 1 on 𝟏, (0 + 1)/M = 1/M on the R-even vectors orthogonal to it, and 0 on the R-odd vectors. The R-odd dimension is the number of 2-cycles of R on {1, …, M−1}, which is ⌊N/2⌋ (the mode a = M/2 is fixed and exists exactly when N is odd), and the remaining R-even dimension beyond the all-ones vector is ⌈N/2⌉ − 1. ∎

Two readings worth keeping. First, the kernel is the chiral-odd sector, so the ⌊N/2⌋ frozen modes at the seed are the antisymmetric combinations of chiral mode pairs, one per pair, with the self-paired middle mode at odd N contributing nothing. That is the same count the site picture calls one frozen mode per balanced pair, and the same seat left empty that [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) counts as the middle. Second, the 1/(N+1) is not a property of the chain but of the transform length: M = N + 1 is the length of the discrete sine transform, and the gap is one over it.

Gate check V5 pins the identity and the whole spectrum with multiplicities, and the deep run carries it to N = 40.

## 7. The floor on the disagreement, and with it the ceiling at every N

Sections 4 and 5 turn the ceiling into one question: does the compression of 𝒦 to LW_ℓ ∩ V₀ with ℓ ≥ 2 have the eigenvalue 1? This section answers it for all N at once, by bounding that compression from below. The operator that supplies the bound is the **other** SU(2): the spin ladder of Section 2, which until here had only carried the off-diagonal floor.

Throughout, v is a vector of LW_ℓ ∩ V₀, and

  **Y_{ab} := X_{ab} v,  X_{ab} := f_{b↓} f_{a↑}**

is what is left of it when one mode is taken off each side. Two abbreviations: **𝔖⁻ := Σ_b f†_{b↓} f_{b̄↑}**, with adjoint **𝔖⁺ = Σ_a f†_{a↑} f_{ā↓}**, is the spin ladder of Lemma 2.3 written in the modes, up to a global sign no norm sees (the staggering Σ is the site-side face of the chiral map a ↦ ā, since v_ā(l) = (−1)^{l+1} v_a(l), so Σ_l (−1)^l c†_{l↑}c_{l↓} = −Σ_a f†_{a↑}f_{ā↓}), and

  **Y⁺_{ab} := Y_{ab} + Y_{b̄ā}**

pairs each index pair with its image under the involution ι: (a, b) ↦ (b̄, ā), whose fixed points are exactly the chiral pairs a + b = M.

**The expansion.** Write D̂ = Σ_l η†_l η_l with η_l := c_{l↓}c_{l↑} the pair annihilator at site l, and expand the site operators in the modes, c_{l↑} = Σ_a v_a(l) f_{a↑} and likewise for ↓. Then η_l = Σ_{a,b} v_a(l)v_b(l) X_{ab}, so for every v

  ⟨v, D̂ v⟩ = Σ_{a,b,c,d} T_{abcd} ⟨Y_{ab}, Y_{cd}⟩,  **T_{abcd} := Σ_l v_a(l)v_b(l)v_c(l)v_d(l)**.  (7.1)

T is the four-mode overlap, and T_{abab} = G_{ab} is the seed-rung matrix of Section 6. Three facts about it are used below. Two are one line each: T is symmetric in all four indices, and T_{abcd} = T_{ab d̄ c̄}, because v_c̄ v_d̄ = v_c v_d. The third is the **selection rule**, and Lemma 7.1 rests on it, so it is worth writing out. Expanding each sine into exponentials of ζ := e^{iπ/M} and summing the resulting geometric series over l gives

  **2M·T_{abcd} = Σ_{σ_b, σ_c, σ_d ∈ {±1}} σ_bσ_cσ_d·[2M | a + σ_b b + σ_c c + σ_d d],**  (7.2)

which is zero whenever a + b + c + d is odd, the summand then being odd under l ↦ M − l. In particular T_{abcd} = 0 unless a ± b ± c ± d ≡ 0 (mod 2M) for at least one choice of the three signs. The same identity is what makes T an integer object; [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md) uses it that way and records the trap that the eight conditions are DIVISIBILITY by 2M and not vanishing of the signed sum.

**What the two constraints do.** V₀ grades: since v has matched Slater energies, Y_{ab} lies entirely in the part of rung ℓ − 1 whose ket-minus-bra energy is ε_b − ε_a, so ⟨Y_{ab}, Y_{cd}⟩ = 0 unless ε_b + ε_c = ε_a + ε_d. Lowest weight gives one linear condition: Σ_a Y_{aa} = Ψv = 0.

**Lemma 7.1 (the grade lemma).** Let 1 ≤ a, b, c, d ≤ N. If T_{abcd} ≠ 0 and ε_b + ε_c = ε_a + ε_d, then a = b and c = d, or (c, d) = (a, b), or (c, d) = (b̄, ā).

*Proof.* Write A := πa/M and likewise B, C, D, so that ε_x = 2cos(πx/M) and the energy condition reads

  cos A + cos D = cos B + cos C.  (★)

By the selection rule (7.2), T_{abcd} = 0 unless a ± b ± c ± d ≡ 0 (mod 2M) for at least one choice of the three signs. Since every index lies in [1, M − 1], the range of each signed sum leaves only the following possibilities, and they exhaust the eight sign patterns:

  (I) a + d = b + c.  (II) a + c = b + d.  (III) a + b = c + d.  (IV) a + b + c + d = 2M.
  (V) one angle is, modulo 2π, the sum of the other three.

Throughout, three range facts are used: 0 < A, B, C, D < π; |x − y| ≤ M − 2 for any two indices, so no half-difference reaches ±π/2; and no two indices sum to 2M.

(I) Sum-to-product on both sides of (★) gives cos((A+D)/2)·cos((A−D)/2) = cos((B+C)/2)·cos((B−C)/2), and a + d = b + c makes the two half-sum factors equal. If that common factor is nonzero, cos((A−D)/2) = cos((B−C)/2) with both arguments in (−π/2, π/2), so |a − d| = |b − c|, which together with a + d = b + c leaves (a, d) = (b, c), that is a = b and d = c, or (a, d) = (c, b), that is (c, d) = (a, b). If the common factor vanishes then (A+D)/2 = π/2, so a + d = M and hence b + c = M, which is c = b̄ and d = ā.

(II) Substitute D = A + C − B. Then (★) becomes cos(A + (C−B)/2)·cos((B−C)/2) = cos((B+C)/2)·cos((B−C)/2), and cos((B−C)/2) ≠ 0, so cos(A + (C−B)/2) = cos((B+C)/2). The two arguments agree up to sign and a multiple of 2π, which gives A = B, or A + C = 0, or one of these shifted, A = B ± 2π and A + C = ±2π. Only the first survives the range facts, and A = B is a = b, which with a + c = b + d makes c = d.

(III) Substitute D = A + B − C. The same manipulation gives cos(A + (B−C)/2) = cos((B+C)/2), whose resolutions are A = C, A + B = 0, A = C ± 2π and A + B = 2π; only A = C survives, and it is a = c and d = b, that is (c, d) = (a, b).

(IV) Here D = 2π − (A+B+C), so cos D = cos(A+B+C) and (★) reads cos((B+C)/2)·cos(A + (B+C)/2) = cos((B+C)/2)·cos((B−C)/2). If cos((B+C)/2) = 0 then b + c = M and a + d = M, which is (c, d) = (b̄, ā). Otherwise cos(A + (B+C)/2) = cos((B−C)/2), whose four resolutions are A + C = 0, A + B = 0, A + C = 2π and A + B = 2π, all excluded.

(V) If D ≡ A + B + C the equation is the one just solved in (IV), forcing b + c = M and then d = a + M or d = a − M, neither in [1, M − 1]. If A ≡ B + C + D, then (★) becomes cos((B+C)/2)·cos((B+C)/2 + D) = cos((B+C)/2)·cos((B−C)/2), and the same two branches force b + c = M, hence a = d + M or a = d − M, again out of range, or C + D = 2πk or B + D = 2πk, excluded. If C ≡ A + B + D, subtract instead: cos A − cos B = cos(A+B+D) − cos D gives sin((A+B)/2)·sin((A−B)/2) = sin((A+B)/2)·sin((A+B)/2 + D), and sin((A+B)/2) ≠ 0, so either B + D = 2πk, excluded, or A + D = π, that is a + d = M, which forces c = b + M or c = b − M, out of range. The remaining pattern B ≡ A + C + D is this one with B and C exchanged, an exchange (★) is invariant under. So (V) has no solutions at all. ∎

**Proposition 7.2 (the exact identity).** For every v in LW_ℓ ∩ V₀,

  **M·⟨v, D̂ v⟩ = ℓ(ℓ + 1)·‖v‖² − ‖𝔖⁻v‖² − ¼·Σ_a ‖Y⁺_{aa}‖² − ⅛·Σ_{a ≠ ā} ‖Y⁺_{a ā}‖².**

*Proof.* By Lemma 7.1 the only terms of (7.1) that survive are those with a = b and c = d, those with (c, d) = (a, b), and those with (c, d) = (b̄, ā). Sorting them, and writing g_a := Y_{aa},

  ⟨v, D̂ v⟩ = Σ_{a,c} G_{ac}⟨g_a, g_c⟩ + Σ_{a ≠ b} G_{ab}‖Y_{ab}‖² + Σ_{a ≠ b, a + b ≠ M} G_{ab}⟨Y_{ab}, Y_{b̄ā}⟩,

where the three sums are disjoint: for a = b the ι-partner (ā, ā) is already in the first sum, and for a + b = M the ι-partner is the pair itself and is already in the second. The coefficients are all G because T_{aacc} = Σ_l v_a²v_c² = G_{ac} for the first sum and T_{abab} = G_{ab} for the second, while for the third the chiral identity gives T_{ab b̄ ā} = T_{abab} = G_{ab}. By Theorem 6.2, M·G = 𝟏𝟏ᵀ + (I + R)/2, and the lowest-weight condition Σ_a g_a = 0 kills the 𝟏𝟏ᵀ term, leaving M·Σ_{a,c}G_{ac}⟨g_a,g_c⟩ = ½Σ_a(‖g_a‖² + ⟨g_a, g_ā⟩) = ¼Σ_a‖Y⁺_{aa}‖². In the remaining two sums M·G_{ab} = 1 off the chiral pairs and 3/2 on them, and reindexing the third sum by ι (which permutes its index set) pairs it with half of the second, so

  M·⟨v, D̂ v⟩ = ¼·Σ_a ‖Y⁺_{aa}‖² + ½·Σ_{a ≠ b, a + b ≠ M} ‖Y⁺_{ab}‖² + ⅜·Σ_{a ≠ ā} ‖Y⁺_{a ā}‖²,  (7.3)

using Y⁺_{a ā} = 2Y_{a ā} on the chiral pairs. It remains to identify the full ι-symmetrised square. Expanding,

  ½·Σ_{a,b} ‖Y⁺_{ab}‖² = Σ_{a,b} ‖Y_{ab}‖² + Σ_{a,b} ⟨Y_{ab}, Y_{b̄ā}⟩ = ℓ² ‖v‖² + ⟨v, Ξ v⟩,  Ξ := Σ_{a,b} X†_{ab} X_{b̄ā},

because Σ_{a,b} X†_{ab}X_{ab} = Σ_{a,b} n_{a↑}n_{b↓} = n_↑n_↓, which is ℓ² on the rung. For Ξ, anticommuting f_{ā↓} past f†_{b↓} in the product 𝔖⁺𝔖⁻ gives 𝔖⁺𝔖⁻ = n_↑ − Ξ, so ⟨v, Ξv⟩ = ℓ‖v‖² − ‖𝔖⁻v‖². Hence ½Σ_{a,b}‖Y⁺_{ab}‖² = ℓ(ℓ+1)‖v‖² − ‖𝔖⁻v‖², and subtracting (7.3) from it leaves exactly ¼Σ_a‖Y⁺_{aa}‖² + ⅛Σ_{a≠ā}‖Y⁺_{aā}‖². ∎

**Corollary 7.3 (the floor, and the ceiling for every N).** For every v in LW_ℓ ∩ V₀,

  **⟨v, 𝒦 v⟩ ≥ (ℓ(N − ℓ)/(N + 1))·‖v‖²,**

so the smallest eigenvalue of 𝒦 compressed to that space is at least ℓ(N − ℓ)/(N + 1), and a vector attaining it is lowest weight for the spin ladder as well as for the η ladder and satisfies Y_{aa} + Y_{āā} = 0 for every a and Y_{a ā} = 0 on every chiral pair. Consequently, for N ≥ 6 and every 2 ≤ ℓ ≤ ⌊N/2⌋ that compression does not have the eigenvalue 1.

*Proof.* Every term subtracted in Proposition 7.2 is a square, so M·⟨v, D̂ v⟩ ≤ ℓ(ℓ+1)‖v‖², and 𝒦 = ℓ − D̂ on the rung gives ⟨v, 𝒦 v⟩ ≥ (ℓ − ℓ(ℓ+1)/M)‖v‖² = (ℓ(N − ℓ)/M)‖v‖². The three equality conditions are the vanishing of the three subtracted terms. On 2 ≤ ℓ ≤ ⌊N/2⌋ the product ℓ(N − ℓ) increases with ℓ, so the binding case is ℓ = 2, where ℓ(N − ℓ) > N + 1 reads 2(N − 2) > N + 1, that is N ≥ 6. ∎

This is the statement the per-rung rank certificate of [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md) was computing one rung at a time, now proved once for all N and all ℓ. What it buys, and what it does not, is worth separating.

**What it buys, exactly.** With Corollary 7.3 the multiplicity of the eigenvalue 1 in the compression of 𝒦 to V₀ on a diagonal band rung can be written down rather than measured. By Corollary 4.2 that spectrum is the union over 0 ≤ ℓ ≤ min(p, N−p): the rung ℓ = 0 carries only 𝒦 = 0; the rung ℓ = 1 carries 𝒦 = 1 with multiplicity exactly ⌊N/2⌋, which is Theorem 6.2's kernel; and for N ≥ 6 no rung ℓ ≥ 2 carries it at all. So the multiplicity is exactly ⌊N/2⌋, and Proposition 5.1 then bounds the frozen space above by ⌊N/2⌋ at large J. Against the floor of Section 3, which holds at every J, that makes **the depth exactly ⌊N/2⌋ in the large-J regime, at every N ≥ 6, with no computation per N.** The threshold in J that "large" means is not made explicit here and may depend on N.

**What it does not buy.** Not every coupling, and not because the argument is weak. Proposition 5.3 carries the bound from large J to all but finitely many J, which is as far as it can be carried: the exceptional couplings exist, and at one of them the depth is ⌊N/2⌋ + 1. So the sharp statement is generic, the exceptions are a finite set of algebraic couplings with no closed form, and what joins there is a spin singlet on the diagonal rungs alone.

Two readings are worth keeping. The first is that **the ceiling belongs to the second SU(2)**. Section 2 introduced the spin ladder to carry the corner's frozen modes sideways, which is a floor; here the same operator returns as the leading exact deficit in the disagreement, and a state maximising the double occupancy must be one it annihilates, the two chiral conditions of Corollary 7.3 holding besides. The band's two ladders turn out to hold the two ends of the same statement, one the floor and one the ceiling. The second is that N = 5 is visible in the bound rather than around it: ℓ(N − ℓ) = N + 1 has, as [ETA_CEILING_REDUCTION](../../experiments/ETA_CEILING_REDUCTION.md) records, the single solution N = 5 at ℓ = 2, so the one N the criterion cannot settle is the one where the floor lands exactly on the frozen value.

Gate checks V10 to V12 pin the section: V10 verifies Lemma 7.1 by exhausting the quadruples, V11 verifies Proposition 7.2 as an operator identity on LW_ℓ ∩ V₀, matrix against matrix rather than spectrum against spectrum, and reads the resulting floor against the measured minimum (the objects it is built from are exact integers, while the basis of the lowest-weight space is numerical, so its residual is machine zero rather than zero), and V12 reads the saturating space, confirming that the bound is attained and that the three equality conditions of Corollary 7.3 hold on all of it.

## 8. What is not proved

**Which couplings are exceptional.** Proposition 5.3 closes the step this list used to open with: the depth is exactly ⌊N/2⌋ for all but finitely many couplings at every N ≥ 6, and the exceptional set is not empty, so nothing stronger is available. What is not proved is the exceptional set itself. Its size is exact per N and has no law (3, 6, 11 at N = 5, 6, 7); that EVERY exceptional mode is an ℓ = 2 spin singlet is proved at N = 5, where the whole band was enumerated at every coupling at once, and is not derived in general, so nothing yet forbids a higher rung from contributing at some larger N; whether the point is semisimple or defective is unread here, where [PROOF_R90_FROZEN_DIVISOR](PROOF_R90_FROZEN_DIVISOR.md) §9 answers defective for its own object, the corner on the R90 locus, whose jump is algebraic where this one is geometric. [THE_EXCEPTIONAL_COUPLINGS](../../experiments/THE_EXCEPTIONAL_COUPLINGS.md) carries the open list.

**The equality in the law.** Corollary 7.3 proves min 𝒦 ≥ ℓ(N − ℓ)/(N + 1), which is the direction the ceiling needs. That the bound is ATTAINED, and by how many states, is answered separately by F145 and F146: the saturating space is the space of rotation-invariant couplings of one spin-1 per chiral pair, of dimension C(⌊N/2⌋, ℓ)·R_ℓ with R the Riordan number, proved in [PROOF_SCALAR_COUNT](PROOF_SCALAR_COUNT.md) for M = N + 1 prime, M = 2p and M = 2^a and measured elsewhere. At ℓ = 2 the maximizers were already known in closed form and are checked per N by gate V9.

**The off-diagonal equality.** Lemma 2.5 derives the floor ⌊N/2⌋ on (2,0) and (0,2). That those blocks carry no MORE than that is the ceiling on the off-diagonal lines, and it stays measured, as [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) states for them.

**The band statement itself.** Section 3 uses |p − q| ∈ {0, 2} as an input. That statement is measured, with a mechanism-unavailability argument for its edge, in the band note; nothing here strengthens it.
