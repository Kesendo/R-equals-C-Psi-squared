# The R₉₀ frozen divisor: watching profiles that pin an eigenvalue for every coupling

**Status:** Theorem (lower bound) Tier 1 derived, and read as an index in Section 3.1 (the bound is two fixed-point counts subtracted, so it cannot depend on the coupling); cofactor closed form + tightness criterion, the latter for γ̄ ≠ 0, Tier 1 derived (Section 6); uniform-endpoint constants in closed form and tightness-for-generic-J at every N and every γ̄ ≠ 0 Tier 1 derived (Section 7, the two boundary clocks); the J-valuation lower bound, total and per pair, Tier 1 derived (Section 8, the distance ladder), with sharpness reduced by the pointed grading to a single nonvanishing (Section 8.1); the exceptional couplings real at N = 3 in closed form, and defective with a single 2×2 block at γ̄ ≠ 0, exact-computed at N = 3, 4, against a 3×3 block on the zero-mean stratum at N = 3 (Section 9)
**Date:** 2026-07-25
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** On the anti-palindromic watching locus (every reflection pair of dephasing rates sums to the same value), the single-excitation corner block of the Liouvillian carries the eigenvalue λ = −4γ̄ with multiplicity at least ⌊N/2⌋, for every Hamiltonian coupling J, and, when the mean rate γ̄ is nonzero, exactly ⌊N/2⌋ for all but finitely many J, at every N on the two chains of Section 1 (Section 7). At γ̄ = 0 the mirror's tax disappears: the same index then bounds the multiplicity below by N rather than ⌊N/2⌋, and N is what it is at every nonzero coupling measured (Section 6). The mechanism is a rank bottleneck of a cell-level mirror, not an invariant subspace and not a spectral symmetry. The nonvanishing cofactor is a single N(N−1)/2 determinant in closed form, (−1)^N(4γ̄)^⌈N/2⌉·det((X P_{O₊} X)|_{V₋}), whose nonvanishing is exactly tightness (Section 6).
**Verification:** [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) (must print "R90 frozen divisor gate: ALL GREEN", 302 checks, ~2 min)
**Depends on:** [PROOF_F91_GAMMA_NINETY_DEGREES](PROOF_F91_GAMMA_NINETY_DEGREES.md) (the R₉₀ reshuffle and its fixed locus), [GAMMA_FOLD_PAIR_OF_MIRRORS](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) (the X^N cross-dock used in the corollary), [PROOF_F139_SEAM_IDENTITY](PROOF_F139_SEAM_IDENTITY.md) (the sibling on the character side)

---

## What this means

Tune the watching so that mirrored sites share their rates evenly around the mean: the first and the last site together, the second and the second-to-last together, every reflection pair carrying the same total. On that one locus, and only there, the chain acquires eigenvalues that the Hamiltonian cannot move. However hard the sites talk to each other, at least ⌊N/2⌋ decay modes, generically exactly that many, sit frozen at one exact rate, one mode per balanced pair. They are not protected by a conserved quantity, and no subspace of states carries them; take any J and the eigenvectors have rearranged, yet the eigenvalue has not. What holds them is a bookkeeping fact: a mirror on the cells has more rooms it leaves alone than rooms it flips, and the surplus has nowhere to go.

The populations are what cuts that surplus down. The watching never touches them, so their rooms sit on the wrong side of the mirror and charge the count instead of feeding it: half the surplus goes, and what is left is one frozen mode per pair rather than one per site. Set the mean rate itself to zero, which this locus permits as soon as the rates are allowed to change sign, and the populations stop charging and start paying, so the count is not merely un-halved: at odd N the middle site's own room comes back with it, and what is frozen is N, one per site, exactly the rooms the mirror leaves alone.

There is a second reading, about how strongly they are held. Switch the coupling off entirely and twice as many modes sit at that rate; switch it on and half of them leave. (With the mean rate set to zero, twice as many leave, and they leave from a taller start, so what stays behind is one per site and not half of anything.) How fast each one leaves is set by nothing but a distance: the mirrored pair whose sites are d apart needs the coupling raised to the power 2d before its mode stirs, because the chain has to walk the excitation from the one site over to the other first. The outermost pair, the two ends of the chain, is the slowest to let go. Section 8 proves the lower half of that law by counting rooms and steps. And on a finite set of couplings one further mode joins the frozen ones, and there the frozen rate stops being a clean eigenvalue altogether: Section 9 pins that down.

F139 taught this arc that a wall can be a divisor instead of a symmetry. This is the same lesson on the home γ axis: a factor of the characteristic polynomial that divides out exactly on a locus, with no symmetry of the spectrum behind it.

## 1. Setting and definitions

Take the open chain of N qubits with an excitation-conserving Hamiltonian H(J) = J·H₁ whose single-excitation matrix h (the N×N matrix ⟨e_a|H₁|e_b⟩, where e_a is the computational state with the single excitation at site a) is real, symmetric, and invariant under the site reversal R: a ↦ N+1−a. The isotropic Heisenberg chain (h = hopping 2 on neighbours plus the ZZ diagonal) and the XY chain (hopping only) both qualify. Site-resolved Z-dephasing acts with rates γ_l (real; positivity is nowhere used in the proof, and the corollary in Section 4 exploits that, as does the zero-mean stratum γ̄ = 0 of Section 6, which is where the tightness statements stop and the count changes), σ = Σ_l γ_l, γ̄ = σ/N.

The Liouvillian L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ) preserves the joint-popcount blocks (popcount of the bra index, popcount of the ket index). The **corner block** is the (1,1) block, spanned by the cells |e_a⟩⟨e_b|, written v_{(a,b)}. On this block

  L_block(J) = J·K − 2Γ,  K v_{(a,b)} = −i Σ_c (h_{ac} v_{(c,b)} − h_{cb} v_{(a,c)}),  Γ v_{(a,b)} = (γ_a + γ_b)·v_{(a,b)} for a ≠ b, 0 for a = b.

Two structures on the cells:

- **The locus.** The R₉₀ reshuffle of [F91](PROOF_F91_GAMMA_NINETY_DEGREES.md) acts on the rate profile; its fixed-point set is the anti-palindromic class **γ_l + γ_{R(l)} = 2γ̄ for every l** (at odd N this forces the middle rate to the mean). All statements below hold on this locus.
- **The mirror.** τQ is the linear involution of cells **(a,b) ↦ (Rb, Ra)** (transpose composed with reversal on both sides). It splits the diagonal cells D = {(a,a)} and the off-diagonal cells O = {(a,b), a≠b} into ±1 eigenspaces D±, O±. Throughout, P_S denotes the orthogonal projector onto the span of the cell set or subspace S (P_D, P_{D₊}, P_{D₋}, P_{O₊}).

Three dimension counts, by inspection: the τQ-fixed cells in O are the anti-diagonal (a, Ra) with a ≠ Ra, so **dim O₊ − dim O₋ = 2⌊N/2⌋**; the diagonal cells pair (a,a) ↔ (Ra,Ra), so **dim D₋ = ⌊N/2⌋**; and Γ vanishes on D.

## 2. The two oddness lemmas

**Lemma 1 (K is τQ-odd): τQ K τQ = −K.**

*Proof.* (τQ K τQ v)_{(a,b)} = (K τQ v)_{(Rb,Ra)} = −i Σ_c [h_{Rb,c}(τQv)_{(c,Ra)} − h_{c,Ra}(τQv)_{(Rb,c)}] = −i Σ_c [h_{Rb,c} v_{(a,Rc)} − h_{c,Ra} v_{(Rc,b)}]. Substituting c → Rc and using h_{Rx,Ry} = h_{xy} (R-invariance) and h symmetric turns this into +i Σ_c [h_{ac} v_{(c,b)} − h_{cb} v_{(a,c)}] = −(Kv)_{(a,b)}. ∎

Consequence: K maps O₊ into D₋ ⊕ O₋. The D₊ component of Kv is zero automatically for v ∈ O₊.

**Lemma 2 (the recentered rate operator is τQ-odd on O, exactly on the locus): τQ (2Γ − 4γ̄)|_O τQ = −(2Γ − 4γ̄)|_O.**

*Proof.* On the cell (a,b), a ≠ b, the conjugated operator reads 2(γ_{Rb} + γ_{Ra}) − 4γ̄, and the locus gives γ_{Ra} + γ_{Rb} = 4γ̄ − γ_a − γ_b, so the entry is −(2(γ_a+γ_b) − 4γ̄). ∎

On the diagonal cells the recentered rate is the constant −4γ̄ on both sides, which is even, not odd; this is why the full-block identity τQ(L_block + 4γ̄)τQ = −(L_block + 4γ̄) + 8γ̄·P_D carries a defect exactly on D, and why the argument below works on O and treats D as a constraint. The defect is not small: the block spectrum is not palindromic about −4γ̄, and no multiset symmetry argument applies.

## 3. The theorem

**Theorem (frozen divisor).** On the R₉₀-fixed locus, for every J,

  **dim ker(L_block(J) + 4γ̄) ≥ ⌊N/2⌋.**

Equivalently: (λ + 4γ̄) raised to the ⌊N/2⌋-th power divides det(λ − L_block(J)) as a polynomial identity in J on the locus. For γ̄ ≠ 0 the bound is an equality for all but finitely many J, at every N and every locus profile, on the two chains whose leading coefficient Section 7 computes: that is the tightness theorem of Sections 6 and 7, which this section does not need. At γ̄ = 0 the bound is not an equality at any J, and the proof below says why: the ⌊N/2⌋ conditions it subtracts are there because the diagonal cells are even under the mirror while everything else is odd, and at γ̄ = 0 the constant that makes them even is zero (Section 6).

*Proof (the pencil argument).* Let W := {v ∈ O₊ : P_{D₋} K v = 0}, a subspace of O₊ cut by dim D₋ = ⌊N/2⌋ linear conditions, so dim W ≥ dim O₊ − ⌊N/2⌋. Take v ∈ W:

1. Kv has no D-part: the D₊ part vanishes by Lemma 1 (oddness), the D₋ part by the definition of W. Hence Kv ∈ O₋.
2. (2Γ − 4γ̄)v ∈ O₋ by Lemma 2 (v is supported on O and τQ-even; an odd diagonal operator maps O₊ to O₋).
3. Therefore (L_block(J) + 4γ̄)v = J·Kv − (2Γ − 4γ̄)v lies in O₋, for every J.

So Φ_J: W → O₋, v ↦ (L_block(J) + 4γ̄)v is a linear map into a space of dimension dim O₋ = dim O₊ − 2⌊N/2⌋, and

  dim ker Φ_J ≥ dim W − dim O₋ ≥ (dim O₊ − ⌊N/2⌋) − (dim O₊ − 2⌊N/2⌋) = ⌊N/2⌋.

Every kernel vector is an exact eigenvector of L_block(J) at −4γ̄. ∎

The multiplicity is a dimension bottleneck of the mirror: the fixed cells of τQ (the anti-diagonal, one per balanced pair plus its transpose) give O₊ a surplus of 2⌊N/2⌋ rooms over O₋, the D₋ constraint taxes away half, and the remainder must freeze. Note what the proof does not use: no invariant subspace (the frozen eigenvectors move with J, only the eigenvalue stands still; Section 10's first bullet), no spectral palindromy (Section 2), no diagonalizability assumption.

Two structural by-products, both verified at machine precision in the gate: the frozen eigenvectors carry **zero weight on the diagonal cells** and have **τQ-even O-part**; the J-dependence of the frozen modes is only the ⌊N/2⌋-dimensional kernel line of Φ_J rotating inside the fixed subspace W.

### 3.1 What the bound is: the mirror's orbit count, subtracted

The paragraph above names the mechanism but not its genre, and the genre is worth naming, because it says in one word why the eigenvalue does not move.

**The surplus is the mirror's trace.** For any involution acting by permuting a basis, the even eigenspace exceeds the odd one by exactly the number of fixed basis vectors: a 2-cycle contributes one even and one odd combination and cancels, and only a basis vector that is **its own image** stands alone, on the even side. So

  dim O₊ − dim O₋ = #{τQ-fixed cells in O} = #{(a, R(a)) with a ≠ R(a)} = 2⌊N/2⌋,

which is the count Section 1 makes by inspection, read as what it is. The fixed cells are the coherences between a site and its own mirror partner, one per site (at odd N the centre site's is a diagonal cell and never enters O, which is the whole difference between N and 2⌊N/2⌋). The tax is the neighbouring number on the other cell set, and it counts the opposite thing: τQ acts on the diagonal cells by (a,a) ↦ (R(a), R(a)), whose only fixed cell is the odd-N centre, so dim D₋ = ⌊N/2⌋ is the count of its **2-cycles** there. Fixed points on one cell set, 2-cycles on the other: the frozen multiplicity is the mirror's own orbit bookkeeping, subtracted, and nothing else enters.

**The bound is an index.** One qualification first, because it is the point: dim W itself is *not* operator-free, since W is cut out by P_{D₋}K v = 0 and K sits in that condition. What is operator-free is the ESTIMATE. W is cut by at most dim D₋ conditions whatever K happens to be, and Φ_J lands in O₋ whatever J happens to be, so the chain dim ker ≥ dim W − dim O₋ ≥ (dim O₊ − dim D₋) − dim O₋ is assembled entirely from three dimensions of the mirror's own eigenspaces. Deform the operator however you like, and as long as it stays odd those three numbers are untouched and the bound cannot move. **That is why the eigenvalue is independent of J**: not a protection, not a symmetry, a difference of dimensions. It also explains, in one stroke, the features Sections 8 and 10 catalogue separately: the eigenvectors may rearrange because only the kernel's dimension is fixed, not its position; no palindromy of the spectrum is needed because an index is a count and not a symmetry; partial balance yields nothing because off the locus the operator is not odd at all, so the hypothesis is absent rather than weakened; and the exceptional couplings of Section 9 can be defective because an index bounds a dimension and says nothing whatever about Jordan structure.

**How little of the physics enters.** The gate makes this concrete by taking the physics away: a **random** τQ-odd matrix, with no Hamiltonian, no chain and no rates in it, still has kernel exactly ⌊N/2⌋ once the even defect on the diagonal cells is added, at N = 3..8. Without that defect the kernel is exactly N, the full fixed-cell count. And an operator that is **not** odd has kernel 0. So oddness is the entire hypothesis **of the counting argument**, which is not to say it comes for free: Lemmas 1 and 2 are exactly the work of buying it, the first from the R-symmetry of h and the second from the locus, and neither is automatic. What the chain supplies beyond that is only the value −4γ̄ and the fine structure of Sections 6 to 9. The freezing itself is older and cheaper than any of it.

**The tax is a lever.** The gap between N and ⌊N/2⌋ is the even defect, and it exists for a plain physical reason: dephasing does not act on populations, so the diagonal cells carry rate zero and their recentered value is the constant −4γ̄, even under the mirror rather than odd. A generator that made the diagonal cells odd too would freeze N modes instead of ⌊N/2⌋, one per site instead of one per pair. The lever has exactly one setting inside the present generator, and it is the degenerate one: at γ̄ = 0 the recentered value is zero, which is odd as well as even, and the count is indeed N (Section 6). The setting is bought either with negative rates or, at its one degenerate point γ ≡ 0, with no watching at all, so what it establishes is that the arithmetic works; it makes no claim about a dissipating channel. What a *channel* that pulls the lever would have to be is not addressed here; it is recorded in Section 12 as the sharpest open question this reading raises.

## 4. Corollary: the antidiagonal corners, by the gamma fold

The block (1, N−1) (single excitation against N−1 excitations) is the image of the corner block under the one-sided X^N bridge, and the [gamma fold](../../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) turns that bridge into algebra: conjugating by right-multiplication with X^N sends L(γ) to L(−γ) − 2σ, writing γ = (γ₁, ..., γ_N) for the whole site profile (profiles carry no arrow here). The theorem applied at the profile −γ (the locus is preserved; the root becomes −4·(−γ̄) = +4γ̄) then gives, for every J,

  **the (1, N−1) block carries the eigenvalue 4γ̄ − 2σ with multiplicity ⌊N/2⌋,**

and the same for (N−1, 1) by transposition. Composing the two one-sided bridges, one on each side, sends the root through r ↦ −r − 2σ twice and lands back on −4γ̄, so **(N−1, N−1) carries −4γ̄ with multiplicity ⌊N/2⌋** as well. Four corner blocks in all, the root selected by the parity of how many one-sided folds separate the block from (1,1).

The two roots −4γ̄ and 4γ̄ − 2σ sum to −2σ: they are partners under the pair of mirrors, and recentered at the palindrome center (x = λ + σ) they sit at x = ±(σ − 4γ̄). At N = 4 the two coincide at the center x = 0. What once looked like two separate root families, including the N = 3 block (1,2) sighting at −2γ̄, is one theorem and its fold image: at N = 3, (1,2) is the antidiagonal corner, and 4γ̄ − 2σ = −2γ̄ there.

## 5. Why only the corners, and on which chain

The proof needs one affine recentering that makes the rate operator odd (Lemma 2). Under reversal a rate −2γ_S (S the disagreement set of the cell) goes to −2(2γ̄|S| − γ_S): the center depends on the size |S|. The corner block has a single off-diagonal size class (|S| = 2), so one center serves; the block (2,2) has classes |S| ∈ {0, 2, 4} and admits no single center. That is an explanation of why the proof stops at the corners, not a proof that nothing else freezes; what closes the gap is measurement. The gate censuses **every** joint-popcount block (p,q) at N = 4, 5, 6 against both roots, and exactly the four corner blocks of Section 4 carry ⌊N/2⌋, each at the root its fold parity selects, while all (N+1)² − 4 remaining blocks carry nothing at either root.

The confinement carries the hypothesis Section 6 isolates, γ̄ ≠ 0, and when it goes the reason the corners were special goes with it. What stops the recentering at the corners is the |S|-dependence of the center, and the center is −2γ̄|S| (the corner's |S| = 2 is where the root −4γ̄ comes from): at γ̄ = 0 it is zero for every size class at once, so no block needs a center of its own and the argument of Section 3 runs on **every** block. Two mirrors serve, one for each diagonal of the block grid. Write R for the site reversal acting on a whole computational index, that is on the bit string rather than on one site (on a single-excitation index this is the R of Section 1, so the corner sees no new map). On (p,p) the cell mirror is τQ itself, (i,j) ↦ (R(j), R(i)), whose fixed cells are the C(N,p) cells with j = R(i); on (p, N−p) it is that map composed with the two-sided X^N bridge, (i,j) ↦ (R(j̄), R(ī)) with the bar for bitwise complement, again with C(N,p) fixed cells. Both are exactly odd there and neither is odd at γ̄ ≠ 0 (entry-wise, as residues mod p, so "odd" means identically zero; both chains, N = 3, 4, 5), so the index of Section 3.1 gives

  **at γ̄ = 0, every block with q = p or q = N − p carries the root with multiplicity at least C(N,p),**

the corner's N being the p = 1 case. The two roots of Section 4 have collapsed onto one by then, since σ = Nγ̄ vanishes with γ̄, so there is only one root left to carry. Measured rather than derived: that these are the ONLY carriers, and that the bound is attained. On the Heisenberg chain both hold at N = 4, 5, 6, where the description gives nine carriers of twenty-five, twelve of thirty-six and thirteen of forty-nine (exact GF(p) ranks). Section 12 records what is left open.

**That census is a Heisenberg statement, which is why the title of this section carries a second clause.** The measurement above runs on the Heisenberg chain. Drop the ZZ diagonal and the confinement is gone: on the XY chain, still at γ̄ ≠ 0, the two roots are carried at the same multiplicity ⌊N/2⌋ by many more blocks (counting a block that carries either: nine of twenty-five at N = 4, twenty-four of thirty-six at N = 5, twenty-one of forty-nine at N = 6, of which nine, twelve and fifteen carry the unfolded root itself), never partially, and every block carrying the unfolded root −4γ̄ has p + q even (necessary, not sufficient: plenty of even blocks carry nothing). The divisor bound of Section 3 and the fold-parity split of Section 4 survive the change of chain; only the confinement does not, and what confines is the ZZ term, though not in the way this paragraph first read it (the end of this section). This is the one place where the two chains this document treats as interchangeable are not, and it was found by the Object Manager witness named at the end of Section 12 rather than by the gate, because the witness reads a chain the gate's census never ran. Gate G2c now pins all three parts of it: the corners still carry, strictly more blocks carry, and every carrier carries the whole ⌊N/2⌋. What the XY chain carries instead of four corners is a diagonal band, and what decides between the two is not one condition but two independent ones: the off-diagonal half of the band needs h to have a bipartite spectrum, while the diagonal half needs only that H be quadratic in the fermions, which the ZZ term destroys and an R-invariant diagonal on h does not. So the confinement to the corners is the ZZ term's after all, though not as a diagonal on h: it is the quarticity, which removes the symmetry that carries the corner's frozen modes up the diagonal. That the Heisenberg band is then empty rather than merely unexplained stays, as above, a measurement. [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) states that law, gives the reason its edge sits at |p − q| ≤ 2, identifies that symmetry, and records which candidate mechanisms are dead.

Away from the R₉₀ locus the whole structure disappears: partial balance (all but one condition satisfied) yields nothing, and the N = 3 closed form shows the defect as an explicit linear factor (M₍₁,₂₎ below is the (1,2) block of L, the N = 3 antidiagonal corner of Section 4),

  det(−2γ₂·I − M₍₁,₂₎) = 512·J⁴·(γ₁+γ₃)²·(γ₁+γ₃−2γ₂)·(4J² + (γ₁+γ₃)(γ₁+γ₃−2γ₂)),

so the balanced root exists exactly on the locus (given J ≠ 0), and the distance of the nearest eigenvalue grows linearly in the defect.

## 6. The cofactor: closed form, tightness, semisimplicity

The pencil argument of Section 3 bounds the multiplicity from below; this section computes what is left of the determinant after the frozen factor divides out, and the answer closes the tightness question along the way.

Write M̃ := L_block(J) + 4γ̄ for the recentered block, and split off its τQ-even part: M̃ = X + 4γ̄·P_D (this X is an operator on the corner block, not the many-body Pauli X^N of Section 4), where

  X = J·K − 2Δ,  Δ v_{(a,b)} = (δ_a + δ_b)·v_{(a,b)} on O, 0 on D,  δ_l := γ_l − γ̄.

X is τQ-odd (Lemmas 1 and 2; on the locus δ_{R(l)} = −δ_l), and X is **γ̄-free**: the mean rate enters M̃ only through the even defect 4γ̄·P_D. Let V₊ = D₊ ⊕ O₊ and V₋ = D₋ ⊕ O₋ be the parity eigenspaces of τQ, and note dim V₋ = ⌊N/2⌋ + dim O₋ = N(N−1)/2.

**Theorem (cofactor).** On the R₉₀-fixed locus, the characteristic polynomial p(ε) = det(εI − M̃) factors as ε^⌊N/2⌋·q(ε) with

  **q(0) = (−1)^N · (4γ̄)^⌈N/2⌉ · det( (X P_{O₊} X)|_{V₋} ),**

a determinant of size N(N−1)/2 whose entries are free of γ̄. Consequently, **for γ̄ ≠ 0**, the frozen multiplicity is exactly ⌊N/2⌋ if and only if this determinant is nonzero, and in that case the frozen eigenvalue is semisimple (algebraic = geometric = ⌊N/2⌋). Section 9 takes up the other case, where the determinant vanishes; the next bullet takes up the prefactor, which is the other way the cofactor can die.

*Proof.* An odd operator exchanges the parity blocks, and an even one preserves them, so in the split V₊ ⊕ V₋ the matrix M̃ − εI has diagonal blocks Λ₊ = diag(4γ̄−ε on D₊, −ε on O₊) and Λ₋ = diag(4γ̄−ε on D₋, −ε on O₋), and off-diagonal blocks X₊₋, X₋₊ (the diagonal rate part −2Δ of X is odd, hence off-diagonal here). For ε ∉ {0, 4γ̄} the Schur complement gives

  det(M̃ − εI) = det(Λ₊) · det( Λ₋ − X Λ₊⁻¹ X |_{V₋} ),  Λ₊⁻¹ = (4γ̄−ε)⁻¹ P_{D₊} − ε⁻¹ P_{O₊}.

Pull the pole at ε = 0 out of the second factor: with A(ε) := Λ₋ − (4γ̄−ε)⁻¹·X P_{D₊} X |_{V₋} (regular at ε = 0),

  det(M̃ − εI) = (4γ̄−ε)^⌈N/2⌉ · (−1)^{dim O₊} · ε^{dim O₊ − dim V₋} · det( X P_{O₊} X |_{V₋} + ε·A(ε) ).

The exponent is dim O₊ − dim V₋ = ⌊N/2⌋ (Section 1's counts), det(XP_{O₊}X + εA(ε)) is regular at ε = 0, and both sides are rational functions equal off finitely many points, hence equal as rational functions; the left side is a polynomial, so the identity extends to ε = 0. This re-proves the divisor bound (order ≥ ⌊N/2⌋ at ε = 0, an independent second proof of the theorem of Section 3), and reading off the coefficient of ε^⌊N/2⌋, with p(ε) = (−1)^{N²}·det(M̃ − εI) (the block is N²×N²), gives q(0) = (−1)^{N² + dim O₊}(4γ̄)^⌈N/2⌉ det(XP_{O₊}X|_{V₋}). The sign: dim O₊ = N(N−1)/2 + ⌊N/2⌋, and N² + N(N−1)/2 + ⌊N/2⌋ ≡ N (mod 2) for every N (check both parities of N). Tightness: q(0) ≠ 0 says the algebraic multiplicity is exactly ⌊N/2⌋; the pencil gives geometric ≥ ⌊N/2⌋, and algebraic ≥ geometric always, so all three coincide. ∎

Four consequences, all pinned in the gate:

- **The γ̄-stratification is exact.** The cofactor is (4γ̄)^⌈N/2⌉ times a polynomial in J and the offset profile δ = (δ₁, ..., δ_N) only: the whole γ̄-dependence of the residual spectrum at the frozen root is the even defect's ⌈N/2⌉ diagonal cells of D₊.
- **The bottom of that stratification is a second regime, not a gap.** The prefactor says where the cofactor dies for a reason that has nothing to do with J, and the locus reaches that place: γ̄ = 0 costs the rates their positivity as soon as any one of them is nonzero, and Section 1 never asked for positivity anyway. There the even defect 4γ̄·P_D vanishes outright, so M̃ = X is τQ-odd on the **whole** corner block, diagonal cells included, and the counting argument of Section 3 runs with no tax at all: an odd operator exchanges the parity eigenspaces, so its kernel is at least dim V₊ − dim V₋ = trace(τQ) = N. That bound is attained: the geometric multiplicity is N on both chains at every nonzero coupling measured (exact GF(p) ranks, N = 3..7, Section 11). The exceptions run parallel to the taxed stratum's rather than disappearing with the cofactor. J = 0 is the familiar one, where the block is diagonal and the kernel is the larger N + 2⌊N/2⌋ (Section 8's Corollary), and there are isolated nonzero couplings besides, where the ALGEBRAIC count rises and the root goes defective: Section 9 gives one in closed form, and its Jordan block is larger than the taxed stratum's. This is the same arithmetic Section 3.1 already ran with the physics taken away, met here inside the theorem's own hypotheses rather than on a random matrix: one frozen mode per **site** instead of one per pair. Note that the two counts do not differ by a factor: the diagonal cells stop charging their ⌊N/2⌋ and, at odd N, contribute the one further room ⌈N/2⌉ − ⌊N/2⌋ themselves, so ⌊N/2⌋ becomes N and not 2⌊N/2⌋. Every "exactly ⌊N/2⌋" below therefore carries γ̄ ≠ 0 as a hypothesis; the lower bound of Section 3 is untouched, since ⌊N/2⌋ ≤ N.
- **Tightness for generic J.** For a fixed locus profile, det(XP_{O₊}X|_{V₋}) is a polynomial in J of degree N(N−1) with leading coefficient det((K P_{O₊} K)|_{V₋}), which is nonzero for every N: Section 7 proves this by computing it in closed form, and the gate scans N = 3..10 on both the Heisenberg and the XY single-excitation matrix as a check. Hence, for γ̄ ≠ 0, the multiplicity is exactly ⌊N/2⌋ for all but finitely many J, at every N and every locus profile.
- **Small N in closed form** (symbolic, in the antisymmetric coordinates δ₁ = γ₁ − γ̄, δ₂ = γ₂ − γ̄):

    N = 3: q(0) = 2¹²·γ̄²·J⁴·(3J² − δ₁²)
    N = 4: q(0) = 2²⁰·γ̄²·J⁸·(8J⁴ − 4J²(3δ₁² + 2δ₁δ₂ + δ₂²) + (δ₁² − δ₂²)²)
    N = 5: q(0) = −2³⁰·γ̄³·J¹²·Q₅(J², δ₁, δ₂) with Q₅ of total degree 8 in (J, δ₁, δ₂), leading term 25J⁸

  (The XY chain differs only in the polynomial coefficients: its J-pure terms are 2J² and 5J⁴ at N = 3, 4.) The J-power in front is **2⌊N/2⌋⌈N/2⌉ = 2⌊N²/4⌋** (4, 8, 12, 18, 24 at N = 3..7), not the 4(N−2) the first three N suggest: the two expressions coincide exactly through N = 5 and split at N = 6 (exact-rational J → 0 valuations 18 and 24 at N = 6, 7; the N = 6 discriminator is gate-pinned). Mode-resolved, the valuation ladder is a distance ladder: the balanced pair (c, Rc) contributes one branch of the D₋ Schur complement with valuation J^{2(N+1−2c)}, twice the pair's site distance, because the Hamiltonian must walk the excitation from c to Rc before that frozen mode moves off its perch. [Section 8](#8-the-valuation-law-the-walk-to-the-anti-diagonal) proves the ladder as a lower bound, total and per pair; the first rung is also available directly, since the order-J² Schur complement on D₋ is −C_anti†C_anti with C_anti = P_anti K P_{D₋} (P_anti the projector onto the anti-diagonal cells) and rank C_anti = 1 at even N, 0 at odd N: only the distance-1 middle pair of even N is one hop from the anti-diagonal. From the theorem of this section only the total degree N(N−1) and the leading coefficient follow.

## 7. The uniform endpoint: the two boundary clocks

At the fully degenerate point of the locus, δ = 0 (uniform watching), X = J·K and the cofactor collapses to a J-monomial:

  q(0)|_{δ=0} = (−1)^N·(4γ̄)^⌈N/2⌉·J^{N(N−1)}·D_N,  D_N := det( (K P_{O₊} K)|_{V₋} ).

This section computes D_N in closed form, for every N. Three lemmas reduce it to a Gram determinant, and the Gram determinant is a pure power of the chain's clock modulus.

**Lemma 3 (the pair basis diagonalizes K²).** Let λ_1 < ... < λ_N and u_1, ..., u_N be the eigenvalues and orthonormal eigenvectors of h, with R-parities π_i (R u_i = π_i u_i; the spectrum is nondegenerate, see Lemma 5). Writing w_{ij} := u_i ⊗ u_j for the cell-space product vectors, τQ w_{ij} = π_iπ_j·w_{ji}, so the vectors φ_{ij} := (w_{ij} − π_iπ_j w_{ji})/√2 over i < j form an orthonormal basis of V₋ (count: N(N−1)/2 ✓), and K² φ_{ij} = −(λ_i−λ_j)²·φ_{ij}.

**Lemma 4 (the Gram reduction).** K = −i·(a real symmetric matrix) is anti-Hermitian, K† = −K, which is what turns −K P_{D₊} K into +C†C below. K is τQ-odd, so P₊ K|_{V₋} = K|_{V₋} (P₊ the projector onto V₊), hence K P_{O₊} K = K² − K P_{D₊} K = K² + C†C on V₋, with C := P_{D₊} K|_{V₋}. By the Weinstein-Aronszajn identity,

  D_N = (−1)^{N(N−1)/2} · Π_{i<j} (λ_i−λ_j)² · det( I − W ),  W := C Λ⁻² C†,  Λ² := diag((λ_i−λ_j)²).

Moreover I − W = V†V with V_{i,k} := ⟨w_{ii}, d_k⟩ (d_k the D₊ basis): K φ_{ij} = −i(λ_i−λ_j)·ψ_{ij} with ψ_{ij} := (w_{ij} + π_iπ_j w_{ji})/√2, so the (λ_i−λ_j)² of Λ⁻² cancels and W_{kl} = Σ_{i<j} ⟨d_k, ψ_{ij}⟩⟨ψ_{ij}, d_l⟩, and {w_{ii}} ∪ {ψ_{ij}} is an orthonormal basis of V₊ containing d_k, so W_{kl} = δ_{kl} − Σ_i ⟨d_k, w_{ii}⟩⟨w_{ii}, d_l⟩. Since the w_{ii} are τQ-even (orthogonal to D₋), the N×N Gram matrix G_{ij} := ⟨w_{ii}, P_D w_{jj}⟩ = Σ_a u_i(a)²u_j(a)² equals VV†, and therefore

  det(I − W) = det(V†V) = pdet(G)  (the product of the nonzero eigenvalues of G, provided rank G = dim D₊ = ⌈N/2⌉, which Lemma 5 gives).

**Lemma 5 (the two boundary clocks).** The single-excitation eigenbasis of the open chain is a cosine angle lattice, with modulus M depending on the chain (not the running clock of [ClockHandLadder](../../compute/RCPsiSquared.Core/Symmetry/ClockHandLadderClaim.cs), whose two hands are the turning ω and the tick γ; here "clock" is the timeless lattice of angles the excitation lives on):

- **Heisenberg (M = N).** λ_k = 4cos(kπ/N) + N − 5 and u_k(a) ∝ cos((2a−1)kπ/(2N)), k = 0..N−1 (the DCT-II basis). The check is two lines: the interior rows force the dispersion, the left end row holds identically for the half-integer cosine (the chain's +2 boundary defect from the ZZ diagonal is exactly what the product-to-sum identity 2cosθ·cos(θ/2) = cos(3θ/2) + cos(θ/2) absorbs), and matching the right end quantizes sin(Nθ) = 0.
- **XY (M = N+1).** λ_k = 4cos(kπ/(N+1)) and u_k(a) ∝ sin(akπ/(N+1)), k = 1..N (the DST-I basis, classical).

Both spectra are nondegenerate (distinct cosines), so Lemma 3 applies. For B_{a,k} := u_k(a)² everything reduces to geometric sums: in the Heisenberg case Σ_{k=0}^{N−1} cos(jkπ/N) = 1 for odd j and N·[2N | j] for even j (the sum's period in j is 2N, not N: at N = 4, j = 4 the sum reads 1 − 1 + 1 − 1 = 0, and the N | j form would fire at every even N with j = N), and, since u_k(a)² contributes the doubled angle cos((2a−1)kπ/N), the products of those cosines split into sums with j = 2(a−b) and j = 2(a+b−1), whose resonances are exactly a = b and b = Ra; the XY case runs the same sums at modulus N+1. The result is one law for both chains:

  **B Bᵀ = (1 − 1/M)·𝟙𝟙ᵀ/N + (I + R)/(2M).**

B is doubly stochastic (B𝟙 = Bᵀ𝟙 = 𝟙), the anti-symmetric space is killed, and on the R-symmetric space (I+R)/2 acts as the identity: the spectrum of BBᵀ, hence of G = BᵀB, is exactly

  {1 (the flat vector), 1/M with multiplicity ⌈N/2⌉ − 1, 0 with multiplicity ⌊N/2⌋},

so rank G = ⌈N/2⌉ and **pdet(G) = M^{−⌊(N−1)/2⌋}**. ∎

**Theorem (uniform-endpoint constants).** For every N ≥ 3, on the Heisenberg (M = N) and XY (M = N+1) open chains,

  **D_N = (−1)^{N(N−1)/2} · Π_{i<j}(λ_i−λ_j)² · M^{−⌊(N−1)/2⌋} ≠ 0.**

Two consequences:

- **Tightness is now a theorem at every N** (upgrading the N = 3..10 scan of Section 6), at every γ̄ ≠ 0: D_N ≠ 0 makes det((XP_{O₊}X)|_{V₋}) a nonzero polynomial in J for every locus profile (leading coefficient J^{N(N−1)}·D_N), so the frozen multiplicity is exactly ⌊N/2⌋ for all but finitely many J; at the uniform point itself the determinant is the monomial J^{N(N−1)}·D_N, so there the multiplicity is exactly ⌊N/2⌋ for **every** J ≠ 0. The γ̄ that D_N is free of is the same γ̄ the cofactor's prefactor carries, and it is the prefactor, not this determinant, that fails at the zero-mean stratum.
- **The J-pure constants of Section 6 are clock numbers.** Π(λ_i−λ_j)² is the discriminant of the (integer) characteristic polynomial of h, and the correction is a pure power of the clock modulus. N = 3 Heisenberg: disc = 2304, M^{−1} = 1/3, D₃ = −768 = −2⁸·3, reproducing the 3J² of Section 6; the gate pins the exact rational assembly at N = 4, 5 (at the gate's point γ̄ = 9/100, J = 4/3, the N = 5 uniform cofactor is exactly −2⁶⁴/2989355625).

The two moduli are the two boundary clocks of the open chain: the XY excitation lives on the Dirichlet lattice sin(akπ/(N+1)), which is the committed SE cosine lattice of modulus N+1 ([NivenRationalityRoot](../../compute/RCPsiSquared.Core/Symmetry/NivenRationalityRootClaim.cs) and the [F65 registry entry](../ANALYTICAL_FORMULAS.md) live there, as do the F129 collision combs and the F139 wall at modulus 11); the Heisenberg excitation lives on the Neumann half-integer lattice cos((2a−1)kπ/(2N)) of modulus N, which had no committed anchor before this proof. The frozen divisor's residual constant reads off which boundary clock the chain carries.

## 8. The valuation law: the walk to the anti-diagonal

Section 6 read the J → 0 valuation of the cofactor determinant off the data and called the per-pair ladder an observation. This section derives it. The single geometric input is that h is tridiagonal (nearest-neighbour hopping plus a diagonal; both chains of Section 1 qualify), and the quantity the answer is written in is the site distance of the balanced pair,

  d_c := N + 1 − 2c,  c = 1, ..., ⌊N/2⌋,  with Σ_c d_c = ⌊N/2⌋·⌈N/2⌉ = ⌊N²/4⌋.

**Theorem (valuation law).** On the R₉₀-fixed locus, for every profile and every N,

  **ord_J det( (X P_{O₊} X)|_{V₋} ) ≥ 2·Σ_c d_c = 2⌊N²/4⌋,**

and if the profile is generic (δ_a + δ_b ≠ 0 for every off-diagonal cell off the anti-diagonal), the ⌊N/2⌋ × ⌊N/2⌋ Schur complement S(J) of that matrix on D₋ satisfies the sharper per-pair form

  **ord_J S_{c,c'} ≥ d_c + d_{c'},  that is  S(J) = Λ·S̃(J)·Λ with Λ := diag(J^{d_c}) and S̃ regular at J = 0.**

Both hold with equality wherever they have been measured: the total at N = 3..10 for both chains, the per-pair orders at N = 3..7, all integer-exact (Section 11).

The proof is a Gram form, a grading, and a counting argument.

**Lemma 6 (the cofactor determinant is a Gram determinant).** In the cell basis h real symmetric makes K complex symmetric (K^T = K), Δ is diagonal, hence X^T = X, and the parity projectors are real. Writing Y := P_{O₊}X|_{V₋} as a matrix from an orthonormal basis of V₋ to one of O₊, of size dim O₊ × dim V₋ with dim O₊ − dim V₋ = ⌊N/2⌋,

  (X P_{O₊} X)|_{V₋} = Y^T Y,  so by Cauchy-Binet  det( (X P_{O₊} X)|_{V₋} ) = Σ_I det(Y_I)²,

the sum running over the maximal square minors of Y (row sets I of size dim V₋). Any other basis of O₊ serves as well: an unnormalized row basis turns the identity into Y^T D Y with D the diagonal of inverse squared norms, and Cauchy-Binet then carries those positive weights along, which is the form Section 11 uses. In particular ord_J of the left side is at least 2·min_I ord_J det Y_I. Cancellation between minors can only raise the left side, so it cannot damage the bound. This is not a positivity statement: Y is complex and the squares are squares of complex numbers. The Gram shape is used only to halve the bookkeeping, so that a count made once on Y pays twice.

**Lemma 7 (gradings, and the level as one of them).** Call an integer function φ on the cells **admissible** if it is τQ-invariant and changes by at most one under a single hop, that is |φ(a±1,b) − φ(a,b)| ≤ 1 and |φ(a,b±1) − φ(a,b)| ≤ 1. An admissible φ is constant on every parity basis vector, so it is defined on the basis vectors of D±, O±, V±, and the entries of Y obey

  Y_{y,x} = 0 whenever |φ(y) − φ(x)| ≥ 2,  ord_J Y_{y,x} ≥ 1 whenever φ(y) ≠ φ(x).

The **level**

  **ℓ(a,b) := |s(a,b)|,  s(a,b) := a + b − (N+1),**

the distance of the cell to the anti-diagonal {b = R(a)}, is admissible.

*Proof.* The entry statements need only admissibility. Δ is diagonal in the cell basis, so it preserves every cell and hence every φ; the J-linear part of X is J·K, which couples (a,b) to (c,b) and to (a,c) with amplitudes h_{ac} and h_{cb}, zero unless the site indices are neighbours or equal, so one power of J moves one site index by one, and an admissible φ by at most one. Since X is affine in J, cells whose φ differ by two or more are uncoupled at every order, and cells with different φ are uncoupled at order J⁰. That the level is admissible: τQ sends (a,b) to (R(b), R(a)) and therefore s to −s, so ℓ is τQ-invariant, and one hop moves s by one, hence |s| by at most one (if s and s' have opposite signs then |s − s'| = |s| + |s'|, so |s − s'| ≤ 1 already forces ||s| − |s'|| ≤ 1). ∎

The Heisenberg ZZ diagonal of h contributes same-cell terms only, so it is level-diagonal and transports nothing. That is why the two chains, whose h differ exactly by that diagonal, carry the same valuation.

**Lemma 8 (the level census).** Write n_j := dim(V₋ ∩ level j) and m_j := dim(O₊ ∩ level j). Then

  n_0 = 0,  m_0 = 2⌊N/2⌋,  and for j ≥ 1:  n_j = N − j,  m_j = N − j − η_j,

where η_j := 1 if j = d_c for some balanced pair c, and 0 otherwise.

*Proof.* Level 0 is exactly the anti-diagonal {(a, R(a))}, N cells, each of them τQ-fixed and hence in V₊: n_0 = 0, and of the N cells the odd-N centre cell is diagonal, leaving m_0 = 2⌊N/2⌋ in O₊. For j ≥ 1 the level consists of the N − j cells with s = +j and the N − j cells with s = −j; τQ maps the first set bijectively onto the second with no fixed cell, so V₊ and V₋ each receive N − j dimensions. A diagonal cell (a,a) has s = 2a − N − 1, so it sits at level j ≥ 1 exactly when j = |2a − N − 1|, that is j = d_c with a = c or a = R(c); those two cells are exchanged by τQ and contribute one dimension to D₊ and one to D₋. Subtracting the D₊ dimension from V₊ ∩ level j leaves m_j. ∎

The census is the whole mechanism in two lines. At each of the ⌊N/2⌋ distances d_c the columns outnumber the rows by exactly one, because the balanced pair puts a diagonal cell there. All the spare rows sit at level 0, on the anti-diagonal, where there are 2⌊N/2⌋ of them and no columns at all. These are the same τQ-fixed rooms whose surplus froze the eigenvalue in Section 3, seen from the other side: there they made the kernel, here they are the only place the excess can go.

**Lemma 9 (the transport bound).** Fix an admissible grading φ. Let A be a set of columns and I a set of rows with |I| = |A|, write Y_{I,A} for the square submatrix they cut out, and put

  F_j(I, A; φ) := #{columns of A with φ ≥ j} − #{rows of I with φ ≥ j}.

Then ord_J det Y_{I,A} ≥ Σ_{j≥1} max(0, F_j(I, A; φ)). Taken with φ = ℓ this pays as follows:

- for A the full column set, ord_J det Y_{I,A} ≥ Σ_{j≥1} Σ_{i≥j} η_i = Σ_i i·η_i = Σ_c d_c = ⌊N²/4⌋;
- for A the O₋ columns together with the single D₋ column of one pair c, ord_J det Y_{I,A} ≥ d_c.

*Proof.* Expand det Y_{I,A} over the bijections π from A onto I (a permutation here, not the R-parities π_i of Lemma 3). By Lemma 7 a term is nonzero only if |φ(π(x)) − φ(x)| ≤ 1 for every column x, and its J-order is then at least the number of columns with φ(π(x)) ≠ φ(x). Fix j ≥ 1 and let E_j := {x ∈ A : φ(x) ≥ j, φ(π(x)) ≤ j−1}. Every column with φ ≥ j outside E_j goes to a row with φ ≥ j, so #{columns of A with φ ≥ j} ≤ #{rows of I with φ ≥ j} + |E_j|, that is |E_j| ≥ F_j(I, A; φ). The one-step constraint forces φ(x) = j exactly for x ∈ E_j, so the sets E_1, E_2, ... are pairwise disjoint, and each of their elements is one mismatch: the order of the term is at least Σ_j |E_j| ≥ Σ_j max(0, F_j(I, A; φ)).

For the two specializations, bound the rows of I at level ≥ j by all the rows there, Σ_{i≥j} m_i, and use Lemma 8. With A the full column set the level counts are the n_i, so F_j ≥ Σ_{i≥j}(n_i − m_i) = Σ_{i≥j} η_i ≥ 0, and summing over j counts each diagonal level i once for every 1 ≤ j ≤ i. With A the O₋ columns plus the D₋ column of the pair c, every diagonal column except that one is dropped, so the count at level i ≥ 1 falls to n_i − η_i + [i = d_c] and F_j ≥ #{i ≥ j : i = d_c}, which sums over j to d_c. ∎

*Proof of the theorem.* The first display is Lemmas 6 and 9. For the second, let A be the O₋ columns together with the D₋ column of the pair c, and B the same with c'. In its rectangular form Cauchy-Binet reads det((Y^TY)_{A,B}) = Σ_I det(Y_{I,A})·det(Y_{I,B}), and Lemma 9's second specialization bounds the two factors by d_c and d_{c'} for every I. Hence the bordered determinant of T := (X P_{O₊} X)|_{V₋} has ord_J ≥ d_c + d_{c'}. For a generic profile the unbordered block T_{O₋,O₋} is a unit at J = 0: its value there is the Gram matrix of the vectors Δw over the O₋ basis, diagonal with entries proportional to (δ_a + δ_b)², nonzero exactly under the stated genericity. So the Schur complement entries S_{c,c'} = det T_{A,B} / det T_{O₋,O₋} inherit the order. Expanding det S over permutations then re-derives the total, every term having order at least Σ_c (d_c + d_{π(c)}) = 2Σ_c d_c. ∎

**The reading.** The frozen mode of the balanced pair (c, R(c)) is carried by a diagonal cell, and a diagonal cell sits at level d_c: as far from the anti-diagonal as the two sites are from each other. The census says that level owes one room, and that the only spare rooms in the whole block are on the anti-diagonal itself. So the debt has to be walked down to level 0, and the walk is the Hamiltonian's: one power of J is one hop, and no hop moves more than one level. The cost of the pair is its site distance; the Gram form charges it twice, once on each side of the determinant. This is the distance-is-t frame inside a determinant: what the excitation would need in time to cross from c to R(c) is what the determinant needs in powers of J before that frozen mode can move.

**Corollary (what the order counts).** At J = 0 the recentered block M̃ is diagonal, carrying −2(δ_a + δ_b) on the cell (a,b) and 4γ̄ on the diagonal cells, so on a generic locus profile with γ̄ ≠ 0 it vanishes exactly on the 2⌊N/2⌋ anti-diagonal cells off the centre: there the frozen root has twice its generic multiplicity (gate, N = 3..6). Exactly ⌊N/2⌋ modes therefore leave the root as the coupling turns on, and since q(0) collects, up to sign, every eigenvalue of M̃ that is not frozen (Section 6), while only the departing ones among them vanish as J → 0, the theorem says the J-orders of the departing modes sum to at least 2Σ_c d_c. Under the measured equality the split is one departing mode per balanced pair, at order J^{2d_c}: the outermost pair, whose two sites are the ends of the chain, is the one that clings to the frozen root longest. On the zero-mean stratum the same count runs one step higher: the diagonal cells carry zero there too, so the J = 0 kernel is N + 2⌊N/2⌋, and 2⌊N/2⌋ modes leave rather than ⌊N/2⌋, which is what has to happen for the N of Section 6 to be what stays.

Two limits of the statement are worth naming. The theorem says only that the mode cannot move sooner; that it does move exactly then is the other half, equivalent to S̃(0) being nonsingular, and it is what the equalities of Section 11 record. Section 8.1 takes that half apart and leaves a single nonvanishing standing. And the only place the chain geometry enters is Lemma 7, through the tridiagonal h: one power of J is one hop. Longer-range hopping crosses several levels per power and the ladder degrades accordingly, so the distance ladder is a statement about the chain, not about the mirror.

### 8.1 Sharpness: the pointed grading, and where each pair can land

The theorem bounds the valuation from below; that it is attained is the other half, and it comes out of running Lemma 9 a second time with a different grading. The level ℓ measures the distance to the anti-diagonal as a whole. Aim instead at one anti-diagonal cell. For a site x with x ≠ R(x) put

  **χ_x(a,b) := |a − x| + |b − R(x)|,**

the walking distance from the cell (a,b) to the single cell (x, R(x)). (Not the basis vectors ψ_{ij} of Lemma 4; this is a function on cells.)

**Lemma 10 (χ_x is admissible, and what it counts).** χ_x is admissible in the sense of Lemma 7: one hop moves one site index by one, so χ_x moves by at most one, and τQ, which sends (a,b) to (R(b), R(a)), exchanges the two summands, so χ_x ∘ τQ = χ_x. Let

  I_x := every O₊ row except the anti-diagonal ones other than (x, R(x)),

a set of the right size, and let A_c be the O₋ columns together with the D₋ column of the pair c, as in the theorem. Then for every j ≥ 1

  **F_j(I_x, A_c; χ_x) = [χ_x(c,c) ≥ j],  hence  ord_J det Y_{I_x, A_c} ≥ χ_x(c,c) = max(d_c, d_x),**

where d_x = |N + 1 − 2x| is the site distance of the pair through x.

*Proof.* The census is trivial in this grading because rows and columns cancel almost everywhere. Every τQ-orbit of off-diagonal cells that is not fixed contributes exactly one O₊ row and one O₋ column, and χ_x, being τQ-invariant, gives them the same value; those contributions cancel in F_j. What remains on the column side is the one D₋ column, at height χ_x(c,c), and on the row side the anti-diagonal cells kept in I_x, of which there is exactly one, namely (x, R(x)) at height 0, which no longer counts once j ≥ 1. So F_j is 1 while j ≤ χ_x(c,c) and 0 after, and Lemma 9 sums it to χ_x(c,c). The last equality is the walk itself: χ_x(c,c) = |x − c| + |R(x) − c| is the sum of the distances from the site c to the two sites of the pair through x, which is d_c when x lies between c and R(c) and d_x otherwise, that is max(d_c, d_x) in both cases. (Check the two readings against each other at N = 7: for the pair c = 2, taking x = 3 gives 1 + 3 = 4 = d_c, while x = 1 gives 1 + 5 = 6 = d_x.) ∎

Two consequences, and they are the sharpness.

**The staircase.** Write ĝ_I(c) for the coefficient of z^{d_c} in det Y_{I,A_c}, z := −iJ as in Section 12's reduction (the two readings differ by the unit i^{d_c}, which moves no valuation and no rank). Lemma 10 says ĝ_{I_x}(c) = 0 whenever d_x > d_c, since then the order max(d_c, d_x) already exceeds d_c. So each pair can only land on the anti-diagonal cells that lie between its own two sites, and because those intervals are nested, so are the supports:

  the pair at distance d reaches exactly the cells (x, R(x)) with c ≤ x ≤ R(c), a set that grows with d.

**Independence, hence tightness.** Order the pairs by decreasing distance and take for each c the row set I_c that keeps the anti-diagonal cell (c, R(c)), the outer end of its own interval. For any tighter pair c′, that cell lies outside the interval, so ĝ_{I_c}(c′) = 0, while ĝ_{I_c}(c) ≠ 0. The matrix ĝ_{I_c}(c′) is therefore triangular, and wherever its diagonal does not vanish the ⌊N/2⌋ vectors ĝ(c) are linearly independent, so by the reduction of Section 12 the leading Schur matrix S̃(0), a positive semidefinite Gram matrix built from them, is nonsingular and the valuation law holds with equality, total and per pair. The next paragraph says what carries that diagonal.

One ingredient of that last paragraph is measured rather than derived: the diagonal entries ĝ_{I_c}(c) ≠ 0, equivalently that the bound of Lemma 10 is attained on a pair's own outer cell. It is exact at N = 3..8 for every pair (Section 11), and the reason it should hold is visible: the walk from (c,c) to (c, R(c)) leaves one index fixed and moves the other straight across, so the cheapest route is unique and there is nothing for it to cancel against. Making that uniqueness an argument is what stands between this section and a two-sided theorem, and it is now a single nonvanishing statement rather than the nonsingularity of a determinant.

## 9. Where tightness fails: the exceptional couplings are defective

Everything here has γ̄ ≠ 0 except the paragraph beginning "The zero-mean stratum has the same phenomenon", which is where that stratum's own exceptional set is taken up, since no cofactor can find it. Sections 6 and 7 say the multiplicity is exactly ⌊N/2⌋ for all but finitely many J, the exceptions being the roots of the cofactor q(0), a polynomial of degree N(N−1) in J. One of those roots is already understood and is not what this section is about: J = 0 is a root of multiplicity 2⌊N²/4⌋, which is precisely the valuation Section 8 computes, and there the block is diagonal and the frozen root carries 2⌊N/2⌋ **semisimple** dimensions (the Corollary of Section 8). **Exceptional coupling** below always means a nonzero root. That leaves a question nobody had asked: are the nonzero roots real, and if a coupling really sits on one, what happens there?

Both halves have an answer, and the first one is already written down. The N = 3 cofactor of Section 6 is q(0) = 2¹²·γ̄²·J⁴·(3J² − δ₁²), so

  **J\* = ±δ₁/√3**

is a genuine real coupling whenever the profile is not uniform, and at it the divisor gains a rung: q(0) is the coefficient of ε^⌊N/2⌋, so q(0) = 0 forces the algebraic multiplicity to at least ⌊N/2⌋ + 1. The exceptions are not an artifact of writing "all but finitely many"; they sit on the real axis, one small coupling away from anywhere.

What happens there is the interesting half. The pencil argument of Section 3 is insensitive to J, so it still delivers at least ⌊N/2⌋ **geometric** dimensions and nothing more, while the algebraic count has gone up by one; whether the geometric count follows it up is exactly what the pencil cannot say. Nothing forces the two to meet, and Section 6's semisimplicity conclusion is exactly the one statement that needed q(0) ≠ 0. Computing the two multiplicities in exact arithmetic (no floating point: a rank test at a tuned degeneracy is the least trustworthy measurement there is) gives, at every exceptional coupling reached so far,

  **algebraic = ⌊N/2⌋ + 1,  geometric = ⌊N/2⌋:  one Jordan block of size exactly two sits at J\*.**

The cleanest way to see it needs no cofactor at all, only kernel dimensions of powers: dim ker M̃^k stabilises at the algebraic multiplicity and its increments give the block sizes. At N = 3 those dimensions are (1, 2, 2) over ℚ(√3), for every profile tried, against (1, 1, 1) at a non-exceptional coupling; at N = 4, on a profile whose exceptional coupling is fully rational (δ = (1, 3, −3, −1)/1000, J\* = 1/1000), they are (2, 3, 3), against (2, 2, 2) at both neighbouring couplings 1/999 and 1/1001. At N = 4 the algebraic count 3 also comes out of the characteristic polynomial itself, so two independent routes agree.

So the frozen divisor is semisimple almost everywhere and defective on a finite set of couplings. That set is the nonzero part of where its own tightness criterion vanishes; the zero part is J = 0 itself, where the criterion also fails but the multiplicity merely doubles and stays semisimple. The criterion cannot tell the two apart, which is why it takes a separate computation to see which kind of failure a given coupling is. Two readings follow. The placement of the next section needs one word of care: the frozen modes and the [Seed](../../compute/MirrorWorld/README.md) count are different objects, and they stay different here, since the Seed count is coupling-independent and vanishes at even N while the N = 4 example below is defective at even N and tuned to one coupling. What is true is weaker and still worth saying: the frozen root, generically a clean eigenvalue, becomes a defective point at J\*, so the divisor produces the kind of object a Seed census is built to count, at isolated couplings rather than generically. And the exceptional couplings are a discrete answer to a discrete condition: the profile is free to be anything on the locus, the coupling at which the extra mode joins is then forced by it, and at N = 3 it is forced in closed form.

**The zero-mean stratum has the same phenomenon, and a bigger block.** At γ̄ = 0 the cofactor is identically zero, so the criterion above says nothing at all and the exceptional couplings have to be found in the characteristic polynomial itself. They are there. On the gate's own zero-mean profile at N = 3, γ = (1/25, 0, −1/25) on the Heisenberg chain, the coefficient of λ³ is

  c₃(J) = (256/625)·J⁴·(75J − 2)(75J + 2),

so J\* = ±2/75 is a real nonzero exceptional coupling, and there the kernel dimensions of M̃, M̃², M̃³ are (3, 4, 5) against (3, 3, 3) at every control coupling tried: geometric N = 3, algebraic 5, **one Jordan block of size three**. That is one larger than anything the taxed stratum has shown. It does not settle Section 12's question, which asks whether the taxed stratum's block stays 2×2 as N grows; what it settles is the wider form of it, that the frozen root's Jordan structure is not universally 2×2 once the profile is free to leave that stratum, and it settles it at the smallest N there is. What the two strata share is the shape: an index that fixes a geometric count at every coupling, and a finite set of couplings where the algebraic count runs ahead of it.

What this section does not claim: that every root of q(0) is real. Only the N = 3 factor is available in closed form, and there the answer is yes. What is settled in general is weaker and still enough to matter: at least one real exceptional coupling exists whenever the profile is not uniform at N = 3, and wherever an exceptional coupling has been reached, the divisor is defective there. How many such couplings a given chain has is the counting question left open in Section 12.

## 10. What it is not (the placement)

- **Not a decoherence-free structure.** The frozen eigenvectors are not in ker(ad_H): the kernel of K intersected with the span of the anti-diagonal cells of O is at most one-dimensional (one at even N, none at odd N; the reading matters, since including the odd-N centre cell, which is diagonal and so not in O, would give one at every N), far below ⌊N/2⌋, and the frozen eigenvectors move with J. Only the eigenvalue stands still.
- **Not the uniform-γ commutant story.** At uniform γ (which sits on the locus as its fully degenerate point) the J-independent spectrum is the committed d_real ladder of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md), explained by weight-sector kernels ([absorption theorem](PROOF_ABSORPTION_THEOREM.md), [F50](PROOF_WEIGHT1_DEGENERACY.md)); those modes have J-independent eigenvectors. The frozen divisor is the site-resolved layer that survives when the profile is generic on the locus, and its mechanism is disjoint from the commutant.
- **Not a defective seed, except where it is.** Away from the exceptional couplings the frozen modes are semisimple (healthy left and right eigenvectors, overlap of order one; the zero-mean stratum behaves the same way at its own count N, semisimple off its own exceptional couplings and defective on them, Section 9), and the Seed count of MirrorWorld concerns defective points while this concerns an eigenvalue pinned across a family: siblings, not the same object. Section 9 is the qualification, and it is not a small one. On the finite set where the cofactor vanishes the frozen root is itself defective, so the two objects do meet, at isolated couplings rather than generically.
- **The F139 kinship is the design lesson.** There the wall factor S₁₀ divides a character polynomial exactly, with no symmetry realizing the reflection; here (λ + 4γ̄)^⌊N/2⌋ divides the corner characteristic polynomial exactly on a locus, with no symmetry of the spectrum. Both walls are divisors.

## 11. Verification

The committed gate [`simulations/r90_frozen_divisor_gate.py`](../../simulations/r90_frozen_divisor_gate.py) checks, and must print "R90 frozen divisor gate: ALL GREEN":

- the block builder against the framework Liouvillian (sub-block equality, exact);
- the mirror identity with defect 8γ̄·P_D at machine zero for N = 4, 5, 6, and that the antilinear variant fails (the mirror is linear);
- the census: corner multiplicity ⌊N/2⌋ at −4γ̄ for N = 3..6 across J ∈ {0.6, 1, 2.3}; the full block census over every (p,q) at N = 4, 5, 6 against both roots, which finds the four corner blocks of Section 4 and nothing else (Section 5); the XY variant (h without the ZZ diagonal) at N = 4, 5 (the corollary's antidiagonal corner is covered by the census above);
- the pencil kernel dimensions (= ⌊N/2⌋), the two eigenvector by-products, and the placement side claim of Section 10 that dim(ker K ∩ span of the anti-diagonal cells of O) is 1 at even N and 0 at odd N (N = 3..8);
- the partial-balance nulls;
- the N = 3 closed form, symbolically exact;
- at N = 6, exact Gaussian-rational arithmetic: the on-locus 36×36 corner determinant is exactly zero, and the transverse vanishing order in the defect is exactly 3 = ⌊N/2⌋;
- the cofactor theorem (Section 6): the closed form against the interpolated exact cofactor in Gaussian-rational arithmetic (N = 4, 5, Heisenberg; float cross-check XY N = 4), the symbolic N = 3 corner cofactor 2¹²γ̄²J⁴(3J² − δ₁²), and the nonvanishing of the leading coefficient det((K P_{O₊} K)|_{V₋}) for N = 3..10, Heisenberg and XY;
- the two boundary clocks (Section 7): the DCT-II / DST-I identification of the SE eigenbasis (machine zero, N = 3..10 both chains), the BBᵀ law with its {1, 1/M, 0} spectrum and pdet(G) = M^{−⌊(N−1)/2⌋}, the D_N closed form against the direct determinant, the exact rational assembly of the uniform cofactor at N = 4, 5 (Heisenberg, sympy discriminant), and uniform tightness at every sampled J (N = 3..6);
- the J-valuation discriminators (Section 6, last consequence): exact-rational cofactor ratios giving ord_J ≈ 18 at N = 6 (against 16) and ≈ 24 at N = 7 (against 20), a two-point estimate on exact arithmetic, since the orders themselves are pinned exactly by the next bullet; and the exact second-order rung rank(P_anti K P_{D₋}) = 1 at even N / 0 at odd N (N = 4..7);
- the counting data behind the open question of Section 12: exact Sturm counts of the distinct real nonzero roots on the polynomial in J², one pair at N = 3 for two profiles, two at N = 4, two at N = 5, and two against four for the two generic N = 6 profiles, with the disagreement itself asserted as a check;
- the exceptional couplings (Section 9), exact: the kernel dimensions of M̃, M̃², M̃³ at N = 3 over ℚ(√3) for three profiles, (1, 2, 2) at J* = δ₁/√3 against (1, 1, 1) at a non-exceptional control; the same three dimensions at N = 4 on a profile with a rational exceptional coupling (δ = (1, 3, −3, −1)/1000, J* = 1/1000), (2, 3, 3) against (2, 2, 2) at both neighbouring couplings 1/999 and 1/1001; and the N = 4 algebraic multiplicity confirmed a second time from the characteristic polynomial, so two independent routes agree;
- the index reading of Section 3.1: trace(τQ) = the fixed-cell count = dim(+) − dim(−) = N (N = 3..10); and the physics taken away, a RANDOM τQ-odd matrix keeping the whole effect (N = 3..8, six draws each): kernel N with no even defect, ⌊N/2⌋ with it, and 0 for an operator that is not odd, which is the control showing oddness is the entire hypothesis;
- the zero-mean stratum γ̄ = 0 (Sections 5, 6, 9), which is the previous item's first regime met on the locus rather than on noise, reached by handing the profile builder the pair total it should use, zero: the corner multiplicity, exact over GF(p) with i a square root of −1 and two primes, ⌊N/2⌋ at γ̄ = 0.09 against N at γ̄ = 0, semisimple at both (the nullity of M̃^k flat in k at the couplings sampled, which are not the exceptional ones), N = 3..7 on both chains at two nonzero couplings; the mechanism without an eigensolver, that M̃ is τQ-odd on the whole corner block at machine zero there, diagonal cells included, and that the attained count is the index dim V₊ − dim V₋; the two block-grid mirrors of Section 5, τQ on (p,p) and its X^N-bridged partner on (p, N−p), exactly odd at γ̄ = 0 and exactly not odd at γ̄ = 0.09 (residues mod p, so "odd" means identically zero), each with C(N,p) fixed cells, N = 3, 4, 5; the census the two of them open, the Heisenberg carriers at γ̄ = 0 being exactly q = p and q = N − p at multiplicity C(N,p), N = 4, 5, 6, and the XY figures Section 12 quotes, thirteen of twenty-five and thirty-six of thirty-six (one prime for the censuses, which is the direction that matters here: a rank mod p can only fall short of the rank over ℚ, so a nullity read this way can only overstate, and overstating is what the index bound already forbids); and the stratum's own exceptional coupling, at N = 3 on its own profile, both ways: the closed form of the λ³ coefficient symbolically, and the kernel dimensions (3, 4, 5) at J\* = 2/75 against (3, 3, 3) at three controls, mod two primes;
- the pointed grading of Section 8.1: that χ_x is admissible in both halves, τQ-invariant (constant on every parity basis vector) and one-hop Lipschitz (N = 3..9), that χ_x(c,c) = max(d_c, d_x) is the L¹ walk length, the case split of Lemma 10's proof exhaustively at N = 3..12, and that the census collapses to F_j = [χ_x(c,c) ≥ j]; then the order law ord_J det Y_(I_x, A_c) = max(d_c, d_x) for every pair and every anti-diagonal cell (N = 4..7), the identification of each support with the interval between that pair's own sites, and the strict nesting that makes the matrix triangular;
- the valuation law (Section 8), in integer and rational arithmetic: the two grading properties of Lemma 7 and the level census of Lemma 8 against their closed forms (N = 3..10, both chains); the transport bound of Lemma 9 tested exhaustively over **all** maximal row sets (N = 3..6 Heisenberg, up to 816 minors), both as the inequality ord_J det Y_I ≥ its transport bound and as the two sharpness statements (the bound never falls below ⌊N²/4⌋, and the measured minimum over row sets attains it); the total ord_J det((X P_{O₊} X)|_{V₋}) = 2⌊N²/4⌋ for N = 3..10 on both chains, by exact interpolation rather than by ratio; the per-pair orders ord_J S_{c,c'} = d_c + d_{c'} for every pair, N = 3..7 Heisenberg; Lemma 6 entry by entry against an independent cell-basis build of the left side, and the same grading and valuation for a third h with non-uniform R-symmetric bonds and an arbitrary R-symmetric diagonal (N = 4..6), which is the generality the section claims; and the positive-weight Gram form of the leading Schur matrix with its positive definiteness (N = 4..6, the reduction recorded in Section 12). The one floating-point item in this group is the departure census of the Corollary (N = 3..6).

## 12. Open

- The upper half of the valuation law, now down to one statement. Section 8.1 derives the staircase that makes the ĝ vectors triangular, so the only thing still measured rather than proved is that each pair reaches its OWN outer anti-diagonal cell, that is ĝ_{I_c}(c) ≠ 0, equivalently that Lemma 10's bound is attained there. The route is the uniqueness of the monotone walk from (c,c) to (c, R(c)): one index stands still, the other crosses, so there is a single cheapest route and nothing to cancel against. Everything downstream is already in place, and is recorded here because it is what that one statement buys: reading Cauchy-Binet at leading order gives

    S̃(0)·det T_{O₋,O₋}(0) = Σ_I w_I·ĝ_I(c)·ĝ_I(c'),  ĝ_I(c) := the coefficient of z^{d_c} in det Y_{I,A_c},

  where A_c is the column set of Section 8 (the O₋ columns plus the D₋ column of the pair c), z := −iJ, and the w_I are strictly positive. In the variable z every entry of Y is real, since X = z·(iK) − 2Δ with iK and Δ real, and z is a unit times J so no valuation moves. The leading matrix is therefore a **real positive semidefinite Gram matrix**, and its determinant cannot cancel: nonsingularity is exactly the linear independence of the ⌊N/2⌋ vectors ĝ(c). The gate confirms the identity and positive definiteness at N = 4..6. Section 8.1 supplies the independence from the staircase; the single nonvanishing above is what it rests on.
- The exceptional couplings of Section 9, as a counting question: how many of the nonzero roots of the cofactor are real (distinct roots, which is what the Sturm count gives), as a function of the profile and N. Exact Sturm counts on the polynomial in J² give the first data, and they already rule out the easy guess: the answer is not a function of N. At N = 3 both nonzero roots are real for every profile tried; at N = 5 only half of them are; and at N = 6 two generic profiles on the same locus give four real couplings and eight, respectively. What decides it is open, and so is whether the block is always 2×2 at larger N on this stratum. Across strata it is not: Section 9's zero-mean example carries a block of size three already at N = 3.
- The lever Section 3.1 names: the tax is the even defect on the diagonal cells, which exists because dephasing leaves populations alone. A generator under which the diagonal cells were ODD too would freeze N modes instead of ⌊N/2⌋, one per site rather than one per pair. Section 6 settles the arithmetic half, and the question that is left is narrower than it first looks. One completely positive generator does reach the count N, and it is the trivial one: γ ≡ 0 lies on the locus with γ̄ = 0, and there the N frozen modes are just the commutant of h in the single-excitation block, its N spectral projectors, which is not freezing so much as not moving. Every OTHER point of the stratum has a negative rate somewhere, since a zero mean and a nonzero profile force one. So the question is whether any generator that actually dissipates can hold N modes at one rate, and the two known ways to the number, the even defect removed and the Hamiltonian left alone, are exactly the two the answer must avoid.
- What the zero-mean stratum leaves open in the block census. The lower bound is settled, and by the same index: Section 5 exhibits the two mirrors, so every block with q = p or q = N − p carries at least C(N,p). What is measured and not derived is the two edges of that, namely that the bound is attained and that no other block carries at all (Heisenberg, N = 4, 5, 6). Whether the family is one object rather than a list is the question behind both, and C(N,p) being the popcount-p sector's dimension is the hint. On the XY chain the same γ̄ = 0 opens the census further and unevenly, thirteen carriers of twenty-five at N = 4 but all thirty-six at N = 5. The C(N,p) floor still holds there, since the two mirrors of Section 5 never used the chain, but it stops being an equality: the middle block (N/2, N/2) at even N carries 8 against C(4,2) = 6 at N = 4 and 24 against C(6,3) = 20 at N = 6; whether those belong to the frozen band of [XY_FROZEN_BAND](../../experiments/XY_FROZEN_BAND.md) is untouched.
- The uniform-γ endpoint: how the frozen divisor's ⌊N/2⌋ modes embed into the enhanced d_real counts when all rate classes collapse. Section 7 supplies the frozen side in closed form; the d_real side of the ledger is still the committed open problem of [DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md).

Two items that stood here have since closed, and are recorded so the list reads as what is left rather than as what was ever asked. The divisor is adopted in MirrorWorld as [`Divisor`](../../compute/MirrorWorld/Divisor.cs), beside Seed, with run mode `divisor N`; its counts hold past the spectral wall, to N = 20. And it is typed in the Object Manager as [`FrozenDivisorClaim`](../../compute/RCPsiSquared.Core/BlockSpectrum/FrozenDivisorClaim.cs), Tier 1 derived under F91 and the joint-popcount decomposition, with the live witness [`FrozenDivisorWitness`](../../compute/RCPsiSquared.Diagnostics/Foundation/FrozenDivisorWitness.cs) at `inspect --root divisor`. What that witness types is the divisor bound and the fold-parity placement of the corners; its multiplicities are exact GF(p) ranks rather than eigenvalue counts, because the departing modes sit at spacing J^{2d} from the root and a floating-point count would go wrong exactly where this document is sharpest. It meets the tightness theorem at every coupling it samples, and it cannot see the defectiveness of Section 9, since a kernel dimension does not move there. Reading the same object on the other chain is what turned Section 5's census into a scoped statement; `inspect --root divisor --N 5 --chain xy` is the invocation that shows it.
