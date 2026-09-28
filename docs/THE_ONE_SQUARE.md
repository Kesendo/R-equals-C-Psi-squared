# The One Square: the diagonals of the repository as the two lines of a square

**Status:** A map of existing results, Tier 1 where it restates them; its own additions (the name anti-diagonal for the R₉₀ locus, the corner's paying moves failing on the diagonal) are exact identities checked in Python. The gate checks the page's algebra with the rates and the Gram matrix written out; it does not read them off the cited objects. The moves, their fixed lines, the disagreement sets, the block grid's two charges and the corner's half-turn law are checked exactly by [`simulations/diagonal_square_gate.py`](../simulations/diagonal_square_gate.py) (must print "diagonal square gate: ALL GREEN", under a second); each check is a finite case of an identity that one line of algebra proves for every N. What the page says about spectra, dynamics, hardware and qudits is cited from the documents named with it.
**Date:** 2026-09-28
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Told in plain words:** [On the One Square](../reflections/ON_THE_ONE_SQUARE.md)

## What the repo already held

The question was Tom's: the word "diagonal" appears everywhere in this repository; is there one system behind it? The sweep went by the primitive, a square of cells with a swap of its two labels and an opposite for each label, not by the word.

- **docs/ANALYTICAL_FORMULAS.md.** F118 (the mirror group), F89c (the cell cost 2γ₀·n_diff and the column flip that supplies the spectral pair-sum), F89d (one leg of the bit-flip Klein group on the coherence-block lattice), F91 (R₉₀ "reflects each pair-sum about 2γ_avg"), F140 (at γ̄ = 0, "two mirrors serving the two diagonals of the block grid"), F143, F158.
- **docs/proofs/.** [The Π factorization](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) §2–4 (the eight mirrors, the two Klein subgroups, the Hamiltonian split, the class that swaps lit and dark letters); [the fold lattice](proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §6–7 (the charges d = p − q̃ and m = p + q̃ − N, the fold parity with kernel {1, t, Klein, t·Klein}); [the absorption theorem](proofs/PROOF_ABSORPTION_THEOREM.md) §4.7; [the frozen divisor](proofs/PROOF_R90_FROZEN_DIVISOR.md) §1, §3.1 (the surplus and the tax), §4 (the four corner blocks), §5 ("two mirrors serve, one for each diagonal of the block grid"), §6–8 and its open list; [F91](proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md); [F103](proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) (the anti-diagonal commutant element); [F158](proofs/PROOF_PALINDROME_TWO_END_COUNT.md); [the uniform law](proofs/PROOF_UNIFORM_LAW.md) (the anti-diagonal cells of block (1, N − 1)); [the qudit palindrome](proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) §6 (the wreath family).
- **experiments/.** [The price-pair flight](../experiments/PRICE_PAIR_HARDWARE_PREDICTION.md); [XOR Space](../experiments/XOR_SPACE.md) (Result 4: the XOR modes purely anti-diagonal); [the gamma fold](../experiments/GAMMA_FOLD_PAIR_OF_MIRRORS.md) (the one-sided X^N bridge of the frozen divisor's §4); [the XY frozen band](../experiments/XY_FROZEN_BAND.md).
- **docs/.** [The Three Diagonals](THE_THREE_DIAGONALS.md), [the qubit necessity](QUBIT_NECESSITY.md), [the glossary](GLOSSARY.md) (the XOR modes, the light-and-lens entry, which is the one-site square; no entry for the word before this page), [CAUGHT_ERRORS](CAUGHT_ERRORS.md) (the Klein-group trap on the rate side).
- **reflections/.** [On the One Diagonal](../reflections/ON_THE_ONE_DIAGONAL.md); [On the Ninety-Degree Gamma](../reflections/ON_THE_NINETY_DEGREE_GAMMA.md) (R₉₀ as the rate-side shadow of the operator quarter-turn); [On Leaving the Circle](../reflections/ON_LEAVING_THE_CIRCLE.md) (the turn and the fade).
- **The typed layer.** `MirrorGroupD4Claim`, `F89CrossFoldSimilarityClaim`, `FrozenDivisorClaim` and `FrozenDivisorWitness` (the "anti-diagonal coherences (a, R(a))"), `PalindromeTwoEndCountClaim` and its witness, `SeedRungGramClaim`, `ThreeDephasingDiagonalsOrbitClaim`, `QuditProductMirrorCap`, `QubitNecessityPi2Inheritance`, `DiagonalWitness` (`inspect --root diagonal`).
- **MirrorWorld.** `Mirror` and `MirrorGroup` (whose header says the group "is Mirror's fold lattice read at the operator level"), `Divisor`, `ParameterKlein`, `GammaFold`, `Restless` (the anti-trace), `Lattice`, `Cat`, `SpookyAction`.
- **OpenArcs.** `one_diagonal_mirror_group`, `sideways_spin_ladder`.
- **fw.Confirmations.** `price_pair_locality_marrakesh_july2026`.
- **hypotheses/.** Nothing beyond generic uses of the word.

Checked for adjacency, pair by pair: the Π factorization against the fold lattice (joined already, in MirrorWorld's `MirrorGroup` and in F89d); the frozen divisor's §5 against F140 (the same sentence in two places); F158 against the XOR-mode entry of the glossary (joined in F158); the qudit palindrome against the qubit necessity (joined in `QuditProductMirrorCap`). Nothing found that joins the squares below on one page, and no document names the R₉₀ locus an anti-diagonal. That name and this page are what is added.

## 1. The square and its moves

A square here is a grid of cells whose rows and columns carry the same labels, with a label involution c ≠ id, the opposite. The **swap** t exchanges row and column; the **two-sided flip** c × c sends both labels to their opposites. The diagonal is Fix(t); the anti-diagonal is Fix(t ∘ (c × c)); with the identity and c × c these are four moves. The **one-sided flip** 1 × c sends only the column label to its opposite; with the swap it generates the eight symmetries of the square (D₄). On paper every square has all eight; the last column of the table below says which of them are known to act on the physics, meaning that they carry the generator, or the one object named, to itself up to sign and shift. A glyph note: on this page R written alone is the one-sided flip ρ ↦ ρ·F of §2–3; the chain's site reversal is written R(a), γ_R(l) and R × R. D is always the transpose here, not the frozen divisor's set of diagonal cells.

| The square | The opposite of a label | Diagonal | Anti-diagonal | What the moves are known to do |
|---|---|---|---|---|
| the cells \|a⟩⟨b\| of a density matrix | every bit flipped, ā | a = b, the populations | b = ā, every site disagreeing | eight mirrors of operator space; the four that swap the lines pay −2Σγ on the dephasing (§2) |
| the blocks (p, q), p and q the excitation counts | p ↦ N − p | p = q | p + q = N | eight index maps of the fold lattice; its folds pair the spectra where its §7 says |
| the mode indices of F143's seed-rung Gram matrix | a ↦ M − a | a = c | a + c = M | eight invariances of that one matrix, Ĝ = 2 + [a = c] + [a + c = M] |
| the single-excitation cells (a, b) of the frozen divisor's corner | the chain reflection a ↦ R(a) | the N populations | the N cells τQ fixes | four; the two that pay, τQ and the half-turn R × R, keep the lines and fail together on the populations (§4) |
| a mirrored pair of rates (γ_l, γ_R(l)), at fixed mean γ̄ | reflection about γ̄ | palindromic profiles (F71) | the R₉₀ locus γ_l + γ_R(l) = 2γ̄ | four parameter maps: F71 and R₉₀, a Klein group |

In the corner, τQ is the frozen divisor's mirror: swap the cell's two sites and reflect both across the chain. On the rate square the swap is F71 and the swap followed by the two-sided flip is R₉₀, γ_l ↦ 2γ̄ − γ_R(l); F91's "reflects each pair-sum about 2γ_avg" is the reflection across that anti-diagonal. ⟨F71, R₉₀⟩ is a Klein group with no element of order four (F91; CAUGHT_ERRORS records the trap of reading a quarter-turn into it). No one-sided flip of this square is known (MirrorWorld's `GammaFold` turns one site's rate about 0, not about γ̄), and the obvious candidate is not one: flipping one pair's rate about γ̄ moves the mean, so it leaves the fixed-mean square.

## 2. The cell square: the half that pays

The swap is the transpose D(ρ) = ρᵀ, the one-sided flip is R(ρ) = ρ·F with F = X^⊗N, the two-sided flip is 𝓕(ρ) = F·ρ·F. Their eight products are the mirror group of the Π factorization, each a permutation of the cells (gate S1, S1b).

| Move | On a cell \|a⟩⟨b\| | Disagreement set a ⊕ b | Fixed cells |
|---|---|---|---|
| 1, 𝓕, D, 𝓕D | (a, b), (ā, b̄), (b, a), (b̄, ā) | kept | all; none; the diagonal; the anti-diagonal |
| R, 𝓕R, Π_Z = R·D, Π_Y | (a, b̄), (ā, b), (b, ā), (b̄, a) | complemented | none; each maps the diagonal onto the anti-diagonal |

A cell decays under the dephasing at −2 times the sum of γ over its disagreeing sites, so the first half keeps the dephasing and the second half turns it into −L_D − 2Σγ (the Π factorization §3; F89c is the same column-flip relation stated as the spectral pair-sum). The paying class is the one the Π factorization §4(f) describes as exchanging the dark letters {I, Z} with the lit letters {X, Y}. The Hamiltonian follows a different split: for the XXZ-type chains here (every term with an even number of Y letters and invariant under F) D turns it over and R keeps it.

The complement law was flown. The price-pair campaign on ibm_marrakesh measured the decay of every coherence pattern on a three-qubit line and tested Γ(D) + Γ(D̄) = Σ_j Γ_j, a pattern's rate plus its complement's equal to the sum of the single-coherence rates (in Lindblad coefficients, 2Σγ). Its readings include local T1 terms. Run 1's headline was withheld for detuning drift, and run 2 showed an unexplained covariance. In the clean run 3 the local-dephasing premise was supported, every covariance within 2σ (fw.Confirmations, `price_pair_locality_marrakesh_july2026`).

## 3. The block grid: two charges

Counting excitations, (p, q) = (popcount a, popcount b), carries each of the eight cell moves to a map of (p, q) alone, distinct and composing as before (gate S2); F89d is the leg R of the Klein subgroup {1, 𝓕, R, 𝓕R}, block-resolved and realized there as an antiunitary similarity, and MirrorWorld's `MirrorGroup` says the same in words. The two lines are the zero sets of the fold lattice's two charges, d = p − q (the excitation difference) and m = p + q − N (the weight of its §6 raising ladder). The moves act on (d, m) by signed permutations:

| Move | (d, m) ↦ |
|---|---|
| 1, D, 𝓕, 𝓕D | (d, m), (−d, m), (−d, −m), (d, −m) |
| R, 𝓕R, Π_Z, Π_Y | (m, d), (−m, −d), (m, −d), (−m, d) |

The free half re-signs the charges and the paying half exchanges them. The fold lattice's §7 finds the same split as the kernel {1, t, Klein, t·Klein} of its fold parity. The frozen divisor puts this grid to work in §4. The diagonal corners (1, 1) and (N − 1, N − 1) carry the frozen root −4γ̄. The anti-diagonal corners (1, N − 1) and (N − 1, 1) carry its fold image 4γ̄ − 2σ. A one-sided X^N bridge is the move between them.

## 4. The corner: the paying moves keep the lines

On the corner the moves known to act are the four. For a single-excitation hopping h that is real, symmetric and invariant under the site reversal (the frozen divisor's §1), the Hamiltonian part, −iK with K(ρ) = hρ − ρh (this page's K; the proof's K carries the −i and the J), commutes with the half-turn R × R and changes sign under the swap t and under τQ, since K(ρᵀ) = −K(ρ)ᵀ (Lemma 1 of the frozen divisor). On the R₉₀ locus, with δ_a := γ_a − γ̄, the recentred rate of a cell, r(a, b) = rate + 4γ̄, is −2(δ_a + δ_b) off the diagonal and 4γ̄ on it, and it is:
- odd under τQ and under the half-turn R × R on every cell off the diagonal (for τQ this is Lemma 2 of the frozen divisor);
- even under both on the populations, where it stays 4γ̄. Since t fixes each population, τQ and R × R act there alike; the evenness is the even defect of the frozen divisor's §3.1, the "tax" (gate S3).

So on this square τQ plays the palindromizer's part, turning both the Hamiltonian and the recentred rates over, and the half-turn R × R plays the part of the one-sided flip on the cell square, keeping the Hamiltonian and turning the rates over. Both keep the two lines in place, and both fail exactly on the diagonal: the populations stand still at rate 0, while the pairing is about the frozen rate −4γ̄. The proof's count runs on τQ: its fixed cells off the diagonal give a surplus of 2⌊N/2⌋, the populations' evenness takes ⌊N/2⌋ of it back, and ⌊N/2⌋ remain (§3.1). At γ̄ = 0 the populations' recentred value is 0, both moves negate every cell, and the lower bound rises to N; the proof's words are that the populations then "stop charging and start paying", and N is measured to be attained. A zero mean with a nonzero profile needs a negative rate somewhere, so this is arithmetic, not a dissipating channel (§3.1). The one-sided reflection (a, b) ↦ (a, R(b)) also keeps the Hamiltonian part, but it sends −2(δ_a + δ_b) to −2(δ_a − δ_b), neither keeping nor turning over the rates (gate S3).

The thesis of §2 is therefore a statement about the cell square, not a law of every square. On the density matrix the half that swaps the diagonal with the anti-diagonal is the half that pays. In the corner the paying moves keep the lines, and the diagonal is where paying fails. Whether a dissipative generator can make the diagonal pay too is the lever the frozen divisor's open list already names.

## 5. The anti-diagonal as a place

- **The anti-trace.** MirrorWorld's `Restless` sums ρ over b = ā: the fastest-dying content in the ordinary world, the conserved one in the anti-watched world. Of `Lattice`'s four worlds of watching, {1, 𝓕, R, 𝓕R}, two keep the immortal set on the diagonal and two put it on the anti-diagonal: two moves from each half of §2.
- **The two ends of the palindrome.** [F158](proofs/PROOF_PALINDROME_TWO_END_COUNT.md): the palindrome holds exactly when dim ker L = dim ker(L + 2σ). With every rate positive, an eigenmode at −2Σγ must disagree at every site and one in the kernel at none: the Hamiltonian part is anti-Hermitian in the Hilbert–Schmidt product and the dephasing is diagonal on cells, so Re λ = ⟨v, L_D v⟩/⟨v, v⟩, which reaches −2Σγ only on the anti-diagonal and 0 only on the diagonal. So the modes at −2Σγ, the XOR modes, live on the anti-diagonal and the kernel on the diagonal; [XOR Space](../experiments/XOR_SPACE.md) (its Result 4) already writes the XOR modes as purely anti-diagonal in the cell basis.
- **The cat.** A single anti-diagonal cell, |0…0⟩⟨1…1|, is the GHZ coherence of MirrorWorld's `Cat`; its two-site version |00⟩⟨11| is `SpookyAction`'s. Under local dephasing it decays at 2Σγ, every site's full share. On hardware the price-pair flight found it outliving the product of its single-site parts. The campaign traced the bulk of that to always-on ZZ, which the GHZ coherence does not feel (XXX commutes with every Z_iZ_j) while the product test-bed's single sites do, and confirmed the ZZ directly in its fourth run; a genuine anti-correlated component of at most about 15 % is left as an upper bound, not a claim.
- **F87.** On F87's diagonal cell (the Z letter) a pair is soft exactly when its Hamiltonian commutes with an operator supported on b = ā ([F103](proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md), the anti-diagonal commutant element).

## 6. In more dimensions

At one site the cell square is 2 × 2: I and Z are diagonal matrices, X and Y anti-diagonal ones. On N sites the cells with a fixed pattern x of disagreeing sites, a ⊕ b = x, are the 2^N parallels of the diagonal in the hypercube's sense, from x = 0 (the diagonal) to x = 1…1 (the anti-diagonal). Grouped by k = popcount(x) they are the levels of Q = N − 2k ([the absorption theorem](proofs/PROOF_ABSORPTION_THEOREM.md) §4.7), the dephasing diagonal, which is the height across the parallels. For the XY and XXZ chains here, whose bonds flip two sites at a time, k changes by 0 or ±2 and its parity is conserved. On the XY block grid the arc `sideways_spin_ladder` places spin multiplets on p + q = N ± 1 and η multiplets on d = ±1.

For a site of dimension d the off-diagonal cells fall into d − 1 wrapped parallels, b = a + s mod d. The palindromizer takes two chiralities along b = a ± 1, and with the transpose generates the wreath product Z_d ≀ Z₂ of order 2d², whose d = 2 column is the group of eight ([the qudit palindrome](proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) §6). `QuditProductMirrorCap`: "at d = 2 the two off-diagonals coincide, the chiralities merge, and the mirror is full: that degeneracy IS the qubit magic". In this page's terms, on a qubit the whole off-diagonal is one line, the anti-diagonal. The count is [the qubit necessity](QUBIT_NECESSITY.md): a complete local exchange of the d diagonal cells with the d² − d others needs d² − 2d = 0. Across sites a palindrome can do more, which is what the qudit palindrome's partial and non-product mirrors are.

## 7. What is not this

- **The dephasing diagonal Q** and its X and Y siblings ([The Three Diagonals](THE_THREE_DIAGONALS.md)) are diagonal operators: the height across the lines, not lines. The whole construction above is written for Z-dephasing; under X- or Y-dephasing the square stands in that basis.
- **The diagonal blocks** of F91 and F63 are an operator's blocks in the eigenspaces of an involution.
- **The diagonal cell** of the F87 table is where a dephasing letter meets its own Klein cell.
- **One coincidence.** D and 𝓕D are the square's two diagonal mirrors and, on Pauli strings, the diagonal sign matrices (−1)^(n_Y) and (−1)^(n_Z) (the Π factorization §2; gate S5).
- **Angles.** On the block grid Π is a quarter-turn. The pairing it makes on the spectrum, λ ↦ −λ − 2Σγ, is an involution. The first statement is about where blocks go, the second about eigenvalues.

## 8. Open

- A dissipative generator under which the corner's diagonal pays too, so that N modes freeze at γ̄ ≠ 0: the lever of the frozen divisor's open list, restated here as the paying moves failing on one line.
- This page's own additions, the name anti-diagonal for the R₉₀ locus and the reading of the tax as the paying moves failing on the diagonal, are checked in Python only; under the repository's witness-first rule they belong in `FrozenDivisorWitness` if they become load-bearing.
