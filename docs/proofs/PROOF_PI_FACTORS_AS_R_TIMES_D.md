# PROOF: the canonical palindromizer factors, Π_Z = R·D, and the mirror inventory closes into one dihedral group D₄

**Status:** Tier 1 derived (signed-permutation identities, exact; group closure, factorization, and palindrome split verified bit-exact, dev 0.00e+00, at N = 1..3 with an N = 5 spot check; the palindrome factorization rows at machine precision ≤ 5.6·10⁻¹⁷ on XXZ with site-dependent γ).
**Date:** 2026-06-10
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [MIRROR_SYMMETRY_PROOF.md](MIRROR_SYMMETRY_PROOF.md): the April palindrome proof and the canonical per-site Π (I → X, X → I, Y → iZ, Z → iY) that this note factors.
- [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2: the windowed-converse spine's two involutions 𝓕 = F ⊗ F and R = I ⊗ F and their sign table.
- The F114 entry ("F114. Commutator Superoperator D-Conjugation Parity") in [ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md): the sign law D·L_σ·D = ε(σ)·L_σ, ε(σ) = (−1)^{n_Y(σ)+1}.
- [reflections/D_PI_Z_EQUALS_PI_Y.md](../../reflections/D_PI_Z_EQUALS_PI_Y.md) and [PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md](PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md): Welle 12's identity D·Π_Z·D = Π_Y, re-read here as the dihedral inversion relation.
- [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md): the truly criterion (n_Y even AND n_Z even) that turns out to be a D₄ character cell.
- [PROOF_CEILING_GOLDEN_ROUTER.md](PROOF_CEILING_GOLDEN_ROUTER.md) (F116) and [PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md): the mirrors that are deliberately NOT in the group (§5).

## Abstract

Two months apart, this repository built two mirror toolkits that never spoke to each other. April 2026 built the canonical palindromizer Π_Z, a per-site signed permutation of Pauli letters, and proved Π_Z·L·Π_Z⁻¹ = −L − 2Σγ ([MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md)). June 2026 built the windowed-converse spine from two plain involutions, the charge conjugation 𝓕 = F ⊗ F and the ket reflection R = I ⊗ F with F = X^⊗N ([PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)). In between, Welle 12 found a third object, the diagonal sign D = diag((−1)^{n_Y}) on the Pauli basis, and noticed it swaps dephase letters: D·Π_Z·D = Π_Y ([D_PI_Z_EQUALS_PI_Y](../../reflections/D_PI_Z_EQUALS_PI_Y.md)).

This note states what they are to each other. Define, on the coherence space of an N-qubit chain:

- **R**, the spine's ket reflection: R(ρ) = ρ·F with F = X^⊗N. In the row-stacking vec convention this is the operator I ⊗ F; it flips the ket (column) index of every coherence |i⟩⟨j|, j ↦ j̄ (all bits flipped).
- **D**, the transpose superoperator: D(ρ) = ρᵀ. On the Pauli basis D is exactly the diagonal sign diag((−1)^{n_Y}), F114's mirror, because transposition flips the sign of each Y letter and nothing else.

**Theorem (factorization).** The canonical palindromizer is their product:

  **Π_Z = R ∘ D,  i.e.  Π_Z(ρ) = ρᵀ · X^⊗N**  (transpose first, then right-multiply by F),

per site σ ↦ σᵀ·X, which reproduces the April rule I → X, X → I, Y → iZ, Z → iY with no extra phase.

**Theorem (the group).** ⟨R, D⟩ is dihedral of order 8, D₄: the rotations are {I, Π_Z, 𝓕 = Π_Z², Π_Y = Π_Z³} and the reflections are {D, 𝓕D, R, 𝓕R}. Every named mirror of the palindrome story sits in this one group: the April Π_Z and Π_Y, the spine's 𝓕 and R (its Klein four-group {I, 𝓕, R, 𝓕R} is a Klein subgroup of D₄), F114's D, and one diagonal mirror the repository had never named, 𝓕D = diag((−1)^{n_Z}). The palindrome identity itself factors along the generators: D flips the Hamiltonian commutator and fixes the dissipator; R fixes the Hamiltonian commutator and reflects the dissipator, carrying the entire −2Σγ shift.

All identities are equalities of signed permutation matrices, verified at dev 0.00e+00 by [`simulations/mirror_inventory_d4.py`](../../simulations/mirror_inventory_d4.py).

## §1 The factorization, with the exact convention

Start at one site, with one concrete matrix. Take the population ρ = |0⟩⟨0|. Transpose does nothing to it (it is diagonal); right-multiplying by X slides the ket: |0⟩⟨0|·X = |0⟩⟨1|. A population became a coherence. That is the palindromizer's signature move (it exchanges what decays slowly with what decays fast), performed by two operations so plain that each one alone looks harmless.

Now the general statement. Work in the row-stacking (C-order) vec convention, |i⟩⟨j| ↦ e_i ⊗ e_j, where kron(A, B) acts as ρ ↦ A·ρ·Bᵀ; this matches `framework.lindblad` and [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §1](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md). With F = X^⊗N (so F² = I and Fᵀ = F), define

  D(ρ) = ρᵀ  (the transpose superoperator; the SWAP permutation on coherence space),
  R(ρ) = ρ·F  (the operator I ⊗ F; it touches only the ket index).

**Π_Z = R ∘ D: apply D first, then R.**

  Π_Z(ρ) = R(D(ρ)) = ρᵀ·F.

The order matters, because R and D do not commute; the opposite order is the inverse (§4a). The check that this is really the April Π is a four-line per-site computation, since both transpose and right-multiplication by X^⊗N factor over sites, leaving σ ↦ σᵀ·X per site:

  Iᵀ·X = X,  Xᵀ·X = I,  Yᵀ·X = −Y·X = iZ,  Zᵀ·X = Z·X = iY.

That is exactly the canonical rule of [MIRROR_SYMMETRY_PROOF](MIRROR_SYMMETRY_PROOF.md), including the essential factors i that the April search had to discover by trial: here they fall out of Y's antisymmetry (Yᵀ = −Y) meeting the Pauli product YX = −iZ. The hard-won phase is the meeting point of the two generators.

On the Pauli basis, D is diagonal: σᵀ = (−1)^{n_Y(σ)}·σ, so D = diag((−1)^{n_Y}), which is precisely the D of [F114](../ANALYTICAL_FORMULAS.md) and Welle 12. The verifier confirms the action Π_Z(ρ) = ρᵀ·F on random ρ at N = 1, 2, 3 (and rejects the wrong-sided alternative F·ρᵀ, which is Π_Y, by a margin of order 1), and the superoperator identity Π_Z = R·D as a matrix equation at dev 0.00e+00.

## §2 The group: eight mirrors, one table

R and D are involutions, and their product has order 4 (Π_Z⁴ = I, with Π_Z² = 𝓕 the charge conjugation). Two involutions whose product has order 4 generate the dihedral group of the square, D₄, with eight elements. Closing ⟨R, D⟩ numerically gives exactly eight signed permutation matrices, and every one of them is a mirror this repository already uses, plus one it never named:

| Element | Acts on ρ as | Order | Pauli-basis form | Where it lives |
|---|---|---|---|---|
| I | ρ | 1 | identity | trivial |
| Π_Z | ρᵀ·F | 4 | signed permutation (I↔X, Y↔iZ, Z↔iY per site) | the April palindromizer (P1) |
| 𝓕 = Π_Z² | F·ρ·F | 2 | diag((−1)^{n_Y+n_Z}) | F1² charge conjugation; the spine's 𝓕 |
| Π_Y = Π_Z³ | F·ρᵀ | 4 | signed permutation | the Y-dephasing palindromizer |
| D | ρᵀ | 2 | diag((−1)^{n_Y}) | F114; Welle 12 |
| 𝓕D | F·ρᵀ·F | 2 | diag((−1)^{n_Z}) | unnamed until today (§4e) |
| R | ρ·F | 2 | signed permutation (I↔X, Y↦−iZ, Z↦iY per site) | the spine's ket reflection |
| 𝓕R | F·ρ | 2 | signed permutation (I↔X, Y↦iZ, Z↦−iY per site) | the spine's bra reflection |

The defining relations are the dihedral ones, all verified as matrix equations:

  R² = D² = I,  Π_Z = R·D,  Π_Z⁴ = I,  **D·Π_Z·D = Π_Z⁻¹.**

The geometry is worth pausing on. D₄ is the symmetry group of a square: one 90° rotation (Π_Z), its powers, and four mirrors in two conjugacy classes, the two diagonal mirrors and the two edge mirrors. Here the two classes have a meaning you can read off the table: **the diagonal mirrors of the square, {D, 𝓕D}, are literally the diagonal matrices** of the Pauli basis (pure ±1 signs, no string moves), while the edge mirrors {R, 𝓕R} are the one-sided F multiplications that move strings (I↔X) and touch one tensor leg each.

D₄ has three subgroups of order 4, and all three are already characters in the story:

- the rotation subgroup {I, Π_Z, 𝓕, Π_Y}: the palindromizer family;
- the spine's Klein four-group **{I, 𝓕, R, 𝓕R}**: the two-sided and one-sided F multiplications, exactly the involution set [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) runs its sign table on;
- the diagonal Klein four-group {I, 𝓕, D, 𝓕D}: the sign gradings, whose joint eigenbasis is the Pauli strings themselves (§4d).

All three intersect in {I, 𝓕}: the charge conjugation 𝓕 = Π_Z² = F1² is the center of D₄. Every mirror in the inventory passes through it.

## §3 The palindrome, factored along the generators

The April proof had three steps: Π flips the XY-weight of the dissipator (producing the −2Σγ shift), Π anti-commutes with the Hamiltonian commutator, combine. The factorization splits those steps cleanly between the two generators. Write L = L_H + L_diss with L_H = −i[H, ·] and L_diss the Z-dephasing dissipator (diagonal on coherences, entry −2·Σ_{l: bit l of i⊕j set} γ_l). Then, verified on the XXZ chain (Δ = 0.7) with site-dependent γ = (0.05, 0.11, 0.07) at N = 3:

  D·L_H·D = −L_H    (D carries the Hamiltonian flip)
  D·L_diss·D = +L_diss    (D fixes the dissipator)
  R·L_H·R = +L_H    (R fixes the Hamiltonian commutator)
  R·L_diss·R = −L_diss − 2Σγ·I    (R reflects the dissipator and carries the entire shift)
  Π_Z·L·Π_Z⁻¹ = −L − 2Σγ·I    (the product is the palindrome)

Each row is one generator doing one job. D flips L_H because transposition reverses commutators ([σ, ρᵀ]ᵀ = −[σᵀ, ρ]) and every XXZ term has n_Y even, so σᵀ = +σ; this is F114's sign law ε(H) = −1 in action. D fixes L_diss because the dephasing entry depends only on i ⊕ j, which is symmetric under i ↔ j. R fixes L_H because conjugating the ket leg by F leaves every XXZ term invariant (F·σ·F = +σ when n_Y + n_Z is even). And R reflects the dissipator because flipping the ket index j ↦ j̄ complements the set of differing bits: the lit sites become the dead sites, the rate −2Σ_{lit} γ_l becomes −2Σγ + 2Σ_{lit} γ_l, and the constant −2Σγ appears. In the recentered language of the spine this last row is exactly the sign-table entry R·Q·R = −Q.

So the 2026-04 palindrome proof and the 2026-06 windowed-converse spine are not two toolkits; they are the same two generators, deployed at different angles. April multiplied them (Π = R·D, one operator, conjugation, the full palindrome). June kept them separate (𝓕 = (RD)², R, two involutions, trace parities, the power-sum spine). The group ⟨R, D⟩ is where both live.

## §4 Corollaries

**(a) Π_Y = Π_Z⁻¹, and Welle 12 is the dihedral inversion relation.** In any dihedral group, conjugating the rotation by a reflection inverts it: s·r·s = r⁻¹. With r = Π_Z and s = D this reads D·Π_Z·D = Π_Z⁻¹, and the verifier confirms Π_Z⁻¹ = Π_Y at dev 0.00e+00. Welle 12's identity D·Π_Z·D = Π_Y ([D_PI_Z_EQUALS_PI_Y](../../reflections/D_PI_Z_EQUALS_PI_Y.md), universal-N proof in [PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N](PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md)) was the dihedral inversion relation in disguise: the Z↔Y dephase-letter swap and "running the palindromizer backwards" are the same operation. Equivalently, the two palindromizers are the two orderings of the same factors:

  Π_Z = R ∘ D (transpose, then reflect: ρᵀ·F),  Π_Y = D ∘ R (reflect, then transpose: F·ρᵀ).

**(b) 𝓕 = Π_Z² = F1².** The charge conjugation X^⊗N, registered as F1² in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md), is the square of the palindromizer and the center of D₄ (§2). That it commutes with everything in the group is no longer an observation but a group-theoretic necessity.

**(c) F114 is D's row of the table.** The sign law D·L_σ·D = ε(σ)·L_σ with ε(σ) = (−1)^{n_Y(σ)+1} becomes transparent once D is recognized as the transpose superoperator: one sign from transposition reversing the commutator, one sign (−1)^{n_Y} from σᵀ = (−1)^{n_Y}·σ. The general-N statement and its derivation are owned by the F114 entry ("F114. Commutator Superoperator D-Conjugation Parity") in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md); the verifier here adds an N = 5 spot check (σ = XIYIY, n_Y = 2, dev 0.00e+00), one size beyond F114's originally verified N ≤ 4.

**(d) The truly criterion is a D₄ character cell.** The diagonal Klein four-group {I, 𝓕, D, 𝓕D} is simultaneously diagonalized by the Pauli strings, and its joint eigenvalues are the pair of parities ((−1)^{n_Y}, (−1)^{n_Z}). The four eigenvalue patterns are the four Klein cells of a Pauli string, and the F85/F87 truly criterion (n_Y even AND n_Z even, [PROOF_F85_KBODY_GENERALIZATION](PROOF_F85_KBODY_GENERALIZATION.md)) is exactly the **joint-fixed cell**: σ is truly iff both diagonal mirrors fix it. At the commutator level the same condition reads through the factorization: Π_Z·L_σ·Π_Z⁻¹ equals −L_σ iff D flips L_σ (n_Y even, F114's ε = −1) and R fixes it (F·σ·F = +σ), which together force n_Z even. The verifier checks this equivalence against the framework classifier `_pauli_tuple_is_truly` for all 63 non-identity strings at N = 3: 63/63.

**(e) The fourth mirror gets a name.** 𝓕D = diag((−1)^{n_Z}), the product of the charge conjugation and the transpose, acting as ρ ↦ F·ρᵀ·F. Its character (−1)^{n_Z} has been doing work in this repository for weeks (it is the second leg of the truly criterion and the grading behind the n_Z-parity arguments of the F85 chain), but the operator itself was never written down or named. It completes the diagonal mirror pair: D grades by Y content, 𝓕D grades by Z content, their product 𝓕 grades by both.

**(f) Why the ceiling cases escape the canonical Π.** The factorization turns "which terms does Π flip" into a two-character question. Per Pauli string σ, conjugation by Π_Z = R·D gives

  Π_Z·L_σ·Π_Z⁻¹ = −L_σ (truly cell: n_Y even, n_Z even),
  Π_Z·L_σ·Π_Z⁻¹ = +L_σ (n_Y odd, n_Z odd),
  Π_Z·L_σ·Π_Z⁻¹ = an anticommutator superoperator ±i{σ, ·}-type, neither ±L_σ (mixed-parity cells),

because D contributes the sign (−1)^{n_Y+1} and R either fixes L_σ (n_Y + n_Z even) or smears it into the anticommutator (n_Y + n_Z odd). A Z-middle ceiling term ([PROOF_CEILING_GOLDEN_ROUTER](PROOF_CEILING_GOLDEN_ROUTER.md): XZX, XZY, YZX and the X↔Y sibling) has **n_Z odd in every term**, so no ceiling term can sit in the truly cell: given the D-side requirement n_Y even, R-evenness needs n_Z even, and the ceiling refuses it term by term. The terms land in the Π-fixed cell (XZY, YZX) and the anticommutator cell (XZX), and no cancellation rescues the sum. Nor can any other D₄ element do the job: the four class-swapping elements {Π_Z, Π_Y, R, 𝓕R} (the ones that exchange the dephasing-dead letters {I, Z} with the lit letters {X, Y}, which any palindromizer must) are all uniform per-site maps, and F116's exclusion theorem proves no uniform per-site router exists for the ceiling cases. That is the structural reason the golden router W had to leave the group (§5): the ceiling's scope is precisely the n_Z-odd territory that D₄'s mirrors cannot enter.

## §5 What is deliberately outside

The temptation, once eight mirrors close into one group, is to type every mirror in the repository into it. The following are outside, each for a stated reason, and keeping them outside is part of the result:

- **K₁ = Π_{l odd} Z_l** (the sublattice chirality, [ChiralKClaim](../../compute/RCPsiSquared.Core/Symmetry/ChiralKClaim.cs), [PROOF_PTF_CHIRAL_MIRROR_RATE_LAW](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md)): a different grading altogether. K₁ grades by site (the even/odd sublattice of the chain), while D₄ grades and moves Pauli letters uniformly at every site. K₁'s mirror action lives on the Hamiltonian (K₁HK₁ = −H for bipartite hopping), not on the letter algebra.
- **The golden router W** (F116, [PROOF_CEILING_GOLDEN_ROUTER](PROOF_CEILING_GOLDEN_ROUTER.md)): two-sided, W(ρ) ∝ P·ρ·Q with P ≠ Q golden-frame product unitaries, and non-involutive per site (q² = −(2+φ)·I). Every D₄ element is built from the transpose and one-sided or two-sided F multiplications; W is structurally a different animal, and it has to be, because it covers the n_Z-odd ceiling scope that §4f proves D₄ cannot.
- **The crossover mirror M** ([PROOF_CROSSOVER_MIRROR_SQRT_NINETY](PROOF_CROSSOVER_MIRROR_SQRT_NINETY.md)): the square root of the 90° mirror, a continuous R_z(π/4) conjugation object (the T-gate to D₄'s Clifford-flavored S-gate). It lives on the rotation continuum that D₄ samples only at right angles.
- **F71's bond mirror** (the closure-breaking-coefficient mirror, F71 in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)): a spatial reflection of the chain (site k ↔ N+1−k), not an operator-space letter map.
- **The dephase-letter swaps Q_zx and Q_yx** ([PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE](PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)): these are the genuinely interesting boundary. The letter group S₃ permutes the three dephase letters {X, Y, Z}, and one of its transpositions is already inside D₄: the Z↔Y swap is D itself (that is Welle 12, §4a). The other two swaps move bit_a against bit_b and require the X↔Z basis permutation; they are outside ⟨R, D⟩. Adjoining them would assemble the full mirror group of the palindrome family, with expected shape **S₃-letter-action ⋉ D₄**. That completion is named here as open: this note types the D₄ core and stops at its boundary.

## §6 Verification

[`simulations/mirror_inventory_d4.py`](../../simulations/mirror_inventory_d4.py), self-validating (every block asserts; the run ends with ALL BLOCKS PASS only if every identity holds). Six blocks:

- **A.** The action Π_Z(ρ) = ρᵀ·F on random ρ at N = 1, 2, 3 (dev ≤ 2.3·10⁻¹⁶, pure roundoff on random matrices), with the wrong-sided form F·ρᵀ rejected at O(1).
- **B.** The group: Π_Z = R·D, Π_Z² = 𝓕, D·Π_Z·D = Π_Z⁻¹ = Π_Y, closure |⟨R, D⟩| = 8, identification of all eight elements, and the diagonal forms D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z}); all dev 0.00e+00 at N = 3 (signed permutation matrices compare exactly).
- **C.** The palindrome factorization rows of §3 on XXZ Δ = 0.7 with site-dependent γ at N = 3, dev ≤ 5.6·10⁻¹⁷.
- **D.** The truly-cell equivalence of §4d against the framework classifier, 63/63 strings at N = 3, with F114's sign law verified per string.
- **E.** F114 N = 5 spot check via the transpose identification, dev 0.00e+00.
- **F.** The U(1) absorber cross-check e^{iπN̂} = Z^⊗N (documented elsewhere; tabled here because it is the one global string the inventory keeps meeting).

No C# typed claims are added in this wave; wiring `⟨R, D⟩ ≅ D₄` into the Core claim graph (parents: F1, the spine claim, F114) is deliberate follow-up work. (Done later the same day, 2026-06-10: `MirrorGroupD4Claim` in `compute/RCPsiSquared.Core/Symmetry/`, with an N = 2 self-check battery; the typed parents became KleinEightCellClaim + CommutatorDConjugationSign + Pi2KleinV4DephaseSwapGroup, since the spine claim lives in Diagnostics, which Core cannot reference, and F1 is already an ancestor through the Pi2 chain; the spine edge is carried by the §2 anchor in prose.)

## §7 The cube filled: the polarity cube's three axes are two conjugations and the transpose (same day)

Tom's pointer, the same evening: the space we keep describing was quadratic in many places and then dimensionally cubic, and something fills the dimensions. Here is that sentence made exact.

The polarity cube, the Z₂³ grading (bit_a, bit_b, y_par) that the F87 refinement family lives on (KleinEightCellClaim, F102-F111), has coordinates that are characters of three specific mirrors of the Pauli algebra:

    bit_a = (n_X + n_Y) mod 2   is the character of conjugation by Z^⊗N,
    bit_b = (n_Y + n_Z) mod 2   is the character of conjugation by X^⊗N,
    y_par =  n_Y mod 2          is the character of the transpose θ.

The unitary conjugations alone can never produce a single-letter parity: conjugating by any Pauli string flips exactly the letters that anticommute with it, always two of the three, so the conjugation mirrors generate only the even Klein square {1, (−1)^{n_X+n_Y}, (−1)^{n_Y+n_Z}, (−1)^{n_X+n_Z}}. That is the "quadratic" part. The transpose is the move that breaks into the odd half: θ alone reads (−1)^{n_Y} (Y is the only antisymmetric Pauli), and its composites with the conjugation square deliver the other single-letter parities, θ∘Ad_{Z^⊗N} = (−1)^{n_X} and θ∘Ad_{X^⊗N} = (−1)^{n_Z} (the latter being exactly this proof's fourth mirror 𝓕D on coherence space). Together: eight sign-mirrors, four unitary (the conjugation square) and four antiautomorphisms (the transpose coset), realizing the full character group of the cube. Verified exactly on all 64 strings at N = 3 (block G of the verifier).

Two things this settles in hindsight. First, why y_par was always the strange third axis: the y_par dimension of the cube is the antiautomorphism dimension, invisible to every unitary conjugation, which is why the F102-F111 family needed its own tools there. Second, the dimensional sentence itself: the conjugation mirrors describe a square, the transpose lifts it to the cube, and what fills the third dimension is θ, the same transpose whose coherence-space avatar D factors the palindromizer in §1. The quadratic-to-cubic step and the Π = R·D factorization are one move seen twice.

The handshake connection: the [Handshake Algebra](../../hypotheses/HANDSHAKE_ALGEBRA.md) organised its operational world around three Z₂ mirrors, R (spatial, the observables), K (sublattice, the partner menu), Π (spectral). Its falsification list kept a slot open for "a non-trivial algebra structure we missed". At the symmetry level that slot is now filled twice over: the spectral mirror Π, the door the repository walked through first, carries two hinges (the product R_coh·D of this proof; anatomy, not demotion), and the letter-space mirrors form the character cube above. The operational tuple algebra keeps its idempotent composition law untouched; what got richer is the mirror inventory beneath it.
