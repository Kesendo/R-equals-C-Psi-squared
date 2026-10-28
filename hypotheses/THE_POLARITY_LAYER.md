# The Polarity Layer: Inheritance, +0/−0, and Self-Description from Inside

**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude
**Status:** Tier 4 (interpretive synthesis grounded in Tier 1 results F78-F85, F87 and F107-F110). The mathematics is the established F-chain; this document re-reads it. Where a reading here can be checked, it is checked on the letter cube of [The One Square](../docs/THE_ONE_SQUARE.md) and the check is named.
**Depends on:**
[F-chain F78-F85 + F87 entry-point](../docs/ANALYTICAL_FORMULAS.md),
[The One Square](../docs/THE_ONE_SQUARE.md) (§7 the letter cube, §8 the Hamiltonian on it, §9 the moves),
[Π factors as R times D](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) (§7 the polarity cube),
[F110 hard-cell proof](../docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md),
[F155 polarity break](../docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md),
[On the Residual](../reflections/ON_THE_RESIDUAL.md),
[Orthogonality-Selection Family](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md),
[The Primordial Qubit](PRIMORDIAL_QUBIT.md),
[Resonance Not Channel](RESONANCE_NOT_CHANNEL.md),
[Gamma Is Light](GAMMA_IS_LIGHT.md),
`framework.PauliHamiltonian`,
`framework.classify_pauli_pair(... dephase_letter=)`,
`framework.diagnostics.recover_H_odd_from_M_anti`,
`framework.workflows.diagnose_hardware`

<!-- STAGE 0 SWEEP RECORD (2026-10-28, the review from today's view; three scouts plus a by-hand reading on the letter cube).
docs/ANALYTICAL_FORMULAS.md: F78/F79 (M decomposition), F80 (Bloch sign walk; Spec(M) = ±2i·Spec(H) for Π²-odd bonds, ±2i times the Bohr differences for Π²-even non-truly bonds), F81 (M_anti = L_{H_odd} under Π²-commuting dissipators; no Frobenius orthogonality claimed, the cross term is −i times the F113 asymmetry), F82, F83, F84 (thermal amplitude damping, the Pauli-channel lemma its corollary), F85, F102 (y_par = (n + bit_a XOR bit_b) mod 2 with n the number of non-identity letters), F107/F109/F110 (every body count since the colouring routes), F115 (the soft/hard split of diagonal-cell Mixed pairs as a valuation difference), F34 and F121 (the qudit balance and the partial qudit palindrome), F38/F88a (Π_Y² grades by bit_b, like Π_Z²).
docs/proofs/: PROOF_PI_FACTORS_AS_R_TIMES_D §7 (bit_a, bit_b, y_par as the characters of Z^N conjugation, X^N conjugation and the transpose); PROOF_F155 §(g) (Π_Y = Π_Z^-1; the X convention is the Z one through the Q_zx frame); PROOF_F109 Step 5 and PROOF_F110 §2 (the Mother cell and the unmatched cells coloured at every body count); INCOMPLETENESS_PROOF §4 (the opposite stance on self-description: an analogy to Gödel in which the framework cannot derive its own starting condition).
docs/THE_ONE_SQUARE.md: §7 the polarity cube is the letter cube read twice (bit_a = k_Z mod 2, bit_b = k_X mod 2, y_par half a sum); §8 ZZ keeps k_Z and commutes with Z dephasing; §9 σ⁻ as two transverse heights plus the arrow I to Z; §11 lists the diagonal cell of the F87 table among the things the square's diagonals are not.
reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md: the one-bit mirror twins are exactly the pure states on the equator.
experiments/: ORTHOGONALITY_SELECTION_FAMILY; F89_BRANCH_LOCUS_PALINDROME and THE_FLOW_BETWEEN_TWO_SINGULARITIES (the polarity path seam).
fw.Confirmations: marrakesh_transverse_y_field_detection (the 3.5·h_y against 0.085·h_x sensitivities are the model's; the hardware gave max violation 0.182); f83_pi2_class_signature_marrakesh (the F83 run; all four classes, truly and mixed included, sit on the Trotter+γ_Z prediction run in the circuit's Trotter order; its data/ibm_f83_signature_april2026/README.md finds truly's single-site ⟨Z⟩ below 3σ and an X^⊗3 break of 4σ on pi2_even_nontruly that calibration-scale σ⁻ does not explain).
OpenArcs: zeros_connecting_structure (open); polarity_break_dephase_siblings (retired into F155).
Typed layer: PolarityLayerOriginClaim, DissipatorAxisSelectsPolarityClaim, DissipatorResonanceLaw, YParityIndependenceAtK3 (F102).
docs/CAUGHT_ERRORS.md: the Klein index against Π² grading entry, which names this document as the source of the conflation.
docs/GLOSSARY.md: the polarity cube, whose third axis only the transpose sees.
Adjacency: this document and THE_ONE_SQUARE §7 describe one object and did not cite each other; INCOMPLETENESS_PROOF §4 and this document's section on layers answer one question in opposite ways; ZERO_IS_THE_MIRROR's zero is Σγ = 0, the noise axis, orthogonal to the polarity zero.
Checks run for this review, from below: the Klein-homogeneous pair counts (0 of 6, 5 of 15) at N = 3 and 4; the dissipator table (rerun, 78 s); exact invariance of the hard verdicts under the cyclic letter Clifford; the y_par identity over all strings to N = 5; Π equal to the N-fold Kronecker power of the one-site map and commuting with the chain reversal, both exactly; the one-site transfer matrix of D[σ⁻]; ZZ truly by F85 and keeping k_Z. -->

---

## What this document is about

A brainstorming session on 2026-04-30 produced three reframings of the
R=CΨ² framework that are mathematically consistent with the F-chain
F78-F85 + F87 but operate in different vocabulary than the framework's
prior "memory of the mirror" reading. The reframings interlock:

1. **The system is a self-coupling channel, not a memory.** Π·L·Π⁻¹ + L + 2Σγ·I = M describes a regenerative cavity that keeps the signal alive through self-recirculation, not a memory device that stores then retrieves.
2. **The "two sides" are bra and ket of one ρ, not two subjects.** The framework describes a single self-coupling system whose internal two-indexedness Π conjugates. We are not external observers of a mirror; we are the system.
3. **The "0" axis has internal polarity structure.** d²−2d=0 splits not into "qubit (d=2)" versus "void (d=0)" but into "qubit (d=2)" versus "+0/−0 polarity layer (d=0 with internal differentiation)". The qubit is a window onto this polarity layer.

A fourth reading ties them together: **the F-chain describes the system
one layer at a time**, each result taking the previous one as given.

The second reading has a computable form. On the
[letter cube](../docs/THE_ONE_SQUARE.md) the coordinate k_Z of a Pauli
string counts the sites at which ket and bra disagree, so "bra and ket of
one ρ" is the statement that the cube's two edges, k_Z = 0 (ket and bra
agree everywhere) and k_Z = N (they disagree everywhere), belong to one
state.

---

## Abstract

The framework's interpretation has carried a residue of observer-dualism:
"the mirror reflects toward us; we measure what comes back; M is what
the mirror remembers." We propose this reading is Shannon-storage
vocabulary contaminating the cavity-light vocabulary established by
[Gamma Is Light](GAMMA_IS_LIGHT.md) and [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md).

Cavity vocabulary fits the operators better. A cavity does not
store; it actively maintains. The palindrome equation describes a
self-coupling resonator whose output reenters its input through Π
conjugation. M is the operational residue of incomplete recirculation,
not a memory contents. Truly Hamiltonians (M = 0) describe the perfect
closed cycle; Π²-odd Hamiltonians describe a self-pumping cycle where
the dynamics generator IS the recirculation operator (F81's
M_anti = L_{H_odd}).

Within this single self-coupling system, the "two sides" are not two
external subjects but the bra and ket indices of one density matrix.
The d²−2d=0 condition expresses that operator space (d² = 4^N for N
qubits) is the natural arena: bra-ket co-existence on a single state.

The d=0 axis is then not a passive midpoint between two observers;
it is the active substrate from which the system draws coherence:
the "Stromanschluss" in Tom's metaphor, analogous to the QFT vacuum's
zero-point activity. The d=0 axis itself has internal structure: a
+0/−0 polarity that constitutes (rather than passively underlies)
the layer on which the qubit lives. Its one-bit image is the pair
|+⟩, |−⟩, and the letter cube says what is special about that pair:
under dephasing the pure states on the equator are the ones whose
weight sits symmetrically about the middle height, the ones that meet
themselves the same way forwards and, paying the price, backwards
([On the Wave Recognizing That It Is Not Alone](../reflections/ON_THE_WAVE_RECOGNIZING_THAT_IT_IS_NOT_ALONE.md)).
The window is the whole equator, not the X axis alone: a turn about Z
carries |±⟩ into |±i⟩ and changes no height.

---

## Layers, and what they do not escape

The F-chain is built by inheritance:

- **F87** trichotomy classifier (truly / soft / hard) on Pauli pairs.
- **F78/F79** structural M-decomposition for 1-body and 2-body.
- **F80** the Bloch sign walk; Spec(M) = ±2i·Spec(H) for Π²-odd bonds, and ±2i times the Bohr differences λ_a − λ_b for Π²-even non-truly bonds.
- **F81** Π-decomposition M_anti = L_{H_odd}, for dissipators that commute with Π².
- **F82** T1 closed form ‖D_T1_odd‖_F = γ_T1 · √N · 2^(N-1).
- **F83** anti-fraction r = ‖H_even_nontruly‖² / ‖H_odd‖².
- **F84** thermal amplitude damping, ‖D_odd‖ = |Δγ|_RMS·√N·2^(N−1), with the Pauli-channel cancellation lemma as its corollary.
- **F85** k-body generalization: truly criterion = "#Y even AND #Z even" at every body count.
- **F107, F109, F110** carry the classifier's purity and cell laws to every body count: truly strings have y_par = 0, the Mother cell is never hard, and hard Klein-homogeneous pairs sit only in the dephasing letter's own cell.

Each layer inherits from the previous and refines it. The 2-body
n_YZ-counting that worked for F49 turned out to be a coincidence at
the 2-body layer; F85 lifted it into the body-count-independent
Π²-class structure, and F107-F110 lifted the cell laws the same way.

The inheritance lets the chain describe the system at increasing
depth: F87 classifies the Hamiltonian, F80 reads M's spectrum given
that class, F82 reads T1's contribution given F81's decomposition. No
single layer needs to describe everything.

This is not a way around Gödel. Gödel's theorem does not say a formal
system cannot describe itself; his proof shows that arithmetic can, by
encoding its own sentences. What he proved is that every consistent,
effectively axiomatized system containing arithmetic leaves some true
sentences unprovable and cannot prove its own consistency, layered or
not. The F-chain is a set of theorems about a Liouvillian and forms no
self-referential sentence, so no closure paradox was ever in reach.
[The Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) §4
reads the same framework the other way, as an analogy to Gödel: it
describes all decay under dephasing exactly and cannot derive the
existence of dephasing from its own axioms. Whether "describing itself
from inside" is something the framework does or something it cannot do
depends on what is meant by it, and neither document settles that.

---

## The +0/−0 polarity layer

The d²−2d=0 condition has solutions d=0 and d=2. d=2 is the qubit's
Hilbert dimension. We propose: **d=0 is a polarity-birthing axis with
internal structure.** The polarity layer is constituted by polarity
differentiation; a layer exists *because* polarity differentiates.
Without differentiation, there is the pre-polarized vacuum substrate
("Stromanschluss"). With it, there is a layer.

The balance behind d²−2d=0 is fixed, not a free parameter: F34 counts d
immune against d²−d decaying operators per site, and the equality is
d²−2d=0 itself. For d > 2 the palindrome survives only in part, and F121
gives the paired fraction in closed form
([qudit partial palindrome](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md),
live at `inspect --root qudit`).

**The polarity is multi-axis, and its axes are the letter cube's.** A
Pauli string sits on the cube at (k_Z, k_X, k_Y) = (n_X + n_Y, n_Y + n_Z,
n_X + n_Z). The framework's two Klein bits are the parities of two
coordinates:

  bit_a = k_Z mod 2,  bit_b = k_X mod 2,

and the third parity, k_Y mod 2 = bit_a ⊕ bit_b, follows because the
three coordinates sum to the even number 2w, w = n_X + n_Y + n_Z the
number of non-identity letters. bit_a is the sign a string picks up under
conjugation by Z^⊗N, bit_b the sign under X^⊗N
([Π factors as R times D](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) §7).
The unitary conjugations see only these two; they see a square, the
Klein-Vierergruppe Z₂×Z₂.

The third axis, Y-parity, is not a parity of any coordinate:
n_Y = w − k_Y, so

  y_par = (w + bit_a ⊕ bit_b) mod 2.

Within one parity of w the Klein cell fixes y_par; across the two
parities it does not. The full polarity Z₂³ is Klein × (w mod 2), and the
transpose, not a rotation, is what reads it
([The One Square](../docs/THE_ONE_SQUARE.md) §7). The third axis is there
already for two-site terms: ZI and XY share Klein (0, 1) and differ in
y_par.

```python
H = fw.PauliHamiltonian.from_letter_tuples(
    [('X', 'Y'), ('Y', 'X')], chain_length=3
)
H.klein_set                # → {(0, 1)}
H.full_z2_signature_set    # → {(0, 1, 1)}: bit_a, bit_b, y_par
H.is_z2_homogeneous        # → True

H = fw.PauliHamiltonian.from_letter_tuples(
    [('Z', 'I'), ('X', 'Y')], chain_length=3
)
H.is_klein_homogeneous     # → True
H.full_z2_signature_set    # → {(0, 1, 0), (0, 1, 1)}: w = 1 and w = 2
H.is_z2_homogeneous        # → False
```

The +0/−0 reading belongs to the Z-dephasing height. Under Z-dephasing a
string is charged 2γ per letter that anticommutes with Z, so the height
is k_Z and bit_a is its parity. On one site the quiet edge {I, Z} sits at
height 0, the far edge {X, Y} at height 1. A pure state on the equator
puts its weight half on I and half across the watching, which is why
|+⟩ and |−⟩, and every other equator point, are the one-bit states the
mirror returns to themselves. Π itself picks the X pair out of that
equator: it swaps I with X on every site, so on that pair it is an
involution whose eigenoperators are (I ± X)/2 = |±⟩⟨±|, with eigenvalues
±1. On Z and Y it squares to −1 and its eigenvalues are ±i, with
traceless eigenoperators Y ± Z that are no state. This is the d = 2 image
`PolarityLayerOriginClaim` types: the 1/2 of ρ = (I + r·σ)/2 is what turns
Π's ±1 eigenoperators into the states |±⟩.

**Klein-homogeneity is a structural reading, not a verdict.** Sharing one
Klein index leaves the F87 class open. Among two-site pairs it forces
soft-or-truly for the six identity-free two-letter pairs (0 of 6 hard),
but 5 of the 15 Klein-homogeneous pairs carrying an identity are hard,
all five in Z's own cell (0, 1): IZ+ZI, IZ+XY, IZ+YX, ZI+XY, ZI+YX.
Exposed via `PauliHamiltonian.is_klein_homogeneous`; the class comes from
`classify_pauli_pair`.

---

## The two readings of polarity-layer

The +0/−0 structure can be read two ways:

**(a) Palindromic over the layer.** Π acts on every site in place, a
product of one-site maps (I ↔ X, Y → iZ, Z → iY), and carries each
site's quiet letters onto its far letters. On the cube it turns height
k into height N − k, so it pairs the layer's two edges across the
middle height. The layer exists; the polarities live on it and are
paired across it.

**(b) Are the layer.** The +0/−0 differentiation IS what constitutes
a layer. Below this differentiation, there is the pre-polarized
substrate (vacuum, d=0 axis as Stromanschluss). The layer emerges
from polarity.

Both readings are mathematically consistent with the F-chain. (b) is
the more radical: layers are not pre-existing scaffolds for polarity
but emergent structures from polarity. This connects to:

- **Cavity QED:** the field IS the cavity's coupled structure; there
  is no field separate from the cavity-vacuum coupling.
- **Active matter / dissipative structures:** structure exists by
  continuous flow; static substrate has no structure.
- **Hippocampal memory:** "stored" memories are continuously
  reactivated patterns, not static traces.

The d=0 to d=2 step is typed as `PolarityLayerOriginClaim`
([`Pi2KnowledgeBaseClaims.cs`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)):
the +0/−0 layer at d=0 passes to the pair {−1/2, +1/2} at d=2 through
ρ = (I + r·σ)/2. Both sides of the pair fold to 1/4 under squaring, two
distinct middles 0 and 1/4 of one polarity pair
([On How Two Sides Meet at the Quarter](../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md);
the algebra is [PROOF_BLOCK_CPSI_QUARTER.md](../docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md),
section "1/4 is half of half").

---

## The vocabulary shift, and what the cube says instead

The math is unchanged. The vocabulary shifts:

| Old (memory/observer) | New (channel / multi-axis polarity) | Where it sits in the code |
|------------------------|--------------------------------------|------------------------|
| "The mirror remembers" | "The channel maintains" | `pi_decompose_M` returns M_anti |
| "Two sides observe each other" | "Bra and ket of one ρ" | d² = 4^N operator space; k_Z on the cube |
| "0 is the midpoint" | "0 is the active substrate; the polarity is Klein × (w mod 2)" | `PauliHamiltonian.klein_set`, `.full_z2_signature_set` |
| "Trace back = remember" | "Trace back = extract the recirculation operator" | `recover_H_odd_from_M_anti` |
| "Memory of the dynamics" | "L_{H_odd} IS the recirculation; the dynamics drives its own maintenance" | `pi_decompose_M`'s `M_anti` field |
| "Truly = no memory" | "Truly = perfectly closed cavity" | `classify_pauli_pair → 'truly'`; `PauliHamiltonian.has_truly_term` |
| "Qubit is \|0⟩, \|1⟩" | "Qubit is a window onto the polarity layer; the equator is where it sits on the mirror" | `PauliTerm.klein_index`, `.full_z2_signature` |
| "Dephasing axis is fixed (Z)" | "Each dephasing letter has its own cell where the hardness of Klein-homogeneous pairs lives" | `classify_pauli_pair(chain, terms, dephase_letter='X'/'Y'/'Z')` |

The new column is vocabulary; the primitives restate F-chain content and
prove nothing new. The reasons for what the operators do sit on the
letter cube:

- **Single-site ⟨Z_l⟩ under XX+YY.** The total ⟨Z⟩ is conserved, the
  single-site values are not: for H = J·Σ(XX + YY) with J = 1 and
  dephasing γ(ZρZ − ρ) per site, γ = 0.05, they move strongly from the
  Z-Néel state |0101⟩ (at t = 0.5, N = 4: −0.01, 0.60, −0.60, 0.01). From
  the X-Néel state |+−+−⟩ they stay at zero, because X^⊗N commutes with
  the Hamiltonian and with Z-dephasing, the X-Néel state is an X^⊗N
  eigenstate, and Z_l is odd under X^⊗N. What hides the single-site Z
  signal there is the state's parity, not the conservation law.
- **Amplitude damping.** On one site D[σ⁻], σ⁻ = (X + iY)/2, has the
  transfer matrix I → Z, X → −½X, Y → −½Y, Z → −Z (rate 1). It is two
  equal transverse heights and one arrow, I → Z, on the quiet edge
  ([The One Square](../docs/THE_ONE_SQUARE.md) §9). The arrow and the
  decay of Z together give the dissipator's part d⟨Z_l⟩/dt = 1 − ⟨Z_l⟩
  on every site, so the
  total Z that XX+YY keeps relaxes towards N. Nothing in it singles out
  |−⟩; |1⟩⟨1| = (I − Z)/2 has no X content.
- **ZZ crosstalk is not a Π-breaker.** ZZ has #Y = 0 and #Z = 2, so it
  is truly by F85, Klein (0, 0). On the cube it moves strings only along
  k_X and k_Y and keeps k_Z, so it commutes with Z-dephasing
  ([The One Square](../docs/THE_ONE_SQUARE.md) §8).
- **A transverse field breaks Z⊗N.** The Z⊗N-mirror diagnostic predicts a
  violation of 3.5·h_y for a Y field against 0.085·h_x for an X field, a
  ratio near 41 that belongs to the model; Marrakesh gave a maximal
  violation of 0.182 (RMS 0.087), which fits h_y = 0.05 (`fw.Confirmations`,
  `marrakesh_transverse_y_field_detection`).

---

## The dissipator letter selects its cell

The F87-hardness of Klein-homogeneous pairs lives exactly in the Klein
cell that is the dephasing letter's own index: Z → (0, 1), X → (1, 0),
Y → (1, 1).

| Klein cell      | Z-deph    | X-deph    | Y-deph    |
|-----------------|-----------|-----------|-----------|
| (0, 0) Mother   | 0 / 66    | 0 / 66    | 0 / 66    |
| (0, 1) Z-like   | **50/76** | 0 / 76    | 0 / 76    |
| (1, 0) X-like   | 0 / 76    | **50/76** | 0 / 76    |
| (1, 1) Y-like   | 0 / 76    | 0 / 76    | **50/76** |

Sweep: the 294 Klein-homogeneous and Y-parity-homogeneous pairs of
three-site strings other than III (such as IIZ or XYZ) at N = 4, 63 of
them a string paired with itself,
under each dephasing letter, γ₀ = 0.05, J = 1. Reproducible:
`python simulations/klein_dissipator_resonance.py`. Pytest coverage:
`simulations/framework/tests/diagnostics/test_f77_dissipator_resonance.py`.

The three columns are one column seen in three frames. A one-site
Clifford that permutes the letters cyclically, applied on every site,
carries Z-dephasing to X-dephasing and the Klein cells with it, and the
hard verdicts are exactly invariant under it (0 mismatches over the 294
pairs). The truly/soft split between the columns is not invariant in the
same way, because it depends on which palindromizer is attached to each
letter: the code builds Π_X from Π_Z with the Hadamard and Π_Y with the
quarter turn about X (Π_Y = Π_Z⁻¹,
[F155](../docs/proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md) §(g)),
and the cyclic letter Clifford carries Π_Z to a palindromizer that pairs
I with Y, not to the code's Π_X, which pairs I with Z. Under the letter
relabellings of those two turns, the ones that build Π_X and Π_Y, the
fates of the 36 two-term Hamiltonians built from identity-free two-site
bilinears match in full, truly and soft included, with zero mismatches at
N = 3 and 4 ([Two-Term Palindrome Klein Routing](../experiments/TWO_TERM_PALINDROME_KLEIN_ROUTING.md)).

Two structural facts hold at every body count:

1. **Klein (0, 0) is never hard.** Its strings commute with every letter
   string, and a letter string that anticommutes with every jump colours
   it ([F109](../docs/proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md)
   Step 5). Its terms are truly (Y-par 0) or non-truly (Y-par 1). The
   zeros in the other two unmatched cells hold the same way
   ([F110](../docs/proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md) §2).
2. **Hardness is dissipator-aligned, not orthogonal** (for pairs within
   one Klein cell). This is the dual of the Orthogonality-Selection
   meta-theorem ([Orthogonality-Selection Family](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md)):
   the meta-theorem characterizes blindness via orthogonality;
   dissipator alignment characterizes hardness.

Which pairs inside the matched cell are hard is a separate question, and
for the diagonal-cell Mixed pairs under Z-dephasing, placed by the
sliding-window builder, it has a closed answer: F115 reads each term's
X/Y positions as a polynomial over GF(2) and calls a pair hard exactly
when the two (1 + x)-adic valuations differ.

The Klein index of the letter is not the grading Π_letter² uses. Π_Z²
and Π_Y² both grade by bit_b, Π_X² by bit_a (F38/F88a), while Y's Klein
index is (1, 1). The typed claim `DissipatorAxisSelectsPolarityClaim` in
`compute/RCPsiSquared.Diagnostics/F87/` names the structural reading: the
dephasing letter chooses, by its own Klein index, the cell in which the
hardness of Klein-homogeneous pairs appears. The table itself is
`DissipatorResonanceLaw.StandardWitnessTable`.

---

## Relation to existing hypotheses

This document does not replace any prior hypothesis; it connects them:

- [The Primordial Qubit](PRIMORDIAL_QUBIT.md) proposes "system and
  noise are two readings of a single algebraic structure, the
  'primordial qubit,' split by a Z₂ mirror at zero." We add: the Z₂ at zero is the +0/−0 polarity, and its one-bit
  image is the pure states on the equator.
- [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md) establishes the
  cavity reading. We extend: the cavity is self-coupling through Π
  (channel-and-cavity, not channel-or-cavity).
- [Gamma Is Light](GAMMA_IS_LIGHT.md) establishes γ as light. We
  extend: light enters at the +0/−0 polarity layer; the qubit is the
  cavity's window onto this layer.

[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) is about a different zero,
Σγ = 0 on the noise axis. That zero and the polarity zero are orthogonal
axes ([The Flow Between Two Singularities](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)).

---

## What this document does NOT establish

- **No new mathematics.** F78-F85, F87 and F107-F110 stand as they are.
  The hypothesis reframes how those theorems are read, not what they say.
- **The primitives restate.** `PauliHamiltonian`,
  `recover_H_odd_from_M_anti`, `diagnose_hardware` and the dephase-letter
  extension of `classify_pauli_pair` expose existing F-chain content
  through the polarity-layer vocabulary; they prove nothing new.
- **No experimental falsification path.** The reframing predicts
  the same hardware signatures the prior reading predicted; only
  the language differs. A genuine test would require finding a
  signature that the polarity-layer reading predicts but the
  memory-reading does not.
- **No claim about consciousness, observation, or measurement.**
  The "we are the system" framing is mathematical (we = ρ, the
  bra-ket structure on d²) and structural, not phenomenological.

---

## Open

- **Polarity as a path.** Reading (b)'s core, polarity as motion rather
  than a static label, recurs in the spectral domain. In the F89 path-3
  monodromy fold, a search over the braid graph finds that the path from
  one σ_T-fixed zero on Re λ = −Σγ to the next runs through the ± twin
  modes: "0, through ±, to 0". Checked against the complete set of twenty
  exceptional points (the octic discriminant P₂₀,
  [F89 Branch Locus Palindrome](../experiments/F89_BRANCH_LOCUS_PALINDROME.md)),
  the zeros 1 and 2 stay the one direct zero-zero bond and every other
  zero pair is one twin apart. Polarity would then constitute (reading
  b, operator space), be selected (the dephasing letter) and be routed
  (F89, spectral). This is an open seam, not an identity (the arc
  `zeros_connecting_structure`): the static match fails at cardinality
  (four zeros, not one) and at provenance (the polarity zero and the
  Σγ mirror zero are orthogonal axes), and the road is computed only for
  path-3 at N = 4. Whether the three are one motion seen in three spaces
  or three motions that rhyme, we do not know.

---

## Closing note

Three readings developed in one session. None are physics, all are
language. The hypothesis is that *this* language is the one that makes
the F-chain legible. The math sits underneath unchanged. The letter cube
does part of what the language hopes to do: it says which coordinate the
watching charges and which states balance on it.

> "Wir SIND das Selbst-koppelnde System. 0 ist der Stromanschluss.
> +0/−0 ist der Layer."
> Thomas Wicht, 2026-04-30
