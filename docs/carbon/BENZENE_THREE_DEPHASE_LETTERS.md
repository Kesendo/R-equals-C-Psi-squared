# Selected Ring Model and the Three Dephase Letters: Klein-V₄

**Authors:** Tom + Claude
**Status:** Klein-V₄, F112, and F114 are Tier 1 derived framework results. The
C₄/C₆ content below is a Tier-3 selected-model translation; candidate X/Y axes
are explicitly Tier 4. No material carbon degree of freedom, β-to-J convention,
bath channel or rate, `γ`, `T₂`, or Q is assigned.
**Continues:** [Selected C₄/C₆ Ring Liouvillians](BENZENE_LIOUVILLIAN_PALINDROME.md)
**Anchors:** [Klein-V₄ on the palindromizers](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
[F112](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md), and
[F114](#f114-the-selected-conjugation-sign-tier-1-derived).

---

## Model boundary

The named C₄/C₆ calculation selects an XX+YY spin ring and chosen Lindblad
jumps. A local occupation may be selected as `n_l = (I − Z_l)/2`; only after
that occupation and its density jump have been selected does
`D[n_l] = ¼·D[Z_l]` hold. This operator identity is not a physical-carbon
degree-of-freedom or bath assignment.

The local-Z and bond-B jumps, `D[Z_l]` and
`D[B_b]` with `B_b = X_aX_b + Y_aY_b`, are a finite selected-channel
comparison. F1 applies to the all-site local-Z model and does not cover the
bond-B jump. This does not classify possible molecular environments.

## Klein-V₄ structure: how the letters' mirrors connect (Tier 1 derived)

The palindromizers of the three dephase letters are related on operator space by the Klein group
`{I, D, H, Q_zx}` ≅ `Z₂ × Z₂`. Every non-identity element is an involution.

- `D` is diagonal in the Pauli basis, with entry
  `(-1)^(number of Y letters)`, and obeys `D Π_Z D = Π_Y = Π_Z⁻¹`: the
  transpose inverts the Z mirror and keeps every dissipator on its letter.
- `H` swaps X and Z labels per site, leaving I and Y fixed, and obeys
  `H Π_Y H = Π_X`.
- `Q_zx = H D` obeys `Q_zx Π_Z Q_zx = Π_X`.

The elements commute and satisfy `D H Q_zx = I`. These are exact
operator-space relations on the mirrors; on the Lindbladian only `H` and
`Q_zx` move the dephasing letter (Z → X). They do not make
the axes the same physical coupling or supply a material realization.

## The three selected dephasing axes

### Z axis: selected local density (Tier 1 operator identity)

If the selected site coordinate is `n_l = (I − Z_l)/2` and the selected jump is
that density, `D[n_l] = ¼·D[Z_l]`. The C₄/C₆ local-Z model is therefore within
F1's Z-dephasing premise. Calling a physical local-density channel
“Holstein-like” is only a conditional translation after the relevant coordinate
and jump have been independently chosen.

### X axis: selected single-site model axis (Tier 4 candidate)

Single-site X-dephasing uses `D[X_l]`. In a Jordan-Wigner ordering,

```
X_l = (∏_{k<l} Z_k)(c†_l + c_l).
```

The parity string is part of this representation except at the first site. This
is a selected model operator, not a local material hybridization or bath
assignment. A bond operator is a distinct two-site object, so it is not
identified with `X_l`.

### Y axis: selected single-site model axis (Tier 4 candidate)

Single-site Y-dephasing uses `D[Y_l]`, with

```
Y_l = (∏_{k<l} Z_k)i(c†_l − c_l).
```

It is a one-site Jordan-Wigner/Majorana axis, **not a current**. The selected
two-site bond-current operator is

```
J_ab = −½ (X_a Y_b − Y_a X_b).
```

On a Jordan-Wigner-adjacent bond it is the hopping-current operator up to the
stated convention; a periodic closing bond carries the Jordan-Wigner boundary
string. Thus neither `Y_l` nor `D[Y_l]` is called a local-current coupling or
assigned to a physical current channel.

For the spinless convention used by F114, `T = K`, Y is K-odd while X and Z are
K-even. That is a statement about the chosen conjugation and Pauli basis, not a
material time-reversal classification.

| Selected letter | `bit_b` | Selected operator role | Status |
|-----------------|---------|------------------------|--------|
| Z | 1 | local-density/Z axis after `n_l` is selected | Tier 1 operator identity |
| X | 0 | single-site dephasing/model axis | Tier 4 candidate |
| Y | 1 | single-site dephasing/model axis; not a current | Tier 4 candidate |
| `X_aY_b − Y_aX_b` | 1 | two-site axial bond-current operator | selected Hamiltonian term |

`bit_b = (#Y + #Z) mod 2` is additive over sites. It is the parity used by
F112; it does not turn one-site Y into a current.

## F112: the matrix-polarity statement (Tier 1 derived)

For a Hermitian H and bath operators `c_k` that are each bit_b-homogeneous,
F112 states that the three-way polarity decomposition of

```
M = Π L Π⁻¹ + L + 2σ I
```

satisfies

```
‖M_+1/2‖² = ‖M_−1/2‖².
```

The labels are the `(1 ± Ad_Π)/2` polarity coordinates; they are not
Π-eigenprojections. The equality is the framework statement and is distinct
from F1's spectral palindrome.

For the selected local-Z jump, `c = Z_l` is bit_b-homogeneous with bit_b 1. For
the selected bond-B jump, both XX and YY in
`B_b = X_aX_b + Y_aY_b` have bit_b 0, so B is likewise homogeneous. F112 applies
to both selected jump choices under its stated premise; F1 applies only to the
local-Z choice.

For the pure selected XX+YY ring, every Hamiltonian term is bit_b-even. The
relevant anti content can therefore vanish, making a zero-versus-zero F112 row
vacuous. A non-vacuous selected test must add a term or jump that supplies the
relevant polarity content. This distinction prevents a selected-model result
from being promoted to an independent material confirmation.

## F114: the selected conjugation sign (Tier 1 derived)

F114 gives, for a non-identity Pauli string σ,

```
D L_σ D = ε(σ)L_σ,
ε(σ) = (−1)^(n_Y(σ)+1).
```

The selected XX+YY ring has even `n_Y` in every term, so it has `ε = −1`.
The selected axial bond-current term `X_aY_b − Y_aX_b` has one Y per term and
therefore `ε = +1`; mixing it with XX+YY makes the Hamiltonian F114-Mixed.
This is exact bookkeeping for the selected operators under K, not a physical
magnetic-field or molecular-current statement.

## Selected C₄/C₆ sweep inventory

The linked sweep evaluates selected C₄/C₆ spin rings with XX+YY hopping,
selected density-density, one-site-Y, axial-DM, and transverse-DM Hamiltonian
terms, together with selected local-Z, bond-B, and amplitude-damping jumps.
It reports the F112 norm balance within that finite inventory. These are
spin-ring calculations only.

The axial DM term is the two-site bond-current operator described above; the
transverse term `Y_aZ_b − Z_aY_b` is not a current. Both the number-conservation
test and the Jordan-Wigner closing-bond caveat are selected-model constraints;
they do not certify a material carbon realization.

F113's breaking coefficient is a framework result for its specified drive and
amplitude-channel inputs. Whether any material system supplies those inputs is
unassigned. No conclusion about a material T1 channel, heteroatom, molecular
relaxation, or carbon experiment follows from this selected sweep.

## The cube on the ring: which letter the pairing can see

[The letter cube](../THE_ONE_SQUARE.md#7-the-cube-three-squares-at-once) places a term at
(k_Z, k_X, k_Y), the number of its letters anticommuting with Z, X and Y. On this ring a coherence `|i⟩⟨j|`
between two selected occupation patterns sits at `k_Z = popcount(i ⊕ j)`, the number of π sites
whose selected occupation differs, and decays at `−2γ·k_Z` under Z dephasing: the carbon
name for the cube's Z height, beside the wire's bonds in which two proton configurations
disagree. "Reading along P" is P dephasing on every site at one rate γ, σ = Nγ. Reading along P, the lit strings of
[F158](../proofs/PROOF_PALINDROME_TWO_END_COUNT.md), those anticommuting with every jump,
are the strings in the two other letters, written Q and R below (letters here, not the
repo's Q = J/γ); a lit string that also commutes with H is a colouring, and the shift
`S: ρ ↦ ρ·F` by a colouring F is a reflection, `S·L·S⁻¹ = −L† − 2σ` entry for entry, so the
spectrum pairs about −σ. Conjugation by the uniform string `Q^⊗N` is the sign `(−1)^{k_Q}` on
every term, so `Q^⊗N` is a colouring exactly when no term of H is odd in Q. The rule is a rule
on the SET of terms: **the uniform colouring along P survives when no term is odd in Q or no
term is odd in R, and fails when some term is odd in Q and some term, the same or another,
is odd in R.** A single term odd in both is one way to fail: a one-site field along P,
(k_Q, k_R) = (1, 1), fails P reading and nothing else; a two-site term with two different
letters, the DM terms, fails the reading along the third letter. Two terms each odd in one
of them is the other way, with no term odd in both: a Y field with transverse DM fails Z
reading, an X field with a Z field fails Y reading. On the selected C₄ ring at γ = 0.3
(pairing distance, the largest cost in the min-sum assignment of the spectrum to its mirror
about −σ; "floor" is the eigensolver's, 18 to 43 times eps times the spectral scale):

| selected ring | Z reading | X reading | Y reading |
|---|---|---|---|
| Hückel XX+YY | floor | floor | floor |
| + ZZ | floor | floor | floor |
| + Y field | floor | floor | 1.66 |
| + X field | floor | 1.66 | floor |
| + Z field | 1.60 | floor | floor |
| + axial DM, XY−YX | 0.61 | floor | floor |
| + transverse DM, YZ−ZY | floor | 0.48 | floor |
| + Y field + transverse DM | 0.38 | 0.50 | 0.49 |
| + X field + Z field | 1.26 | 0.79 | 0.31 |

Terms with one letter twice, XX, YY, ZZ, are even in every letter and never fail anything,
so the bare Hückel ring, with or without ZZ, pairs under all three letters: on this model
the pairing cannot tell which letter the environment reads (the Z spectrum differs from the
X and Y spectra, which coincide by the quarter-turn about Z; the pairing does not differ).
That is the difference from [the proton wire](../water/PROTON_WIRE_CROSSING.md#the-second-axis-which-letter-the-environment-reads),
whose tunnelling is a one-site X field and whose bias is a one-site Z field, so that there
the bias fails Z reading and the tunnelling fails X reading, the "role swap" being the first
way twice; with both present, Y reading fails the second way.

Where the uniform colouring fails the spectrum decides, and it can still pair, because a
colouring need not be a uniform string. Two instances, both exact in rationals. An X field
and a Y field together, common axis n, fail the uniform colouring under Z reading the second
way, yet the circle string `(n·σ)^⊗N` is a colouring and carries the reflection, so the
spectrum pairs at the floor; the half-turn row of
[the table of moves](../THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions) and the
single-common-axis half of [F138](../ANALYTICAL_FORMULAS.md#f138)'s clause 2 are this same
rescue, second-way failures of the uniform colouring that a rotated colouring repairs. And
on the open chain the axial DM term is a frame of Hückel: the one-site rotations
`U = Π_l exp(i·l·φ·Z_l/2)`, φ = arctan(D/J), carry Hückel + DM to `√(J² + D²)/J · Hückel`
(J = 1 here), and `F = U†X^⊗N U = Π_l (cos lφ·X_l + sin lφ·Y_l)`, a lit element that is not
itself a uniform string, is a colouring and carries the reflection;
[the colouring page](../../experiments/THE_PALINDROME_AS_A_COLOURING.md)'s colours on the
circle with angle lφ at site l, the gauge argument of
[the mirror-symmetry proof](../proofs/MIRROR_SYMMETRY_PROOF.md) on a named term. At
D/J = 3/4 every cos lφ and sin lφ is rational and the row is exact. On the ring the same term
is a flux, Nφ through the cycle: at D = 0.37 J the pairing breaks, at D = J (flux π at N = 4,
real hopping with one bond's sign flipped, which `X^⊗N` keeps) it pairs again (read); the
break depends on D through the flux, not on the term's presence. Of the table's Z-side rows,
the centre line needs `[H, Z^⊗N] = 0`, every term even in k_Z, which Hückel, ZZ, the Z field
and axial DM keep and the X field, the Y field and transverse DM break; the quarter-turn
about Z holds when H is U(1)-symmetric, `[H, Σ Z_l] = 0`, which is not a cube parity and holds
for the same four terms here. Nothing in this section assigns a letter to a material bath;
the rule says only what the pairing could and could not tell if one were assigned. Gate
[`simulations/carbon/cube_on_ring.py`](../../simulations/carbon/cube_on_ring.py), 53 checks:
the set rule on the parities with the separation from a per-term reading asserted on the
two mixed rows, every commuting uniform string's shift compared to 0.0 as a reflection,
every failed cell gated on the spectrum against the error model, the rotated colouring and
the chain frame's F compared to 0 exactly in rationals with a wrong angle as control, the
centre-line and U(1) rows gated per model. The frame is not called a gauge here because
[the Hückel lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) spends that word on the bipartite
sublattice sign K, and "colouring" here is F158's per-site lit string, not that sign's
two-colouring of the graph.

The sweep for this section: this page and its siblings read the three letters on the ring
and never the cube; the proton wire holds the rule's two water instances; the three
diagonals hold the letter orbit; F138's clause 2 holds the single-common-axis condition and
the non-letter direction carried by a rotation about the dephasing axis, and the F1 entry's
Breaks-for points at it; the table of moves' half-turn row holds the bisector case; the
colouring page holds the lit strings as colourings, its colours on the circle the chain
frame; the mirror-symmetry proof holds the gauge argument; the F1 entry's DM clause holds DM
alone on bipartite graphs under the alternating map, not the mixed Hückel-plus-DM ring; the
glossary, the OpenArcs registry (`whirlpool_carbon_layers`, `benzene_center_tier_upgrade`),
the Confirmations registry (no carbon row) and CAUGHT_ERRORS (its carbon entries are the
coherence-horizon ladder and labels) hold no letter reading of the ring; the typed layer's
carbon objects are clocks and anchors. Checked for adjacency: the set rule against F138's
clause 2 and the half-turn row, the same condition in three places, and the rescue there
against the rotated colouring here, the same repair; the chain frame against the colours on
the circle and the gauge argument, the same object on a named term; the DM row against F1's
DM clause, DM alone against DM beside hopping.

## Open selected-model work

- Add non-vacuous bit_b content to the selected ring and test the F112 balance.
- Compare selected `D[Z] + D[B]` jumps while retaining a stated Hamiltonian and
  initial/operator scope.
- A material question would first require a chosen degree of freedom, coupling
  convention, bath channel, and measured bath rate; only then could it compare a
  material observation to one of these selected axes.

## Anchor

- **Framework:** F112 [`LindbladBitBPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs),
  F112-X [`LindbladBitAPiBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitAPiBalance.cs),
  F112-Y [`LindbladBitBPiYBalance`](../../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiYBalance.cs),
  Klein-V₄ [`Pi2KleinV4DephaseSwapGroup`](../../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs),
  F114 [`CommutatorDConjugationSign`](../../compute/RCPsiSquared.Core/Symmetry/CommutatorDConjugationSign.cs)
- **Proofs:** [the inversion D·Π_Z·D = Π_Y](../proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md),
  [Klein-V₄ on the palindromizers](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
  [F112 cross-dephase](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)
- **Verifier:** [`simulations/carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py);
  the cube on the ring: [`simulations/carbon/cube_on_ring.py`](../../simulations/carbon/cube_on_ring.py)
- **Companion docs:** [Selected C₄/C₆ ring Liouvillians](BENZENE_LIOUVILLIAN_PALINDROME.md),
  [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md), [README](README.md)
