# Selected Ring Model and the Three Dephase Letters: Klein-V₄

**Authors:** Tom + Claude
**Status:** Klein-V₄, F112, and F114 are Tier 1 derived framework results. The
C₄/C₆ content below is a Tier-3 selected-model translation; candidate X/Y axes
are explicitly Tier 4. No material carbon degree of freedom, β-to-J convention,
bath channel or rate, `γ`, `T₂`, or Q is assigned.
**Continues:** [Selected C₄/C₆ Ring Liouvillians](BENZENE_LIOUVILLIAN_PALINDROME.md)
**Anchors:** [Klein-V₄ dephase swaps](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
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

## Klein-V₄ structure: how the letters connect (Tier 1 derived)

The three dephase letters are related on operator space by the Klein group
`{I, D, H, Q_zx}` ≅ `Z₂ × Z₂`. Every non-identity element is an involution.

- `D` is diagonal in the Pauli basis, with entry
  `(-1)^(number of Y letters)`, and obeys `D Π_Z D = Π_Y`.
- `H` swaps X and Z labels per site, leaving I and Y fixed, and obeys
  `H Π_Y H = Π_X`.
- `Q_zx = H D` obeys `Q_zx Π_Z Q_zx = Π_X`.

The elements commute and satisfy `D H Q_zx = I`. These are exact
operator-space relations. They relate selected dephasing axes; they do not make
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
- **Proofs:** [Z↔Y dephase-letter swap](../proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md),
  [Klein-V₄ dephase swaps](../proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md),
  [F112 cross-dephase](../proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md)
- **Verifier:** [`simulations/carbon_realistic_sweep.py`](../../simulations/carbon_realistic_sweep.py)
- **Companion docs:** [Selected C₄/C₆ ring Liouvillians](BENZENE_LIOUVILLIAN_PALINDROME.md),
  [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md), [README](README.md)
