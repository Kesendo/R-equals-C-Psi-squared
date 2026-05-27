# PROOF: Klein-V₄ Dephase-Letter Swaps on Operator Space (Tier1Derived)

**Status:** Tier1Derived universal N via per-site factorization. Welle 12 Task 2, 2026-05-27.
**Date:** 2026-05-27 (Welle 12)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Verifier:** [`simulations/_klein_dephase_swap_explore.py`](../../simulations/_klein_dephase_swap_explore.py) (numpy bit-exact + sympy-symbolic at N = 1, 2, 3, 4).
**Builds on:** [PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md](PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md) (Welle 12 Task 1, the Z↔Y diagonal involution).
**Connects:** [PiOperator.ActOnLetter](../../compute/RCPsiSquared.Core/Symmetry/PiOperator.cs), [reflection D_PI_Z_EQUALS_PI_Y.md](../../reflections/D_PI_Z_EQUALS_PI_Y.md).

## Statement

The Klein four-group V₄ = Z₂ × Z₂ acts on the set of dephase letters {I, X, Y, Z} by Klein-product (additive on the (bit_a, bit_b) labels). Its three non-trivial elements correspond to the three letter-swaps {Z↔Y, Z↔X, Y↔X}. We claim these swaps lift to operator-space involutions on the 4^N Pauli basis, forming a faithful Klein-V₄ subgroup.

Specifically, define on the single-site Pauli basis ordered (I, X, Z, Y):

- **d_l := diag(1, 1, 1, −1)**  (the diagonal involution from PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md)
- **h := X↔Z basis-permutation matrix**, the 4×4 matrix

  ```
            col=I  col=X  col=Z  col=Y
  row=I   [   1     0      0      0   ]
  row=X   [   0     0      1      0   ]
  row=Z   [   0     1      0      0   ]
  row=Y   [   0     0      0      1   ]
  ```

  i.e. h fixes the I and Y basis vectors and swaps the X and Z basis vectors per site. (Note: h does NOT swap the physical Pauli operators X ↔ Z. It permutes their *basis indices* in the Pauli-string enumeration. See the convention note at the end.)

Per-site operators per Klein swap:

- **q_zy := d_l**  (Z↔Y swap, from PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md)
- **q_zx := h · d_l**  (Z↔X swap)
- **q_yx := h**  (Y↔X swap, pure permutation, no signs)

N-site lifts by tensor power: D := ⊗_l d_l (= q_zy^⊗N), H := ⊗_l h, Q_zx := ⊗_l q_zx = H · D, Q_yx := ⊗_l q_yx = H.

**Conjugation identities (Tier1Derived universal N):**

```
    D · Π_Z · D⁻¹    = Π_Y           (Welle 12 Task 1)
    Q_zx · Π_Z · Q_zx⁻¹ = Π_X        (this proof)
    Q_yx · Π_Y · Q_yx⁻¹ = Π_X        (this proof)
```

**Klein-V₄ group structure on operator space (Tier1Derived universal N):**

```
    D² = Q_zx² = Q_yx² = I           (all three are involutions)
    D · Q_zx · Q_yx = I              (Klein closure)
    [D, Q_zx] = [D, Q_yx] = [Q_zx, Q_yx] = 0   (pairwise commutativity; Klein V₄ is abelian)
```

Together: {I, D, Q_zx, Q_yx} is a Klein-V₄ subgroup of U(4^N), and the assignment

```
    (Klein letter v) ↦ Π_v   for v ∈ {Z, X, Y}
```

is V₄-equivariant under the conjugation action by {I, D, Q_zx, Q_yx}.

## Proof

The proof follows the same per-site factorization structure as PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md: reduce to single 4×4 identities, then lift to universal N by the mixed-product property of the Kronecker product.

### Step 1: Per-site π_local matrices

From `PiOperator.ActOnLetter` (cross-checked against `framework.symmetry.build_pi_full(N=1, dephase_letter=...)` in the verifier, max diff = 0.000e+00):

```
                            col=I  col=X  col=Z  col=Y
π_Z_local : I↔X·1, Z↔Y·i      [   0     1      0      0   ]
                                [   1     0      0      0   ]
                                [   0     0      0      i   ]
                                [   0     0      i      0   ]

π_X_local : I↔Z·1, X↔Y·−i     [   0     0      1      0   ]
                                [   0     0      0     −i   ]
                                [   1     0      0      0   ]
                                [   0    −i      0      0   ]

π_Y_local : I↔X·1, Z↔Y·−i     [   0     1      0      0   ]
                                [   1     0      0      0   ]
                                [   0     0      0     −i   ]
                                [   0     0     −i      0   ]
```

Each is a 4×4 signed permutation (column-as-input convention).

### Step 2: Tensor-product factorization

Π_Z, Π_X, Π_Y each factorize as N-fold Kronecker products of their per-site versions:

```
    Π_Z = ⊗_l π_Z_local,   Π_X = ⊗_l π_X_local,   Π_Y = ⊗_l π_Y_local.
```

This follows from `PiOperator.BuildFullUncached` constructing Π by independent per-site action plus N-site flat-index reassembly under the `a + 2·b` packing convention (Step 1 of PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md).

Hence by the mixed-product property of the Kronecker product, the N-site conjugation identities

```
    (Q_zx) · Π_Z · (Q_zx)⁻¹ = Π_X,       (Q_yx) · Π_Y · (Q_yx)⁻¹ = Π_X
```

reduce to per-site identities

```
    (q_zx) · π_Z_local · (q_zx)⁻¹ = π_X_local                              (Z↔X)
    (q_yx) · π_Y_local · (q_yx)⁻¹ = π_X_local                              (Y↔X)
```

It suffices to verify these on a single 4×4 block.

### Step 3: Per-site Z↔X identity (q_zx · π_Z · q_zx⁻¹ = π_X)

With q_zx := h · d_l:

Step 3a: Compute q_zx · π_Z_local. Left-multiplication by d_l negates row Y of π_Z_local (entries with row index 3). Then left-multiplication by h swaps rows X and Z (indices 1, 2) of the result:

```
                                       col=I  col=X  col=Z  col=Y
d_l · π_Z_local            row=I   [   0     1      0      0   ]
                            row=X   [   1     0      0      0   ]
                            row=Z   [   0     0      0      i   ]
                            row=Y   [   0     0     −i      0   ]

h · d_l · π_Z_local        row=I   [   0     1      0      0   ]   (unchanged)
                            row=X   [   0     0      0      i   ]   (was row Z)
                            row=Z   [   1     0      0      0   ]   (was row X)
                            row=Y   [   0     0     −i      0   ]   (unchanged)
```

Step 3b: q_zx⁻¹ = (h · d_l)⁻¹ = d_l⁻¹ · h⁻¹ = d_l · h (since d_l² = h² = I). Right-multiplication by d_l negates col Y. Then right-multiplication by h swaps cols X and Z:

```
q_zx · π_Z_local · d_l    col=I  col=X  col=Z  col=Y
   row=I   [   0     1      0      0   ]   (col Y negated, still 0)
   row=X   [   0     0      0     −i   ]
   row=Z   [   1     0      0      0   ]
   row=Y   [   0     0     −i      0   ]

q_zx · π_Z_local · d_l · h   col=I  col=X  col=Z  col=Y
                            (swap cols 1 and 2)
   row=I   [   0     0      1      0   ]
   row=X   [   0     0      0     −i   ]
   row=Z   [   1     0      0      0   ]
   row=Y   [   0    −i      0      0   ]
```

This equals π_X_local from Step 1. ∎ (Z↔X per-site)

### Step 4: Per-site Y↔X identity (q_yx · π_Y · q_yx⁻¹ = π_X)

With q_yx := h:

Step 4a: h · π_Y_local. Left-multiplication by h swaps rows X and Z of π_Y_local:

```
                            col=I  col=X  col=Z  col=Y
h · π_Y_local      row=I  [   0     1      0      0   ]   (unchanged)
                    row=X  [   0     0      0     −i   ]   (was row Z)
                    row=Z  [   1     0      0      0   ]   (was row X)
                    row=Y  [   0     0     −i      0   ]   (unchanged)
```

Step 4b: q_yx⁻¹ = h⁻¹ = h. Right-multiplication by h swaps cols X and Z:

```
h · π_Y_local · h  col=I  col=X  col=Z  col=Y
                  (swap cols 1 and 2)
   row=I   [   0     0      1      0   ]
   row=X   [   0     0      0     −i   ]
   row=Z   [   1     0      0      0   ]
   row=Y   [   0    −i      0      0   ]
```

This equals π_X_local from Step 1. ∎ (Y↔X per-site)

### Step 5: Klein-V₄ group structure (per-site)

Per-site involutivity:

- d_l² = diag(1, 1, 1, −1)² = diag(1, 1, 1, 1) = I.
- h² = I (transposition swap is self-inverse).
- q_zx² = (h · d_l) · (h · d_l). Compute: d_l · h has the form
  ```
  d_l · h = diag(1,1,1,−1) · h, which scales the rows of h by (1,1,1,−1):
  row=I [ 1 0 0 0 ],  row=X [ 0 0 1 0 ],  row=Z [ 0 1 0 0 ],  row=Y [ 0 0 0 −1 ].
  ```
  And h · d_l has the form
  ```
  h · d_l = h · diag(1,1,1,−1), which scales the cols of h by (1,1,1,−1):
  row=I [ 1 0 0 0 ],  row=X [ 0 0 1 0 ],  row=Z [ 0 1 0 0 ],  row=Y [ 0 0 0 −1 ].
  ```
  So h · d_l = d_l · h (they commute), and q_zx² = (h · d_l)² = (h · d_l) · (h · d_l) = h · (d_l · h) · d_l = h · (h · d_l) · d_l = h² · d_l² = I · I = I.

So all three per-site operators d_l, h, q_zx are involutions of order 2.

Klein-V₄ closure (per-site): need d_l · q_zx · q_yx = I.

```
d_l · q_zx · q_yx = d_l · (h · d_l) · h = d_l · h · d_l · h.
```

Using d_l · h = h · d_l from above, this becomes h · d_l · d_l · h = h · I · h = h² = I. ∎

Pairwise commutativity (per-site):

- [d_l, h] = 0: shown above (d_l · h = h · d_l).
- [d_l, q_zx] = d_l · (h · d_l) − (h · d_l) · d_l = (d_l · h) · d_l − h · d_l² = h · d_l² − h · I = h − h = 0.
- [q_zx, q_yx] = (h · d_l) · h − h · (h · d_l) = h · d_l · h − h² · d_l = h · (h · d_l) − I · d_l = h² · d_l − d_l = d_l − d_l = 0.

So {I, d_l, q_zx, q_yx} is a Klein-V₄ subgroup of U(4) per-site, abelian, all elements order 2.

### Step 6: Lift to universal N by tensor power

By the mixed-product property of the Kronecker product, every per-site identity in Steps 3, 4, 5 lifts to the N-site tensor power:

- D² = (⊗ d_l)² = ⊗ d_l² = ⊗ I = I (and similarly Q_zx² = Q_yx² = I).
- D · Q_zx · Q_yx = (⊗ d_l) · (⊗ q_zx) · (⊗ q_yx) = ⊗ (d_l · q_zx · q_yx) = ⊗ I = I.
- All commutators [D, Q_zx], [D, Q_yx], [Q_zx, Q_yx] vanish because their per-site factors all vanish (the per-site versions commute, so their tensor products also commute).
- The conjugation identities Q_zx · Π_Z · Q_zx⁻¹ = Π_X and Q_yx · Π_Y · Q_yx⁻¹ = Π_X follow from the per-site Steps 3, 4 by mixed-product.

Hence {I, D, Q_zx, Q_yx} is a Klein-V₄ subgroup of U(4^N) for every N, and the dephase-letter assignment v ↦ Π_v (v ∈ {Z, X, Y}) is V₄-equivariant under this group. ∎

## Verification

**Numerical (numpy double precision):** [`simulations/_klein_dephase_swap_explore.py`](../../simulations/_klein_dephase_swap_explore.py) computes both sides of every identity at N = 1, 2, 3, 4 and reports residual = 0.000e+00 in every case:

- per-site π matches `build_pi_full(N=1, ...)`: max diff = 0.000e+00
- per-site `q_zx · π_Z · q_zx⁻¹ = π_X`: residual = 0
- per-site `q_yx · π_Y · q_yx⁻¹ = π_X`: residual = 0
- N-site `Q · Π_Z · Q⁻¹ = Π_X`: residual = 0 at N = 1, 2, 3, 4
- N-site `Q' · Π_Y · Q'⁻¹ = Π_X`: residual = 0 at N = 1, 2, 3, 4
- N-site Klein closure `D · Q_zx · Q_yx = I`: residual = 0 at N = 1, 2, 3, 4
- N-site pairwise commutators: all = 0 at N = 1, 2, 3, 4

**Symbolic (sympy exact rationals + I):** the same script's Step 2 verifies the Z↔X per-site identity `q_zx · π_Z · q_zx⁻¹ = π_X` symbolically (entry-wise zero residual matrix in sympy). Combined with the tensor-product lift in Step 6, this closes the universal-N case in finite symbolic computation.

The brute-force per-site search (Steps 1, 5 of the verifier) enumerated all 16 × 16 = 256 sign-vector combinations S_L · h · S_R and found 64 valid Z↔X solutions and 64 valid Y↔X solutions. The canonical involutions (q_zx = h · d_l with order 2, q_yx = h with order 2) were selected as the simplest representatives; the other matches are order-4 cyclic variants that conjugate Π_Z to Π_X correctly but do not form an order-2 group.

## Implications

- **F1 palindrome family is fully Klein-V₄ equivariant.** Any F1 diagnostic (residual norm, Frobenius inner product, spectrum, Master Lemma identity) computed under one dephase letter automatically transfers to the other two via unitary conjugation by the appropriate element of {I, D, Q_zx, Q_yx}. The Welle 11 F112 universal-N closure under Z-dephasing implies the same under X- and Y-dephasing without re-proof.
- **The controller's order-4 conjecture is falsified.** Pre-dispatch analysis suggested Z↔X would be order 4 (cyclic) because Π_Z and Π_X act on different Klein axes (bit_a vs bit_b). The actual structure is order 2: the basis permutation h plus the sign-flip diagonal d_l together intertwine the two axes, and the composite (h · d_l)² = h² · d_l² = I (since h and d_l commute). This is the structural reason the Klein-V₄ lift works.
- **Composability of Klein swaps.** The closure D · Q_zx · Q_yx = I means the three swaps are not independent: any two determine the third. Operationally, the Z↔Y diagonal D (the simplest swap, no permutation) combined with the Y↔X pure permutation H = Q_yx gives the Z↔X swap as their product (D · H = H · D = Q_zx since D and H commute).
- **F108 cross-dephasing.** F108 Parts 1, 2, 3 (Pi5bilinear under Z-, X-, Y-dephasing respectively) are conjectured equivalent by this Klein-V₄ equivariance; an F108-style proof under Z-dephasing implies the X- and Y-dephasing versions automatically. This sharpens the F108 picture: the three parts are not independent constructions but Klein-V₄ images of a single structural identity.

## Convention note: the basis permutation h

h does NOT swap the physical Pauli operators σ_X ↔ σ_Z. Rather, h is a permutation of the **basis indices** in the Pauli-string enumeration `(I, X, Z, Y)` at positions (0, 1, 2, 3). At the level of the Pauli-basis coordinate representation (a 4^N-dim vector or 4^N × 4^N matrix indexed by Pauli-string flat indices), h acts as a column-swap of indices 1 and 2 (X and Z). Conjugation by h on a Pauli-basis matrix M permutes the rows and columns of M corresponding to letter-Z and letter-X entries.

The reason h works for the Z↔X swap on Π is structural: Π_Z permutes the I↔X and Z↔Y Pauli-string orbits (bit_a flip), while Π_X permutes the I↔Z and X↔Y orbits (bit_b flip). The X↔Z basis-index permutation maps the bit_a orbit structure into the bit_b orbit structure exactly, because the Klein V₄ on (bit_a, bit_b) labels treats X = (1, 0) and Z = (0, 1) as the two single-axis generators: swapping them is precisely the Klein-V₄ involution (a, b) ↔ (b, a) that maps the bit_a-axis to the bit_b-axis.

This is the conceptual content of the lift: the Klein-V₄ acting on dephase letters via Klein-product is the same Klein-V₄ acting on basis indices via index-permutation, intertwined through the Pauli-string enumeration. The diagonal sign d_l handles the phase difference between Π_Z (phase i) and Π_Y (phase −i) on the Z↔Y orbit. Together they realize the full Klein-V₄ on operator space.
