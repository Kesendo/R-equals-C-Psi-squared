# PROOF F103: F87 Trichotomy Z₂³ Refinement at k=3 (N=4 Empirical Anchor)

**Status:** Tier 1 derived (empirical anchor; closed-form derivation pending)
**Date:** 2026-05-24
**Anchor:** N=4, k_body=3, 294 Z₂³-homogeneous + Y-par-homogeneous Pauli pairs (pair count is N-independent at fixed k; the empirical anchor is N=4)
**Regenerate:** `simulations/f87_z2cubed_split_n4_k3.py` (~60s)

## 1. Context

The F87 trichotomy classifies a 2-letter Pauli pair as truly, soft, or hard
under a chosen dephase letter, based on whether the bilinear's M-residual
vanishes (truly), spectrally pairs (soft), or remains unpaired (hard).

F102 (`YParityIndependenceAtK3`) established that Y-parity y_par = (#Y) mod 2
is independent of the Klein (bit_a, bit_b) signature at k_body≥3 (canonical
counterexample: XYZ vs III, same Klein (0,0), different y_par). This is the
operator-algebra premise that makes a Z₂³ refinement of any k_body-stable
classification meaningful at k_body≥3.

F103 asks: does the F87 trichotomy actually refine into y_par sub-cells at
k_body=3, or is the trichotomy y_par-blind? Method: enumerate the 294
Z₂³-homogeneous + Y-par-homogeneous k_body=3 Pauli pairs (this pair count is
alphabet-only, N-independent; classification is then carried out at the
N=4 chain), classify each under Z / X / Y dephasing, bucket by (Klein ×
dephase letter × y_par × trichotomy class).

### Notation (shared across F103, F105, F106)

Pauli letters carry a two-bit Klein-Vierergruppe index
([`framework/pauli.py`](../../simulations/framework/pauli.py)):

| letter | (bit_a, bit_b) |
|--------|----------------|
| I      | (0, 0)         |
| X      | (1, 0)         |
| Z      | (0, 1)         |
| Y      | (1, 1)         |

For a Pauli string σ with letter counts n_X, n_Y, n_Z,

    bit_a(σ) = (n_X + n_Y) mod 2   (F61 axis: Π²_X = Z⊗N parity)
    bit_b(σ) = (n_Y + n_Z) mod 2   (F63 axis: Π²_Z = X⊗N parity)
    y_par(σ) = n_Y mod 2

The Klein signature (bit_a, bit_b) and y_par together form the Z₂³ axis
decomposition (one bit per axis). Π² is the squared conjugation operator;
under Z-dephasing it acts on Pauli strings as Π²·σ_α = (−1)^bit_b(α)·σ_α
(F81 Step 1; the X-deph and Y-deph variants swap to bit_a / bit_b
respectively per F108). The **diagonal Klein cell** for dephase letter D
is the cell whose (bit_a, bit_b) equals D's: Z → (0,1), X → (1,0),
Y → (1,1). F105 and F106 refer back to this notation block.

## 2. Method

The Python framework's `classify_pauli_pair(chain, terms, dephase_letter=...)`
returns one of {'truly', 'soft', 'hard'} for a Pauli pair under single-letter
dephasing at chain N=4. For k_body=3 pairs the framework dispatches to a
k-body builder (sliding-window over chain sites).

Enumeration constraints:
- Both terms have k_body=3 (no identity-padded letters)
- Both terms share the same Klein index (bit_a, bit_b)
- Both terms share the same y_par (#Y mod 2)
- Pair is unordered (deduplicated)

Result: 294 pairs partitioned across 4 Klein cells × 2 y_par values:

```
Klein    y_par=0  y_par=1  total
(0, 0)        45       21     66
(0, 1)        55       21     76
(1, 0)        55       21     76
(1, 1)        21       55     76
                             294
```

The (1,1) inversion (21 / 55 instead of 55 / 21 like (0,1) / (1,0)) reflects
that at k=3 Klein (1,1) admits 6 y_par=0 letter-triples and 10 y_par=1
letter-triples, inverted relative to (0,1) and (1,0)'s 10 / 6 split (Y is the
unique Klein-(1,1) letter, so a Klein-(1,1) k=3 string with even #Y must use
exactly one Y plus two non-Y Klein-(1,1)-summing letters; the count
mechanically inverts vs the off-axes).

## 3. Five Observed Patterns

### 3.1 Truly is y_par=0-pure

Across all 12 (Klein × dephase) cells, every truly classification has y_par=0.
Total truly classifications across the grid: 300. y_par=1 truly count: 0.

### 3.2 Hard in diagonal cells splits 42:8 with Y-inversion

Hard appears only when the Klein cell of the pair matches the Klein cell of
the dephase letter (Z → (0,1), X → (1,0), Y → (1,1)). In these 3 diagonal
cells:

```
Klein (0,1) Z-deph hard = (42, 8)   total 50
Klein (1,0) X-deph hard = (42, 8)   total 50
Klein (1,1) Y-deph hard = ( 8, 42)  total 50   ← Y-inversion
```

The Y-dephase swap reflects that Y itself carries y_par=1, so the "y_par
favored by the dephase letter" inverts.

### 3.3 Same diagonal cells contain a soft 13:13 split

In addition to the hard 42:8, the diagonal cells contain a y_par-symmetric
soft 13:13 split:

```
Klein (0,1) Z-deph soft = (13, 13)   total 26
Klein (1,0) X-deph soft = (13, 13)   total 26
Klein (1,1) Y-deph soft = (13, 13)   total 26
```

Unlike hard's 42:8 asymmetry with Y-inversion, soft in these cells is
y_par-symmetric and independent of which Klein cell is on the diagonal.

### 3.4 Mother sector (0,0) soft is y_par=1-pure

For Klein (0,0) (the Mother sector), soft cells under any dephase letter are
y_par=1-pure:

```
Z-deph: (0, 21)   X-deph: (0, 21)   Y-deph: (0, 21)
```

Zero y_par=0 soft pairs, 21 y_par=1 soft pairs per letter.

### 3.5 Off-diagonal soft cells split into Pattern B + Pattern C

The 6 off-diagonal soft cells (Klein non-mother, Klein ≠ dephase Klein) split
into two sub-patterns:

```
Pattern B (proportional to (Klein, y_par) enumeration breakdown):
Klein (0,1) Y-deph soft = (55, 21)   matches (0,1) enum split
Klein (1,1) Z-deph soft = (21, 55)   matches (1,1) enum split (inverted)
Klein (1,1) X-deph soft = (21, 55)   matches (1,1) enum split (inverted)

Pattern C (y_par=1-pure):
Klein (0,1) X-deph soft = ( 0, 21)
Klein (1,0) Z-deph soft = ( 0, 21)
Klein (1,0) Y-deph soft = ( 0, 21)
```

The exact rule connecting (pair Klein, dephase letter) to which sub-pattern
fires is observed but not yet algebraically closed.

## 4. Full Count Tables

### Truly classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)          45   0   45    45   0   45    45   0   45
(0, 1)           0   0    0    55   0   55     0   0    0
(1, 0)          55   0   55     0   0    0    55   0   55
(1, 1)           0   0    0     0   0    0     0   0    0
```

### Soft classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0  21   21     0  21   21     0  21   21
(0, 1)          13  13   26     0  21   21    55  21   76
(1, 0)           0  21   21    13  13   26     0  21   21
(1, 1)          21  55   76    21  55   76    13  13   26
```

### Hard classifications by (Klein × dephase × y_par)

```
                  Z-deph         X-deph         Y-deph
Klein           y0  y1  tot    y0  y1  tot    y0  y1  tot
(0, 0)           0   0    0     0   0    0     0   0    0
(0, 1)          42   8   50     0   0    0     0   0    0
(1, 0)           0   0    0    42   8   50     0   0    0
(1, 1)           0   0    0     0   0    0     8  42   50
```

## 5. Open Questions

1. **Closed-form derivation of 42:8.** The 50-pair hard count itself is
   already F87-derived (see [F87 entry in ANALYTICAL_FORMULAS.md](../ANALYTICAL_FORMULAS.md)
   and [PROOF_F85_KBODY_GENERALIZATION.md](PROOF_F85_KBODY_GENERALIZATION.md):
   4-cell × 3-letter trichotomy table at N=4 k=3); the 42:8 split under
   y_par would follow from Pauli-letter Klein arithmetic + Y-as-y_par-1
   weighting, but the precise derivation is open.

2. **N>4 and k>3 universality.** The (42, 8, 50) numbers are N=4 k=3
   specific. Does the structural pattern (asymmetric hard split + Y-inversion)
   carry to other (N, k)? Cheap enumeration to test: run the same script with
   N∈{5, 6}, k∈{3, 4} and observe.

3. **Pattern B vs Pattern C selection rule for off-diagonal soft.** Six cells
   partition into B (proportional) and C (y_par=1-pure); the (pair Klein,
   dephase letter) → pattern mapping is observed but the algebraic rule is
   not yet stated.

4. **Hardware confirmation.** No k≥3 F87 confirmations exist; all 5 Marrakesh
   F87 confirmations (palindrome trichotomy, π-protected XIZ/YZZY, Lebensader
   skeleton/trace, d_zero sector trichotomy, F83 Π²-class signature) are
   k=2. A k=3 QPU run targeting the diagonal-cell 42:8 prediction would be
   the natural next hardware probe.
