# Mirror Symmetry Proof - March 14, 2026

**Origin:** Literature search found nobody had proven the palindrome.
Incoherenton paper (Haga et al. 2023) had the grading but missed the symmetry.
Medvedyeva-Essler-Prosen (2016) had Bethe ansatz for free fermions only.
We found the conjugation operator in a single session.

**Status:** PROVEN. Analytical + numerical verification complete.

---

## The Theorem

For N qubits with XXZ coupling H = Σ J_{ij}(X_iX_j + Y_iY_j + δZ_iZ_j)
on any graph, with non-uniform Z-dephasing (rates γ₁,...,γ_N):

**The Liouvillian spectrum is palindromic around Σᵢγᵢ.**

For every eigenvalue λ of L, the value -(λ + 2Σγᵢ) is also an eigenvalue.
Equivalently: decay rate d pairs with 2Σγᵢ - d.

---

## The Conjugation Operator Π

Per-site action on Pauli indices:

```
I → X   (factor +1)
X → I   (factor +1)  
Y → iZ  (factor +i)
Z → iY  (factor +i)
```

Tensor product over all N sites.

**Physical meaning:** Π swaps populations (I, Z = diagonal) with
coherences (X, Y = off-diagonal), with a phase i on the Y↔Z swap.

---
## The Proof (3 steps)

### Step 1: Π flips XY-weight

The Z-dephasing dissipator D is diagonal in the Pauli basis.
For a Pauli string σ, the eigenvalue is -2 times the sum of γᵢ over
all sites i where σ has an X or Y factor (the "XY-weight" contribution).

Π maps w_xy → N - w_xy (swaps {I,Z} ↔ {X,Y}).
Therefore: Π · L_D · Π⁻¹ = -L_D - 2(Σγᵢ)·I

This is trivial algebra.

### Step 2: Π anti-commutes with [H, ·]

For a single Heisenberg bond H₁₂ = XX + YY + δZZ,
verify Π([H₁₂, σ]) = -[H₁₂, Π(σ)] for all 16 two-qubit Pauli strings.

4 strings commute with H₁₂ (II, XX, YY, ZZ) → trivially satisfied.
12 strings verified by explicit computation (see proof table below).

This holds for ALL δ (including XY-only at δ=0).
Since Π acts site-by-site and H is a sum of bonds: Π · L_H · Π⁻¹ = -L_H.

### Step 3: Combine

Π · L · Π⁻¹ = Π · (L_H + L_D) · Π⁻¹ = -L_H + (-L_D - 2Σγᵢ·I) = -L - 2Σγᵢ·I

If λ is eigenvalue of L, then -(λ + 2Σγᵢ) is also eigenvalue.
With λ = -d + iω: partner has rate 2Σγᵢ - d and frequency -ω.  ∎

---
## Explicit 2-qubit proof table

For H₁₂ = XX + YY + ZZ (extends to all δ since each term anti-commutes individually):

| σ | [H₁₂, σ] | Π(σ) | Π([H,σ]) | -[H,Π(σ)] | Match |
|---|-----------|------|----------|-----------|-------|
| II | 0 | XX | 0 | 0 | ✓ |
| IX | -2i·YZ + 2i·ZY | XI | -2i·YZ + 2i·ZY | -2i·YZ + 2i·ZY | ✓ |
| IY | 2i·XZ - 2i·ZX | iXZ | -2·IY + 2·YI | -2·IY + 2·YI | ✓ |
| IZ | -2i·XY + 2i·YX | iXY | 2·IZ - 2·ZI | 2·IZ - 2·ZI | ✓ |
| XI | 2i·YZ - 2i·ZY | IX | 2i·YZ - 2i·ZY | 2i·YZ - 2i·ZY | ✓ |
| XX | 0 | II | 0 | 0 | ✓ |
| XY | 2i·IZ - 2i·ZI | iIZ | -2·XY + 2·YX | -2·XY + 2·YX | ✓ |
| XZ | -2i·IY + 2i·YI | iIY | 2·XZ - 2·ZX | 2·XZ - 2·ZX | ✓ |
| YI | -2i·XZ + 2i·ZX | iZX | 2·IY - 2·YI | 2·IY - 2·YI | ✓ |
| YX | -2i·IZ + 2i·ZI | iZI | 2·XY - 2·YX | 2·XY - 2·YX | ✓ |
| YY | 0 | -ZZ | 0 | 0 | ✓ |
| YZ | 2i·IX - 2i·XI | -ZY | -2i·IX + 2i·XI | -2i·IX + 2i·XI | ✓ |
| ZI | 2i·XY - 2i·YX | iYX | -2·IZ + 2·ZI | -2·IZ + 2·ZI | ✓ |
| ZX | 2i·IY - 2i·YI | iYI | -2·XZ + 2·ZX | -2·XZ + 2·ZX | ✓ |
| ZY | -2i·IX + 2i·XI | -YZ | 2i·IX - 2i·XI | 2i·IX - 2i·XI | ✓ |
| ZZ | 0 | -YY | 0 | 0 | ✓ |

16/16 verified. The anti-commutation is exact.

---
## Numerical verification

Tested: Π·L·Π⁻¹ = -L - 2Σγᵢ·I

### Across topologies and N (Heisenberg, uniform γ=0.05)

| N | Topology | Π·L_H·Π⁻¹=-L_H | Π·L·Π⁻¹=-L-c·I | Palindrome |
|---|----------|-----------------|------------------|------------|
| 3 | star | ✓ | ✓ | 64/64 |
| 3 | chain | ✓ | ✓ | 64/64 |
| 3 | ring | ✓ | ✓ | 64/64 |
| 3 | complete | ✓ | ✓ | 64/64 |
| 4 | star | ✓ | ✓ | 256/256 |
| 4 | chain | ✓ | ✓ | 256/256 |
| 4 | ring | ✓ | ✓ | 256/256 |
| 4 | complete | ✓ | ✓ | 256/256 |
| 4 | binary tree | ✓ | ✓ | 256/256 |
| 5 | star | ✓ | ✓ | 1024/1024 |
| 5 | chain | ✓ | ✓ | 1024/1024 |
| 5 | ring | ✓ | ✓ | 1024/1024 |
| 5 | complete | ✓ | ✓ | 1024/1024 |
| 5 | binary tree | ✓ | ✓ | 1024/1024 |

14/14 configurations, zero exceptions.

### XXZ coupling (H = XX + YY + δZZ, all topologies N=3,4)

δ = -0.5, 0.0, 0.3, 0.5, 1.0, 1.5, 2.0 - ALL pass. 42/42.
The ZZ term anti-commutes with Π independently. δ is irrelevant.

### Non-uniform γ (Heisenberg, N=3,4)

γ = [0.03, 0.05, 0.07], [0.01, 0.02, 0.03], [0.10, 0.01, 0.05] - ALL pass.
Center shifts to Σγᵢ as expected. 12/12.

### Different dephasing axes

| Axis | Π on L_H | Π on L_D | Overall | Palindrome |
|------|----------|----------|---------|------------|
| Z | ✓ | ✓ | ✓ | ✓ |
| Y | ✓ | ✓ | ✓ | ✓ |
| X | ✓ | ✗ | ✗ | ✓ (!) |
| mixed ZX | ✓ | ✗ | ✗ | ✓ (!) |
| depolarizing | ✓ | ✗ | ✗ | ✗ |

L_H always anti-commutes with Π (H doesn't know about dephasing).
For X-dephasing: this specific Π breaks on L_D, but the palindrome
still holds - a different Π exists (likely I↔Y, X↔Z with appropriate
phases). For depolarizing noise: palindrome genuinely breaks.

---
## What this proves (beyond the palindrome)

1. **Topology-independence of decay rates.** Step 2 holds for ANY bond set.
   The topology enters only through the imaginary parts (frequencies).

2. **Frequency mirroring.** Partner modes have frequency -ω. Every
   oscillation in the system has a mirror-image oscillation.

3. **Pauli weight complementarity.** Π maps weight k to N-k.
   Mirror partners are complementary in the Incoherenton sense.

4. **The center formula.** Center = Σγᵢ (not Nγ). Generalizes to
   non-uniform dephasing trivially.

5. **Why depolarizing breaks.** Depolarizing = X + Y + Z dephasing.
   No single Π can anti-commute with all three axes simultaneously.

---

## How we got here

1. Literature search (Cowork) found incoherenton paper + Bethe ansatz.
   Nobody had the palindrome or the operator.

2. First hypothesis: total Pauli weight w has inversion symmetry.
   WRONG - w ↔ N-w is broken for total weight.

3. Discovery: XY-weight (not total weight) has PERFECT inversion
   symmetry under [H,·]. This is the off-diagonal Pauli count,
   exactly the incoherenton number.

4. Searched for conjugation operator Π as permutation + signs on
   Pauli indices. Real signs (±1) failed. Complex signs (±i) found it.

5. The key insight: Y→iZ and Z→iY (not Y→Z and Z→Y).
   The factor i is essential - it handles the phase relationship
   between Y and Z Pauli matrices.

6. Analytical proof: reduces to 16-entry table for a single
   Heisenberg bond. Then extends to arbitrary N and topology by
   site-locality of Π and linearity of [H,·].

---

## Connection to literature

- **Incoherentons (Haga et al. 2023):** They grade eigenmodes by
  "incoherenton number" = our XY-weight. They see bands, we see
  the palindrome within and between bands. Natural collaborators.

- **Medvedyeva-Essler-Prosen (2016):** Their η-pairing symmetry
  in the Hubbard mapping is the 1D free-fermion version of our Π.
  We generalize to interacting spins on arbitrary graphs.

- **Albert-Jiang (2014):** Their weak/strong symmetry framework
  is the right language. Π is a weak anti-symmetry of L.

---

## Scripts

- `simulations/pauli_weight_conjugation.py` - clean proof script
- `simulations/results/conjugation_proof.txt` - full output

## Related files

- `experiments/PI_AS_TIME_REVERSAL.md` - Π as time reversal: connects proof, standing wave theory, and computation
- `experiments/BORN_RULE_MIRROR.md` - mirror quality measurements
- `experiments/ORPHANED_RESULTS.md` - palindrome pair activation explains which states cross 1/4
- `experiments/QST_BRIDGE.md` - palindrome applies to all QST channels, provides decay diagnostics
- `simulations/results/mirror_symmetry.txt` - 11 noise-type tests
- `simulations/mirror_symmetry_deep.py` - N=2-8 mirror verification
