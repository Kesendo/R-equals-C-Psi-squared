# Proof of F102: Y-Parity Term-Level Z₂ Independence at k_body≥3

**Status:** Tier 1 derived, mechanically verified (9 k_body=2 bilinears + k_body=3 counterexample).
**Date:** 2026-05-24
**Anchors:** `docs/ANALYTICAL_FORMULAS.md` F102; `compute/RCPsiSquared.Core/Symmetry/YParityIndependenceAtK3.cs`.

## Abstract

The Klein signature (bit_a, bit_b) classifies every Pauli letter at single-site level: I = (0,0), X = (1,0), Z = (0,1), Y = (1,1). For a Pauli string of more than one letter, the bits add componentwise modulo 2. At first glance, Y-parity (the number of Y letters mod 2) looks like it should just be one of those bits, or a fixed function of them. F102 asks: is it?

The answer depends on body count. At k_body = 2 (Pauli bilinears like XX, XY, YZ), Y-parity is exactly bit_a XOR bit_b. The two are locked together, and asking about Y-parity adds nothing beyond what Klein already knows. At k_body = 3 the lock breaks. The same Klein signature can carry either Y-parity, and the two become independent classifiers. By k_body = 4 the lock returns; the general pattern is: even k_body, they agree; odd k_body, they disagree by exactly the parity of k_body.

The canonical witness is the comparison between the all-identity string III at k_body = 0 and the all-different string XYZ at k_body = 3. Both have Klein signature (0, 0), giving Y-parity zero by Klein arithmetic, but XYZ has actual Y-parity 1 because it contains exactly one Y. The mismatch is the sign that Klein and Y-parity are different axes once we leave the body-count-even sector.

The diagnostic upshot is that the polarity cube of Pauli classification has three axes once higher body counts come into play: bit_a, bit_b, AND y_par. The third axis is what the F87Z₂³ refinement family (F103, F105, F106) explores empirically, finding which trichotomy classes split further along the y_par axis and which stay y_par-blind. F102 is the precondition: y_par is a real third axis, not a derived quantity, once the body count goes up.

## Statement

For any Pauli string σ = ⊗_l σ_α_l on N qubits, define the term-level Y-parity:

    y_par(σ) = (Σ_l [σ_α_l = Y]) mod 2 = n_Y(σ) mod 2

Define the Klein signature (bit_a, bit_b) per the F61 / F63 conventions:

    bit_a(σ) = (n_X + n_Y) mod 2     (F61 axis: Π²_X = Z⊗N parity)
    bit_b(σ) = (n_Y + n_Z) mod 2     (F63 axis: Π²_Z = X⊗N parity)

**Claim:** At k_body=2 (n_X + n_Y + n_Z = 2), y_par = bit_a XOR bit_b for every Pauli bilinear. The identity equivalently holds iff k_body is even. At k_body=3 (and odd k_body in general) the identity strictly fails by y_par = (k_body + (bit_a XOR bit_b)) mod 2; in particular Y-parity is independent of the Klein signature (bit_a, bit_b) once Pauli strings of different k_body values are admitted at the same Klein signature.

## Per-letter table

| Letter | n_X | n_Y | n_Z | bit_a | bit_b | y_par |
|--------|-----|-----|-----|-------|-------|-------|
| I      | 0   | 0   | 0   | 0     | 0     | 0     |
| X      | 1   | 0   | 0   | 1     | 0     | 0     |
| Y      | 0   | 1   | 0   | 1     | 1     | 1     |
| Z      | 0   | 0   | 1   | 0     | 1     | 0     |

## Algebraic derivation

```
bit_a XOR bit_b = (n_X + n_Y) + (n_Y + n_Z) mod 2
                = n_X + 2·n_Y + n_Z mod 2
                = (n_X + n_Z) mod 2
```

So bit_a XOR bit_b counts X's plus Z's, modulo 2. Independently:

```
y_par = n_Y mod 2
```

The identity y_par = bit_a XOR bit_b becomes:

```
n_Y mod 2 = (n_X + n_Z) mod 2
```

equivalently:

```
(n_X + n_Y + n_Z) mod 2 = 0
```

i.e., the total k_body weight is even.

## k_body=2 verification (9 bilinears)

The k_body counts non-identity letters: k_body = #(non-I letters). At k_body=2, all 9 Pauli bilinears (n_X + n_Y + n_Z = 2):

- XX: n_X=2, n_Y=0, n_Z=0. bit_a=0, bit_b=0. y_par=0. bit_a XOR bit_b = 0. Match.
- XY: n_X=1, n_Y=1, n_Z=0. bit_a=0, bit_b=1. y_par=1. bit_a XOR bit_b = 1. Match.
- XZ: n_X=1, n_Y=0, n_Z=1. bit_a=1, bit_b=1. y_par=0. bit_a XOR bit_b = 0. Match.
- YX: same as XY. Match.
- YY: n_X=0, n_Y=2, n_Z=0. bit_a=0, bit_b=0. y_par=0. bit_a XOR bit_b = 0. Match.
- YZ: n_X=0, n_Y=1, n_Z=1. bit_a=1, bit_b=0. y_par=1. bit_a XOR bit_b = 1. Match.
- ZX: same as XZ. Match.
- ZY: same as YZ. Match.
- ZZ: n_X=0, n_Y=0, n_Z=2. bit_a=0, bit_b=0. y_par=0. bit_a XOR bit_b = 0. Match.

All 9 k_body=2 bilinears satisfy y_par = bit_a XOR bit_b. The k_body=0 case II (any identity-only tensor) also satisfies trivially (0 = 0).

## k_body parity rule

For general k_body, the identity becomes:

```
y_par = (k_body + (bit_a XOR bit_b)) mod 2
```

Reason: y_par = n_Y mod 2, bit_a XOR bit_b = (n_X + n_Z) mod 2, k_body = n_X + n_Y + n_Z. So n_Y + (n_X + n_Z) = k_body and y_par + (bit_a XOR bit_b) = k_body mod 2.

- **Even k_body** (0, 2, 4, ...): y_par = bit_a XOR bit_b. Identity holds.
- **Odd k_body** (1, 3, 5, ...): y_par = 1 + (bit_a XOR bit_b) mod 2. Identity strictly fails (off by exactly 1).

In particular the 6 k_body=1 cases at 2 tensor positions (IX, IY, IZ, XI, YI, ZI) all violate the identity (each has Klein signature equal to the non-identity letter's row in the per-letter table, k_body=1, hence y_par − (bit_a XOR bit_b) = 1):

- IX: bit_a=1, bit_b=0, y_par=0. 0 ≠ 1. Fail.
- IY: bit_a=1, bit_b=1, y_par=1. 1 ≠ 0. Fail.
- IZ: bit_a=0, bit_b=1, y_par=0. 0 ≠ 1. Fail.
- XI: bit_a=1, bit_b=0, y_par=0. 0 ≠ 1. Fail.
- YI: bit_a=1, bit_b=1, y_par=1. 1 ≠ 0. Fail.
- ZI: bit_a=0, bit_b=1, y_par=0. 0 ≠ 1. Fail.

## k_body=3 counterexample

At k_body=3 we want two Pauli strings with identical Klein signature (bit_a, bit_b) but different y_par. Per the parity rule above, every k_body=3 string has y_par = 1 + (bit_a XOR bit_b) mod 2, i.e. y_par is determined by Klein **within** k_body=3. So independence does not surface among k_body=3 strings alone.

The k_body=3 independence surfaces by mixing k_body values at the same Klein signature. The canonical example is:

- σ = XYZ at k_body=3: n_X=1, n_Y=1, n_Z=1. bit_a=0, bit_b=0, y_par=1. Klein (0, 0).
- σ' = II...I (identity, k_body=0): n_X=0, n_Y=0, n_Z=0. bit_a=0, bit_b=0, y_par=0. Klein (0, 0).

Both at the same N (≥3) tensor positions; both Klein (0, 0); different y_par.

**Search for k_body=3 Klein (0, 0).** We need (n_X + n_Y) mod 2 = 0 AND (n_Y + n_Z) mod 2 = 0 with n_X + n_Y + n_Z = 3. The first two give n_X ≡ n_Y mod 2 and n_Y ≡ n_Z mod 2, so all three have the same parity. With sum 3 (odd), all three must be odd, so (n_X, n_Y, n_Z) = (1, 1, 1). The k_body=3 strings with Klein (0, 0) are therefore the 6 letter permutations of XYZ (XYZ, XZY, YXZ, YZX, ZXY, ZYX = S_3 on {X, Y, Z}), all with y_par=1.

This is exactly the structural independence at k_body≥3: the Klein signature does not pin y_par once strings of different k_body values are admitted. The witness in `compute/RCPsiSquared.Core.Tests/Pauli/PauliHamiltonianKleinHelpersTests.cs` (XYZ_AtK3_IsKleinHomogeneousButZ2HomogeneityRefinesViaYParity) demonstrates this in the `PauliHamiltonian.IsZ2Homogeneous` refinement of `IsKleinHomogeneous` at k_body=3.

## Summary

The k_body=2 collapse holds at every even k_body (trivially at k_body=0, by the per-letter table at k_body=2, and by induction at higher even k_body). At odd k_body the identity strictly fails by y_par = (k_body + (bit_a XOR bit_b)) mod 2. The structural Z₂³ Y-parity axis is independent in the sense that the Klein signature does not pin y_par across the full Pauli lattice once k_body varies; in particular mixing k_body=0 with k_body=3 at the same Klein signature yields different y_par.

This is the term-level Z₂ classifier whose F-number is F102, registered as the YParity slot of `PolarityCubeMap`.
