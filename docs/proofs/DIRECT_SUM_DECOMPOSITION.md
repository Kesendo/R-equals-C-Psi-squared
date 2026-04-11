# Direct-Sum Decomposition of the Liouvillian

**Status:** PROVEN. Follows algebraically from the Mirror Symmetry Proof
and the Parity Selection Rule. No new numerical verification required;
it is a corollary of two established results.

**Date:** April 11, 2026

**Motivation:** The paper "A new understanding of Einstein-Rosen bridges"
(Gaztañaga, Kumar, Marto; Class. Quantum Grav. 43, 015023, 2026;
[arXiv:2512.20691](https://arxiv.org/abs/2512.20691)) proposes that physical systems admit a direct-sum
quantum theory: two Hilbert-space sectors of equal dimension, related
by a discrete transformation, carrying opposite time orientations.
The present proof shows that the Liouvillian of the Heisenberg chain
under Z-dephasing exhibits exactly this structure, with the conjugation
operator Π as the discrete transformation. The two frameworks were
developed independently.

---

## What this document proves, in plain language

The full space of quantum operators for N qubits splits into two
halves of equal size (the "even sector" and the "odd sector"). The
Liouvillian never mixes the two halves: dynamics in one half stays
in that half. This is the block-diagonal structure proven in the
[Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md).

This proof adds one observation: for systems with an **odd** number
of qubits, the mirror operator Π maps one half onto the other,
with reversed dynamics. One half decays, the other is its mirror
image. The entire system is a direct sum of two time-reversed copies,
joined by Π. Knowing the spectrum of one half gives you the other
for free.

For an **even** number of qubits, Π acts within each half. Each
half is independently palindromic: it contains its own mirror. No
sector exchange occurs.

---

## Theorem

**The Liouvillian of the Heisenberg/XXZ chain under local Z-dephasing
admits the direct-sum decomposition L = L_even ⊕ L_odd, where:**

1. **Both sectors have equal dimension:** dim(V_even) = dim(V_odd) = 2^(2N−1).

2. **For odd N, Π exchanges sectors:**
   Π: V_even → V_odd and Π: V_odd → V_even.
   The sector spectra are mutual palindromes:
   L_odd = −Π · L_even · Π⁻¹ − 2Σγ · I.

3. **For even N, Π preserves sectors:**
   Π: V_even → V_even and Π: V_odd → V_odd.
   Each sector is independently palindromic.

4. **The superselection charge** is the parity operator
   P_XY = (−1)^{n_XY}, which commutes with L and has eigenvalues ±1
   defining the two sectors.

---

## Definitions

All definitions (Pauli strings, n_XY weight, Liouvillian L = L_H + L_D)
are as in the [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md)
and the [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md).

**V_even** = span{P : n_XY(P) is even}, the "population-like" sector.
Contains the identity operator, all pure {I,Z} strings, and all
operators with an even number of X/Y factors.

**V_odd** = span{P : n_XY(P) is odd}, the "coherence-like" sector.
Contains all operators with an odd number of X/Y factors.

**P_XY** is the linear operator on the 4^N-dimensional operator space
defined by its action on each Pauli basis element:

    P_XY(σ) = (−1)^{n_XY(σ)} · σ

and extended by linearity. Equivalently, P_XY = Σ_σ (−1)^{n_XY(σ)} |σ⟩⟨σ|
in the Hilbert-Schmidt inner product, where {|σ⟩} is the orthonormal
Pauli basis. P_XY is diagonal with eigenvalues +1 (on V_even) and
−1 (on V_odd).

---

## Proof

### Step 1: Block-diagonal structure (established)

The Parity Selection Rule proves:
- L_D preserves n_XY exactly.
- L_H preserves n_XY parity (changes n_XY by 0 or ±2).
- Therefore L preserves n_XY parity: L(V_even) ⊆ V_even, L(V_odd) ⊆ V_odd.

This gives L = L_even ⊕ L_odd as a direct sum of linear maps. ∎

### Step 2: Equal dimensions

For an N-qubit system, the Pauli basis has 4^N elements. To count
dim(V_k) for each n_XY weight k:

- Choose which k of N sites carry X or Y: C(N,k) ways.
- At each of those k sites, choose X or Y: 2^k options.
- At each of the remaining N−k sites, choose I or Z: 2^(N−k) options.

Therefore:

    dim(V_k) = C(N,k) · 2^k · 2^(N−k) = C(N,k) · 2^N

Summing over even k:

    dim(V_even) = 2^N · Σ_{k even} C(N,k) = 2^N · 2^(N−1) = 2^(2N−1)

where we used the identity Σ_{k even} C(N,k) = 2^(N−1) (half the
binomial coefficients, valid for all N ≥ 1).

By the same argument (or by subtraction from 4^N):

    dim(V_odd) = 2^(2N−1)

Both sectors have exactly half the total dimension. ∎

**Concrete values:**

| N | 4^N (total) | 2^(2N−1) (per sector) |
|---|-------------|----------------------|
| 2 | 16 | 8 |
| 3 | 64 | 32 |
| 4 | 256 | 128 |
| 5 | 1,024 | 512 |
| 6 | 4,096 | 2,048 |
| 7 | 16,384 | 8,192 |
| 8 | 65,536 | 32,768 |

### Step 3: Π is invertible (Π⁴ = I)

The per-site action of Π is (from the Mirror Symmetry Proof):

    Π₁: I → X → I       (period 2 on this pair)
    Π₁: Y → iZ → i·(iY) = −Y → −iZ → −i·(iY) = Y   (period 4)

Explicitly, applying Π₁ repeatedly to each single-site Pauli:

| Apply | I | X | Y | Z |
|-------|---|---|---|---|
| Π₁    | X | I | iZ | iY |
| Π₁²   | I | X | −Y | −Z |
| Π₁³   | X | I | −iZ | −iY |
| Π₁⁴   | I | X | Y | Z |

So Π₁⁴ = id on every single-site Pauli. Since the N-qubit Π is the
tensor product Π = Π₁^{⊗N}, and each factor satisfies Π₁⁴ = id:

    Π⁴ = (Π₁^{⊗N})⁴ = (Π₁⁴)^{⊗N} = id^{⊗N} = I

Therefore Π is invertible with Π⁻¹ = Π³. ∎

**Note:** Π² ≠ I in general. From the table: Π₁²(Y) = −Y and Π₁²(Z) = −Z,
so Π² acts as (−1)^{n_YZ} on Pauli strings, where n_YZ counts sites
with Y or Z. This is consistent with the PT Symmetry Analysis, which
identifies Π as order 4, not order 2.

### Step 4: Π maps between or within sectors depending on N parity

Recall from the Mirror Symmetry Proof, Step 1: Π maps a Pauli string
with XY-weight w to one with XY-weight N − w. (This follows from the
per-site action: each {I,Z} site becomes {X,Y} and vice versa, so the
count of {X,Y} sites changes from w to N − w.)

We need to determine whether N − w has the same parity as w. Since
(N − w) − w = N − 2w, which has the same parity as N:

- **N even:** parity(N − w) = parity(w) for all w.
- **N odd:** parity(N − w) = opposite of parity(w) for all w.

Therefore:

**N odd:** Π maps every even-weight string to an odd-weight string and
vice versa. Hence **Π: V_even → V_odd** and **Π: V_odd → V_even**.

**N even:** Π maps every even-weight string to an even-weight string and
every odd-weight string to an odd-weight string. Hence
**Π: V_even → V_even** and **Π: V_odd → V_odd**. ∎

### Step 5: The palindrome restricted to sectors

From the Mirror Symmetry Proof:

    Π · L · Π⁻¹ = −L − 2Σγ · I

**Case N odd.** Let v ∈ V_even be an eigenvector of L with eigenvalue λ.
Since L preserves V_even (Step 1), v is an eigenvector of L_even.

Starting from L · v = λv, multiply both sides on the left by Π:

    Π · L · v = λ · Πv                              ... (i)

We want to rewrite the left side using the palindrome. Insert
Π⁻¹ · Π = I between L and v (this is valid because Π is invertible
by Step 3):

    Π · L · v = Π · L · (Π⁻¹ · Π) · v
              = (Π · L · Π⁻¹) · (Πv)                ... (ii)

Combining (i) and (ii):

    (Π · L · Π⁻¹) · Πv = λ · Πv

Substituting the palindrome relation Π · L · Π⁻¹ = −L − 2Σγ · I:

    (−L − 2Σγ · I) · Πv = λ · Πv
    L · Πv = (−λ − 2Σγ) · Πv

Since Πv ∈ V_odd (Step 4, N odd), this shows Πv is an eigenvector
of L_odd with eigenvalue −λ − 2Σγ.

**This map is a bijection.** Π is invertible (Step 3: Π⁴ = I implies
Π⁻¹ = Π³ exists), so the map v ↦ Πv is injective: Πv = Πw implies
Π⁻¹Πv = Π⁻¹Πw, hence v = w. Since dim(V_even) = dim(V_odd) (Step 2),
an injective linear map between equal-dimensional spaces is a bijection.
Every eigenvector of L_even maps to a unique eigenvector of L_odd.

**L_odd is entirely determined by L_even.** Since L preserves both
sectors (Step 1), the palindrome relation Π · L · Π⁻¹ = −L − 2Σγ · I
restricted to the action on V_even reads:

    Π · L_even · Π⁻¹ = −L_odd − 2Σγ · I

(The left side takes V_even → V_even → V_odd via L_even then Π.
The right side takes V_even → V_odd via Π, then acts with −L_odd − 2Σγ·I
on V_odd.) Rearranging:

    L_odd = −Π · L_even · Π⁻¹ − 2Σγ · I

Every eigenvalue of L_odd is the palindromic partner of an eigenvalue
of L_even. The two sectors carry mirror-image dynamics. ∎

**Case N even.** Π maps V_even to itself and V_odd to itself (Step 4).
Let v ∈ V_even be an eigenvector of L_even with eigenvalue λ. The same
three-line argument as in the odd case (multiply L · v = λv on the left
by Π, insert Π⁻¹ · Π between L and v, substitute the palindrome relation)
gives:

    L · Πv = (−λ − 2Σγ) · Πv

Since Πv ∈ V_even (Step 4, N even), Πv is an eigenvector of L_even with
eigenvalue −λ − 2Σγ. The palindromic partner lives in the same sector.
Restricted to V_even, the palindrome relation reads:

    Π · L_even · Π⁻¹ = −L_even − 2Σγ · I

So L_even is self-palindromic: its own spectrum is a palindrome. The
identical argument applied to a vector w ∈ V_odd gives Πw ∈ V_odd and
Π · L_odd · Π⁻¹ = −L_odd − 2Σγ · I, so L_odd is independently
self-palindromic. Each sector is a self-contained palindromic system. ∎

### Step 6: The superselection charge

The parity operator P_XY was defined in the Definitions section as a
linear operator diagonal in the Pauli basis, with eigenvalues +1 on
V_even and −1 on V_odd.

**Claim:** \[P_XY, L\] = 0.

**Proof.** Since P_XY is diagonal in the Pauli basis with eigenvalue
(−1)^{n_XY(σ)} on basis element σ, and V_even, V_odd are the eigenspaces
of P_XY with eigenvalues +1, −1 respectively, it suffices to show that
L and P_XY have the same eigenspaces, i.e., that L maps each eigenspace
of P_XY to itself.

This is exactly what Step 1 proves: L(V_even) ⊆ V_even and
L(V_odd) ⊆ V_odd. Therefore, for any ρ in V_even:

    P_XY(L(ρ)) = +1 · L(ρ) = L(ρ) = L(+1 · ρ) = L(P_XY(ρ))

and for any ρ in V_odd:

    P_XY(L(ρ)) = −1 · L(ρ) = L(−1 · ρ) = L(P_XY(ρ))

(using L(ρ) ∈ V_odd in the first equality, and linearity of L in the
last). Since V_even ⊕ V_odd spans the full space and P_XY · L = L · P_XY
on both subspaces, \[P_XY, L\] = 0 on the entire operator space. ∎

P_XY is the conserved charge of the direct-sum decomposition. In the
language of superselection: no physical process generated by L can change
the P_XY eigenvalue of a state. The two sectors are dynamically
disconnected.

---

## What this means

### The direct-sum structure for odd N

For odd N, the system is literally composed of two mirror-image halves:

- **V_even** (the "population" sector): contains the identity, all
  {I,Z}-only strings, and all operators with even numbers of coherence
  factors. Physically: the sector accessible from single-excitation
  initial states (proven in the Parity Selection Rule, Part 3).

- **V_odd** (the "coherence" sector): contains all operators with odd
  numbers of X/Y factors. Physically: the sector that single-excitation
  states can never reach.

Π is the bridge between them. It maps every mode in V_even to its
palindromic partner in V_odd, reversing the time direction (L ↦ −L − 2Σγ·I).
This is the structure of a direct-sum quantum theory with two sectors
carrying opposite arrows of time, connected by a discrete transformation.

### Self-duality for even N

For even N, the mirror is internal to each sector. The palindrome still
holds, but every mode and its partner live in the same sector. There is
no sector exchange; each half is a self-contained palindromic system.

### Connection to Gaztañaga et al.

The direct-sum quantum theory in "A new understanding of Einstein-Rosen
bridges" postulates:

1. Two Hilbert-space sectors of equal dimension → **Step 2: dim = 2^(2N−1) each**
2. A discrete transformation exchanging sectors → **Step 4 (odd N): Π: V_even ↔ V_odd**
3. Opposite time orientation in each sector → **Step 5: L_odd = −Π L_even Π⁻¹ − 2Σγ I**
4. A superselection rule preventing transitions → **Step 6: \[P_XY, L\] = 0**

All four conditions are satisfied for odd N. The two frameworks were
developed without knowledge of each other.

The even-N case has no direct analogue in Gaztañaga's framework (the
self-dual palindrome is not a two-sector theory).

### The N-parity distinction

The odd/even N split has a simple physical origin: Π flips n_XY weight
to its complement (w → N−w). When N is odd, complements have opposite
parity; when N is even, they have the same parity. The crossover between
"two-sector direct sum" and "self-dual" is controlled by a single bit.

This also explains the fragile bridge observation that odd chain lengths
are less stable than even ones (N=3 is 35× less stable than N=2; N=4
recovers 2.3× stability). The gain-loss coupling in the fragile bridge
behaves differently depending on whether Π exchanges sectors or not.

---

## Computational implication (prediction)

**Status:** This section is a prediction, not a measurement. The current
compute engine does not yet exploit the direct-sum structure. The numbers
below follow from matrix size and the cubic scaling of dense
eigendecomposition; they describe what a sector-aware engine *should*
achieve, not what has been observed. A C# engine refactored around the
analytical formulas and the per-sector decomposition is open work.

The block-diagonal structure means eigendecomposition can in principle
be performed on each sector independently. For odd N, only one sector
needs to be diagonalized (the other follows from Π).

| N | Full matrix | Per-sector matrix | RAM savings (predicted) |
|---|-------------|-------------------|-------------------------|
| 7 | 16,384² | 8,192² | 4× |
| 8 | 65,536² | 32,768² | 4× |

The 4× RAM factor follows directly from halving each matrix dimension
(2² = 4 for the dense storage). Wall-clock predictions follow from cubic
scaling of dense eigendecomposition (8× speedup per sector, halved again
for odd N where only one sector must be diagonalized). Concrete RAM and
time numbers for specific N will be added once the sector-aware engine
exists and produces verified output.

---

## Scope and limitations

**Valid for:** all conditions under which the Mirror Symmetry Proof and
Parity Selection Rule hold: Heisenberg/XXZ/XY coupling on any graph,
any site-dependent Z-dephasing profile γ_k, any N.

**The even/odd distinction applies to N specifically,** not to the
graph structure. An N=5 ring and an N=5 chain both have direct-sum
structure (odd N). An N=6 ring and an N=6 chain are both self-dual
(even N).

**Breaks for:** the same conditions that break the Parity Selection Rule
(transverse fields with odd n_XY terms, amplitude damping). See the
[Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md) scope section.

---

## Proofs this builds on

- [Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md): the palindrome Π·L·Π⁻¹ = −L − 2Σγ·I
- [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md): L preserves n_XY parity

## Related results

- [Fragile Bridge](../../hypotheses/FRAGILE_BRIDGE.md): gain-loss coupled system; odd/even N stability asymmetry
- [Standing Wave Theory](../STANDING_WAVE_THEORY.md): palindromic mode pairing
- [PT Symmetry Analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md): AIII chiral classification

## External reference

- Gaztañaga, Kumar, Marto: ["A new understanding of Einstein-Rosen bridges."](https://arxiv.org/abs/2512.20691)
  Class. Quantum Grav. 43, 015023 (2026).
  Independent development of direct-sum quantum theory with geometric
  superselection sectors.
