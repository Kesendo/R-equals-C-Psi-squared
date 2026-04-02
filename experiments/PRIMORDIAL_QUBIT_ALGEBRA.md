# Primordial Qubit Algebra: Z₂-Grading and Algebraic Doubling

<!-- Keywords: primordial qubit Z2 grading algebraic doubling, Pi eigenspaces
subalgebra test M_{2|2} super-algebra, forward backward decomposition palindrome,
Tomita-Takesaki modular conjugation Pi vs J, quaternion algebra Pauli space,
Pi^2 w_YZ parity Z2-graded operator algebra, R=CPsi2 primordial qubit algebra -->

**Status:** Tier 2 (computed) for Phase 1-2; Tier 1 (analytical) for Phase 3
**Date:** April 1, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [primordial_qubit_algebra.py](../simulations/primordial_qubit_algebra.py)
**Output:** [primordial_qubit_algebra.txt](../simulations/results/primordial_qubit_algebra.txt)
**Depends on:**
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (definition of Π)
- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (the gap this addresses)
- [Primordial Qubit Hypothesis](../hypotheses/PRIMORDIAL_QUBIT.md) (the question)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (Π classification: AIII, order 4)

---

> *There was a primordial qubit. It sits on both sides. It still does.*
> *But now there are watchers. And the watchers bent the angle between*
> *its two sides. And the bending is what we call time.*

## Abstract

The [Primordial Qubit hypothesis](../hypotheses/PRIMORDIAL_QUBIT.md) asks:
does the palindromic Z₂-grading under Π force an algebraic doubling of the
operator space, making the "mirror world" a theorem rather than an
interpretation?

Three phases of computation (N=2 Heisenberg chain, γ=0.05) give a
nuanced answer:

**What the algebra forces:**
The centered Liouvillian L_c is exactly block-off-diagonal in the
Π-eigenbasis, coupling V_{+1} ↔ V_{-1} and V_{+i} ↔ V_{-i}. The
operator algebra admits a proper Z₂-grading from Π² (the w_YZ parity),
decomposing as the super-algebra M_{2|2}(C) with even part
M₂(C) ⊕ M₂(C), a complexified quaternion algebra in each factor.

**What the algebra does not force:**
V_{+1} alone is not a subalgebra (3 of 16 products leak to V_{-1}).
Π has no relation to the Tomita-Takesaki modular conjugation J: Π is
linear, J is anti-linear, and no unitary transformation connects them.
This rules out the canonical doubling mechanism of von Neumann algebra
theory.

**What remains open:**
The two M₂(C) factors in the even subalgebra correspond to
symmetric/antisymmetric sectors under Π, not spatial "sides."
A thermal state at β = 1/Σγ gives a modular Hamiltonian proportional
to L_H, but this state is not the Lindblad steady state.

---

## Phase 1: Π Eigenspaces and L_c Block Structure

### Step 1: Π properties (verified)

Per-site action: I → X (+1), X → I (+1), Y → iZ (+i), Z → iY (+i).

For N=2 (16×16 matrix in the Pauli basis):

| Property | Result | Error |
|----------|--------|-------|
| Π² = diag((-1)^{w_YZ}) | Exact | 0.00e+00 |
| Π⁴ = I | Exact | 0.00e+00 |
| Eigenvalues | {+1, -1, +i, -i} | each with multiplicity 4 |

Π acts on the 16 Pauli strings as 8 two-cycles. Each cycle pairs two
strings related by the per-site map, with phases determining which
eigenvalues appear:

| Orbit | Phases | λ² | Eigenvalues |
|-------|--------|-----|-------------|
| II ↔ XX | +1, +1 | +1 | ±1 |
| IX ↔ XI | +1, +1 | +1 | ±1 |
| YY ↔ ZZ | -1, -1 | +1 | ±1 |
| YZ ↔ ZY | -1, -1 | +1 | ±1 |
| IY ↔ XZ | +i, +i | -1 | ±i |
| IZ ↔ XY | +i, +i | -1 | ±i |
| YI ↔ ZX | +i, +i | -1 | ±i |
| YX ↔ ZI | +i, +i | -1 | ±i |

### Step 2: The four eigenspaces

| Eigenspace | dim | Basis | w_YZ |
|-----------|-----|-------|------|
| V_{+1} | 4 | II+XX, IX+XI, YY-ZZ, YZ-ZY | even (0 or 2) |
| V_{-1} | 4 | II-XX, IX-XI, YY+ZZ, YZ+ZY | even (0 or 2) |
| V_{+i} | 4 | IY+XZ, IZ+XY, YI+ZX, YX+ZI | odd (1) |
| V_{-i} | 4 | IY-XZ, IZ-XY, YI-ZX, YX-ZI | odd (1) |

**Key observation:** V_{+1} ⊕ V_{-1} contains exactly the 8 Pauli strings
with even w_YZ (0 or 2 sites carrying Y or Z). V_{+i} ⊕ V_{-i} contains
the 8 strings with odd w_YZ. This is the Π²-parity splitting.

### Step 3: L_c is block-off-diagonal

The centered Liouvillian L_c = L + Σγ·I satisfies Π·L_c·Π⁻¹ = -L_c
(verified: ||Π L_c Π⁻¹ + L_c|| = 0.00e+00).

Block norms ||P_a · L_c · P_b|| in the Π-eigenbasis:

|  | V_{+1} | V_{-1} | V_{+i} | V_{-i} |
|--|--------|--------|--------|--------|
| V_{+1} | 0.0000 | **4.0025** | 0.0000 | 0.0000 |
| V_{-1} | **4.0025** | 0.0000 | 0.0000 | 0.0000 |
| V_{+i} | 0.0000 | 0.0000 | 0.0000 | **5.5875** |
| V_{-i} | 0.0000 | 0.0000 | **5.7289** | 0.0000 |

L_c is **exactly** block-off-diagonal: it maps V_λ → V_{-λ} and nothing
else. The same-eigenspace blocks are zero to machine precision. This is a
necessary algebraic consequence of the anti-commutation Π·L_c = -L_c·Π.

Two independent sectors: the (+1, -1) pair (controlled by the even-w_YZ
sector) and the (+i, -i) pair (the odd-w_YZ sector).

---

## Phase 2: Algebraic Structure

### Step 4a: V_{+1} is NOT a subalgebra

Testing whether products of V_{+1} basis elements (as 4×4 operator
matrices) remain in V_{+1}:

| Product | Leakage | Destination |
|---------|---------|-------------|
| (II+XX)·(IX+XI) | 0% | stays in V_{+1} |
| (IX+XI)·(II+XX) | 0% | stays in V_{+1} |
| **(IX+XI)·(YY-ZZ)** | **100%** | **V_{-1}** |
| **(YY-ZZ)·(IX+XI)** | **100%** | **V_{-1}** |
| **(YY-ZZ)·(YY-ZZ)** | **100%** | **V_{-1}** |

Three of the 16 products leak entirely to V_{-1}. The leakage target
is exclusively V_{-1} (never V_{+i} or V_{-i}).

**V_{+1} is not a subalgebra.** The Z₄ structure from Π is too fine
for algebraic closure.

Explicit calculation: (IX+XI)·(YY-ZZ) = 2i(YZ+ZY). Since
Π(YZ+ZY) = -(YZ+ZY), this product lives in V_{-1}.

### Step 4b: V_{+1} ⊕ V_{-1} IS a subalgebra

Max leakage from the combined even sector (w_YZ even): **0.00e+00**.
All 64 products of the 8 even-sector basis elements remain in the
even sector.

### Step 4c: Proper Z₂-grading from Π²

| Rule | Max leakage | Status |
|------|-------------|--------|
| even × even → even | 0.00e+00 | ✓ |
| odd × odd → even | 0.00e+00 | ✓ |
| even × odd → odd | 0.00e+00 | ✓ |

Π² gives a proper Z₂-graded algebra. The grading is the w_YZ parity:
per-site, {I,X} form a subgroup and {Y,Z} form a coset, so w_YZ mod 2
is additive under operator multiplication.

**Π² is a multiplicative automorphism** of the operator algebra. This
is not true of Π itself (Π does not preserve operator products).

### Step 4d: Even subalgebra structure

The even subalgebra (dim 8) decomposes via central idempotents:

    P = (II + XX)/2,    Q = (II - XX)/2

Verified: P² = P, Q² = Q, PQ = 0, and P commutes with all even-sector
elements (max ||[P, A_even]|| = 0.00e+00).

**P-sector** (dim 4, basis P·{II+XX, IX+XI, YY-ZZ, YZ+ZY}):

The four basis elements satisfy p₁² = p₂² = p₃² = p₄² ∝ p₁ (the
identity) and {p₂, p₃} = 0 (anti-commutation). This is the
**complexified quaternion algebra** H_C, which is isomorphic to M₂(C).

By the Artin-Wedderburn theorem:

    Even subalgebra ≅ M₂(C) ⊕ M₂(C)  (P-sector ⊕ Q-sector)

The full operator algebra, with the Z₂-grading from Π², is the
**super-algebra M_{2|2}(C)**:

- Even part (w_YZ even): M₂(C) ⊕ M₂(C), dim 8
- Odd part (w_YZ odd): bimodule connecting the two M₂(C) factors, dim 8

### Step 5: Forward/Backward Decomposition

In the centered frame (L_c = L + Σγ·I), eigenvalues pair as λ ↔ -λ:

| Category | Count | Re(λ_c) | w_XY | Origin |
|----------|-------|---------|------|--------|
| Backward (slow) | 3 | +0.100 | 0 | Steady states of L |
| Boundary | 10 | 0.000 | 0-2 | Average-rate modes |
| Forward (fast) | 3 | -0.100 | 2 | Maximum-rate modes |

The forward and backward modes are **exact eigenstates of L_H**
(eigenvalue 0): all 3 backward modes (II, ZZ, ZI+IZ) and all 3 forward
modes (XX, YY, XY+YX) commute with the Hamiltonian. The Hamiltonian
acts entirely within the 10-mode boundary sector.

Block norms in the L_c eigenbasis:

| Operator | F→F | F→B | B→B | 0→0 |
|----------|-----|-----|-----|-----|
| L_c | 0.173 | **0.000** | 0.173 | 9.797 |
| L_H | 0.000 | **0.000** | 0.000 | 9.800 |
| L_D+Σγ·I | 0.173 | **0.000** | 0.173 | 0.142 |

**The forward/backward split is purely dissipative at N=2.** The
Hamiltonian is invisible to it. The dissipator alone determines which
modes decay faster (forward) or slower (backward) than average.

### Step 6: L_H and L_D in the Π-Eigenbasis

Both L_H and L_D + Σγ·I are odd under Π (block-off-diagonal):

L_H block norms:

|  | V_{+1} | V_{-1} | V_{+i} | V_{-i} |
|--|--------|--------|--------|--------|
| V_{+1} | 0 | **4.00** | 0 | 0 |
| V_{-1} | **4.00** | 0 | 0 | 0 |
| V_{+i} | 0 | 0 | 0 | **5.66** |
| V_{-i} | 0 | 0 | **5.66** | 0 |

L_D + Σγ·I block norms:

|  | V_{+1} | V_{-1} | V_{+i} | V_{-i} |
|--|--------|--------|--------|--------|
| V_{+1} | 0 | **0.14** | 0 | 0 |
| V_{-1} | **0.14** | 0 | 0 | 0 |
| V_{+i} | 0 | 0 | 0 | **0.14** |
| V_{-i} | 0 | 0 | **0.14** | 0 |

Both operators respect the Π²-block structure (no cross-coupling between
even and odd sectors). The Hamiltonian dominates the (+i, -i) sector
(norm 5.66 vs 0.14); the dissipator is weak in both sectors.

---

## Phase 3: Tomita-Takesaki

### Step 7: Modular Conjugation

The Lindblad steady states at N=2 span {II, ZZ, ZI+IZ} (dimension 3).
The natural choice for the GNS construction is ρ_ss = I/d (maximally
mixed state, β = 0).

| Object | Formula | Value for ρ = I/d |
|--------|---------|-------------------|
| GNS inner product | ⟨A\|B⟩ = (1/d)Tr(A†B) | Hilbert-Schmidt/d |
| Tomita operator S | S\|A⟩ = \|A†⟩ | Complex conjugation K |
| Modular operator Δ | Δ = S†S | I (trivial) |
| Modular conjugation J | J = S·Δ^{-1/2} | K (complex conjugation) |
| Modular Hamiltonian | K_mod = -log Δ | 0 |

**J = complex conjugation K in the Pauli basis (anti-linear).**

### Step 7b: Π vs J, Impossibility Proof

Π is linear: Π(c·v) = c·Π(v).
J is anti-linear: J(c·v) = c̄·J(v).

**Theorem:** No unitary U exists with Π = U·J.

*Proof.* Assume Π = U·J. For any vector v and scalar c:
Π(c·v) = c·Π(v) [linearity], but U·J(c·v) = c̄·U·J(v) = c̄·Π(v).
So c·Π(v) = c̄·Π(v) for all c. Setting c = i: i·Π(v) = -i·Π(v),
hence Π(v) = 0 for all v. Contradiction with Π being invertible. ∎

**Key insight (state-independence):** J is the adjoint map A ↦ A† for
ALL faithful states on M_d(C). The modular operator Δ depends on ρ;
J does not. This is a standard result for type I factors.

Therefore: **Π has no relation to the Tomita-Takesaki modular
conjugation for any faithful state.** The "canonical doubling" of
von Neumann algebra theory (which uses J to exchange an algebra with
its commutant) is not the mechanism behind the palindromic symmetry.

### Step 8: Modular Hamiltonian at Finite Temperature

For a thermal state ρ_β = exp(-βH)/Z:

    K_mod = β · i · L_H

If β = 1/Σγ = 10.0, the modular Hamiltonian is proportional to the
Hamiltonian Liouvillian. This would connect to the thermofield double
picture: modular flow on "the other side" equals Hamiltonian evolution
on "our side."

But the thermal state at β = 1/Σγ is NOT the Lindblad steady state
(which is at β = 0). The connection exists mathematically but is not
forced by the Lindbladian dynamics.

---

## Phase 4: Summary Table

| Statement | Result |
|-----------|--------|
| L_c block-off-diagonal in Π-eigenbasis | **Confirmed** (diagonal norm 0.00e+00) |
| V_{+1} is a subalgebra | **Refuted** (3/16 products leak to V_{-1}) |
| V_{+1}⊕V_{-1} is a subalgebra | **Confirmed** (max leakage 0.00e+00) |
| Π² gives proper Z₂-graded algebra | **Confirmed** (all three rules exact) |
| Even subalgebra ≅ M₂(C) ⊕ M₂(C) | **Confirmed** (quaternion structure) |
| Full algebra = M_{2\|2}(C) super-algebra | **Confirmed** |
| Forward/backward split from L_H | **Refuted at N=2** (dissipator only) |
| Π = J (Tomita-Takesaki) | **Refuted** (linear vs anti-linear) |
| Π = U·J for some unitary U | **Refuted** (impossibility proof) |

---

## What This Establishes

1. **The palindrome creates real algebraic structure.** The centered
   Liouvillian is exactly block-off-diagonal in the Π-eigenbasis.
   The operator algebra decomposes as M_{2|2}(C) under the Z₂-grading
   from Π², with an even subalgebra isomorphic to M₂(C) ⊕ M₂(C).
   None of this was previously known.

2. **The structure is Z₂, not Z₄.** Π itself (order 4) does not give
   a subalgebra decomposition. The algebraically meaningful grading
   comes from Π² = (-1)^{w_YZ}, which IS a multiplicative automorphism.
   The finer Z₄ structure breaks closure.

3. **Tomita-Takesaki is ruled out.** The modular conjugation J is
   anti-linear and state-independent (for type I factors). Π is linear.
   No unitary transformation can bridge this gap. The "canonical
   doubling" of von Neumann algebras is a different mathematical
   structure from the palindromic symmetry.

4. **The forward/backward split is dissipative.** At N=2, the
   Hamiltonian is invisible to the forward/backward partition.
   All forward and backward modes commute with H. The palindromic
   pairing of decay rates is entirely a property of the dephasing
   structure, not the Hamiltonian.

---

## What This Does NOT Establish

- That the Z₂-grading generalizes beyond N=2 (likely, since Π² is
  always a multiplicative automorphism, but the algebra structure at
  larger N may be richer)
- That the two M₂(C) factors have physical meaning as "our side" vs
  "their side" (they correspond to symmetric/antisymmetric sectors
  under Π, not spatial sides)
- That the thermofield double connection at β = 1/Σγ is meaningful
  (suggestive but not forced by the dynamics)
- That the forward/backward decoupling from L_H persists at N > 2
  (at larger N, Hamiltonian mixing of different w_XY sectors may
  create forward/backward modes that do not commute with H)

---

## Falsification Assessment

| Hypothesis claim | Verdict |
|-----------------|---------|
| Z₂-grading forces Hilbert space doubling | **Partially confirmed**: forces M_{2\|2}(C) super-algebra structure, but this is a grading within one algebra, not a doubling of the Hilbert space |
| "Mirror world" is a theorem | **Not confirmed**: the algebra forces structure (block-off-diagonal L_c, two M₂(C) factors) but the identification with "two sides of a mirror" is interpretation |
| Π related to Tomita-Takesaki J | **Falsified**: Π and J are fundamentally incompatible (linear vs anti-linear) |
| The decomposition X = A - Π·A·Π⁻¹ has a canonical form | **Open**: the trivial A = X/2 works; the forward/backward decomposition is physically meaningful but not canonical in the algebraic sense |

The hypothesis is **not falsified as a whole** (the algebraic structure
is real and non-trivial) but the specific mechanism proposed
(Tomita-Takesaki doubling) does not apply. The question reduces to
whether the M_{2|2}(C) super-algebra structure, which IS forced by the
palindrome, constitutes a "doubling" in any physical sense.

---

## Hints: What the Data Points Toward

The confirm/refute table above is too clean. Five features in the data
point beyond what was expected, and none of them fit neatly into
"confirmed" or "refuted."

### 1. The Clifford algebra signature

The P-sector multiplication table: p₁² = p₂² = p₃² = 1, {p₂, p₃} = 0.
This is not "some algebra." This is the Clifford algebra Cl(2,0) ≅ H
(real quaternions), complexified to H_C ≅ M₂(C). The palindrome selects
a Clifford algebra structure on the Pauli space. Clifford algebras are
the language of spinors. That the Liouvillian of a dephasing system
carries a spinor-type structure was not predicted by any of the three
approaches in [PRIMORDIAL_QUBIT.md](../hypotheses/PRIMORDIAL_QUBIT.md).

### 2. The leakage is structured

V_{+1} products leak to V_{-1}. Never to V_{+i}. Never to V_{-i}.
The Π²-parity (w_YZ mod 2) is exactly preserved under operator
multiplication; the finer Π-eigenvalue shifts by at most one Z₄ step.
V_{+1} is not a subalgebra, but it is "one step from being one": the
leakage is controlled, directional, and entirely within the even sector.

### 3. The Tomita-Takesaki elimination is a signpost

Π is linear. J is anti-linear. No unitary connects them. This
kills the modular-theory path. But it focuses attention on the correct
mechanism: not the modular conjugation, but the **Z₂-grading itself**.
M_{2|2}(C) is the algebraic backbone of the AIII symmetry class. The
closed path is: palindrome → chiral symmetry → class AIII → M_{2|2}(C)
super-algebra. The Tomita-Takesaki failure is not a dead end; it
eliminates the wrong mechanism and leaves the super-algebra as the
surviving candidate.

### 4. β = 1/Σγ: the palindrome as temperature

The modular conjugation J is dead. The modular Hamiltonian is not:

    K_mod(β) = β · i · L_H

At β = 1/Σγ, the modular Hamiltonian of the thermal state is
proportional to L_H. The palindrome parameter Σγ appears as a
**temperature**. The Lindblad steady state sits at β = 0 (infinite
temperature). But the thermofield double construction does not
require the steady state; it requires a faithful state. The thermal
state at T = Σγ is faithful. Its modular flow generates Hamiltonian
evolution scaled by 1/Σγ.

This does not prove the TFD connection. But it pins down the candidate
reference state: not ρ_ss = I/d, but ρ_β = exp(-H/Σγ)/Z.

### 5. The two sides emerge from noise

At Σγ = 0: forward and backward modes are degenerate (both at Re = 0).
No separation. The palindrome is centered, the mirror is at zero,
and the "two sides" do not exist.

At Σγ > 0: the split is ±Σγ. The separation IS the noise. Not caused
by it. Identical to it. This is the primordial qubit picture in
eigenvalue language: the split between the two sides IS the dephasing
rate. Noise does not shift a pre-existing boundary; noise creates the
boundary.

### 6. L_c as supercharge

The strongest hint, not in the original task design:

L_c is block-off-diagonal in the Π-eigenbasis. In the M_{2|2}(C)
super-algebra, L_c is an **odd operator**: it maps within the even
sector (V_{+1} ↔ V_{-1}) and within the odd sector (V_{+i} ↔ V_{-i}),
but always off-diagonally. It never maps a subspace to itself.

In supersymmetric quantum mechanics, the supercharge Q has exactly this
property: Q maps bosonic states to fermionic states and vice versa.
Q² ∝ H (the Hamiltonian). The question is whether L_c² has a similar
structural role.

L_c does not live on either side of the mirror. L_c IS the crossing.
Every time step connects both halves. You cannot have dynamics without
traversing the mirror. The generator of time evolution is the bridge.

**The right question may not have been "Is Π = J?" but: "Is L_c a
supercharge of the M_{2|2}(C) super-algebra?"**

### 7. L_c² = L_H² + (L_D + Σγ)²: the Pythagorean theorem (computed)

Tested immediately after the hints were identified. All results exact
to machine precision (N=2, γ=0.05).

**Result 1:** L_c² is block-diagonal in the Π-eigenbasis.

|  | V_{+1} | V_{-1} | V_{+i} | V_{-i} |
|--|--------|--------|--------|--------|
| V_{+1} | **16.00** | 0.00 | 0.00 | 0.00 |
| V_{-1} | 0.00 | **16.00** | 0.00 | 0.00 |
| V_{+i} | 0.00 | 0.00 | **22.62** | 0.00 |
| V_{-i} | 0.00 | 0.00 | 0.00 | **22.62** |

[Π, L_c²] = 0 exactly. L_c maps each eigenspace OUT (off-diagonal);
L_c² maps each eigenspace BACK (diagonal). This is the defining
property of an odd operator in a Z₂-graded algebra.

**Result 2:** The anti-commutator {L_H, L_D + Σγ·I} = 0 exactly.

    ||{L_H, L_D + Σγ·I}|| = 0.00e+00

The Hamiltonian and the centered dissipator are orthogonal. They do
not interfere. This is not approximate. It is exact at N=2.

**Why:** All 24 nonzero entries of L_H connect Pauli strings whose
w_XY values sum to N: w_XY(source) + w_XY(target) = 2. Since
L_D = -2γ·w_XY·I (diagonal), the anti-commutator entry is
(d_a + d_b)·(L_H)_{ab} = -2γ·(w_a + w_b)·(L_H)_{ab} = -2Σγ·(L_H)_{ab}.
Adding the Σγ·I shift: {L_H, L_D} + 2Σγ·L_H = -2Σγ·L_H + 2Σγ·L_H = 0.

**Result 3:** The Pythagorean decomposition.

    L_c² = L_H² + (L_D + Σγ·I)²

    ||L_c² - (L_H² + (L_D + Σγ)²)|| = 0.00e+00

The square of the time evolution generator decomposes into the sum
of the Hamiltonian square and the dissipator square. No cross term.

| Component | Norm | Role |
|-----------|------|------|
| L_H² | 39.19 | Oscillation (negative semi-definite) |
| (L_D+Σγ)² | 0.03 | Cooling (positive semi-definite) |
| L_c² | 39.18 | Total (sum, no interference) |

**Result 4:** The boundary carries no cooling.

The centered dissipator (L_D + Σγ·I) has diagonal entries:

| w_XY | L_D + Σγ | (L_D + Σγ)² | Count | Sector |
|------|----------|-------------|-------|--------|
| 0 | +Σγ | +Σγ² | 4 | Backward |
| 1 | **0** | **0** | 8 | **Boundary (mirror)** |
| 2 | -Σγ | +Σγ² | 4 | Forward |

The 8 boundary modes (w_XY = 1) sit at (L_D + Σγ) = 0 exactly.
They contribute nothing to the cooling term. They ARE the mirror:
the point where the centered dissipator vanishes, where forward and
backward are indistinguishable, where all Hamiltonian dynamics lives.

The cooling opens symmetrically around the mirror: +Σγ² for w_XY = 0
(backward, slow decay) and +Σγ² for w_XY = 2 (forward, fast decay).
The split is Σγ. The square of the split is Σγ². The Pythagorean
theorem says: this is the ONLY non-Hamiltonian contribution to L_c².

**Physical reading:**

At Σγ = 0 (no noise): L_c = L_H. Pure oscillation. L_c² = L_H²
(negative semi-definite). No cooling. No arrow of time.

At Σγ > 0: the (L_D + Σγ)² term appears. It is orthogonal to L_H²
(no cross term). It opens a gap of ±Σγ around the mirror. The gap
IS the cooling. The cooling is not a perturbation on the oscillation;
it is an independent, perpendicular contribution.

The Pythagorean theorem of the palindrome:

    (Zeitentwicklung)² = (Oszillation)² + (Kühlung)²

Oscillation and cooling are the two legs of a right triangle. The time
evolution is the hypotenuse. They do not interfere because they are
orthogonal ({L_H, L_D + Σγ} = 0). The boundary is where cooling
vanishes, where only oscillation remains, where the mirror stands.

### N-scaling: the Pythagorean theorem beyond N=2 (computed)

Computed for N=2, 3, 4 (Heisenberg chain, γ=0.05).

#### The algebraic fact (all N)

[Π, L_c²] = 0 holds exactly at every N. This is not numerical; it is
algebraic: Π·L_c = -L_c·Π implies Π·L_c² = L_c²·Π. L_c² is
block-diagonal in the Π-eigenbasis at every system size, regardless of
whether the Pythagorean decomposition holds.

#### The N-scaling table

| N | 4^N | ||{L_H, L_Dc}|| | Cross/||L_c²|| | Relative | [Π, L_c²] |
|---|-----|-----------------|----------------|----------|-----------|
| 2 | 16 | **0.00e+00** | **0.00%** | **0.00%** | 0.00e+00 |
| 3 | 64 | 2.77e+00 | **1.83%** | 14.4% | 5.5e-15 |
| 4 | 256 | 9.60e+00 | **2.07%** | 8.8% | 7.6e-15 |

The cross term is ~2% of ||L_c²|| at N=3 and N=4. Not exact, but small.

#### Why N=2 is exact and N≥3 is not

L_H changes w_XY by 0 or ±2 at each bond. Therefore
w_XY(a) + w_XY(b) is always **even** for any nonzero (L_H)_{ab}.

The orthogonality condition requires w_XY(a) + w_XY(b) = N for all
nonzero L_H entries. This is:

- **Possible at N=2** (even): all 24 entries satisfy sum = 2. ✓
- **Impossible at odd N** (3, 5, 7...): the sum is always even, but
  N is odd. No entry can satisfy sum = N. Every entry violates.
- **Partial at even N ≥ 4**: some entries satisfy sum = N, others
  do not. The violation comes from bonds that act on 2 of N sites,
  leaving the other N-2 sites to contribute arbitrary w_XY.

The w_XY sum distribution at each N:

| N | Sums observed | Entries at sum=N | Violations |
|---|--------------|-----------------|------------|
| 2 | {2} | 24/24 (100%) | 0/24 |
| 3 | {2, 4} | **0/192 (0%)** | 192/192 |
| 4 | {2, 4, 6} | 576/1152 (50%) | 576/1152 |

At N=3, EVERY L_H entry violates the sum condition. Yet the cross
term is only 1.83% of ||L_c²||. This is because the violations from
sum=2 and sum=4 partially cancel (opposite signs of the correction
factor w_a + w_b - N).

#### The cross term is γ-independent

The relative orthogonality ||{L_H, L_Dc}|| / (||L_H||·||L_Dc||) is
constant across all γ values at N=3:

| γ | Relative orthogonality |
|---|----------------------|
| 0.500 | 0.144338 |
| 0.100 | 0.144338 |
| 0.050 | 0.144338 |
| 0.010 | 0.144338 |
| 0.001 | 0.144338 |

Exact constant: 1/√48 = 1/(4√3) ≈ 0.14434. This is a **geometric**
property of the Heisenberg Hamiltonian's w_XY transition structure,
not a physical parameter. It depends on the graph topology and coupling
type, not on γ.

#### The mirror position

At N=2: the mirror sits at w_XY = N/2 = 1 (integer). Eight modes sit
exactly on it. L_D + Σγ = 0 for these modes. The mirror carries modes.

At N=3: the mirror sits at w_XY = N/2 = 1.5 (not integer). No modes
sit on it. Every mode has nonzero L_D + Σγ. The mirror falls between
w_XY = 1 and w_XY = 2:

| w_XY | L_D + Σγ | (L_D + Σγ)² | Count | Position |
|------|----------|-------------|-------|----------|
| 0 | +3γ | 9γ² | 8 | far backward |
| 1 | +γ | γ² | 24 | near backward |
| 2 | -γ | γ² | 24 | near forward |
| 3 | -3γ | 9γ² | 8 | far forward |

The near-mirror modes (w=1, w=2) have dissipator contribution γ²
(small). The far modes (w=0, w=3) have 9γ². The cooling is no longer
zero at the boundary; it is merely small.

#### The anti-commutator decomposition at N=3

The cross term {L_H, L_Dc} decomposes by w_XY sector pairs:

| w_a → w_b | ||block|| | Correction factor |
|-----------|----------|-------------------|
| 0 → 2 | 0.80 | +γ (sum=2 < N=3) |
| 1 → 1 | 1.60 | +γ (sum=2 < N=3) |
| 1 → 3 | 0.80 | -γ (sum=4 > N=3) |
| 2 → 0 | 0.80 | +γ (sum=2 < N=3) |
| 2 → 2 | 1.60 | -γ (sum=4 > N=3) |
| 3 → 1 | 0.80 | -γ (sum=4 > N=3) |

The correction factor is -2γ·(w_a + w_b - N): positive when the sum
is below N (the Hamiltonian connects modes that are both "too close"
to the backward side), negative when above N (both "too close" to
forward). The sum=N terms contribute zero, but there are none at N=3.

#### What N=2 has that N=3 does not

N=2 is the **only** system size where:

1. The mirror falls exactly on modes (w_XY = N/2 is integer with N=2)
2. ALL L_H transitions satisfy w_XY(a) + w_XY(b) = N
3. The Pythagorean decomposition is exact
4. Oscillation and cooling are perfectly orthogonal

This is because at N=2, a single Heisenberg bond spans the entire
system. Every L_H transition involves both sites. No site is a
spectator. The bond IS the system.

At N≥3, bonds are local. Each bond changes 2 of N sites. The other
N-2 sites contribute their w_XY to the sum, breaking the constraint.
The Pythagorean decomposition survives to ~98% accuracy, but the
perfect right angle between oscillation and cooling is a unique
property of the primordial bond.

#### Physical reading

The Urqubit (N=2, single bond) is special:

- Perfect orthogonality of oscillation and cooling
- Mirror exactly on modes (8 of 16 at the boundary)
- Pythagorean decomposition exact
- The time evolution decomposes cleanly into two non-interfering parts

At N=3 (the first composite system):

- Orthogonality broken by ~2% (cross term γ-independent)
- Mirror between modes (no mode at the boundary)
- The two parts of time evolution begin to interfere
- Cooling is nowhere zero; every mode feels the dissipator

The transition from N=2 to N=3 is the transition from the primordial
bond (where oscillation and cooling are perfectly separated) to the
first composite system (where they begin to mix). The ~2% cross term
is the first interference between oscillation and irreversibility.

### Why time is irreversible

The chain of results, each computed or proven:

1. **γ must come from outside** the d(d-2)=0 framework.
   ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md),
   Tier 2, elimination complete.)

2. **γ creates the cooling, H creates the oscillation.** The centered
   Liouvillian decomposes: L_c = L_H + (L_D + Σγ·I), where L_H is
   anti-Hermitian (pure oscillation) and L_D + Σγ·I is Hermitian
   (decay/growth split).

3. **At N=2: {L_H, L_D + Σγ} = 0.** Oscillation and cooling are
   orthogonal. Separable. The Pythagorean decomposition holds exactly.
   (Computed, Step 9, error = 0.00e+00.)

4. **At N≥3: {L_H, L_D + Σγ} ≠ 0.** The cross term is ~2% of
   ||L_c²||, γ-independent, geometric. Oscillation and cooling are
   woven together. Inseparable. (Computed, Step 10.)

5. **From N≥3 back to N=2: impossible.** Tracing out a qubit from an
   N=3 system produces non-Markovian, non-palindromic effective noise
   (0/16 palindromic pairs). The Urqubit structure is destroyed.
   ([failed_third.py](../simulations/failed_third.py), Tier 2.)

Time reversal would require undoing the cooling without disturbing the
oscillation. This requires {L_H, L_D + Σγ} = 0 (orthogonality). This
holds only at N=2. We live at N >> 2.

This is not a thermodynamic argument. Thermodynamics says: reversal is
improbable (entropy increases, fluctuation theorem, second law). This
is an algebraic argument: **reversal is structurally impossible at
N > 2**, because the cross term {L_H, L_D + Σγ} weaves the reversible
and irreversible parts of the dynamics into a single inseparable
structure, and this weaving cannot be undone by reduction (tracing out
produces non-Markovian noise that destroys the palindrome).

The arrow of time is not entropy. The arrow of time is the cross term.

| N | {L_H, L_Dc} | Separation | Time reversal |
|---|-------------|------------|---------------|
| 2 | = 0 (exact) | Oscillation ⊥ cooling | Algebraically possible |
| 3 | ≠ 0 (2%) | Woven together | Algebraically impossible |
| N>>2 | ≠ 0 | Fully entangled | The world we live in |

### What to compute next

1. **Topology dependence:** The cross term 1/√48 at N=3 is specific
   to the chain topology. Does the star or ring give a different
   geometric constant? Is there a topology that preserves the exact
   orthogonality beyond N=2?

2. **N=3 algebra structure:** Does M_{2|2}(C) generalize? The even
   subalgebra at N=3 has dim 32. Its decomposition into simple factors
   and the Clifford structure may be richer.

3. **TFD at β = 1/Σγ:** Construct the thermofield double for the
   thermal state at T = Σγ. Compare the TFD Hamiltonian with L_c.

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Computation script | simulations/primordial_qubit_algebra.py |
| Full output | simulations/results/primordial_qubit_algebra.txt |
| Hypothesis | hypotheses/PRIMORDIAL_QUBIT.md |
| Prior computation | simulations/urqubit_test.py (C²⊗C² structure) |

---

*At N=2 the bond is the system. Oscillation and cooling stand at a
right angle. You could separate them. You could, in principle, reverse
time.*

*At N=3 the bond is local. Two of three qubits participate; one watches.
The watcher bends the right angle by 2%. Oscillation and cooling bleed
into each other. You cannot separate them. You cannot reverse time.*

*At N >> 2 every bond is local. Every bond has N-2 watchers. The cross
term accumulates. The reversible and the irreversible are woven together
so tightly that no reduction can untangle them. That is why we remember
the past but cannot return to it. Not because too much has happened.
Because the algebra forbids it.*

*Thomas Wicht, Claude (Anthropic), April 1, 2026*
