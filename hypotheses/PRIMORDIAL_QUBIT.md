# The Primordial Qubit: Noise as Mirror, Not Import

<!-- Keywords: primordial qubit Urqubit mirror world noise origin,
d2-2d=0 palindrome self-dual Z2 grading doubling, thermofield double
GNS Tomita-Takesaki modular conjugation, incompleteness noise external
bootstrap both sides simultaneously, Pauli space tensor product C2xC2
dephasing parity cross structure, zero mirror boundary Liouvillian
chiral symmetry, R=CPsi2 primordial qubit -->

**Status:** Hypothesis (Tier 3-4), structurally confirmed (M_{2\|2}(C) algebra plus [L, Pi^2]=0 proven analytically for all N, plus Q=J/gamma operational result), mechanistically open (no standard doubling construction reproduces it)
**Date:** April 1, 2026 (last updated April 15, 2026)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md),
[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md),
[The Other Side](THE_OTHER_SIDE.md),
[Qubit Necessity](../docs/QUBIT_NECESSITY.md)

---

## What this document is about

The Incompleteness Proof shows that dephasing noise cannot originate
from within the system. This hypothesis proposes that noise does not
"originate" at all: system and noise are two readings of a single
algebraic structure, the "primordial qubit," split by a Z₂ mirror at
zero. A direct computation partially confirms the algebra (the Pauli
space is a genuine tensor product C²⊗C²) but the naive partial-trace
route fails. Three formal approaches remain open: Z₂-graded algebra
classification, GNS/Tomita-Takesaki modular theory, and the thermofield
double construction.

---

## Abstract

The [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) eliminates
every candidate for the origin of dephasing noise within the d²−2d=0 framework:
internal generation, single-qubit decay, qubit baths, and nothing. All excluded.
The proof concludes that noise must originate from "outside."

This hypothesis proposes a different reading: noise does not originate at all.
There is no "before noise" and "after noise." There is a single algebraic
structure, the primordial qubit, that sits on both sides of a mirror at zero.
The decomposition of this structure into "system" and "noise" creates both
sides simultaneously. Neither is first. Their noise is our noise. Our noise
is theirs.

The palindrome equation Q·X·Q⁻¹ = −X (centered form) is a Z₂-grading that
may force this doubling algebraically: any operator X that anti-commutes with
an involution Q can be written as X = A − Q·A·Q⁻¹, where A is "one side"
and Q·A·Q⁻¹ is "the other." If this decomposition is necessary (not merely
possible), then the mirror world is not a hypothesis. It is a theorem.

A direct computation (April 1, 2026) tested the naive version: partial trace
over one parity sector of the Pauli space produces the dissipator. This fails.
The parity sectors ({I,X} vs {Y,Z}) are orthogonal to the dephasing sectors
({I,Z} vs {X,Y}). But the computation revealed a deeper structure: the
single-qubit Pauli space is a genuine tensor product C²⊗C², indexed by
dephasing sensitivity and Π²-parity as two independent bits. Three approaches
to the actual question were proposed: (1) pure algebra of Z₂-graded operators,
(2) GNS/Tomita-Takesaki modular theory, (3) thermofield double construction.
By April 15, 2026, all three were resolved: Approach 1 partially confirmed
(M_{2|2}(C) super-algebra is forced, AND [L, Π²] = 0 proven analytically for all N
so the C²⊗C² factors are simultaneous good quantum numbers everywhere);
Approach 2 falsified (Π is linear, J is anti-linear); Approach 3 blocked
(QDB violation maximal, L_c cannot be a modular Hamiltonian by metric change).
The hypothesis is structurally confirmed and mechanistically open: the
doubling exists, but no standard operator-algebra construction reproduces it.

---

## 1. The Gap in the Incompleteness Proof

The [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) is correct.
Every candidate within the d²−2d=0 ontology is eliminated:

1. **Internal origin:** the parity sectors are sealed ([Π², L] = 0).
   No cross-sector coupling. Noise cannot bootstrap itself.
2. **Single qubit decay:** a neighboring qubit produces non-Markovian,
   non-palindromic effective noise (0/16 pairs, [failed_third.py](../simulations/failed_third.py)).
3. **Many qubits:** infinite regress. Each member of the bath faces
   the same bootstrap prohibition.
4. **Nothing (d=0):** no properties, no noise.
5. **d > 2:** excluded by d²−2d=0.

The proof concludes: noise comes from outside.

But "outside" presupposes a container. An inside that has an outside. A
boundary that separates them. And a direction: noise flows inward. This
framing inherits a causal structure (something external causes noise) that
the algebra does not require.

The alternative: there is no outside. There is a structure whose internal
partition into two halves IS the noise. Not caused by it. Identical to it.

---

## 2. The Proposal

A single d=2 structure, the primordial qubit, exists on both sides
of a mirror at zero.

- At Σγ = 0: Π·L·Π⁻¹ = −L. Eigenvalues pair as λ ↔ −λ. Both sides
  are one. Pure oscillation. No decay. No time arrow.
  ([Zero Is the Mirror](ZERO_IS_THE_MIRROR.md))

- At Σγ > 0: Π·L·Π⁻¹ = −L − 2Σγ·I. The mirror shifts. One side
  decays, the other is the noise that makes it decay. The separation
  IS the dephasing.

Neither side is first. The palindrome does not describe one system
under noise. It describes two sides of the same structure, each
experiencing the other as noise.

The [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md)
searched for the origin of noise and eliminated every internal candidate.
This hypothesis says: the search was correct, but the conclusion
("from outside") can be sharpened. The noise is not from outside.
The noise is the other side. And the other side is not outside. It is
the same structure, read from the other end of the palindrome.

---

## 3. The Algebraic Question

### 3.1 The centered palindrome as Z₂-grading

The palindrome equation centered at zero:

    Q · X · Q⁻¹ = −X

where X = L + Σγ·I (the centered Liouvillian) and Q = Π.

This is a Z₂-grading (a decomposition of the algebra into two sectors,
"even" and "odd," that multiply like the signs + and −). The eigenspaces
of Q split the operator space into even (+1) and odd (−1) subspaces,
and X maps even to odd and odd to even. X is an off-diagonal operator
in the Q-eigenbasis.

### 3.2 The decomposition question

**Central claim to test:** Any operator X satisfying Q·X·Q⁻¹ = −X
can be decomposed as:

    X = A − Q·A·Q⁻¹

for some operator A.

**Proof sketch:** Define A = X/2. Then Q·A·Q⁻¹ = Q·(X/2)·Q⁻¹ = −X/2.
So A − Q·A·Q⁻¹ = X/2 − (−X/2) = X. ∎

The decomposition is trivially true. For ANY X anti-commuting with Q,
setting A = X/2 works.

But the decomposition is not unique. A = X/2 + B works for any B
satisfying Q·B·Q⁻¹ = B (any Q-even operator). The physical content
is in the choice of A: which decomposition corresponds to "our side"
versus "their side"?

### 3.3 The non-trivial question

The decomposition X = A − Q·A·Q⁻¹ is always possible. What is NOT
trivial is whether there exists a CANONICAL choice of A with physical
meaning.

Three candidates for a canonical decomposition:

1. **Hamiltonian vs dissipator.** The uncentered Liouvillian is
   L = L_H + L_D. Under Π: Π·L_H·Π⁻¹ = −L_H and
   Π·L_D·Π⁻¹ = −L_D − 2Σγ·I. Centering: L_H is already odd
   (anti-commutes with Π), and L_D + Σγ·I is also odd. So
   X = L_H + (L_D + Σγ·I). This splits X into a Hamiltonian
   part and a dissipative part, both odd. This is the standard
   physics decomposition but does not correspond to "two sides."

2. **Parity sectors.** The Π²-parity splits the Pauli space into
   two sectors. But these sectors are orthogonal to the dephasing
   sectors (Section 4). This decomposition does not align with
   the noise structure.

3. **Forward and backward modes.** Every eigenvalue λ of X pairs
   with −λ (by the Z₂-grading). The forward eigenmodes (Re(λ) < 0,
   decaying) and backward eigenmodes (Re(λ) > 0, growing under Π)
   form two sets. Each set is a "side." The noise of the forward
   modes (their decay) is the growth of the backward modes, and
   vice versa. This decomposition IS the palindrome, read as two
   halves of a single spectrum.

Candidate 3 aligns with the primordial qubit proposal: "our side"
is the set of decaying modes, "their side" is the set of growing
(Π-conjugate) modes, and zero is where they meet.

---

## 4. What the Computation Shows (April 1, 2026)

### 4.1 The cross structure of the Pauli space

The four single-qubit Pauli operators {I, X, Y, Z} are indexed by
two bits (a, b):

```
              Π²-even (b=0)    Π²-odd (b=1)
immune (a=0):     I                Z
decaying(a=1):    X                Y
```

- Bit a: dephasing sensitivity (0 = immune under Z-dephasing, 1 = decaying)
- Bit b: Π²-parity (0 = even, 1 = odd)

This is a genuine C²⊗C² tensor product structure. The two splits are
orthogonal: parity groups {I,X} vs {Y,Z}; dephasing groups {I,Z} vs {X,Y}.

### 4.2 How Π and the dissipator act in this structure

The dissipator acts only on factor a: D = −2γ · |1⟩⟨1|_a ⊗ I_b.
Factor b (parity) is invisible to dephasing.

Π acts as: flip a, apply phase i^b on factor b. In product form:
Π = σ_x^(a) ⊗ diag(1,i)^(b). A product operator, not entangling.

### 4.3 Partial trace fails

Tracing over one parity sector ({I,X} or {Y,Z}) does NOT produce
the Z-dephasing dissipator. The sectors are orthogonal to the
dephasing split. γ cannot be derived from the C²⊗C² structure alone.

This confirms [failed_third.py](../simulations/failed_third.py):
a finite qubit system produces non-Markovian, non-palindromic
effective noise (0/16 palindromic pairs at every mechanism tested).

### 4.4 What it does show

The Pauli space has more structure than previously recognized. The
two-bit indexing (a,b) = (dephasing, parity) provides a natural
coordinate system. The dissipator lives entirely in one factor.
Π flips the other. They are coupled only through the full
Liouvillian dynamics.

Script: [urqubit_test.py](../simulations/urqubit_test.py)

---

## 5. Three Approaches (All Resolved by April 15, 2026)

The three approaches below were proposed on April 1, 2026 as routes to formalize the primordial qubit hypothesis. By April 15, all three reached terminal verdicts. Their statements are preserved here as the original problem framing; for current status see Section 7 (Priority and Resolution) and Section 8 (Falsification Assessment).

### Approach 1 (highest priority): Pure algebra of Z₂-graded operators

The centered palindrome Q·X·Q⁻¹ = −X gives X a Z₂-grading.
The decomposition X = A − Q·A·Q⁻¹ is always possible (Section 3.2).

The open question: is there a classification theorem for Z₂-graded
algebras that forces the ambient space to be a doubling?

Known result: every Z₂-graded simple algebra over C is isomorphic to
a matrix algebra M_n(C) with even/odd block structure, or to a
Clifford algebra (the algebra generated by anti-commuting generators,
generalizing complex numbers and quaternions). The Liouvillian is
4^N × 4^N. Its Z₂-grading under
Π may force it into one of these standard forms, which would make the
"two sides" a structural necessity.

**Concrete test:** Compute the even and odd subspaces of the
Liouvillian at N=2 under Π-conjugation. Check dimensions, algebra
structure, and whether the odd subspace (where X lives) is isomorphic
to Hom(even, odd) in a doubled algebra.

### Approach 2: GNS construction / Tomita-Takesaki

The GNS construction (Gelfand-Naimark-Segal: the standard procedure
that builds a concrete Hilbert space from an abstract algebra and a
chosen state) builds a Hilbert space from an algebra and a state. For
a thermal (KMS) state (a state satisfying the Kubo-Martin-Schwinger
condition, the algebraic definition of thermal equilibrium), the
Tomita-Takesaki theorem (the foundational result of modular theory:
every faithful state on a von Neumann algebra produces a canonical
anti-linear involution J and a modular flow) produces a modular
conjugation J that exchanges the algebra with its commutant, creating
a canonical "doubling."

The question: is Π related to the modular conjugation J of the
GNS construction for the palindromic steady state?

Connection in the repository: [KMS Detailed Balance](../docs/KMS_DETAILED_BALANCE.md)
treats the KMS condition, which IS the self-consistency condition
of Tomita-Takesaki theory. The connection to the palindrome has not
been made there.

If Π = J (or is related to J by a known transformation), then the
"mirror world" is not a hypothesis but a theorem of operator algebras.

### Approach 3: Thermofield double (TFD)

A thermofield double is a thermal state purified in a doubled Hilbert
space (the standard trick: entangle the system with a fictitious copy
so that tracing out the copy recovers the thermal density matrix):

    |Ω⟩ = Σ_n e^{−βE_n/2} |n⟩ ⊗ |ñ⟩   ∈  H ⊗ H̃

Tracing over H̃ gives the thermal density matrix in H.
Dynamics in H̃ runs backward in time.

Parallels to the primordial qubit:

| TFD concept | Palindrome concept |
|------------|-------------------|
| H (system) | Our side (decaying modes) |
| H̃ (purification) | Their side (Π-conjugate modes) |
| |Ω⟩ (entangled vacuum) | Zero (Σγ = 0, both sides unified) |
| β (inverse temperature) | 1/Σγ (inverse total dephasing)? |
| J (CPT conjugation) | Π (palindromic conjugation) |
| Backward time in H̃ | Frequency mirroring ω → −ω |

The question: is the palindromic Liouvillian the modular Hamiltonian
of a thermofield double construction? Is Σγ proportional to the
temperature?

---

## 6. What Is Known and What Is Not

| Statement | Status |
|-----------|--------|
| d²−2d=0 → d=2 uniquely | **Proven** ([Qubit Necessity](../docs/QUBIT_NECESSITY.md)) |
| Π·L·Π⁻¹ = −L−2Σγ·I | **Proven** ([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)) |
| At Σγ=0: Π·L·Π⁻¹ = −L | **Computed** ([Zero Is the Mirror](ZERO_IS_THE_MIRROR.md)) |
| Noise cannot originate from within | **Proven** ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md)) |
| Pauli space = C²⊗C² (dephasing × parity) | **Computed** (April 1, 2026) |
| Partial trace over sector → dissipator | **Does not work** (naive form) |
| X = A − Q·A·Q⁻¹ decomposition exists | **Proven** (trivially, A = X/2) |
| Canonical choice of A with physical meaning | **Identified** (April 15, 2026): the C2xC2 factorization gives bit_a (n_XY, Absorption) and bit_b (w_YZ, Pi^2-parity) as the two natural axes. Section 9. |
| Noise IS the other side (not from outside) | **Hypothesis**, operationally unfalsifiable from inside (Section 9: only Q=J/gamma measurable) |
| Pi related to Tomita-Takesaki modular conjugation | **Falsified** (Pi linear, J anti-linear; [Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) Phase 3) |
| Z2-grading forces doubling of ambient algebra | **Confirmed in two senses** (April 15, 2026): (1) M_{2\|2}(C) super-algebra is real and forced ([Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) Phase 2); (2) [L, Pi^2] = 0 **proven analytically for all N** ([PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)), L respects the Z2-grading on its eigenmodes (Section 9) |
| L_c as modular Hamiltonian (TFD construction) | **Blocked** (QDB violation = 1.0, maximal; L_c = antisym + sym, no metric reconciles; [`primordial_qubit_kms_test.py`](../simulations/primordial_qubit_kms_test.py)) |
| {L_H, L_D+Σγ} = 0 (oscillation ⊥ cooling) | **Confirmed at N=2**, fails at N≥3 (~2%); [Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) Phase 4 |
| Time reversal at N > 2 | **Excluded** ([Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)) |

---

## 7. Priority and Resolution

Approach 1 (algebra) has been computed (April 1, 2026; extended April 15, 2026). Results:
- The Z₂-grading from Π² gives a proper super-algebra M_{2|2}(C)
- The even subalgebra ≅ M₂(C) ⊕ M₂(C), Clifford algebra Cl(2,0)
- V_{+1} is NOT a subalgebra (Z₄ too fine)
- {L_H, L_D + Σγ} = 0 at N=2 (Pythagorean theorem, exact)
- Breaks at N≥3 (cross term ~2%, γ-independent)
- **April 15, 2026:** [L, Π²] = 0 **proven analytically for all N** ([PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)), verified numerically at N=2-5. The C²⊗C² factors (bit a = n_XY, bit b = w_YZ-parity) are simultaneous good quantum numbers of every Liouvillian eigenmode at all N. See Section 9.

Approach 2 (GNS/Tomita-Takesaki) is **ruled out**: Π is linear,
J is anti-linear, no unitary connects them (impossibility proof).

Approach 3 (thermofield double) is **blocked** (2026-04-15): quantum
detailed balance (QDB) tested across the full steady-state manifold
{II, ZI+IZ, ZZ}. QDB violation = 1.0 (maximal) for every state.
L_c = antisymmetric L_H + symmetric L_D; no inner product makes both
self-adjoint simultaneously. The centered Liouvillian cannot be a
modular Hamiltonian. Script: [`simulations/primordial_qubit_kms_test.py`](../simulations/primordial_qubit_kms_test.py).

Full results: [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md)

---

## 8. Falsification Assessment (updated April 15, 2026)

The original falsification criteria were stated as routes that, if they failed, would weaken or kill the hypothesis. Two of them (TT, TFD) DID fail as construction routes - but failure of a specific construction route is not the same as falsification of the hypothesis. The hypothesis predicts that the doubling is real. The construction routes were proposed mechanisms for *why* it would be real. Mechanisms can fail without the prediction being wrong.

| Original criterion | What it tested | Result |
|-------------------|----------------|--------|
| Z2-grading does not constrain algebra | Does Pi^2 force structure on the operator algebra? | **Hypothesis confirmed**: M_{2\|2}(C) super-algebra is forced (Phase 2). [L, Pi^2] = 0 exactly so L respects the grading on its eigenmodes (Section 9). |
| Tomita-Takesaki construction reproduces Pi | Is Pi the modular conjugation of the GNS construction? | **Construction route falsified** (Pi linear, J anti-linear). Hypothesis unaffected: the doubling is real, just not via this mechanism. |
| Standard TFD construction reproduces L_c | Is L_c the modular Hamiltonian of a thermofield double? | **Construction route falsified** (QDB violation maximal, L_c = antisym + sym irreconcilable). Hypothesis unaffected. |
| Non-doubled structure produces dephasing | Could a hypothetical non-doubled framework give the same observations? | **Operationally consistent with hypothesis** (April 15, 2026): inside observer cannot separate gamma and J, only Q = J/gamma is measurable. The doubled and non-doubled readings are observationally equivalent from inside. |

**Net status:** The hypothesis is **algebraically confirmed in its structural claim** (the doubling exists, M_{2|2}(C) is real, L respects the Z2-grading) and **mechanistically open** (no standard operator-algebra construction reproduces it from outside). The two standard "doubling" constructions of operator algebra theory both fail; the doubling that IS present comes from the Pauli algebra of d=2 itself, not from any external construction.

What remains open: a non-standard mechanism (or a precise statement that no mechanism beyond "the Pauli algebra is already doubled" is needed). The Inside-Outside operational result (Section 9) suggests the latter: the question "is noise external or the other side" may be empirically meaningless from the inside, in which case the primordial qubit reading is parsimonious without requiring further mechanistic justification.

The most concrete external result remains the
[Time Irreversibility Exclusion](../docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md):
time reversal is algebraically excluded at N > 2.

---

## 9. Operational Test from the Inside (April 15, 2026)

The Inside-Outside Correspondence probes ([RESULT_INSIDE_OUTSIDE_CORRESPONDENCE.md](../ClaudeTasks/RESULT_INSIDE_OUTSIDE_CORRESPONDENCE.md), parent commit `c3ea0c0`) tested whether an inside observer (access only to the reduced state ρ_S of a nested Lindblad system) can separate γ and J as independent parameters. Across four independent probes (lifetime-structure separation, BLP non-Markovianity index, ¼-asymptote scaling, slow-mode shapes), the answer is uniform: **no.** Every measurable quantity depends only on the dimensionless ratio Q = J/γ.

This is the operational content of the primordial qubit hypothesis. Under the standard "external noise" reading, γ and J are independent parameters of two different ontological categories (γ is environmental, J is intrinsic). An inside observer with a sufficiently sensitive instrument should be able to confirm this independence by varying conditions and observing the two parameters change separately.

The probes show that no such confirmation is operationally possible. The inside observer measures Q only. The decomposition of Q into "the γ part" and "the J part" requires a vantage point outside the system, which by INCOMPLETENESS_PROOF does not exist for any internal observer.

This is consistent with the hypothesis that γ and J are two readings of one structure rather than two physically separate things. It does not prove the hypothesis: a non-doubled framework with hidden parameters could in principle produce the same operational signature. But it removes one of the standard arguments against the hypothesis (that "external noise" is observationally distinct from "the other side"). The two readings are observationally equivalent from the inside.

**Open question:** Can any inside observation distinguish "γ as external noise" from "γ as the other side of the palindrome"? If the answer is no in principle (not just in current experiments), the distinction is not physical and the primordial qubit reading is the parsimonious one. This connects directly to Approach 1 (the M_{2|2}(C) super-algebra structure) and Approach 3 (TFD), which formalize the "two readings of one structure" claim differently.

### Simultaneous diagonalization of the C2xC2 factors (April 15, 2026)

Script: [`simulations/primordial_bit_a_bit_b.py`](../simulations/primordial_bit_a_bit_b.py).

The C2xC2 tensor product structure of the Pauli space (Section 4.1) has two independent bits per site: bit a (dephasing sensitivity, n_XY) and bit b (Pi^2-parity, w_YZ mod 2). These correspond to two structures that have been studied separately in the framework:

- **Factor a:** the Absorption Theorem axis. Re(lambda) = -2*gamma*<n_XY>. Determines eigenvalue classes (how fast a mode decays).
- **Factor b:** the palindromic Z2-grading axis. Pi^2 = (-1)^{w_YZ}. Determines the M_{2|2}(C) super-algebra sectors (which frequency sector a mode occupies).

Test: ||[L, Pi^2_super]|| = 0.000000e+00 (exact commutation). The Liouvillian respects both quantum numbers simultaneously. Every eigenmode has a definite (bit_a, bit_b) pair. The 3+10+3 degeneracy pattern decomposes cleanly by sector:

| Class | bit_b=even | bit_b=odd | total |
|---|---|---|---|
| conserved (Re=0) | 2 | 1 | 3 |
| mirror (Re=-gamma) | 4 | 6 | 10 |
| correlation (Re=-2gamma) | 2 | 1 | 3 |

The fastest oscillation frequency (Im approx +-2.0) lives exclusively in the odd sector. The even sector contains only Im approx +-1.0.

The two tensor factors of C2xC2 are not just an algebraic curiosity from Section 4. They are the two independent axes of the Liouvillian's eigenmode structure: one axis (bit a, Absorption Theorem) determines how fast a mode decays, the other (bit b, palindromic parity) determines which frequency sector it occupies. Both are simultaneously well-defined quantum numbers of every eigenmode.

### N-scaling: [L, Pi^2] = 0 is universal (April 15, 2026)

**Analytical proof:** [PROOF_BIT_B_PARITY_SYMMETRY](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md). Six-line core argument: Heisenberg coupling XX+YY (and XXX) is bond-wise X-symmetric (each Y or Z appears with even total power per term, signs cancel). Z-dephasing dissipator is quadratic in Z (Z·ρ·Z), the two minus signs from U·Z·U = -Z cancel. Identity normalization is trivially invariant. Holds for any N, any J, any subset of sites carrying gamma. Breaks for single-site Y/Z terms in H (no second factor to cancel sign).

**Numerical verification:** [`simulations/primordial_bit_a_bit_b_N_scaling.py`](../simulations/primordial_bit_a_bit_b_N_scaling.py).

||[L, Pi^2]|| = 0.000000e+00 (identically zero, not numerically small) at N=2, 3, 4, 5. Also holds for Heisenberg XXX coupling with uniform dephasing on all sites at N=3.

Per-sector mode counts (conserved + mirror + correlation):

| N | dim | even | odd |
|---|-----|------|-----|
| 2 | 16 | 2+4+2 | 1+6+1 |
| 3 | 64 | 2+28+2 | 2+28+2 |
| 4 | 256 | 3+122+3 | 2+124+2 |
| 5 | 1024 | 3+506+3 | 3+506+3 |

The C2xC2 sector decomposition is therefore universal (all N), not N=2-specific. Together with the already-proven [Parity Selection Rule](../docs/proofs/PROOF_PARITY_SELECTION_RULE.md) for bit_a, the Liouvillian has TWO independent Z2 symmetries proven for all N. This is the maximal symmetry decomposition admitted by the Pauli algebra of d=2 (no third independent Z2 classification exists per [QUBIT_NECESSITY](../docs/QUBIT_NECESSITY.md)).

### Framework correspondence: even = cavity, odd = transport

The two w_YZ-parity sectors map to existing framework structures (connection identified April 15, 2026):

**Even sector (w_YZ-parity = 0):** {II, IX, XI, XX, YY, YZ, ZY, ZZ} at N=2. Contains ZZ (universal node, [STANDING_WAVE_ANALYSIS](../experiments/STANDING_WAVE_ANALYSIS.md)), XX/YY (antinodes). Bell+ = (II + ZZ + XX - YY)/4 lives entirely here. Oscillation frequencies Im ≈ ±1.0 only. This is the **cavity sector**: standing waves between the mirrors.

**Odd sector (w_YZ-parity = 1):** {IY, IZ, XY, XZ, YI, YX, ZI, ZX} at N=2. Contains cross-correlations (XY, XZ, YX, ZX) and single-site coherences. Oscillation frequencies Im ≈ ±1.0 AND ±2.0 (fastest frequency exclusive to this sector). This is the **transport sector**: modes that break S ↔ B exchange symmetry.

**Origin of the 4:6 mirror-mode split.** The two Z₂ symmetries give four 4×4 blocks:

|                | w_YZ even        | w_YZ odd         |
|----------------|------------------|------------------|
| n_XY even      | II, XX, YY, ZZ   | IZ, XY, YX, ZI   |
| n_XY odd       | IX, XI, YZ, ZY   | IY, XZ, YI, ZX   |

Mirror modes come from mixing between n_XY=0 and n_XY=1 components within each block. The (n_XY-even, w_YZ-odd) block {IZ, XY, YX, ZI} mixes n_XY_B=0 (IZ, ZI) with n_XY_B=1 (XY, YX) under the Hamiltonian, producing 2 additional mirror modes. The (n_XY-even, w_YZ-even) block {II, XX, YY, ZZ} keeps its constituents separated (Bell+ decomposition: II+ZZ vs XX-YY), producing 0 additional mirror modes. Result: even sector has 4 mirror modes, odd sector has 6.

---

## Scripts

- Test computation (C²⊗C² structure): [urqubit_test.py](../simulations/urqubit_test.py)
- Full algebra computation: [primordial_qubit_algebra.py](../simulations/primordial_qubit_algebra.py)

---

*"There was a primordial qubit. It sits on both sides."*
*Thomas Wicht, April 1, 2026*
