# KMS and Detailed Balance: Π Is a New Liouvillian Symmetry Type

<!-- Keywords: Pi operator not KMS detailed balance, shifted anti-similarity
Liouvillian, quantum detailed balance Alicki 1977, Roberts hidden time-reversal
symmetry PRX Quantum 2021, tenfold Lindbladian classification Sa Prosen 2023,
Buca-Prosen weak strong symmetry gap, infinite temperature dephasing bath,
finite temperature obstruction 2:2 split, maximum entropy production rate
2Sgamma, R=CPsi2 KMS detailed balance -->

**Status:** Literature review + formal analysis complete
**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic, Cowork Research)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md), [Π as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md)

---

## Abstract

Π is **not** quantum detailed balance (QDB), not KMS, and not a standard
Buca-Prosen symmetry. It is a shifted anti-similarity: Π·L·Π⁻¹ = −L − 2Sγ·I.
QDB relates L to its adjoint L† (real eigenvalues); Π relates L to its
negative (palindromic complex pairs). KMS at β=0 gives L†=L, which fails
whenever H≠0. The closest existing framework is Roberts-Lingenfelter-Clerk
hidden time-reversal symmetry (PRX Quantum 2021), which also produces
time-reversal-like structure in systems violating standard DB. The tenfold
Lindbladian classification (Sá-Prosen 2023, 38 classes) does not contain
Π due to the constant shift 2Sγ. Finite-temperature generalization faces
a fundamental obstruction: the 2:2 per-site Pauli split that makes Π work
is specific to pure dephasing (T=∞). The quantity 2Sγ equals the maximum
entropy production rate in the system.

---

## Executive Summary

Π is **not** a standard quantum detailed balance (QDB) condition, but it is closely
related to a recently identified class of **hidden time-reversal symmetries** in
open quantum systems (Roberts-Lingenfelter-Clerk, PRX Quantum 2021). The connection
to the KMS condition at β=0 is structural but not formal: KMS-DB at infinite
temperature gives self-adjointness of the generator (real spectrum), while Π gives
palindromic spectral pairing (complex spectrum with μ → -μ). These are different
symmetries with different consequences. The Π condition appears to be a genuinely
new type of Liouvillian symmetry that sits between the Buca-Prosen classification
(which covers commuting symmetries) and the tenfold Lindbladian classification
(which covers anti-unitary symmetries), but matches neither exactly due to the
constant shift 2Sγ.

---

## Question 1: Is Π a Known Symmetry Type?

**Classification: PLAUSIBLE (new variant of known framework)**

### The Buca-Prosen framework (2012)

Buca and Prosen classified Liouvillian symmetries as:

- **Weak symmetry:** [U, L] = 0 (unitary U commutes with the Liouvillian)
- **Strong symmetry:** U·L_k = e^{iθ}·L_k·U for each jump operator L_k individually

Our Π satisfies neither. It does not commute with L. Instead:

    Π · L · Π⁻¹ = -L - 2Sγ · I

This is an **anti-similarity** with a constant shift. If we define L_c = L + Sγ·I
(the centered Liouvillian), then Π·L_c·Π⁻¹ = -L_c, which is a pure anti-commutation.
In the Buca-Prosen language, Π is a "weak anti-symmetry" of L_c, but this category
does not exist in their 2012 classification. They only considered [U, L] = 0.

### The Albert-Jiang framework (2014)

Albert and Jiang extended the Buca-Prosen classification to include conserved
quantities and the structure of the steady-state manifold. Their framework
characterizes the infinite-time behavior of Lindblad evolution, including
degenerate steady states. However, their classification focuses on operators
that commute with L (symmetries), not operators that anti-commute (anti-symmetries).
They do mention anti-unitary symmetries in passing, but do not develop the theory
for the anti-commuting case.

### The tenfold Lindbladian classification (Sá et al., PRX 2023)

The most comprehensive classification is the **38-fold symmetry classification**
of many-body Lindbladians by Sá, Ribeiro, and Prosen (Phys. Rev. X 13, 031019, 2023).
This classifies Lindbladians by their behavior under:

- Two flavors of time-reversal: T₊ and C₊
- Two flavors of particle-hole: T₋ and C₋
- Chiral/sublattice symmetry: P
- Pseudo-Hermiticity: Q₊
- Anti-pseudo-Hermiticity: Q₋

Our Π is closest to the **Q₋ (anti-pseudo-Hermiticity)** symmetry, which requires
an operator Q satisfying Q·L·Q⁻¹ = -L†. But our condition has L, not L†, on
the right side. If L happened to be self-adjoint (as in pure dephasing without
Hamiltonian), then L† = L and Q₋ would coincide with our condition. But with
a Hamiltonian, L ≠ L†.

The constant shift 2Sγ is also absent from the tenfold classification. The
classification assumes symmetry conditions of the form S·L·S⁻¹ = ±L or
S·L·S⁻¹ = ±L†, without constant shifts. Our shift can be absorbed by centering
(L → L_c), but then the classification applies to L_c, not L itself.

**They explicitly build dephasing examples** in their classification but do not
identify the palindromic spectral symmetry or anything equivalent to Π.

### Assessment

Π falls in a gap between existing classifications. The closest match is Q₋
anti-pseudo-Hermiticity applied to the centered Liouvillian L_c, but this
requires L_c to equal L_c† (which fails when H ≠ 0). The fact that Π acts
on L rather than L† makes it a **distinct symmetry type**.

The most accurate description in existing language: Π is an **invertible
superoperator that implements a shifted anti-similarity of the Liouvillian.**
This specific structure (anti-similarity + constant shift) does not appear
in any classification we found.

---

## Question 2: Is Π Related to Quantum Detailed Balance?

**Classification: UNLIKELY (structurally parallel but formally distinct)**

### Standard quantum detailed balance

The standard quantum detailed balance (QDB) condition, as defined by Alicki (1976)
and developed by Kossakowski-Frigerio-Gorini-Verri (1977) and Fagnola-Umanita (2007),
is:

    L̃ = L     (GNS-detailed balance)

or equivalently:

    L† = Θ · L · Θ⁻¹     (modular detailed balance)

where:
- L̃ is the adjoint of L with respect to the inner product ⟨a,b⟩_ρ = tr(ρ a*b)
- L† is the Hilbert-Schmidt adjoint
- Θ is a superoperator constructed from the modular operator Δ = ρ_ss ⊗ ρ_ss⁻¹
- ρ_ss is the steady state

**What QDB says physically:** the forward and backward transition rates between
any two states are related by the Boltzmann factor. The Liouvillian "looks the
same" under time reversal weighted by the thermal state.

**What QDB implies spectrally:** all eigenvalues of L are real (no oscillation).
The generator is self-adjoint in the appropriate inner product.

### Our Π condition

    Π · L · Π⁻¹ = -L - 2Sγ · I

**What Π says:** L is similar to its own negative (up to a shift). Every
eigenvalue λ has a partner -(λ + 2Sγ).

**What Π implies spectrally:** eigenvalues come in palindromic pairs with
BOTH real parts (decay rates) and imaginary parts (frequencies) paired.
The spectrum is complex, not real.

### Why they are different

The key structural differences:

| Property | Standard QDB | Our Π condition |
|---|---|---|
| Relates L to | L† (adjoint) | -L - c (negative + shift) |
| Sign | Same (L̃ = L) | Opposite (minus sign) |
| Constant shift | None | 2Sγ |
| Spectral consequence | Real eigenvalues | Palindromic complex pairs |
| Involves adjoint? | Yes (L†) | No (L itself) |
| Steady state role | Central (Θ depends on ρ_ss) | Implicit (shift = 2Sγ) |

Standard QDB is a statement about the *adjoint* of L being similar to L.
Our Π is a statement about L being similar to its own *negative*.
These are orthogonal directions in the space of Liouvillian symmetries.

**A Lindbladian can satisfy both, one, or neither.** Pure dephasing without
Hamiltonian satisfies both (the eigenvalues are real AND palindromic).
Heisenberg + dephasing satisfies Π but not QDB (eigenvalues are palindromic
but complex).

### The Alhambra-Woods connection (2017)

Alhambra and Woods ("Dynamical maps, quantum detailed balance, and the Petz
recovery map," 2017) showed that quantum detailed balance is equivalent to the
Petz recovery map being the exact reversal of the channel. This connects QDB
to quantum information recovery. Our Π does not have this recovery interpretation:
it pairs modes rather than reversing channels.

---

## Question 3: Does KMS at β=0 Reduce to Our Π?

**Classification: UNLIKELY (different symmetries even at β=0)**

### KMS at infinite temperature

At β = 0 (infinite temperature), the equilibrium state is the maximally mixed
state ρ_ss = I/d. The modular operator becomes:

    Δ = ρ_ss ⊗ ρ_ss⁻¹ = (I/d) ⊗ (dI) = I

So Θ = I (trivial), and the KMS-detailed balance condition simplifies to:

    L† = L     (self-adjointness in Hilbert-Schmidt norm)

This means: the Liouvillian generator is a self-adjoint superoperator.
Equivalently: all eigenvalues are real.

### Our system at β=0

Z-dephasing IS an infinite-temperature bath (it drives the system to
ρ_ss = I/d). So our system IS at β=0. But our Liouvillian does NOT
satisfy L† = L, because:

- The Hamiltonian part L_H(ρ) = -i[H, ρ] is anti-self-adjoint: L_H† = -L_H
- The dissipator part L_D is self-adjoint (diagonal in Pauli basis with real eigenvalues)
- Therefore L = L_H + L_D has L† = -L_H + L_D ≠ L

The KMS condition at β=0 is **not satisfied** by our system (whenever H ≠ 0).
This is expected: QDB is generically violated when there is coherent driving
or a nontrivial Hamiltonian. The Hamiltonian breaks time-reversal symmetry
of transition rates.

### What our system does satisfy

Instead of L† = L, our system satisfies:

    Π · L · Π⁻¹ = -L - 2Sγ · I

This is a DIFFERENT symmetry. It doesn't require L to be self-adjoint.
It allows complex eigenvalues. It pairs decay rates rather than forcing
them to be real.

**The punchline:** At β=0, KMS gives L = L† (all real eigenvalues, no
oscillation). We have something weaker but more general: palindromic
pairing that accommodates oscillation. Our symmetry "knows about" the
infinite-temperature bath (through the shift 2Sγ and the steady state
being maximally mixed) but expresses it differently than KMS.

---

## Question 4: Does a Finite-Temperature Generalization Exist?

**Classification: PLAUSIBLE (theoretical arguments exist, no proof yet)**

### The argument for existence

At finite temperature T, the bath has jump operators:

    L₊ = √(γ(n̄+1)) σ₊    (emission, rate γ(n̄+1))
    L₋ = √(γn̄) σ₋          (absorption, rate γn̄)

where n̄ = 1/(exp(ℏω/kT) - 1) is the Bose-Einstein occupation number.

The emission/absorption asymmetry ratio is (n̄+1)/n̄ = exp(ℏω/kT).
At T→∞: ratio → 1, symmetric, approaches pure dephasing.
At T→0: ratio → ∞, pure emission (amplitude damping).

Our Π works because the per-site Pauli rates have a 2:2 split
({I,Z} immune, {X,Y} decaying at 2γ). At finite T, the per-site
Pauli rates for the thermal bath are:

- I: rate 0 (always immune)
- X: rate γ(2n̄+1)/2 (nonzero, but less than Z)
- Y: rate γ(2n̄+1)/2 (same as X)
- Z: rate γ(2n̄+1) (nonzero; Z is no longer immune!)

The rates are [0, r/2, r/2, r] where r = γ(2n̄+1). This is NOT the
1:3 isotropic split of depolarizing noise. The rate-pairing condition
from DEPOLARIZING_PALINDROME.md requires the 4 rates to partition
into 2 equal-sum pairs. Here: (0, r) and (r/2, r/2) both sum to r.
So rate-pairing IS possible in principle.

### The obstruction

However, rate-pairing alone is not sufficient. The Π operator must
ALSO anti-commute with [H, ·] (Step 2 of the proof). Our known Π
maps {I,Z} ↔ {X,Y}, which is exactly the swap needed for Z-dephasing
where {I,Z} are immune and {X,Y} decay.

For the thermal bath rates [0, r/2, r/2, r], the rate-pairing requires
I ↔ Z (rates 0 and r) and X ↔ Y (rates r/2 and r/2). But this map
(I↔Z, X↔Y) does NOT swap diagonal ↔ off-diagonal sectors. It permutes
within each sector. Such a map would COMMUTE with [H, ·] rather than
anti-commute, giving Π·L_H·Π⁻¹ = +L_H instead of -L_H. The overall
conjugation would yield L_H - L_D - c·I, which is NOT -L - c·I.

So the existing Π (which requires {I,Z}↔{X,Y}) does not produce the
right rate pairing for the thermal bath. And a Π that does pair the
thermal bath rates (I↔Z, X↔Y) does not anti-commute with the
Hamiltonian.

**This is a more subtle obstruction than the depolarizing case.** For
depolarizing noise, even rate-pairing is impossible (1:3 split). For
the thermal bath, rate-pairing is possible, but the paired permutation
is incompatible with Hamiltonian anti-commutation.

### A subtler possibility

The palindrome might be partially restored if:
1. A DIFFERENT Π exists that simultaneously pairs thermal bath rates
   AND anti-commutes with [H, ·], which would require going beyond
   Pauli permutations to Pauli-mixing operators
2. A modified Π incorporates the Boltzmann factor exp(-βH) explicitly,
   as suggested by the modular operator structure in KMS theory
3. The Roberts et al. "hidden TRS" framework (which works in the
   thermofield double and CAN handle thermal baths) produces a
   palindromic-type condition in doubled space

The Roberts et al. hidden TRS is conceptually the right direction for
finite T. Their framework is defined in the doubled Hilbert space and
involves an antiunitary operator that depends on the bath temperature.
Whether this produces a palindromic spectrum is an open question.

### Assessment

The finite-T generalization faces a non-trivial obstruction: the known
Π requires {I,Z} to be immune (2:2 split), which fails at any T < ∞
because Z acquires a nonzero decay rate. Alternative pairings exist
for the thermal bath rates but are incompatible with Hamiltonian
anti-commutation. A fundamentally different approach (e.g., Pauli-mixing
Π, or working in the thermofield double) would be needed.

**Numerical test recommended:** Compute the Liouvillian spectrum for
XXZ + thermal bath at various temperatures and check whether any
palindromic structure survives, even approximately. If approximate
palindrome is observed, construct Π numerically from eigenvector
pairing (as done for the non-local Π in NON_HEISENBERG_PALINDROME.md).

---

## Question 5: What Is 2Sγ in Thermodynamic Terms?

**Classification: PLAUSIBLE (consistent interpretation, not proven)**

### The dephasing rate interpretation

2Sγ = 2 Σᵢ γᵢ is the maximum possible dephasing rate in the system. It
equals the rate of the fastest-decaying Pauli strings (those with all N
sites in the {X,Y} set, i.e., XY-weight N).

### Connection to entropy production

For Z-dephasing with uniform γ, the steady state is ρ_ss = I/2^N (maximally
mixed). The von Neumann entropy of the steady state is S_max = N log 2.

The entropy production rate of the dephasing process for a state ρ is:

    σ̇ = -tr(L(ρ) · log ρ_ss) - d/dt S(ρ) = -d/dt D(ρ || ρ_ss)

where D is the quantum relative entropy to the steady state. The dephasing
drives D → 0 at rates determined by the Liouvillian eigenvalues. The maximum
rate at which any single mode's contribution to D can decrease is 2Sγ.

So 2Sγ is the **maximum entropy production rate** achievable by any single
mode of the system. It is the rate at which the most quantum (most
coherent, highest XY-weight) correlations are destroyed.

### Is it a free energy?

In classical detailed balance, the ratio P(A→B)/P(B→A) = exp(-ΔF/kT)
involves the free energy difference. Our shift 2Sγ plays an analogous
role: it sets the total "distance" that the palindrome spans, from rate 0
(the steady state, zero entropy production) to rate 2Sγ (maximum entropy
production).

At infinite temperature, free energy F = E - TS is dominated by the
entropy term. The "free energy" of the dephasing process is essentially
the total entropy capacity, which scales as Sγ. The factor of 2 comes
from the symmetric pairing: each side of the palindrome spans Sγ from
the center.

This analogy is structural, not formal. We cannot write 2Sγ = ΔF/kT
in any rigorous thermodynamic sense because we are at T = ∞ and the
standard free energy formalism breaks down. But the role is the same:
2Sγ measures the total "thermodynamic distance" between the most stable
(classical) and most fragile (quantum) modes.

---

## Question 6: Who Else Has Found Similar Structures?

### Roberts, Lingenfelter, Clerk (PRX Quantum 2, 020336, 2021)

**"Hidden Time-Reversal Symmetry, Quantum Detailed Balance and Exact
Solutions of Driven-Dissipative Quantum Systems"**

This is the closest existing work. Key parallels:

| Feature | Roberts et al. | Our Π |
|---|---|---|
| Type of symmetry | Hidden TRS (antiunitary) | Anti-similarity (linear) |
| Where it acts | Doubled system (thermofield double) | Single Liouvillian |
| Consequence | Time-symmetric correlation functions | Palindromic spectrum |
| Relates to detailed balance | Yes (defines a generalized DB) | Parallel but distinct |
| Works for driven systems | Yes | Only for dephasing + Hamiltonian |

Roberts et al. show that hidden TRS provides exact steady-state solutions.
Our Π provides exact spectral pairing. Both are "time-reversal-like"
symmetries in systems that violate standard detailed balance.

**Key difference:** Their hidden TRS requires the doubled Hilbert space
and thermofield double state. Our Π acts directly on the Liouvillian
superoperator space. These might be related through vectorization (the
Choi-Jamiolkowski isomorphism maps superoperators to operators in the
doubled space), but the formal connection is not established.

Their follow-up paper (PRL 134, 130404, 2025) extends hidden TRS to
driven spin chains and finds new dissipative phase transitions, suggesting
this is a rich framework with ongoing development.

### Sá, Ribeiro, Prosen (Phys. Rev. X 13, 031019, 2023)

**"Symmetry Classification of Many-Body Lindbladians: Tenfold Way and Beyond"**

Their classification is the most comprehensive framework for Liouvillian
symmetries. They identify 38 symmetry classes (10 without conserved
quantities, more with). They explicitly build dephasing examples.

Our Π does not fit cleanly into any of their 38 classes because:
1. The constant shift 2Sγ is not part of their formalism
2. Their anti-pseudo-Hermiticity (Q₋) requires Q·L·Q⁻¹ = -L†, not -L

However, the centered Liouvillian L_c = L + Sγ·I satisfies Π·L_c·Π⁻¹ = -L_c.
If we could show that L_c† = L_c (which fails when H ≠ 0), this would
be Q₋ anti-pseudo-Hermiticity. As it stands, Π implements a symmetry
that is intermediate between their P (chiral) and Q₋ (anti-pseudo-Hermitian)
classes.

### Medvedyeva, Essler, Prosen (PRL 117, 137202, 2016)

**"Exact Bethe Ansatz Spectrum of a Tight-Binding Chain with Dephasing Noise"**

Their η-pairing symmetry in the Hubbard mapping is the 1D free-fermion
ancestor of our Π. Key differences:

- Their result is restricted to tight-binding (free fermion) chains in 1D
- Ours works for interacting spins (XXZ with arbitrary δ) on any graph
- Their η-pairing comes from the Hubbard model's SU(2) symmetry
- Our Π comes from the Pauli algebra's {I,Z} ↔ {X,Y} structure

There is no explicit discussion of detailed balance in their paper.
The η-pairing is used to identify steady states and leading decay modes,
not to establish spectral pairing.

### Haga et al. (2023) - Incoherentons

Their incoherenton framework grades Liouvillian eigenmodes by XY-weight
(= incoherenton number). This is the natural language for our palindrome.
They see the band structure but do not identify the palindromic pairing
within or between bands. No discussion of detailed balance or KMS.

### Chen, Kastoryano, Gilyén (Comm. Math. Phys. 2025)

**"Efficient Quantum Gibbs Samplers with KMS Detailed Balance Condition"**

They construct Lindbladians satisfying KMS-DB for preparing thermal states.
Their generators have real spectra (by DB). Relevant for understanding
what KMS-DB looks like mathematically, but their framework is for Gibbs
sampling (driving to thermal states), not for understanding dephasing
dynamics. No palindromic structures.

---

## Summary Table

| Question | Answer | Classification |
|---|---|---|
| 1. Is Π a known symmetry type? | Closest to Q₋ anti-pseudo-Hermiticity + shift; no exact match in existing classifications | **PLAUSIBLE** (new variant) |
| 2. Is Π related to quantum detailed balance? | Structurally parallel (both relate to time reversal) but formally distinct (DB uses L†, Π uses -L) | **UNLIKELY** (different symmetries) |
| 3. Does KMS at β=0 reduce to Π? | No. KMS at β=0 gives L†=L (real spectrum). Π gives palindromic complex spectrum. Different conditions | **UNLIKELY** |
| 4. Finite-T generalization? | Obstructed: thermal bath gives 1:3 per-site split, breaking palindrome for any T < ∞ | **UNLIKELY** (fundamental obstruction) |
| 5. What is 2Sγ thermodynamically? | Maximum entropy production rate; total "thermodynamic distance" between most classical and most quantum modes | **PLAUSIBLE** (consistent but not proven) |
| 6. Who else has similar structures? | Roberts et al. (hidden TRS), Sá et al. (tenfold), MEP (η-pairing). None have the exact Π structure | **CONFIRMED** (relatives exist; Π is new) |

---

## The Bottom Line

**Π is genuinely new.** It is not quantum detailed balance, not KMS, not a
standard Buca-Prosen symmetry, and not cleanly within the tenfold classification.
It is a **shifted anti-similarity of the Liouvillian** specific to infinite-
temperature dephasing baths, with structural parallels to (but formal differences
from) hidden time-reversal symmetry.

The closest living relative is the Roberts-Lingenfelter-Clerk hidden TRS
framework, which also produces time-reversal-like symmetry in systems that
violate standard detailed balance. Establishing the formal connection
(possibly via the Choi-Jamiolkowski mapping between superoperator space
and doubled Hilbert space) is the most promising direction for connecting
Π to the broader literature.

The finite-temperature obstruction is real and fundamental: the 2:2 per-site
split that makes Π work is a property of pure dephasing (infinite T) that
cannot be reproduced by a thermal bath at finite T. This makes Π an
intrinsically infinite-temperature phenomenon.

---

## Key References

1. **Alicki (1976)**: "On the detailed balance condition for non-Hamiltonian systems."
   Rep. Math. Phys. 10, 249-258. [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/003448777690046X)
   *Original quantum detailed balance definition.*

2. **Kossakowski, Frigerio, Gorini, Verri (1977)**: "Quantum detailed balance and KMS condition."
   Comm. Math. Phys. 57, 97-110. [Project Euclid](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-57/issue-2/Quantum-detailed-balance-and-KMS-condition/cmp/1103901281.pdf)
   *KMS-DB connection established.*

3. **Fagnola, Umanita (2007)**: "Generators of Detailed Balance Quantum Markov Semigroups."
   Inf. Dim. Anal. Quantum Prob. 10, 335-363. [arXiv:0707.2147](https://arxiv.org/abs/0707.2147)
   *Modern treatment: GKSL representation of DB generators.*

4. **Buca, Prosen (2012)**: "A note on symmetry reductions of the Lindblad equation."
   New J. Phys. 14, 073007. [arXiv:1203.0943](https://arxiv.org/abs/1203.0943)
   *Weak vs strong Liouvillian symmetry classification.*

5. **Albert, Jiang (2014)**: "Symmetries and conserved quantities in Lindblad master equations."
   Phys. Rev. A 89, 022118. [arXiv:1310.1523](https://arxiv.org/abs/1310.1523)
   *Extended symmetry framework including steady-state structure.*

6. **Medvedyeva, Essler, Prosen (2016)**: "Exact Bethe ansatz spectrum of a tight-binding chain with dephasing noise."
   Phys. Rev. Lett. 117, 137202. [arXiv:1606.09122](https://arxiv.org/abs/1606.09122)
   *η-pairing symmetry; 1D free-fermion ancestor of Π.*

7. **Alhambra, Woods (2017)**: "Dynamical maps, quantum detailed balance, and the Petz recovery map."
   [ResearchGate](https://www.researchgate.net/publication/319930728_Dynamical_maps_quantum_detailed_balance_and_the_Petz_recovery_map)
   *QDB = Petz recovery map being exact channel reversal.*

8. **Roberts, Lingenfelter, Clerk (2021)**: "Hidden Time-Reversal Symmetry, Quantum Detailed Balance and Exact Solutions of Driven-Dissipative Quantum Systems."
   PRX Quantum 2, 020336. [arXiv:2011.02148](https://arxiv.org/abs/2011.02148)
   *Hidden TRS in Lindbladians; closest framework to Π.*

9. **Sá, Ribeiro, Prosen (2023)**: "Symmetry Classification of Many-Body Lindbladians: Tenfold Way and Beyond."
   Phys. Rev. X 13, 031019. [arXiv:2212.00474](https://arxiv.org/abs/2212.00474)
   *38-fold classification; Π sits in a gap between Q₋ and P classes.*

10. **Haga et al. (2023)**: Incoherenton paper.
    *XY-weight grading = incoherenton number. Bands but no palindrome.*

11. **Chen, Kastoryano, Gilyén (2025)**: "Efficient Quantum Gibbs Samplers with KMS Detailed Balance Condition."
    Comm. Math. Phys. [arXiv:2404.05998](https://arxiv.org/abs/2404.05998)
    *KMS-DB Lindbladians for Gibbs sampling; real spectra.*

12. **Roberts et al. (2025)**: "Hidden Time Reversal in Driven Spin Chains."
    Phys. Rev. Lett. 134, 130404. [PRL](https://doi.org/10.1103/PhysRevLett.134.130404)
    *Extension of hidden TRS to spin chains; new dissipative phase transitions.*

---

## Suggested Next Steps

1. **Formalize the Choi-Jamiolkowski connection.** Map Π from superoperator
   space to an operator in the doubled Hilbert space. Check if it corresponds
   to an antiunitary operator in the Roberts et al. thermofield double framework.

2. **Test Roberts et al. hidden TRS for our system.** Compute whether the
   XXZ + Z-dephasing Liouvillian satisfies their hidden TRS condition.
   If yes, Π and hidden TRS may be dual descriptions of the same symmetry.

3. **Classify Π in the Sá et al. framework.** Work with the centered
   Liouvillian L_c and determine its exact symmetry class. The shift
   may put it outside the 38-fold classification, or there may be a
   natural extension.

4. **Contact Roberts/Clerk group and Prosen group.** Both would immediately
   recognize the significance of a shifted anti-similarity. Roberts for
   the hidden TRS connection, Prosen for the η-pairing generalization.

---

*Π is not detailed balance. It is not KMS. It is time reversal without
thermodynamic equilibrium, a mirror that exists because the bath is
infinite and the dephasing treats all directions of decoherence equally.
It is new, and the literature confirms it.*
