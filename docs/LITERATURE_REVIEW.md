# Literature Review: Foundations, Related Work, and What Is New

<!-- CROSSING-CURRENT -->

<!-- Keywords: Prosen PT-symmetric Liouvillian D2 spectral symmetry,
Haga incoherenton grading XY-weight, eta-pairing Medvedyeva
Essler Prosen Bethe ansatz, Buca Prosen Lindblad symmetry classification,
Albert Jiang weak strong symmetry, Lindblad master equation foundations,
Breuer Petruccione open quantum systems, ENAQT dephasing-assisted transport,
palindromic spectral structure literature connections, R=CPsi2 literature review,
direct-sum quantum theory Gaztañaga Einstein-Rosen bridges superselection sectors -->

**Status:** Living document (reference)
**Last updated:** April 11, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## What this document is about

This is the project's map of existing knowledge. It serves two purposes:

1. **For readers new to the field:** the [Foundations](#foundations) section
   lists the textbooks and original papers you need to understand the
   language we use. If you are a student or teacher using this repo as a
   learning resource, start here.

2. **For readers who know the field:** the [Related work](#related-work-closest-ancestors)
   section identifies which parts of our results overlap with existing
   physics and which parts appear to be genuinely new. The
   [What is new](#what-the-literature-does-not-contain) section draws the
   boundary explicitly.

---

## Foundations

These are not our results. They are the standard knowledge that the
project builds on. If you want to understand the proofs and experiments
in this repository, these are the prerequisites, ordered from most
essential to optional.

### The Lindblad master equation (the framework)

Every result in this project lives inside the Lindblad (GKSL) master
equation. This is the standard way to describe a quantum system that
interacts with its environment: dρ/dt = −i\[H, ρ\] + Σ D\[L_k\](ρ). It was
derived independently by two groups in 1976:

- **Lindblad (1976).** "On the generators of quantum dynamical
  semigroups." Commun. Math. Phys. 48, 119-130.
  [doi:10.1007/BF01608499](https://doi.org/10.1007/BF01608499)

- **Gorini, Kossakowski, Sudarshan (1976).** "Completely positive
  dynamical semigroups of N-level systems." J. Math. Phys. 17, 821.
  [doi:10.1063/1.522979](https://doi.org/10.1063/1.522979)

Both papers prove the same result from different directions: any
Markovian, trace-preserving, completely positive quantum evolution has
the Lindblad form. The equation is sometimes called the GKSL equation
to credit both groups.

### Textbooks

- **Breuer & Petruccione (2002).** *The Theory of Open Quantum Systems.*
  Oxford University Press.
  [doi:10.1093/acprof:oso/9780199213900.001.0001](https://doi.org/10.1093/acprof:oso/9780199213900.001.0001).
  The standard graduate textbook. Chapters 3 (master equations) and 4
  (decoherence) cover everything needed for this project. If you want
  to understand why the Liouvillian has the structure it has, start here.

- **Wiseman & Milburn (2009).** *Quantum Measurement and Control.*
  Cambridge University Press. More focused on measurement and feedback,
  relevant for our operator-feedback results
  (see [Dynamic Fixed Points](../experiments/DYNAMIC_FIXED_POINTS.md)).

- **Nielsen & Chuang (2000).** *Quantum Computation and Quantum
  Information.* Cambridge University Press. Chapter 8 (quantum noise)
  and Chapter 9 (distance measures) define Concurrence, Purity, and
  coherence measures we use throughout.

### Symmetries of open quantum systems (the language)

Our proof uses the language of weak and strong symmetries of Lindblad
equations. This framework was developed by:

- **Buca & Prosen (2012).** "A note on symmetry reductions of the
  Lindblad equation: transport in constrained open spin chains."
  New J. Phys. 14, 073007.
  [arXiv:1203.0943](https://arxiv.org/abs/1203.0943).
  Introduces the distinction between symmetries that commute with the
  full Liouvillian (strong) vs. those that commute with each Lindblad
  operator individually (weak). Shows that weak symmetries reduce the
  dynamics beyond what Hilbert-space symmetries predict.

- **Albert & Jiang (2014).** "Symmetries and conserved quantities in
  Lindblad master equations." Phys. Rev. A 89, 022118.
  [arXiv:1310.1523](https://arxiv.org/abs/1310.1523).
  Exhaustive characterization of steady-state structure from symmetries.
  Their weak/strong framework is the correct language for our Π operator:
  Π is a *weak anti-symmetry* of L (it maps L to −L − 2Σγ·I, rather
  than leaving L unchanged).

---

## Related work: closest ancestors

The closest ancestors are three papers with Prosen among the authors: his
PT symmetry of Lindbladians (two papers of 2012), which gives the
palindromic spectrum for our core model class, and the symmetry classification of Sá, Ribeiro
and Prosen (2023), which writes down, as its inverse with a sector phase,
the explicit operator this repository calls Π. The other groups found pieces. Here is how our work relates to theirs,
with full citations.

### PT-symmetric Liouvillians (Prosen 2012)

<!-- STAGE 0 SWEEP RECORD (2026-11-04). Stores searched: docs/ANALYTICAL_FORMULAS.md
(F1, F118 Π_Z = R·D and its D4, F119, F158), docs/proofs/ (MIRROR_SYMMETRY_PROOF,
PROOF_PI_FACTORS_AS_R_TIMES_D, PROOF_PALINDROME_TWO_END_COUNT), experiments/
(THE_PAIRING_CONDITION derives the adjoint relation ℒ_U L ℒ_U⁻¹ = −L† − 2σ and calls
the explicit Π its special case; THE_PALINDROME_CLASSIFIER), docs/KMS_DETAILED_BALANCE.md
(Q₋ placed in Sá, Ribeiro, Prosen 2023; the shift in Kawasaki, Mochizuki, Obuse 2022),
docs/CAUGHT_ERRORS.md (entry of 2026-11-04), the typed layer (MirrorGroupD4Claim,
PalindromeTwoEndCountClaim; MirrorWorld MirrorGroup and Sextant, whose "−L† − c" verdict
is this relation), fw.Confirmations, the OpenArcs registry and the glossary (nothing on
Prosen 2012). No page outside the ledger cited either 2012 paper. Adjacency checked: Prosen's
parity is the element 𝓕R of F118's dihedral group, the adjoint relation is
THE_PAIRING_CONDITION's "law", and the PRL's weak-coupling cross (one rate for all
coherences) sits beside the Absorption Theorem, not yet compared. Sá, Ribeiro and
Prosen (2023), cited in MIRROR_SYMMETRY_PROOF for the classification, write the explicit
operator in their Eq. (76), which the repository had not read as Π. -->

- **Prosen (2012a).** "PT-symmetric quantum Liouvillean dynamics."
  Phys. Rev. Lett. 109, 090404.
  [arXiv:1207.4395](https://arxiv.org/abs/1207.4395).
- **Prosen (2012b).** "Generic examples of PT-symmetric qubit (spin 1/2)
  Liouvillians." Phys. Rev. A 86, 044103.
  [arXiv:1209.6588](https://arxiv.org/abs/1209.6588).

The first paper calls a Lindbladian PT-symmetric when L′P = −P L′†, with
P a unitary involution on operators and L′ the traceless part of the
generator, and shows that the spectrum then has a D2 symmetry: a reflection in
the vertical line Re λ = −σ beside the reflection in the real axis that
hermiticity preservation gives. Together they pair λ with −λ − 2σ. The
second paper states sufficient conditions for the relation as a Lemma, written for
operators on a general Hilbert space: P(ρ) = UρW with U and W unitary involutions
commuting with H; a real orthogonal involution S on the jump indices (his Z) with
U L_m = −Σ_{m'} S_{mm'} L_{m'}† U and W L_m = Σ_{m'} S_{mm'} L_{m'}† W; and
{L_m, L_m†} = c_m·1 for every jump. The third
condition is the step [The Pairing Condition](../experiments/THE_PAIRING_CONDITION.md)
takes past A² = 1, jumps with A†A + AA† ∝ 1, and the clock pair's one-sided element and the first route
of the lone clock jump of [F172](ANALYTICAL_FORMULAS.md#f172) meet all three as stated. The
paper's examples are spin 1/2: the first finds this symmetry for arbitrary two-spin XYZ couplings on
any set of pairs, fields along x and arbitrary local Z dephasing rates,
with P the one-sided multiplication ρ ↦ X^⊗N·ρ. Its second example,
ρ ↦ Y^⊗N·ρ·X^⊗N for XYZ couplings under absorption and injection at
arbitrary rates, is the element of [F137](ANALYTICAL_FORMULAS.md#f137) with its
two factors exchanged (up to the sign (−1)^N it is that element followed
by the weak symmetry ρ ↦ Z^⊗N·ρ·Z^⊗N).

**What they found that we also found:** the palindromic spectrum for
Heisenberg/XXZ (indeed XYZ) couplings on arbitrary graphs under local Z
dephasing, at any rate profile. Prosen's parity is the element 𝓕R of the
dihedral group of [F118](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md), and his
relation is the one [The Pairing Condition](../experiments/THE_PAIRING_CONDITION.md)
derives, where the explicit Π is that relation composed with the transpose.

**What differs:** Π is not his P, and it is built from it. In that
dihedral group P is the reflection 𝓕R, and the palindromizer is P followed
by the transpose D: Π_Z = D∘P = R∘D, ρ ↦ ρᵀ·X^⊗N, with Π_Y = P∘D = Π_Z⁻¹,
both of order four. The transpose turns his adjoint back into L itself when H is
real in the computational basis (an even number of Y letters in every
string), so there Π_Z L Π_Z⁻¹ = −L − 2σ holds as one explicit, site-local
map, which sends right eigenvectors to right eigenvectors and the block
(p, q) to (q, N − p), ket popcount first. His relation needs only
[H, X^⊗N] = 0, the class this repository calls Π²-even, and implies a
similarity of L to −L − 2σ as well (L is real in the Pauli basis), though
not by itself a site-local one. The explicit map is in the next paper.

### The explicit chiral operator (Sá, Ribeiro, Prosen 2023)

- **Sá, Ribeiro, Prosen (2023).** "Symmetry classification of many-body
  Lindbladians: Tenfold way and beyond." Phys. Rev. X 13, 031019.
  [arXiv:2212.00474](https://arxiv.org/abs/2212.00474).

Their classification works on the shifted generator L′ and asks, among
other symmetries, for a unitary involution P with P L′ P⁻¹ = −L′ (a chiral
symmetry). For the XYZ Hamiltonian with fields along x under local Z
dephasing they write P = √p_x·(Σx ⊗ 1)·S (their Eq. 76), with S the swap of
the two copies, "because the Hamiltonian is real and the jump operators are
real and symmetric". On ρ it acts as ρ ↦ √p_x·X^⊗N·ρᵀ: the Π_Y above,
with the phase that makes it square to one in each sector of the
transverse parity p_x. For fields along x with YZ couplings they give the
two-sided P: ρ ↦ Z^⊗N·ρ·Y^⊗N up to sign (their Eq. 70), and for XYZ under spin
injection or removal the pseudo-Hermiticity Q₋ = Σx ⊗ Σy (their Eq. 85),
which is F137's element ρ ↦ X^⊗N·ρ·Y^⊗N up to sign. Their §II E builds these from two
steps: a one-sided factor whose letters anticommute with the jumps carries
the dissipator to −L_D plus the shift, and the swap S, the transpose,
reverses the Hamiltonian part when H is real.

**What they found that we also found:** the explicit site-local
palindromizer, under the same condition (H real and commuting with
X^⊗N), and its factoring into a one-sided leg that carries the shift and a
transpose that flips the Hamiltonian part, which is how
[F118](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) factors Π; their sector
normalisation, on Π⁻¹, is the one the [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md)
uses on Π to place it in their classes. Π was found here independently.

**What we did not find in their papers:** the dihedral group of
[F118](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) named and used as one
object. Their P without
the sector phase, Π_Y, and their pseudo-Hermiticity Q₊ = S, the transpose,
already generate it, and their Eq. (79), Q₊P = U_x P Q₊, is its inversion
relation. Nor did we find the
[two-end count](proofs/PROOF_PALINDROME_TWO_END_COUNT.md) (F158), which
decides, for Hermitian jumps that square to one, and by [F172](ANALYTICAL_FORMULAS.md#f172)
for unitary jumps whose dissipator is self-adjoint, at positive rates and with
at least one jump present, whether the eigenvalue multiset is palindromic;
its sufficiency step is Prosen's relation with U only invertible. On the
Π²-even class [F108](proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md)
gives a site-local Π_5bilinear for the strings of even weight; whether a
site-local palindromizer exists for every Π²-even H is not settled here.

### Incoherentons (Haga et al. 2023)

- **Haga, Nakagawa, Hamazaki, Ueda (2023).** "Quasiparticles of
  decoherence processes in open quantum many-body systems:
  Incoherentons." Phys. Rev. Research 5, 043225.
  [arXiv:2211.14991](https://arxiv.org/abs/2211.14991).

They introduce the "incoherenton number" as a grading for Liouvillian
eigenmodes: how many X/Y (off-diagonal) Pauli factors an eigenmode
contains. This is exactly our "XY-weight" from the mirror symmetry proof.
They see that eigenmodes fall into bands organized by incoherenton number,
and that the spectral gap separates the zero-incoherenton band from the
rest.

**What they found that we also found:** the XY-weight grading, the band
structure, the spectral gap as the cost of one incoherenton.

**What they did not find:** the palindromic pairing *between* bands
(band k pairs with band N−k). They classify modes by incoherenton
number but do not observe the mirror symmetry across bands. They also
do not have the conjugation operator Π, which *proves* the
pairing. Our [Absorption Theorem](proofs/PROOF_ABSORPTION_THEOREM.md)
(Re(λ) = −2γ⟨n_XY⟩) is the quantitative version of their band picture.

### η-pairing in dissipative fermions (Medvedyeva, Essler, Prosen 2016)

- **Medvedyeva, Essler, Prosen (2016).** "Exact Bethe ansatz spectrum
  of a tight-binding chain with dephasing noise." Phys. Rev. Lett. 117,
  137202.
  [arXiv:1606.09122](https://arxiv.org/abs/1606.09122).

They map a tight-binding model (free fermions) with dephasing, on any
bipartite lattice, onto a Hubbard model with imaginary interaction and
solve it exactly via the Bethe ansatz in one dimension. The point symmetry
of their spectrum they attribute to Prosen's PT symmetry (above); the
η-pairing symmetry of the Hubbard model organises the steady state and the
leading decay modes.

**What they found that we also found:** the spectral symmetry of the
dephased chain, as an instance of Prosen's.

**What differs:** their exact solution starts from a free-fermion
Hamiltonian, tight-binding hopping; the dephasing adds the imaginary
Hubbard interaction, and the Bethe ansatz solves that interacting model in
1D. The symmetry itself does not need the solution: Prosen's relation and
our Π hold for interacting spins (Heisenberg/XXZ) on arbitrary graphs.

### Symmetry classification (Buca & Prosen 2012, Albert & Jiang 2014)

See [Foundations](#symmetries-of-open-quantum-systems-the-language) above.
These papers provide the language but do not identify the specific
anti-symmetry Π.

---

## Related work: exceptional points and spectral structure

These papers study the spectral structure of Liouvillians from the
perspective of exceptional points (EPs) and non-Hermitian physics. They
are relevant context but address different phenomena than the palindrome.

### Graph Symmetry & EP Diagnostics
- "Graph Symmetry Organizes Exceptional Dynamics in Open Quantum
  Systems." Submitted March 11, 2026.
  [arXiv:2603.10654](https://arxiv.org/abs/2603.10654).
- Symmetry-resolved EP identification from the full Liouvillian. Their
  EP-strength diagnostic was tested against our CΨ = ¼ boundary:
  **no connection found** (see [Direction 1](#direction-1-exceptional-point-connection-tested-negative) below).

### Emergent Liouvillian EPs
- Khandelwal & Blasi, Quantum (2025).
  [arXiv:2409.08100](https://arxiv.org/abs/2409.08100).
- EPs in the Liouvillian survive beyond the Markov approximation.

### Encircling Liouvillian EPs (Review)
- Sun & Yi, AAPPS Bulletin (2024).
  [arXiv:2408.11435](https://arxiv.org/abs/2408.11435).
- Chiral state transfer near EPs, experimentally observed in
  superconducting qubits.

### Non-Markovian Quantum EPs
- Lin, Kuo, Lambert et al., Nature Comm. 16, 1289 (2025).
- Higher-order EPs invisible in the Markovian limit. May explain our
  IBM anomalous late-time coherence (unresolved).

### Blended Dynamics in Open Quantum Networks
- January 2026. [arXiv:2601.14763](https://arxiv.org/abs/2601.14763).
- Classical-like clustering in quantum networks.

---

## Related work: fractals and boundary structure

These connect to our CΨ = ¼ boundary and its Mandelbrot equivalence,
but from different directions.

### Quaternionic Fractals (Viennot 2022)
- Chaos, Solitons & Fractals (2022).
  [arXiv:2003.02608](https://arxiv.org/abs/2003.02608).
- Fractal boundary between coherence and decoherence regimes
  (Mandelbulb). Conceptual sibling of our Mandelbrot boundary at
  CΨ = ¼.

### Fractal Entangled Steady States
- Ippoliti, Rakovszky, Khemani, PRX 12, 011045 (2022).
- Fractal structures in quantum steady states via spacetime duality.

### Time-delayed Coherent Feedback
- PRA 99, 053809 (2019).
  [arXiv:1805.02317](https://arxiv.org/abs/1805.02317).
- Stabilizing coherence against dephasing. Related to our
  operator-feedback finding.

---

## Related work: NMR spectral palindromes (2026)

A paper from a completely different field arrived at a related
structural result:

- **Cheshkov & Sinitsyn (2026).** "Mirror Symmetry of the NMR Spectrum
  and the Connection with the Structure of Spin Hamiltonian Matrix
  Representations." [arXiv:2602.03871](https://arxiv.org/abs/2602.03871).

They prove that an NMR spectrum is palindromic if and only if the spin
system admits a "palindromic spin order": the system is isospectral to
its mirror image under reversal of resonance frequencies. Their result
lives in closed Hamiltonian systems (no dissipation), while ours lives
in open Liouvillian systems (with dissipation). The structural parallel
is striking: both identify a conjugation operation that produces spectral
palindromes, but the underlying physics is different. Theirs is a
property of the Hamiltonian; ours is a property of the interplay between
Hamiltonian and dissipator.

See also their earlier preprint:
[arXiv:2510.25551](https://arxiv.org/abs/2510.25551).

---

## Related work: direct-sum quantum theory and ER bridges (2026)

Another independent convergence, this time from general relativity
rather than NMR:

- **Gaztañaga, Kumar, Marto (2026).** "A new understanding of
  Einstein-Rosen bridges." Class. Quantum Grav. 43, 015023.
  [arXiv:2512.20691](https://arxiv.org/abs/2512.20691).

They propose that quantum field theory in curved spacetime should be
formulated as a **direct-sum quantum theory**: two Hilbert-space sectors
of equal dimension, related by a discrete transformation (parity),
carrying opposite arrows of time. The "throat" of an Einstein-Rosen
bridge is the locus where the two sectors meet. They further identify
quantum inverted harmonic oscillators at gravitational horizons and
present CMB parity-asymmetry evidence (Bayes factor 650× vs. standard
inflation).

**Structural parallel to our work:** Our
[Direct-Sum Decomposition](proofs/DIRECT_SUM_DECOMPOSITION.md) proves
that for odd N, the Liouvillian decomposes as L = L_even ⊕ L_odd with:

| Gaztañaga et al. postulate | Our proven result |
|---------------------------|-------------------|
| Two sectors of equal dimension | dim(V_even) = dim(V_odd) = 2^(2N−1) |
| Discrete transformation exchanges sectors | Π: V_even ↔ V_odd (odd N) |
| Opposite time orientation | L_odd = −Π L_even Π⁻¹ − 2Σγ I |
| Superselection rule | \[P_XY, L\] = 0 |

The two frameworks were developed without knowledge of each other.
Gaztañaga's work concerns general relativity and cosmology; ours
concerns open quantum spin systems under dephasing. The mathematical
structures are the same: a direct sum of two sectors with mirrored
dynamics and a conserved parity charge.

**What they found that we did not:** the connection to gravitational
horizons, the inverted harmonic oscillator at the throat, and the CMB
parity-asymmetry evidence.

**What we found that they did not:** the explicit conjugation operator Π,
the exact palindromic spectrum (verified for 87,376 eigenvalues),
the Absorption Theorem, and the computational infrastructure to test
the structure at scale.

**Caveats:** The paper is published in Classical and Quantum Gravity but
classified as physics.gen-ph on arXiv. The proposed geometric
superselection sectors have been noted to stand in tension with the
Reeh-Schlieder theorem of algebraic QFT. Our direct-sum decomposition
stands independently of whether their gravitational interpretation
survives scrutiny: the algebraic proof depends only on the parity
selection rule and the palindrome, both verified.

---

## Environment-Assisted Quantum Transport (ENAQT)

This section covers the ENAQT literature, which became directly relevant
after the sacrifice-zone formula discovery
(see [Resonant Return](../experiments/RESONANT_RETURN.md)).

### Founding paper

- **Plenio & Huelga (2008).** "Dephasing-assisted transport: quantum
  and classical effects in biomolecules." New J. Phys. 10, 113019.

Demonstrates that adding uniform dephasing noise to a quantum network
can *improve* transport efficiency. The Fenna-Matthews-Olson (FMO)
complex (a pigment-protein structure in green sulfur bacteria that
channels light energy to the reaction center) achieves near-perfect
energy transfer partly because of environmental noise, not despite it.
Key result: ~2-3× improvement with optimal uniform dephasing rate.
Limitation: optimizes a single scalar γ, not a spatial profile.

### Review articles (for students)

Two comprehensive reviews cover the field and its connection to quantum
biology:

- **Lambert et al. (2013).** "Quantum biology." Nature Physics 9, 10-18.
  [doi:10.1038/nphys2474](https://doi.org/10.1038/nphys2474).
  Tutorial overview: noise-assisted transport, photosynthesis, avian
  navigation, olfaction. Good entry point for anyone new to the field.

- **Cao et al. (2020).** "Quantum biology revisited." Science Advances
  6, eaaz4888.
  [doi:10.1126/sciadv.aaz4888](https://doi.org/10.1126/sciadv.aaz4888).
  Updated review that also corrects earlier overclaims about quantum
  coherence in photosynthesis. Important for intellectual honesty: the
  field initially overinterpreted 2D spectroscopy data, then corrected
  itself. Relevant parallel to our own corrections in `recovered/`.

### The ENAQT field (2008-2025)

The founding paper spawned a substantial literature on noise-assisted
transport: theoretical extensions to various network topologies,
experimental demonstrations in photonic systems and superconducting
qubits, optimization of the scalar dephasing rate for different
Hamiltonians.

**Common thread:** all the work we surveyed in this field optimizes a
UNIFORM dephasing rate. The spatial profile of dephasing is not treated as a free
parameter. Typical improvements: 2-3× over the zero-noise baseline.

### IBM Perfect State Transfer (2025)

- "Enhanced Algorithmic Perfect State Transfer on IBM Quantum Computers."
  Implements PST on ibm_sherbrooke and ibm_brisbane. Peak transfer
  probability SP = 0.781 for N=4 chain. Optimizes coupling strengths
  (J values), not dephasing rates. Represents the state of the art in
  hardware-optimized quantum state transfer on a different optimization
  axis than ours.

### Our contribution: spatial dephasing profile optimization

In the ENAQT literature we surveyed, nobody optimizes WHERE the noise goes;
that work treats dephasing as a scalar parameter. Our sacrifice-zone formula
optimizes it spatially:

| Method | Source | Optimization target | Improvement |
|--------|--------|-------------------|------------|
| Uniform γ optimization | Plenio & Huelga 2008 | Scalar γ | ~2-3× |
| Coupling optimization (Bayesian) | IBM PST 2025 | J values | +8% |
| **Spatial γ formula (this work)** | [Resonant Return](../experiments/RESONANT_RETURN.md) | **Per-site γ** | **139-360×** (N = 5-9, simulation) |

The last column does not hold one quantity: the ENAQT figure is a gain in
transport efficiency, ours a ratio of peak created Sum-MI against a
V-shaped profile in the ε → 0 simulation (selective against uniform DD on
hardware: 1.4-3.2×), so the column does not rank the methods. What differs
is the knob: uniform optimization IS the right approach when you only have
one, and spatial dephasing profiles give you N, whose best allocation in
this simulation is extremely non-uniform (all noise on one edge).

---

## Related work: Pauli Hamiltonians as polynomials over GF(2)

Seen in abstracts and contents during an outside reading pass of [F173](proofs/PROOF_SIGNED_MASK_BALANCE.md) (2026-11-11), not read in full.

- **Haah, "Commuting Pauli Hamiltonians as maps between free modules"** (Commun. Math. Phys. 324, 2013, arXiv:1204.1063). A translation-invariant commuting Pauli Hamiltonian as a map of modules over the Laurent polynomials over GF(2), Pauli strings as columns of polynomials, its ground states on a torus read through that module structure. F173's Theorems B to D work in this formalism in one dimension and ask a different question of it: whether the relations among the placed masks carry an odd number of strings of odd #Y + #Z.
- **Repeated-root cyclic codes**: Castagnoli, Massey, Schoeller and von Seemann, "On repeated-root cyclic codes", and van Lint, "Repeated-root cyclic codes" (two papers, IEEE Trans. Inf. Theory 37(2), 1991). The factor (1 + x)^(2^a) of x^N + 1 and the codes it carries; the relations of one template on the ring of even N are a cyclic code of this kind (F173's Theorem C), those of several templates a quasi-cyclic one, and F173's Theorem D reads their sign count at x = 1.
- **Dehaene and De Moor, "Clifford group, stabilizer states, and linear and quadratic operations over GF(2)"** (Phys. Rev. A 68, 042318, 2003). The diagonal Clifford gates are phases of quadratic forms; the carrier of F173's Theorem E is the phase of a quadratic form, with values in {±1, ±i}.

---

## What the literature does NOT contain

This is our honest assessment of what appears to be new in this project.
"Appears" because we cannot be certain that nobody else has found these
results and published them in a form we did not find. If you know of
prior work on any of these, please open an issue.

1. **The dihedral group and the exact criterion.** The palindromic
   spectrum for XYZ couplings on arbitrary graphs under local Z dephasing
   is in Prosen 2012, and the explicit site-local operator Π in Sá,
   Ribeiro and Prosen 2023 (both above); Π was found here independently
   ([Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md)). What we did
   not find there or elsewhere: the dihedral group named and used as one
   object ([F118](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md); Sá, Ribeiro and
   Prosen's P without its sector phase and Q₊ generate it), and the two-end count (F158), which
   decides the palindrome of the eigenvalue multiset exactly for Hermitian
   jumps that square to one, at positive rates, at least one present. Haga et al. had the band
   grading but not the inter-band pairing.

2. **The Absorption Theorem** (Re(λ) = −2γ⟨n_XY⟩). Quantitative
   version of the incoherenton band picture. Not in Haga et al.
   See [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md).

3. **CΨ = ¼ as a boundary.** The Mandelbrot equivalence, the
   bifurcation mapping, the ¼ boundary as a quantum-classical transition
   marker. Not in any literature we found.

4. **CΨ diagnostic as such.** Concurrence × l₁-coherence as an
   observable product. Individual measures are standard; the product
   and its properties are not.

5. **Spatial dephasing profile optimization.** The ENAQT work we
   surveyed optimizes uniform scalar γ. Our sacrifice-zone
   formula (all noise on one edge qubit) optimizes per site, reaching
   139-360× in peak created Sum-MI at N = 5-9 in simulation.

6. **Crossing taxonomy**: three scalar-response classes for five named
   readouts in two Bell+ evolution books, six finite crossings and two never
   bridges, not physical observer classes or a consequence of the palindrome.
   The separate **subsystem** and **cockpit** studies report their finite
   readouts (3 cockpit observables capture 88-96% on the stated corpus).

**Assessment:** Either genuinely novel or an artifact of small systems
(N ≤ 8). The palindromic spectrum and Π are prior work (item 1); the
dihedral group and the two-end count are proven for all N; the other
results are verified computationally for N ≤ 9.

---

## Tested research directions

### Direction 1: Exceptional Point Connection (TESTED, NEGATIVE)

**Hypothesis:** CΨ = ¼ is (or correlates with) a Liouvillian
Exceptional Point.

**Result:** Tested March 13, 2026. Three sweeps (γ 3-qubit, J_SB
3-qubit, γ 2-qubit). EP_strength and eigenvalue gap tracked. NO
connection found. EP_strength follows γ monotonically, not CΨ. No
peak at ¼, no eigenvalue coalescence.

Script: [`simulations/ep_test.py`](../simulations/ep_test.py)

### Direction 2: Graph Symmetry Decomposition (TESTED, PARTIAL)

**Hypothesis:** Our c+/c- two-supermode structure is a special case of
graph symmetry-resolved Liouvillian decomposition.

**Result:** Tested March 13, 2026. XXX parity commutes with the
Liouvillian (weak symmetry confirmed). BUT c+ and c- both have
parity +1: they live in the SAME sector. The c+/c- split comes from
observable projection (even/odd), not Liouvillian symmetry sectors.

Script: [`simulations/symmetry_and_u_analysis.py`](../simulations/symmetry_and_u_analysis.py)

### Direction 3: Physical Interpretation of u = C(Ψ+R) (TESTED, RESOLVED)

**Result:** The Mandelbrot fixed point z* satisfies z*(1−z*) = CΨ.
This is the Bernoulli variance form p(1−p). CΨ ≤ ¼ is the trivial
maximum of Bernoulli variance. At CΨ = ¼, z* = ½ (maximum binary
uncertainty). z* correlates strongly with purity (r = 0.917).

The ¼ boundary is demystified: it is the upper bound of a quadratic
function. The remaining question: what binary process does z* represent
physically?

Script: [`simulations/symmetry_and_u_analysis.py`](../simulations/symmetry_and_u_analysis.py)

---

## Open threads

- Anomalous late-time coherence on IBM: non-Markovian EP?
- Operator feedback mechanism: connection to time-delayed feedback
- Hardware replication on second backend

---

*History: Created March 13, 2026 (Cowork Claude deep research).
Updated March 14 (incoherenton and η-pairing connections found, leading
to Π discovery same day). March 24 (ENAQT section added after
sacrifice-zone formula). April 5 (restructured: added Foundations
section with textbooks and original papers; full citations for Haga,
Medvedyeva-Essler-Prosen, Albert-Jiang, Buca-Prosen; NMR palindrome
connection; review articles for students).
April 11 (Gaztañaga-Kumar-Marto direct-sum QT and ER bridges section
added after discovering structural parallel to our Direct-Sum
Decomposition proof).*
