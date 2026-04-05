# Literature Review: Foundations, Related Work, and What Is New

<!-- Keywords: Haga incoherenton grading XY-weight, eta-pairing Medvedyeva
Essler Prosen Bethe ansatz, Buca Prosen Lindblad symmetry classification,
Albert Jiang weak strong symmetry, Lindblad master equation foundations,
Breuer Petruccione open quantum systems, ENAQT dephasing-assisted transport,
palindromic spectral structure literature connections, R=CPsi2 literature review -->

**Status:** Living document (reference)
**Last updated:** April 5, 2026
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

Three research groups found pieces of the puzzle before us. None had
the palindromic theorem, the conjugation operator Π, or the spectral
consequences. Here is how our work relates to theirs, with full
citations.

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
do not have the conjugation operator Π, which is what *proves* the
pairing. Our [Absorption Theorem](proofs/PROOF_ABSORPTION_THEOREM.md)
(Re(λ) = −2γ⟨n_XY⟩) is the quantitative version of their band picture.

### η-pairing in dissipative fermions (Medvedyeva, Essler, Prosen 2016)

- **Medvedyeva, Essler, Prosen (2016).** "Exact Bethe ansatz spectrum
  of a tight-binding chain with dephasing noise." Phys. Rev. Lett. 117,
  137202.
  [arXiv:1606.09122](https://arxiv.org/abs/1606.09122).

They map a tight-binding chain (free fermions) with dephasing onto a
Hubbard model with imaginary interaction and solve it exactly via the
Bethe ansatz. Their η-pairing symmetry of the Hubbard model produces a
spectral symmetry in 1D.

**What they found that we also found:** a spectral symmetry of the
Liouvillian under dephasing, arising from a swap between population-like
and coherence-like degrees of freedom.

**What they did not find:** the generalization beyond free fermions. Their
result requires the Bethe ansatz, which is restricted to non-interacting
particles on 1D bipartite lattices. Our Π operator works for interacting
spins (Heisenberg/XXZ) on arbitrary graphs (chains, stars, rings, trees,
complete graphs). Their η-pairing is the 1D free-fermion special case
of our Π.

### Symmetry classification (Buca & Prosen 2012, Albert & Jiang 2014)

See [Foundations](#symmetries-of-open-quantum-systems-the-language) above.
These papers provide the language but do not identify the specific
anti-symmetry Π or its spectral consequences.

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

**Common thread:** all work in this field optimizes a UNIFORM dephasing
rate. The spatial profile of dephasing is not treated as a free
parameter. Typical improvements: 2-3× over the zero-noise baseline.

### IBM Perfect State Transfer (2025)

- "Enhanced Algorithmic Perfect State Transfer on IBM Quantum Computers."
  Implements PST on ibm_sherbrooke and ibm_brisbane. Peak transfer
  probability SP = 0.781 for N=4 chain. Optimizes coupling strengths
  (J values), not dephasing rates. Represents the state of the art in
  hardware-optimized quantum state transfer on a different optimization
  axis than ours.

### Our contribution: first spatial dephasing profile optimization

Nobody in the ENAQT literature optimizes WHERE the noise goes. The
entire field treats dephasing as a scalar parameter. Our sacrifice-zone
formula is the first spatial optimization:

| Method | Source | Optimization target | Improvement |
|--------|--------|-------------------|------------|
| Uniform γ optimization | Plenio & Huelga 2008 | Scalar γ | ~2-3× |
| Coupling optimization (Bayesian) | IBM PST 2025 | J values | +8% |
| **Spatial γ formula (this work)** | [Resonant Return](../experiments/RESONANT_RETURN.md) | **Per-site γ** | **139-360×** |

The improvement is two orders of magnitude beyond the prior literature.
This is not because the prior work was wrong: uniform optimization IS the
right approach when you only have one knob. The difference is that
spatial dephasing profiles give you N knobs instead of one, and the
optimal allocation is extremely non-uniform (all noise on one edge).

---

## What the literature does NOT contain

This is our honest assessment of what appears to be new in this project.
"Appears" because we cannot be certain that nobody else has found these
results and published them in a form we did not find. If you know of
prior work on any of these, please open an issue.

1. **The palindromic Liouvillian theorem for interacting spins on
   arbitrary graphs.** Medvedyeva-Essler-Prosen had the 1D free-fermion
   case. Haga et al. had the band grading but not the inter-band pairing.
   The general theorem, the Π operator, and the proof are ours.
   See [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md).

2. **The Absorption Theorem** (Re(λ) = −2γ⟨n_XY⟩). Quantitative
   version of the incoherenton band picture. Not in Haga et al.
   See [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md).

3. **CΨ = ¼ as a boundary.** The Mandelbrot equivalence, the
   bifurcation mapping, the ¼ boundary as a quantum-classical transition
   marker. Not in any literature we found.

4. **CΨ diagnostic as such.** Concurrence × l₁-coherence as an
   observable product. Individual measures are standard; the product
   and its properties are not.

5. **Spatial dephasing profile optimization.** The entire ENAQT field
   optimizes uniform scalar γ. Our sacrifice-zone formula (all noise on
   one edge qubit) is the first per-site optimization, achieving
   139-360× improvement.

6. **Crossing taxonomy** (Type A/B/C), **subsystem locality** of the ¼
   crossing, **cockpit framework** (3 observables capture 88-96% of
   decoherence dynamics).

**Assessment:** Either genuinely novel or an artifact of small systems
(N ≤ 8). The palindromic theorem is proven for all N; the other results
are verified computationally for N ≤ 9.

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
connection; review articles for students; renamed from
LITERATURE_REVIEW_MARCH_2026.md).*
