# Literature Review - March 2026

**Source:** Cowork Claude (Deep Research), March 13, 2026
**Purpose:** Literature connections, next research directions

---

## Key Papers

### MAIN HIT: Paper 4 - Graph Symmetry & EP Diagnostics
- **Title:** "Graph Symmetry Organizes Exceptional Dynamics in Open Quantum Systems"
- **Submitted:** March 11, 2026 (!) - [arXiv:2603.10654](https://arxiv.org/abs/2603.10654)
- **What it does:** Introduces symmetry-resolved approach to identify Exceptional Points
  from the full Liouvillian. Defines "EP-Strength" diagnostic E based on eigenvector
  conditioning. Applied to tight-binding models with correlated dephasing.
- **Why it matters:** Our two-supermode result (c+/c-) could be a special case of their
  graph symmetry decomposition. Their E diagnostic is directly applicable to test whether
  CΨ = 1/4 corresponds to an Exceptional Point.

### Paper 1: Emergent Liouvillian EPs from Exact Principles
- Khandelwal & Blasi, Quantum (2025) - [arXiv:2409.08100](https://arxiv.org/abs/2409.08100)
- EPs in Liouvillian survive beyond Markov approximation

### Paper 2: Encircling Liouvillian EPs (Review)
- Sun & Yi, AAPPS Bulletin (2024) - [arXiv:2408.11435](https://arxiv.org/abs/2408.11435)
- Chiral state transfer near EPs, experimentally observed in superconducting qubits

### Paper 3: Non-Markovian Quantum EPs
- Lin, Kuo, Lambert et al., Nature Comm. 16, 1289 (2025)
- Higher-order EPs invisible in Markovian limit. May explain our IBM anomalous late-time coherence.

### Paper 5: Blended Dynamics in Open Quantum Networks
- January 2026 - [arXiv:2601.14763](https://arxiv.org/abs/2601.14763)
- Classical-like clustering in quantum networks. Matches our signal processing view.

### Paper 6: Quaternionic Fractals (Viennot)
- Chaos, Solitons & Fractals (2022) - [arXiv:2003.02608](https://arxiv.org/abs/2003.02608)
- Fractal boundary between coherence and decoherence regimes (Mandelbulb).
  Conceptual sibling of our Mandelbrot boundary at CΨ = 1/4.

### Paper 7: Fractal Entangled Steady States
- Ippoliti, Rakovszky, Khemani, PRX 12, 011045 (2022)
- Fractal structures in quantum steady states via spacetime duality

### Paper 8: Time-delayed Coherent Feedback
- PRA 99, 053809 (2019) - [arXiv:1805.02317](https://arxiv.org/abs/1805.02317)
- Stabilizing coherence against dephasing. Related to our operator-feedback finding.

---

## What the literature does NOT contain

1. CΨ = 1/4 as a specific boundary (our Mandelbrot equivalence is new)
2. CΨ diagnostic as such (concurrence x l1-coherence as observable)
3. Crossing taxonomy (Type A/B/C)
4. Subsystem locality of crossing
5. Five-regulator framework

Assessment: Either genuinely novel or an artifact of small systems.
The EP connection could clarify this.

---

## Three prioritized research directions

### Direction 1: Exceptional Point Connection (TESTED - NEGATIVE)

**Hypothesis:** CΨ = 1/4 is (or correlates with) a Liouvillian Exceptional Point.

**Result:** Tested March 13, 2026. Three sweeps (gamma 3-qubit, J_SB 3-qubit,
gamma 2-qubit). EP_strength and eigenvalue gap tracked. NO connection found.
EP_strength follows gamma monotonically, not CΨ. No peak at 1/4, no eigenvalue
coalescence. The 2-qubit system showed CΨ_max = 1/3 for all gamma values.

CΨ = 1/4 remains algebraically exact but has no detected EP connection.
Script: simulations/ep_test.py

### Direction 2: Graph Symmetry Decomposition (TESTED - PARTIAL)

**Hypothesis:** Our c+/c- two-supermode structure is a special case of the
symmetry-resolved Liouvillian decomposition from Paper 4.

**Result:** Tested March 13, 2026. XXX parity commutes with the Liouvillian
(weak symmetry confirmed). BUT c+ and c- both have parity +1 - they live
in the SAME sector. All oscillatory eigenmodes are 50/50 mixed between sectors.
The c+/c- split comes from observable projection (even/odd), not Liouvillian
symmetry sectors. Graph symmetry applies but does NOT explain the two channels.
Script: simulations/symmetry_and_u_analysis.py

### Direction 3: Physical Interpretation of u = C(Psi+R) (TESTED - RESOLVED)

**Hypothesis:** u has an information-theoretic interpretation.

**Result:** Tested March 13, 2026. The Mandelbrot fixed point z* satisfies
z*(1-z*) = CPsi. This is the Bernoulli variance form p(1-p).
CPsi <= 1/4 is the trivial maximum of Bernoulli variance.
At CPsi = 1/4, z* = 1/2 (maximum binary uncertainty).
z* correlates strongly with purity (r=0.917) and anti-correlates with
von Neumann entropy (r=-0.838).

The 1/4 boundary is demystified: it is the upper bound of a quadratic function.
The remaining question: what binary process does z* represent physically?
Script: simulations/symmetry_and_u_analysis.py

---

## Secondary threads (for later)

- Anomalous late-time coherence on IBM: non-Markovian EP? (Paper 3)
- Operator feedback mechanism: connection to time-delayed feedback (Paper 8)
- Hardware replication on second backend

---

*Source: Cowork Claude deep research, March 13, 2026.*
*Integrated into R=CΨ² repo same day.*


---

## Update: March 14, 2026

A second, more focused literature search was conducted targeting 5 specific
directions (mirror symmetry, rate counting, CΨ=1/4, band structure, rational rates).
Full findings in `ClaudeTasks/LITERATURE_FINDINGS.md`.

Key outcome: the search identified the incoherenton framework (Haga et al. 2023)
and η-pairing symmetry (Medvedyeva-Essler-Prosen 2016) as closest prior work.
Neither had found the palindromic symmetry or the conjugation operator.

This directly led to the discovery and proof of Π on the same day.
See `docs/MIRROR_SYMMETRY_PROOF.md`.
