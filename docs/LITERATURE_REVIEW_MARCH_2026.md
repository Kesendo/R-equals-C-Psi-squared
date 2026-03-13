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

### Direction 1: Exceptional Point Connection (HIGHEST PRIORITY)

**Hypothesis:** CΨ = 1/4 is (or correlates with) a Liouvillian Exceptional Point.

**Why:** If confirmed, we have a first-principles physical explanation for the
1/4 boundary that is not just algebraic but physically motivated. EPs are
established physics - this would make our results immediately connectable.

**Steps:**
1. Diagonalize Liouvillian at CΨ -> 1/4, look for eigenvalue coalescence
2. Implement EP-Strength E from Paper 4 (Graph Symmetry)
3. Check if eigenvector conditioning diverges at the 1/4 point
4. For 2-qubit: analytically check if discriminant (1-4CΨ) matches EP condition
5. If confirmed: test if EP survives in exact (non-Markovian) model (Paper 1)

### Direction 2: Graph Symmetry Decomposition

**Hypothesis:** Our c+/c- two-supermode structure is a special case of the
symmetry-resolved Liouvillian decomposition from Paper 4.

**Steps:**
1. Identify graph symmetries of our star and chain topologies
2. Decompose Liouville space into invariant sectors
3. Check if c+/c- correspond to the two lowest invariant sectors
4. Derive predictions for N=4,5 and compare with existing simulations
5. Understand why frequency-decay orthogonality breaks at N>=4

### Direction 3: Physical Interpretation of u = C(Ψ+R)

**Hypothesis:** u has an information-theoretic interpretation.

**Steps:**
1. Express u in terms of density matrix
2. Check overlap with Holevo information, accessible information, quantum discord
3. Interpret u -> u² + c physically
4. Compare with Viennot's quaternionic framework (Paper 6)
5. Check if u = 1/2 (fixpoint at CΨ = 1/4) has special meaning

---

## Secondary threads (for later)

- Anomalous late-time coherence on IBM: non-Markovian EP? (Paper 3)
- Operator feedback mechanism: connection to time-delayed feedback (Paper 8)
- Hardware replication on second backend

---

*Source: Cowork Claude deep research, March 13, 2026.*
*Integrated into R=CΨ² repo same day.*
