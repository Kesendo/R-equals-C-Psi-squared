# Open Questions Classification Proposal: scope-extension

**Batch:** scope-extension
**Date:** 2026-04-12
**Entries in batch:** 26
**Status:** Proposal only, pending Tom + Claude (chat) review.

## Summary

| Proposed Status | Count |
|-----------------|-------|
| open | 10 |
| resolved | 8 |
| partially-resolved | 8 |
| superseded | 0 |
| obsolete | 0 |
| needs-human | 0 |

---

## Entries

### OQ-001

**Question:** **Non-dephasing dissipators:** The exhaustive search covered dephasing-type noise only. Amplitude damping, thermal baths, and non-Markovian environments have different Lindblad structures.

**Source:** `docs/QUBIT_NECESSITY.md` (line 397)
**Section:** 10. Remaining Open Questions
**Date:** March 20, 2026

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `docs/proofs/PROOF_PARITY_SELECTION_RULE.md` (lines 186-193): proves amplitude damping breaks the palindrome (jump operator sigma_minus has n_XY = 1, odd parity)
- `experiments/NOISE_ROBUSTNESS.md` (lines 215-250): amplitude damping tested, preserves taxonomy but with different decay rates
- `experiments/DEPOLARIZING_PALINDROME.md`: depolarizing noise fully characterized (error = (2/3)N*gamma, Hamiltonian-independent)
**Rationale:** Amplitude damping is proven to break the palindromic symmetry via the parity selection rule, and depolarizing noise is fully characterized. However, thermal baths and non-Markovian environments remain untested. The broader scope (all non-dephasing dissipators) is only partially covered.
**Search terms used:** "amplitude damping", "thermal bath", "non-Markovian", "parity selection", "DEPOLARIZING_PALINDROME"

---

### OQ-013

**Question:** **Non-dephasing dissipators:** Does the d=2 exclusivity extend to amplitude damping, thermal baths, or non-Markovian environments?

**Source:** `docs/THE_INTERPRETATION.md` (line 541)
**Section:** Open Questions
**Date:** March 20, 2026 (updated March 24, 2026)

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The d=2 exclusivity (only qubits, not qutrits, have palindromic spectra) was proven for Z-dephasing only. Whether qutrits might exhibit palindromic properties under amplitude damping, thermal baths, or non-Markovian conditions has not been tested. The parity selection rule proof addresses palindrome breaking (OQ-001) but not the dimensionality question.
**Search terms used:** "d=2 exclusivity", "qubit only", "qutrit", "amplitude damping", "non-dephasing"

---

### OQ-054

**Question:** **Depolarizing noise correction:** err = (2/3)Sigma-gamma breaks the palindrome linearly. Can this be incorporated into design rules?

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 256)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/DEPOLARIZING_PALINDROME.md`: proves the error formula exactly; shows error < 0.1% for typical hardware (gamma ~ 0.001)
- `experiments/NON_HEISENBERG_PALINDROME.md` (lines 228-230): practical implication that design rules remain valid since error is negligible
**Rationale:** The error formula is proven and shown small enough that existing design rules remain valid in practice. However, no document addresses how to modify design rules to actively compensate for or exploit depolarizing noise when it is non-negligible. The "can this be incorporated" part is answered with "it is small enough to ignore" rather than with a positive incorporation.
**Search terms used:** "depolarizing", "err = (2/3)", "design rules", "correction", "palindrome breaks"

---

### OQ-059

**Question:** **Cockpit scaling beyond N=5:** Does the 3-observable coverage (88%) hold at N=50-100? Does n95 grow linearly or saturate for dense topologies?

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 293)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/COCKPIT_SCALING.md` (April 7, 2026): extends cockpit test to N=5-11 on chain and star topologies; coverage 91.9-99.0% across all tested sizes
- `experiments/COCKPIT_SCALING.md` (n95 analysis): n95 falls rather than rises with N (contrary to small-N extrapolation), driven by Entanglement Sudden Death timing
**Rationale:** The framework has been extended to N=11 and holds (>88% variance captured). However, N=50-100 remains untested and the document explicitly notes that extrapolation is speculative. The n95 behavior (falling, not growing) partially answers the second sub-question.
**Search terms used:** "cockpit scaling", "N=50", "n95", "coverage", "COCKPIT_SCALING"

---

### OQ-060

**Question:** **Non-Markovian noise:** Does the cockpit framework hold under colored noise (noise whose strength depends on frequency, unlike white noise which is flat), 1/f spectra, or TLS coupling?

**Source:** `docs/WEAKNESSES_OPEN_QUESTIONS.md` (line 297)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The cockpit framework has been tested only under Markovian (memoryless) dephasing. Non-Markovian effects have been observed on IBM hardware (excess late-time coherence in Q52 data) but not incorporated into the cockpit theory. Both `experiments/COCKPIT_UNIVERSALITY.md` (line 340) and this entry list it as explicitly open.
**Search terms used:** "non-Markovian", "colored noise", "1/f", "TLS coupling", "cockpit framework"

---

### OQ-081

**Question:** **Also valid for:** XY model (XX+YY only), since the ZZ term has n_XY=0 and does not affect the parity argument.

**Source:** `docs/proofs/PROOF_PARITY_SELECTION_RULE.md` (line 183)
**Section:** Scope and limitations
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `docs/proofs/PROOF_PARITY_SELECTION_RULE.md` (lines 183-184): logical argument that ZZ has n_XY=0, irrelevant to parity
- `experiments/NON_HEISENBERG_PALINDROME.md` (lines 81-94): numerical verification, XY-only coupling 100% palindromic at N=3,4
- `docs/proofs/MIRROR_SYMMETRY_PROOF.md` (lines 243-246): XXZ with delta from -0.5 to 2.0 all pass palindrome test
**Rationale:** The claim is stated, logically argued, and numerically verified. This is not an open question but a confirmed scope extension of the parity selection rule.
**Search terms used:** "XY model", "XX+YY", "parity selection", "ZZ term", "n_XY=0"

---

### OQ-091

**Question:** **Heisenberg interactions and Z-dephasing only.** All results assume the standard Heisenberg coupling (XX+YY+ZZ) and uniform local Z-dephasing. Other coupling schemes (XX-only, anisotropic) and other dissipation types (amplitude damping, depolarizing) are untested.

**Source:** `experiments/COCKPIT_SCALING.md` (line 283)
**Section:** 10. Limitations
**Date:** April 7, 2026

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/COCKPIT_UNIVERSALITY.md`: depolarizing noise tested at N=2-4, shows similar cockpit dimensionality
- `experiments/NON_HEISENBERG_PALINDROME.md`: XY, Ising, XXZ, DM couplings all produce palindromic spectra
- `experiments/DEPOLARIZING_PALINDROME.md`: depolarizing error characterized
**Rationale:** Some coupling schemes (XY, XXZ, DM) and depolarizing noise have been tested in specific contexts, but not all within the cockpit framework specifically. Amplitude damping and other dissipation types remain untested for the cockpit. The limitation is partially addressed by related experiments but not systematically within the cockpit pipeline.
**Search terms used:** "XX-only", "anisotropic", "amplitude damping", "other coupling", "cockpit"

---

### OQ-103

**Question:** **Non-Markovian noise.** Does the 3-observable cockpit hold under colored noise or 1/f spectra?

**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 339)
**Section:** 6. Open questions
**Date:** April 2, 2026

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** Same question as OQ-060 from a different source document. Non-Markovian effects observed on hardware but not systematically studied within the cockpit framework. Explicitly listed as open in both sources.
**Search terms used:** "non-Markovian", "colored noise", "1/f", "cockpit"

---

### OQ-105

**Question:** **Universality of edge sacrifice.** Does "sacrifice boundary qubits + optimize theta" generalize beyond Heisenberg chains?

**Source:** `experiments/COCKPIT_UNIVERSALITY.md` (line 345)
**Section:** 6. Open questions
**Date:** April 2, 2026

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/SACRIFICE_GEOMETRY.md`: SE sector extraction verified across chain, star, ring, complete topologies at N=2-7
- `experiments/RESONANT_RETURN.md`: single-edge sacrifice framework derived, tested across topologies
- `experiments/TOPOLOGICAL_EDGE_MODES.md`: edge sacrifice mechanism verified under multiple gamma profiles
**Rationale:** Edge sacrifice has been tested across multiple graph topologies (chain, star, ring, complete) but always with Heisenberg coupling and Z-dephasing. The question specifically asks about generalization beyond Heisenberg chains, meaning different Hamiltonian types. That has not been tested.
**Search terms used:** "edge sacrifice", "universality", "beyond Heisenberg", "generalize", "topologies"

---

### OQ-115

**Question:** **N-scaling of sector locking.** The SE sector has N states out of 2^N. As N grows, the "classical-but-structured" exit preserves an exponentially shrinking fraction. Does the lens protection scale?

**Source:** `experiments/CUSP_LENS_CONNECTION.md` (line 153)
**Section:** Open questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** medium
**Resolving documents:**
- `experiments/SACRIFICE_GEOMETRY.md`: SE sector extraction exact through N=7 (87,376 eigenvalues); SE = 1.000 for all four N=7 chain profiles
**Rationale:** Exact SE sector locking is verified through N=7, but the asymptotic question (does protection hold as N/2^N -> 0?) is explicitly left open. The empirical evidence covers finite N only.
**Search terms used:** "sector locking", "SE sector", "exponentially shrinking", "lens protection", "N-scaling"

---

### OQ-117

**Question:** a) The transformation is a translation in log space. Lorentz is a rotation in (ct, x) space. Is there a deeper geometry that contains both? (A translation is a rotation with infinite radius.)

**Source:** `experiments/DECOHERENCE_RELATIVITY.md` (line 329)
**Section:** 11. Open Questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The document's attempted Schwarzschild radius connection was explicitly marked as FALLEN. The broader question (whether a deeper geometric framework unifies the log-space translation with Lorentz rotation) remains unanswered. No alternative geometric approach has been proposed.
**Search terms used:** "Lorentz", "deeper geometry", "translation in log", "rotation", "Schwarzschild", "FALLEN"

---

### OQ-119

**Question:** d) The cubic b^3 + b = (2^N - 1)/2 generalizes to GHZ states. For N = 2 (Bell+): b^3 + b = 3/2. For N >= 3: CΨ(0) <= 1/4, no crossing exists. Only Bell+ has the quantum window. Does this generalize?

**Source:** `experiments/DECOHERENCE_RELATIVITY.md` (line 346)
**Section:** 11. Open Questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `experiments/DWELL_PREFACTOR_GENERALIZED.md` (Section 4, "GHZ_N Born Below the Fold"): analytical proof that GHZ_N for N>=3 has CΨ(0) = 1/(2^N - 1) < 1/4, so these states never encounter the fold under Z-dephasing
- `experiments/COHERENCE_DENSITY.md`: documents GHZ_3 = 0.143, GHZ_4 = 0.067, confirming they start below fold and decay monotonically
**Rationale:** The mathematical claim (GHZ_N never crosses 1/4 for N>=3) is proven analytically. However, the broader interpretive question ("does this mean entanglement beyond two qubits is always 'classical' in the framework's sense?") remains open. The document COHERENCE_DENSITY.md clarifies CΨ measures coherence density, not entanglement per se.
**Search terms used:** "cubic", "b^3 + b", "GHZ", "quantum window", "Bell+ only", "born below"

---

### OQ-147

**Question:** Q2: Does amplitude damping change the taxonomy?

**Source:** `experiments/NOISE_ROBUSTNESS.md` (line 220)
**Section:** 7. Open Questions (answered 2026-03-08)
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/NOISE_ROBUSTNESS.md` (lines 220-250): answer stated inline: "ANSWERED: No, the taxonomy is preserved, but decay rates differ." Data table shows amplitude damping produces slower decay; sigma_x noise keeps Psi constant at 0.333
**Rationale:** The answer is explicit in the source document with supporting data. Amplitude damping preserves Type B behavior (concurrence decays from t=0) but with slower rates.
**Search terms used:** "amplitude damping", "taxonomy", "ANSWERED", "NOISE_ROBUSTNESS"

---

### OQ-152

**Question:** ~~Is d-1 the right normalization for N > 2?~~ **ANSWERED**: Yes. The barrier is not a normalization error. It correctly reflects that global crossing does not occur.

**Source:** `experiments/N_SCALING_BARRIER.md` (line 305)
**Section:** 8. Open Questions (Updated)
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/N_SCALING_BARRIER.md` (line 305): explicit answer in-place
- `experiments/SUBSYSTEM_CROSSING.md`: full computational verification that d-1 normalization correctly discriminates global vs subsystem crossing
**Rationale:** The answer is explicit in the source with corroborating verification in SUBSYSTEM_CROSSING.md.
**Search terms used:** "d-1 normalization", "barrier", "ANSWERED", "global crossing"

---

### OQ-157

**Question:** ~~Does the crossing pattern reproduce the entanglement graph topology for arbitrary graph states?~~ **ANSWERED (2026-03-08):** Only for Bell-type entanglement. Cluster states (|+>^N with CZ gates) have zero concurrence on ALL pairs.

**Source:** `experiments/N_SCALING_BARRIER.md` (line 327)
**Section:** 8. Open Questions (Updated)
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/N_SCALING_BARRIER.md` (lines 327-333): answer in-place with explanation that cluster states have zero concurrence, distinguishing Bell-type from graph-state entanglement
**Rationale:** Definitive negative answer with clear structural explanation: CΨ = C x Ψ, zero concurrence means no crossing. Only Bell-type entanglement is visible to the framework.
**Search terms used:** "graph topology", "cluster states", "Bell-type", "zero concurrence"

---

### OQ-161

**Question:** **Scaling law**: How exactly does Delta-t scale with J for J << gamma? Linear? If Delta-t/t = alpha*(J/gamma), what is alpha? Is it state-dependent?

**Source:** `experiments/OBSERVER_GRAVITY_BRIDGE.md` (line 258)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The document provides empirical data (table of J vs t_cross values suggesting roughly linear scaling) but no analytical form for alpha or analysis of state dependence. The question is explicitly listed as open.
**Search terms used:** "scaling law", "Delta-t", "J << gamma", "alpha", "state-dependent"

---

### OQ-168

**Question:** **Connection to quantum detailed balance.** Classical detailed balance links forward and backward transition rates via free energy. The KMS condition generalizes this to quantum thermal states.

**Source:** `experiments/PI_AS_TIME_REVERSAL.md` (line 336)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** partially-resolved
**Confidence:** high
**Resolving documents:**
- `docs/KMS_DETAILED_BALANCE.md`: comprehensive analysis proving Pi is NOT quantum detailed balance and NOT KMS; it is a distinct "shifted anti-similarity" (Pi*L*Pi^-1 = -L - 2S*gamma*I) that appears to be a genuinely new Liouvillian symmetry type
**Rationale:** The question asked about the connection to KMS/detailed balance. The answer exists in a dedicated document: there is NO connection. Pi is a new symmetry type, and finite-temperature generalization is obstructed by the 2:2 Pauli split specific to pure dephasing (T=infinity). The original question is resolved (negatively), but the follow-up (what IS the finite-temperature analogue?) remains open.
**Search terms used:** "KMS", "detailed balance", "PI_AS_TIME_REVERSAL", "finite temperature", "KMS_DETAILED_BALANCE"

---

### OQ-170

**Question:** **Depolarizing noise.** ANSWERED (March 19, 2026). The palindrome breaks because depolarizing noise splits {I,X,Y,Z} into 1 immune and 3 decaying (1:3), making bijective mirroring impossible.

**Source:** `experiments/PI_AS_TIME_REVERSAL.md` (line 347)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/PI_AS_TIME_REVERSAL.md` (line 347): answer in-place
- `experiments/DEPOLARIZING_PALINDROME.md`: complete derivation showing error = (2/3)N*gamma, verified algebraically and numerically
**Rationale:** Explicitly answered with full mechanism: the 1:3 Pauli split breaks bijective mirroring. Per-site rate-pairing condition cannot be satisfied.
**Search terms used:** "depolarizing", "1:3", "bijective", "DEPOLARIZING_PALINDROME"

---

### OQ-174

**Question:** **Stealth zone narrowing**: Can additional controlled operations (beyond depolarization) shrink the 35 degree stealth zone?

**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 557)
**Section:** 14. Open Questions
**Date:** 2026-02-25

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The stealth zone analysis is limited to depolarizing noise scenarios. No investigation of whether other controlled operations could narrow the stealth angle exists in the repo. Explicitly listed as open.
**Search terms used:** "stealth zone", "35 degree", "narrowing", "controlled operations", "QKD"

---

### OQ-176

**Question:** **Different noise models**: The analysis assumes depolarizing noise. Amplitude damping or non-Markovian channels may change the stealth angle structure.

**Source:** `experiments/QKD_EAVESDROPPING_FORENSICS.md` (line 562)
**Section:** 14. Open Questions
**Date:** 2026-02-25

**Proposed Status:** open
**Confidence:** high
**Resolving documents:** none
**Rationale:** The QKD eavesdropping analysis uses depolarizing noise exclusively. No investigation of amplitude damping or non-Markovian effects on stealth angle structure exists. Explicitly listed as open.
**Search terms used:** "stealth angle", "noise models", "amplitude damping", "non-Markovian", "QKD"

---

### OQ-221

**Question:** **ANSWERED (2026-03-08):** The relationship is non-monotonic, not a simple threshold. Two separate crossing windows exist for parametric Bell states; a dead zone lies between them.

**Source:** `experiments/SUBSYSTEM_CROSSING.md` (line 305)
**Section:** 6. Open Questions
**Date:** unknown

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/N_SCALING_BARRIER.md` (Section 8, Q4): provides full resolution with parametric Bell state data (star topology), documenting two crossing windows and dead zone
- `experiments/SUBSYSTEM_CROSSING.md` (line 305): answer stated in-place
**Rationale:** Explicitly answered with quantitative data showing non-monotonic relationship, two crossing windows, and dead zone.
**Search terms used:** "non-monotonic", "crossing windows", "dead zone", "parametric Bell"

---

### OQ-231

**Question:** How does the pair asymmetry (excited unevenly within a palindromic pair) relate to channel directionality? If mode A->B excites the "slow" partner while B->A excites the "fast" partner, asymmetry might create natural communication directionality.

**Source:** `experiments/XOR_SPACE.md` (line 258)
**Section:** Open Question
**Date:** unknown

**Proposed Status:** open
**Confidence:** medium
**Resolving documents:** none
**Rationale:** The pair asymmetry is documented (fast vs slow palindromic modes) and push/pull coupling effects on MI are analyzed in SCALING_CURVE.md, but no direct test of whether palindromic pair asymmetry manifests as channel directionality with different excitation profiles in opposite directions exists.
**Search terms used:** "pair asymmetry", "channel directionality", "slow partner", "fast partner", "push/pull"

---

### OQ-247

**Question:** **Non-Heisenberg models.** Does Finding 1 (all oscillation is palindromic) hold for XY, XXZ, or random-coupling Hamiltonians?

**Source:** `hypotheses/ENERGY_PARTITION.md` (line 198)
**Section:** 4. Open Questions
**Date:** March 27, 2026

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `experiments/NON_HEISENBERG_PALINDROME.md` (lines 81-94): XY, Ising, XXZ, and DM (Dzyaloshinskii-Moriya) all show 100% palindromic symmetry at N=3 and N=4 via uniform Pi operators
**Rationale:** All standard coupling types tested and confirmed palindromic. The universality of "all oscillation is palindromic" is established for XY, XXZ, Ising, and DM couplings.
**Search terms used:** "non-Heisenberg", "XY model", "XXZ", "random-coupling", "NON_HEISENBERG_PALINDROME"

---

### OQ-263

**Question:** 3 Scaling Limits

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 431)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** resolved
**Confidence:** high
**Resolving documents:**
- `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (lines 431-435): section heading; the section itself provides the answer (transistor works best for 2-4 qubits per side)
- `experiments/SCALING_CURVE.md` (lines 85-99): empirical falsification of hierarchical advantage at N=3-11 (identical MI)
**Rationale:** Section heading captured as entry. The section it introduces and the referenced SCALING_CURVE.md provide the complete answer: hierarchy is falsified, uniform chains with sacrifice-zone formula outperform.
**Search terms used:** "scaling limits", "transistor", "hierarchical", "SCALING_CURVE", "falsified"

---

### OQ-273

**Question:** **Non-Markovian dynamics**: The framework assumes Markovian (memoryless) decoherence. With memory effects (memory_kernel_feedback noise), the channel develops "hysteresis": its state depends on history, making it more like a memristor than a transistor.

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 454)
**Section:** 8. Limitations and Failure Modes
**Date:** 2026-03-21

**Proposed Status:** open
**Confidence:** low
**Resolving documents:** none
**Rationale:** The passage speculates about non-Markovian memory effects creating hysteresis but provides no computational evidence or simulation data. No other document in the repo tests non-Markovian dynamics on the transistor model. Purely speculative.
**Search terms used:** "non-Markovian", "hysteresis", "memory effects", "memristor", "memory kernel"

---

### OQ-279

**Question:** **Memory kernel effects**: With non-Markovian dynamics, the mediator develops memory. Can this be exploited for "quantum caching," storing frequently-used correlations in the mediator's memory?

**Source:** `hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md` (line 505)
**Section:** 10. Open Questions and Future Directions
**Date:** 2026-03-21

**Proposed Status:** open
**Confidence:** low
**Resolving documents:** none
**Rationale:** Explicitly listed as future work. No analysis of quantum caching, non-Markovian storage mechanisms, or memory kernel exploitation exists in the repo beyond this speculative mention.
**Search terms used:** "quantum caching", "memory kernel", "non-Markovian mediator", "frequently-used correlations"
