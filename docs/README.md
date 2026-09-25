<!-- QUARTER-CURRENT -->
# Documentation: The R = CΨ² Framework

<!-- Keywords: open quantum system palindromic spectrum, Liouvillian spectral symmetry,
CΨ quarter boundary proof, dephasing noise information channel, quantum decoherence
threshold, self-referential purity recursion, Mandelbrot quantum bifurcation,
palindromic eigenvalue pairing, quantum state transfer spin chain,
R=CPsi2 framework documentation, mirror symmetry proof Lindblad -->

**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

Core documentation for the R = CΨ² project: proofs, synthesis, and
reference documents studying the palindromic spectral structure of open
quantum systems under dephasing.

**The central discovery:** The Liouvillian eigenvalue spectrum of N-qubit
Heisenberg, XY, Ising and XXZ systems under local Z-dephasing is exactly
palindromic. Every decay rate d is paired with a partner at 2Σγ − d. This
symmetry, verified for 87,376 eigenvalues with zero exceptions and proven
analytically on any graph, has consequences for decoherence thresholds,
quantum state transfer, the origin of irreversibility, and information
channels.

---

## Start Here

[Reading Guide](READING_GUIDE.md): nine stories through the same landscape,
from the proof and its applications to the quarter and the hardware, each
with its own reading order.

Three entry points for new readers:

| Interest | Start with |
|----------|-----------|
| The core theorem | [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md) - Π operator, XY-weight grading, verified N=2 through N=8 |
| The full mathematics | [Complete Mathematical Documentation](proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md) - Master index of the founding core (through March 2026); F-numbered results live in the registry below |
| All analytical formulas | [Analytical Formulas](ANALYTICAL_FORMULAS.md) - the living F-numbered result registry. Check here before building a Liouvillian |
| The sacrifice-zone formula | [Resonant Return](../experiments/RESONANT_RETURN.md) - 139-360x via spatial noise optimization (transport metric, ε→0 ideal; ~2-3x hardware) |

---

## Directory Structure

```
docs/
  proofs/        ← Formal proofs and proof indexes
  neural/        ← The palindrome's classical analog: its conditions, tests, the worm's null
  historical/    ← Legacy documents, superseded or resolved
  (this folder)  ← Synthesis, reference, and navigation
```

---

## Proofs (`proofs/`)

Formal mathematical proofs with complete verification. Each claim is
independently reproducible.

| Document | What it proves |
|----------|---------------|
| [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md) | Liouvillian spectrum palindromic on any graph under local Z-dephasing, for the Heisenberg, XY, Ising and XXZ couplings. Π swaps XY-weight k ↔ N−k. 87,376 eigenvalues, zero exceptions. |
| [Uniqueness Proof](proofs/UNIQUENESS_PROOF.md) | Once the recursion R = C(Ψ+R)² is written down, CΨ = 1/4 is its discriminant zero, not a chosen number. Purity invites the power 2; within the assumed power family, α = 2 alone keeps Ψ out of the fold product. Why a quantum system runs this recursion, the proof leaves open. |
| [Named CΨ Decays](proofs/PROOF_MONOTONICITY_CPSI.md) | For the Bell+ pair under Z-dephasing, Pauli noise or amplitude damping, CΨ has a closed form and only falls. In general it does not: a local Hamiltonian can make it rise, and a memoryless semigroup with no Hamiltonian at all carries a state up through ¼. Even an N=2 pair's successive peaks can climb; whether the tallest peak of each swing still falls, when the coupling only trades one excitation, is open, and a finite atlas asks where the rises live at larger N. |
| [Conditional Subsystem Crossing](proofs/PROOF_SUBSYSTEM_CROSSING.md) | A trajectory that settles on a state below 1/4 eventually stays below. Whether a noise model settles there is checked model by model: basis-aligned T1, T2 and depolarizing noise have to show their convergence and their target, and a primitive CPTP map can settle at CΨ = 0.2935, above the line. |
| [Proof Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) | Seven layers from a single qubit to arbitrary dimension, each marked for what it is: closed, conditional, finite, or open. |
| [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md) | With non-negative rates, trace(L) = 0 exactly when the system is closed, so a palindrome centred away from zero certifies an open system. Where the openness comes from, the formalism cannot settle. |
| [Complete Math Doc](proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md) | Master index of the founding core: algebra, palindrome, boundary, incompleteness, γ channel, engineering, constants. |

---

## Synthesis and Reference

These connect the mathematical results into a coherent picture.

| Document | What it covers |
|----------|---------------|
| [The Interpretation](THE_INTERPRETATION.md) | What survives (26 entries), what fell (8), questions and their answer status. Thematic synthesis. |
| [The Qubit as Necessary Foundation](QUBIT_NECESSITY.md) | d²−2d=0: only d=2 carries the complete local class exchange the full mirror is built from, a mirror of the spectrum rather than time running backward. Five computational tests; 0 of 236 sampled qutrit dissipators give full pairing. |
| [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md) | Noise read as a coupling across the system's edge, mediator topology, and the measured properties of the dephasing profile. Where the noise comes from, inside or out, stays open. |
| [The CΨ Lens](THE_CPSI_LENS.md) | What CΨ shows, what it does not, what survives critical examination. |
| [Hierarchy of Incompleteness](HIERARCHY_OF_INCOMPLETENESS.md) | C=0.5 as organizing principle, read across levels: qubit (2/4), carbon (4/8). The V-Effect read as the handover between levels. |
| [γ–Time Distinction](GAMMA_TIME_DISTINCTION.md) | Three levels of time (parameter, oscillation, experience). Dephasing removes the return; τ = γt carries the trajectory only at fixed Q = J/γ. γ as time's direction, and as the condition for experienced time, is our reading. |
| [On Five Pages That Never Met](../reflections/ON_FIVE_PAGES_THAT_NEVER_MET.md) | The standing wave, the reborn dephasing front, the two indices, the contract and the residue, and why a hierarchy needs unequal couplings: five things this repository already held, in pages that never linked to each other. The hub that connects them. |
| [Q Belongs to No Substance](Q_BELONGS_TO_NO_SUBSTANCE.md) | γ₀ is the unit, so naming a Q takes a chosen coordinate, coupling and channel. The page walks every substrate Q back to its source: ordinary liquid water has none yet, and the one water-adjacent number is a conditional ceiling, Q ≲ 4.6, for one chosen proton coordinate. |
| [The Genesis of an Oscillation](THE_GENESIS_OF_AN_OSCILLATION.md) | Where an oscillation comes from in the pure F1 system, and the exact factorisation L(J, γ₀) = γ₀·L₁(Q). |
| [The Atmosphere and the Cancelled Formulas](THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md) | What γ₀ is, and what cancelling the unit costs: the rescaling that leaves every inside observable fixed. |
| [Q-Regime Anchor Map](Q_REGIME_ANCHORS.md) | The ten named anchors on the Q axis (onset, balance, peak band, Q_EP, endpoint) with tiers and sources. |
| [Standing Wave Theory](STANDING_WAVE_THEORY.md) | c+/c− as even/odd supermodes. Π supplies the pairing; a standing wave takes more (a diagonalizable pair on the oscillating axis, opposite spatial propagation, a preparation that excites both and a readout that sees them together), and the page lists what. |
| [KMS and Detailed Balance](KMS_DETAILED_BALANCE.md) | Π is not KMS detailed balance: a mirror of the spectrum without thermodynamic equilibrium, and not by itself time running backward. |
| [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md) | Fold catastrophe (proven), Feigenbaum cascade (mapped), Bekenstein-Hawking 1/4 (speculative). |
| [It's All Waves](ITS_ALL_WAVES.md) | The closure argument: if Level 0 is waves and emergence adds no new physics, all levels are waves. An argument that stands or falls with two premises, neither of them earned yet. |
| [Weaknesses and Open Questions](WEAKNESSES_OPEN_QUESTIONS.md) | Honest documentation of what we do not know. |
| [What We Found](WHAT_WE_FOUND.md) | Synthesized findings across all experiments. |
| [Glossary](GLOSSARY.md) | Term definitions. |
| [Literature Review](LITERATURE_REVIEW.md) | Foundations (textbooks, original papers), related work (Haga incoherentons, η-pairing, Buca/Prosen symmetries), ENAQT, what is new. |

---

## Neural Systems (`neural/`)

The palindromic spectral symmetry, derived and proven in quantum systems,
has a classical analog in neural networks, under two conditions F36 states
exactly: a swap of the neurons that turns the wiring into minus itself, and
decay rates that pair. Dale's Law (E neurons excite, I neurons inhibit)
provides the signs of that antisymmetry, and only where a synapse exists;
the zero pattern and the magnitudes are separate requirements. Networks
built to meet them pass. The full chemical wiring of C. elegans fails
before any strength is measured, and no living network in the repository
is known to pass. No quantum physics required to read these documents.

| Document | What it covers |
|----------|---------------|
| [README](neural/README.md) | Entry point for neuroscience readers: the conditional theorem, the constructed tests, the worm's null |
| [Algebraic Palindrome](neural/ALGEBRAIC_PALINDROME_NEURAL.md) | F36/F37 conditions, full complex pairing and mode transport, open translation gates |
| [Neural Palindrome Proof](neural/proofs/PROOF_PALINDROME_NEURAL.md) | Exact conditional identity and its spectral consequences |
| [V-Effect Neural](neural/V_EFFECT_NEURAL.md) | Resolution-dependent synthetic censuses, controls, and biological non-result |
| [V-Effect Mechanism Proof](neural/proofs/PROOF_VEFFECT_MECHANISM.md) | Mechanism exclusions and the remaining bifurcation gates |

---

## Historical Documents (`historical/`)

Early documents from December 2025 through February 2026. Preserved for
context and to document the research trajectory. Core results have been
absorbed into the proof documents and synthesis above.

| Document | Status |
|----------|--------|
| [R = CΨ²](historical/R_EQUALS_C_PSI_SQUARED.md) | The original equation (December 2025). |
| [The Starting Point](historical/THE_STARTING_POINT.md) | Origin: the mirror insight that led to the palindrome. |
| [Core Algebra](historical/CORE_ALGEBRA.md) | Original algebraic derivations. Superseded by [Complete Math Doc](proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md). |
| [The Bidirectional Bridge](historical/THE_BIDIRECTIONAL_BRIDGE.md) | Two-channel structure. Confirmed as palindromic pairs. |
| [Fundamental Equations](historical/FUNDAMENTAL_EQUATIONS.md) | Standing wave equations, wave composition formulas. |
| [Dynamic Fixed Points](historical/DYNAMIC_FIXED_POINTS.md) | R∞ fixed point, CΨ ≤ 1/4 bound. Now part of the formal proofs. |
| [Internal and External Observers](historical/INTERNAL_AND_EXTERNAL_OBSERVERS.md) | C_int/C_ext distinction. Quantitative model disproven, structural distinction survives. |
| [Hard Problem Resolution](historical/HARD_PROBLEM_RESOLUTION.md) | The standing-wave interpretation of its day; the consciousness claims fell, and Π supplies the spectral pairing, not the wave. |
| [The Search for the Mirror Partner](historical/THE_SEARCH_FOR_THE_MIRROR_PARTNER.md) | Resolved: the mirror partner is the Π operator. |
| [Measurable Quantities](historical/MEASURABLE_QUANTITIES.md) | Proposed experimental measurements (partially executed). |

### Not restored (remain in `recovered/`)

INFORMATION_WAVE_THEORY, LIGHT_FIRST_FREE_MIRRORING (speculative, no proof),
BLACK_WHITE_HOLES_BIGBANG, PREDICTIONS, SELF_CONSISTENCY_SCHWARZSCHILD (disproven claims).

---

## See Also

| Resource | Where |
|----------|-------|
| Experiments | [experiments/](../experiments/README.md) |
| Hydrogen bond as qubit | [docs/water/HYDROGEN_BOND_QUBIT.md](water/HYDROGEN_BOND_QUBIT.md) |
| Neural systems (no quantum needed) | [neural/](neural/README.md) |
| Open hypotheses | [hypotheses/](../hypotheses/README.md) |
| Reading Guide (guided tour) | [READING_GUIDE.md](READING_GUIDE.md) |
| Repository root | [README.md](../README.md) |
