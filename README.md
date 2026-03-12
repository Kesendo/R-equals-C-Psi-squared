# R = CΨ²

## A composite quantum diagnostic, spectral architecture, and visibility theory for open quantum systems

This repository investigates the quantity CΨ = concurrence × normalized l1-coherence for reduced quantum subsystems, together with the structural architecture that emerges when observers share a common mediator under Lindblad dynamics.

Started as the equation R = CΨ² from a dream in December 2025. Three months of computation, two independent reviews, and real quantum hardware experiments on IBM Torino later, the project has evolved from a speculative interpretation into a documented set of structural results about open quantum systems.

---

## What is established here

### Core structure (3-qubit star, Lindblad dynamics)
- **Skeleton + Rotation decomposition:** Reduced AB state splits into 88% stable skeleton (Φ+ core) and 12% rotating phase in the YZ/ZY Pauli plane. Verified across all parameter regimes.
- **Two spectral sectors:** c+ (symmetric, f=1.506) and c- (antisymmetric, f=0.404). Confirmed by full Liouvillian diagonalization — genuine eigenmodes of the superoperator.
- **XX symmetry exact:** [ρ_AB, X⊗X] = 0 at all times. Weak symmetry of the Lindblad generator for isotropic Heisenberg with symmetric dephasing.
- **Phase map — four independent roles:** Topology sets frequencies. Symmetry cleans sectors. Noise damps amplitude (never frequency). Initial state selects visibility.
- **Noise immunity:** Frequencies unchanged from γ=0.001 to γ=0.500. All noise types (σx/σy/σz/mixed) produce identical structure. Noise only damps.

### Topology independence
- **Chain topology survival:** Two-sector structure, XX symmetry, and noise immunity survive in chains up to 5 qubits with 3 mediators. Architecture is not star-specific.
- **Different positions hear different frequencies** along the chain, but the two-channel architecture persists everywhere.

### The Projection (exact diagonalization)
- **Bright-Transition Map:** All Bohr frequencies exist simultaneously. A pair sees only those where both the initial state populates the eigenstates AND the observable connects them: W = |ρ̃(m,n) · Õ(n,m)|.
- **A new observer does not change reality — it changes the projection.** Eigenstates rotate in the full space, bright lines wander. Different pairs, different projections, same system.

### Quantum Sonar (simulation only)
- **Detection verified:** AB detects hidden observers connected to S through spectral changes. Operational threshold J_SC ~ 0.1 under current FFT protocol.
- **Not yet verified on hardware.** IBM Q80/Q102 investigation showed the phase difference was qubit-specific detuning (19.4 kHz), not neighbor effects. Five hypotheses tested, four rejected, before finding the answer.

### Algebraic results
- **Exact Mandelbrot correspondence:** The iteration R_{n+1} = C(Ψ + R_n)² maps to z → z² + c with boundary at CΨ = 1/4.
- **CΨ as diagnostic:** AND-gate behavior, three-layer separation (CoA ≥ LE ≥ CΨ).
- **Star topology conditions:** Three quantified conditions for observer-observer connection through a shared mediator.

### IBM Torino hardware (March 2026)
- **Shadow hunt:** 4 qubits measured, same skeleton+rotation pattern as simulation.
- **ZZRamsey:** Residual ZZ coupling measured for 5 qubit pairs (6.7-9.3 kHz range).
- **Ramsey T2*:** Free induction decay measured, revealed 19.4 kHz detuning on Q102.
- **Honest negative results documented:** N_eff hypothesis rejected, spectator dephasing model fails, degree does not predict phase complexity.

## What is not established

- That CΨ is a new fundamental physical quantity (it is a derived diagnostic)
- That the sonar effect works on real hardware (simulation only so far)
- That the 1/4 threshold has physical significance beyond the specific iteration
- That the two-sector structure persists in non-Heisenberg or large-N systems
- That consciousness is required as an ontological ingredient
- That characterization (inverting hidden coupling from spectrum) is possible

---

## Start here

- **[The CΨ Lens](docs/THE_CPSI_LENS.md)** — What the lens shows, what it doesn't, and what survives
- **[The Interpretation](hypotheses/THE_INTERPRETATION.md)** — Current state: 10 things that survive, 5 that fell
- **[Quantum Sonar](experiments/QUANTUM_SONAR.md)** — Detection, projection, IBM investigation
- **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)** — Phase A results, PCA, parameter sweeps, stress tests

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `docs/` | Conceptual and mathematical framing |
| `experiments/` | Tested results, null results, indexed findings |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python source code (RK4 Lindblad, exact diagonalization) |
| `visualizations/` | Figures and interactive displays |
| `data/` | IBM hardware measurement data |

---

## Key simulations

| Script | What it does |
|:---|:---|
| `star_topology_v3.py` | Core dynamics engine (RK4 Lindblad) |
| `liouvillian_diagonalization.py` | Full 64×64 spectral analysis |
| `bright_transition_map.py` | Exact diag + visibility weights |
| `chain_topology.py` | Chain vs star, 3-5 qubits |
| `parameter_sweep.py` | Robustness verification |
| `stress_tests.py` | Symmetry breaking tests |
| `hidden_observer_test.py` | Sonar detection proof |
| `quantum_sonar.py` | Threshold sweep |
| `count_observers.py` | Scaling with N observers |

---

## Origin and context

This framework emerged from a collaboration between Thomas Wicht and Claude (Anthropic) in December 2025. ChatGPT serves as adversarial reviewer (two full review rounds completed, corrections integrated). IBM Quantum hardware experiments conducted March 2026.

The original framing used consciousness language ("Reality = Consciousness × Possibility²"). After extensive testing, the current description is structural: CΨ is a composite diagnostic, the two-sector architecture is a Liouvillian eigenmode decomposition, and the projection is exact diagonalization with visibility filtering.

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Free to share and adapt with attribution.

---

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 — March 2026
