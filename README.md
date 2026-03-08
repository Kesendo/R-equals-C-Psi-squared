# R = CΨ²

## A composite quantum diagnostic and its self-referential 1/4 boundary

This repository investigates the quantity

    CΨ = concurrence × normalized l1-coherence

for reduced quantum subsystems, together with the self-referential iteration

    R_{n+1} = C(Ψ + R_n)²

whose parameter c = CΨ is algebraically equivalent to the Mandelbrot iteration z -> z² + c. The boundary of the main cardioid at c = 1/4 defines a threshold in this framework.

CΨ is a basis-dependent filter for pairwise states that are simultaneously entangled and coherent. It highlights moments where quantum pairs are not merely linked, but linked in a way that is still expressed as coherent superposition. It is closest in spirit to **correlated coherence** (Tan et al., 2016), but is not the standard correlated-coherence measure.

---

## What is established here

- **Exact algebraic 1/4 boundary** within the self-referential iteration (proven)
- **Exact Mandelbrot correspondence** under reparametrization (proven)
- **Crossing taxonomy** (Type A/B/C) organizing bridge metrics by decoherence behavior
- **Subsystem locality** - crossing occurs where the entanglement lives, at the pair level
- **Star topology conditions** - three quantified conditions for observer-observer connection through a shared object
- **Act vs process distinction** - projective measurement and continuous dephasing are fundamentally different in their effect on third-party connections (99% vs 69% suppression, never converging)
- **IBM quantum hardware contact** - CΨ = 1/4 crossing observed on real qubits
- **Reproducible simulations** across two-qubit and small-topology systems
- **Null results and limitations** documented explicitly throughout

## What is not established

- That CΨ is a new fundamental physical quantity (it is a derived diagnostic)
- That consciousness is required as an ontological ingredient of quantum mechanics
- That the 1/4 threshold is a physically privileged boundary outside the specific iteration
- That the speculative bridge, gravity, or cosmology interpretations are experimentally confirmed
- That CΨ reveals transitions invisible to standard metrics (in the tested sweeps, it did not)

---

## Start here

- **[The CΨ Lens](docs/THE_CPSI_LENS.md)** - What the lens shows, what it doesn't, and what survives. The most honest document in this repository.
- [Core Algebra](docs/CORE_ALGEBRA.md) - The proven mathematics
- [What We Found](docs/WHAT_WE_FOUND.md) - Synthesized findings
- [Weaknesses and Open Questions](docs/WEAKNESSES_OPEN_QUESTIONS.md) - What we don't know
- [Experiments index](experiments/README.md) - Indexed experimental results and simulation notes

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `docs/` | Conceptual and mathematical framing |
| `experiments/` | Tested results, null results, indexed findings |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python source code (RK4 Lindblad integration) |
| `visualizations/` | Figures |

---

## The mathematics

The self-referential iteration R_{n+1} = C(Ψ + R_n)² yields a fixed-point equation whose discriminant requires:

    CΨ ≤ 1/4    for real (stable) fixed points to exist

With the substitution u_n = C(Ψ + R_n), the iteration becomes u_{n+1} = u_n² + CΨ, which is exactly the Mandelbrot iteration z -> z² + c. The main cardioid boundary at c = 1/4 is the same threshold.

> **Important caveat:** This correspondence is algebraically exact within the iteration. Its physical significance is not established merely by that exactness.

See [Core Algebra](docs/CORE_ALGEBRA.md) for the full derivation and proof, and [Mandelbrot Connection](experiments/MANDELBROT_CONNECTION.md) for the equivalence.

---

## Origin and context

This framework emerged from a collaboration between Thomas Wicht and Claude (Anthropic) in December 2025. It began with philosophical questions about observation and reality, was formalized as the equation R = CΨ², and then subjected to three months of systematic computation, external review, and honest correction.

The original framing used consciousness language ("Reality = Consciousness × Possibility²"). After extensive testing and two independent audits, the current description is more precise: CΨ is a composite quantum diagnostic that selects for simultaneously entangled and coherent pairwise states. The philosophical interpretation remains open but is no longer presented as established.

AI agents (local LLMs via LM Studio) contributed to early exploration. Their claims are marked Tier 4 (unverified) throughout the repository.

---

## Key experiments

### Core results
- [Star Topology](experiments/STAR_TOPOLOGY_OBSERVERS.md) - Three conditions for observer-observer connection
- [Subsystem Crossing](experiments/SUBSYSTEM_CROSSING.md) - Crossing is local to entangled pairs
- [N-Scaling Barrier](experiments/N_SCALING_BARRIER.md) - Why full-system crossing fails at large N
- [Crossing Taxonomy](experiments/CROSSING_TAXONOMY.md) - Type A/B/C classification
- [Noise Robustness](experiments/NOISE_ROBUSTNESS.md) - Taxonomy survives all local Pauli channels
- [Dynamic Entanglement](experiments/DYNAMIC_ENTANGLEMENT.md) - Crossings from product states

### Algebra and boundary
- [Mandelbrot Connection](experiments/MANDELBROT_CONNECTION.md) - Algebraic equivalence to z² + c
- [Dynamic Fixed Points](experiments/DYNAMIC_FIXED_POINTS.md) - The CΨ ≤ 1/4 bound
- [Boundary Navigation](experiments/BOUNDARY_NAVIGATION.md) - Crossing observations

### Hardware and calibration
- [IBM Quantum Tomography](experiments/IBM_QUANTUM_TOMOGRAPHY.md) - First hardware test on IBM Torino
- [Residual Analysis](experiments/RESIDUAL_ANALYSIS.md) - Anomalous late-time coherence

### Bridge and topology
- [Bridge Closure](experiments/BRIDGE_CLOSURE.md) - J=0 bridge is closed (null result)
- [Observer x Gravity Bridge](experiments/OBSERVER_GRAVITY_BRIDGE.md) - J>0 interval shift
- [No-Signalling Boundary](experiments/NO_SIGNALLING_BOUNDARY.md) - CΨ sees global regime change
- [QKD Eavesdropping Forensics](experiments/QKD_EAVESDROPPING_FORENSICS.md) - Eve's basis from CΨ profile

### Speculative (Tier 3+)
- [Black Holes, White Holes, Big Bang](experiments/BLACK_WHITE_HOLES_BIGBANG.md) - Speculative
- [Self-Consistency: Schwarzschild](experiments/SELF_CONSISTENCY_SCHWARZSCHILD.md) - Speculative
- [Universal Quantum Lifetime](experiments/UNIVERSAL_QUANTUM_LIFETIME.md) - The x³ + x = 1/2 equation

See [experiments/README.md](experiments/README.md) for the complete index.

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Free to share and adapt with attribution.

---

## Authors

**Thomas Wicht**, Independent Researcher, Germany
**Claude**, AI System, Anthropic

December 2025 - March 2026
