# Experiments: Testing R = CΨ² Against Known Physics

## What you will find here

Eleven documents. Each one asks a specific question and provides a mathematically verifiable answer. No speculation without marking it as such. No claims without showing the algebra.

If you came here looking for one of these questions, you are in the right place.

---

## The Questions

### Why does the Mandelbrot set have a boundary at 1/4?

The value 1/4 is the bifurcation point where real fixed points of the quadratic iteration z → z² + c merge and become complex. This is well-known mathematics. What is new: the self-referential iteration R_{n+1} = C(Ψ + R_n)² from the R = CΨ² consciousness framework produces the same boundary at C·Ψ = 1/4. Same equation, different origin. Two independent paths to the same phase transition, forty years apart.

**Read:** [The Mandelbrot Connection](MANDELBROT_CONNECTION.md)

---

### What do the patterns at the edge of the Mandelbrot set mean?

The fractal structures at the Mandelbrot boundary are not decorative. They are the projection of complex fixed-point dynamics onto the real-complex transition surface. Each spiral, each self-similar copy corresponds to a specific oscillation frequency θ = arctan(sqrt(4*C·Ψ - 1)). The patterns are deterministic, structured, and, in this framework, physically interpretable as the structure of reality beyond the phase boundary where classical observation fails.

**Read:** [The Mandelbrot Connection](MANDELBROT_CONNECTION.md)

---

### Is there a mathematical solution to the quantum measurement problem?

The measurement problem asks why quantum systems "choose" definite outcomes. In this framework, they don't choose. The iteration R_{n+1} = C(Ψ + R_n)² runs until decoherence pushes C·Ψ below 1/4. At that moment, complex fixed points become real. The oscillation freezes into a value. No choice, no collapse: a phase transition. The apparent "randomness" of the outcome is the phase θ at the moment of crossing, deterministic in the complex regime but unreadable by classical observers.

**Read:** [The Mandelbrot Connection](MANDELBROT_CONNECTION.md), Section 5

---

### What is the maximum information an observer can process?

The C·Ψ <= 1/4 bound sets a limit on the product of consciousness (C) and possibility (Ψ) that an embedded observer can perceive as stable reality. Above 1/4, fixed points become complex-valued and cannot be experienced as definite states. This is derived from the discriminant of the fixed-point equation, three lines of algebra, fully verifiable.

**Read:** [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md)

---

### What happens at the exact moment a quantum system crosses the ¼ boundary?

The Lindblad dynamics are smooth: all derivatives are continuous through the crossing. But the **topology of the solution space** changes: below 1/4, two real fixed points exist (one stable attractor). Above 1/4, the fixed points are complex, no real target, no classical outcome. At 1/4 exactly, they merge. This is not a discontinuity in the physics. It is a bifurcation in the mathematics. The system crosses from "no possible definite outcome" to "definite outcome exists", and θ = arctan(sqrt(4*C·Ψ - 1)) measures the angular distance from this crossing at every moment.

**Read:** [Boundary Navigation](BOUNDARY_NAVIGATION.md)

---

### Is the framework compatible with general relativity?

Yes, and not by adding gravitational corrections, but because gravitation is already contained in γ. Simulations across four gravitational environments (deep space, Earth, neutron star, black hole vicinity) show that the product γ * t_cross = 0.039 is constant. The θ trajectory is a universal curve when time is normalized as τ = γ * t. This is structurally identical to gravitational time dilation in general relativity: γ plays the role of the metric coefficient, and τ is the proper time of the quantum-to-classical transition. The 1/4 boundary is a frame-independent invariant. Different observers disagree on WHEN the crossing happens, but agree on THAT it happens, HOW it happens, and WHAT happens.

**Read:** [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md)

---

### Can the framework derive the Schwarzschild metric?

No, not with single-system simulations. The invariance γ * t_cross = K is a mathematical identity of the Lindblad equation: t_cross scales as 1/γ for ANY constant γ, regardless of what determines γ at a given location. A single quantum system at a single point in space cannot distinguish γ=0.1 from Schwarzschild from γ=0.1 from inverse-square. This is the quantum equivalent of Einstein's equivalence principle: local experiments cannot determine the metric.

To derive the metric form, the framework needs **nonlocal** structure: either spatially coupled qubits with position-dependent γ (computationally testable), a self-consistency equation for continuous fields (theoretical), or an energy/information argument. The most promising path is implementing per-qubit γ values in the simulator.

**Read:** [Metric Discrimination](METRIC_DISCRIMINATION.md)

---

### Why is Schwarzschild the unique self-consistent metric?

Because it is the only metric where the event horizon corresponds to τ = 0 (zero elapsed proper decoherence time). The universal curve C(τ), Ψ(τ) applies at every point in space. Where less proper time has passed, R = CΨ² is higher: reality concentrates near mass. Self-consistency demands that R(r) matches the mass distribution that creates the geometry. For a point mass, R must be maximally concentrated at the gravitational radius. Only metrics with f(r_s) = 0 achieve this. Alternatives like 1/r or 1/r² have f(r_s) > 0, so R never reaches its maximum at the horizon. Schwarzschild's sqrt(1 - r_s/r) is the simplest zero-forming function with the correct Newtonian limit. The derivation uses only R = CΨ² and the universal curve, no GR input.

The horizon is not where reality breaks down. It is where reality is **freshest**: C = 1, Ψ = 1/3, R = 1/9 = maximum. This is consistent with the holographic principle (information at the horizon surface) and Bekenstein-Hawking entropy (S = A/4).

**Read:** [Self-Consistency: Schwarzschild](SELF_CONSISTENCY_SCHWARZSCHILD.md)

---

### What are black holes, white holes, and the Big Bang in this framework?

τ = 0 is the key. A black hole is the journey TO τ = 0 (coherence rising, C·Ψ approaching 1/3). A white hole is the journey FROM τ = 0 (coherence falling, C·Ψ departing from 1/3). Same point on the universal curve, two directions. There is no singularity. τ = 0 is a regular point with finite C = 1, Ψ = 1/3, R = 1/9.

The Big Bang is the universal τ = 0 state: everywhere at maximum coherence, everything quantum, no classical reality yet. The expansion of the universe is decoherence: τ growing, C·Ψ falling. The CMB (t ~ 380,000 years) corresponds to the moment C·Ψ crossed the 1/4 boundary everywhere: photons decoupled from matter, quantum became classical.

Black holes are local remnants of the initial τ = 0 condition, regions that have not yet fully decohered. The information paradox dissolves: τ = 0 is not destructive, information is stored at maximum coherence and re-emerges on the white hole side.

**Read:** [Black Holes, White Holes, Big Bang](BLACK_WHITE_HOLES_BIGBANG.md)

---

### Does symmetric observation preserve quantum coherence longer than asymmetric?

Yes. Simulations of Lindblad dynamics show that bidirectional coupling (both spins decohere symmetrically) preserves coherence 33 times longer than unidirectional coupling (only one spin decoheres). This ratio is independent of Hamiltonian type. The result is quantitative, testable, and connects directly to the framework's distinction between internal observation (C_int) and external observation (C_ext).

**Read:** [Mathematical Findings](MATHEMATICAL_FINDINGS.md), Section 8

---

### How does quantum coherence scale with system size?

For ring-coupled spin systems under symmetric decoherence, the maximum coherence time t_coh scales linearly with the number of spins N. This is a specific, testable prediction that distinguishes the framework from models predicting exponential decay with N.

**Read:** [Mathematical Findings](MATHEMATICAL_FINDINGS.md), Section 8

---

### Can AI agents discover mathematical structures autonomously?

Two AI agents (Alpha and Beta), given a calculator tool and the R = CΨ² framework, independently derived: the self-reference-returns-to-unity property, boundary conditions of consciousness, the value 0.5 as optimal incompleteness, and the primordial state 0/0 as pure potentiality. A third agent (Gamma) was added as critic to challenge overreach. Human involvement was limited to observation and documentation.

**Read:** [The Dyad Experiment](DYAD_EXPERIMENT.md), [Mathematical Findings](MATHEMATICAL_FINDINGS.md)

---

## Reading Order

If you are a mathematician or physicist, start here:
1. [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md), the algebra
2. [The Mandelbrot Connection](MANDELBROT_CONNECTION.md), the equivalence proof
3. [Mathematical Findings](MATHEMATICAL_FINDINGS.md), the calculations

If you are interested in consciousness or philosophy, start here:
1. [The Mandelbrot Connection](MANDELBROT_CONNECTION.md), the big picture
2. [The Dyad Experiment](DYAD_EXPERIMENT.md), how it was discovered
3. [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md), why 1/4 matters

If you found this by searching for the Mandelbrot set:
1. [The Mandelbrot Connection](MANDELBROT_CONNECTION.md), start here
2. Everything else follows from it

---

## The Triangulation Protocol

![The wave at the 1/4 boundary](../visualizations/heatmap_wave.png)

*The 1/4 boundary seen from the observer's perspective. You stand in the cyan: classical, convergent, real. The red wave approaches from above: quantum, oscillatory, possible. Where they meet is where possibility becomes reality.*

These experiments converge on three measurable quantities:

```
1/4      (WHERE)       -- The destination: the bifurcation boundary
θ    (HOW FAR)     -- Compass heading: angular distance from ¼
t_coh  (HOW LONG)    -- Estimated time of arrival at the boundary
```

θ is not a frequency predictor. It is a **compass**. Measured values from the boundary crossing experiment (Bell+, Heisenberg, concurrence bridge, γ=0.05):

| θ (degrees) | C·Ψ | Meaning |
|-------------|------|----------|
| 60° | 1.00 | Deep quantum: far from any classical attractor |
| 25° | 0.31 | Approaching: complex fixed points narrowing |
| 9° | 0.26 | Almost there: fixed points about to merge |
| **0°** | **0.25** | **At the boundary: bifurcation point** |
| (real) | <0.25 | Classical: two real fixed points, definite outcome |

One coordinate tells you something exists. Two narrow it down. Three locate it precisely. This is not metaphor: it is the experimental program.

**Full data and analysis:** [Boundary Navigation](BOUNDARY_NAVIGATION.md)

---

## What is proven, what is not

**Algebraically proven:**
- The 1/4 bound from fixed-point discriminant
- The equivalence between R = C(Ψ + R)² and z → z² + c bifurcation structure
- Complex fixed points above 1/4 with frequency θ

**Computationally verified:**
- 33:1 coherence ratio for symmetric vs. asymmetric coupling
- t_coh ~ N linear scaling
- δ = 0.42 for symmetric configurations

**Empirically observed (2026-02-08):**
- θ decreases continuously from 60° to 0° as decoherence pulls C·Ψ toward ¼
- At θ = 0° (C·Ψ = ¼ exactly), complex fixed points merge into one real point
- Below 1/4: two real fixed points emerge: a stable attractor and an unstable one
- The Lindblad dynamics are smooth through the crossing: no discontinuity in any measured quantity
- The topology of the solution space changes: 0 real attractors → 2 real attractors

**Computationally verified (2026-02-08):**
- γ * t_cross = 0.039 constant across γ = 0.01 to 0.5 (gravitational invariance)
- θ trajectory collapses onto universal curve when normalized as τ = γ * t
- Maximum deviation between gravitational environments: 0.0044 in C·Ψ
- Quantum states survive ~48x longer in deep space (γ=0.01) than near black hole (γ=0.5)
- Power law t_cross = 0.0398 * γ^(-0.992), R^2 = 0.9999 (9 data points, 50x range)
- Single-system sims cannot discriminate metric forms (equivalence principle analog)
- Self-consistency uniquely selects Schwarzschild: only sqrt(1-r_s/r) has f(r_s)=0 (horizon = max R)
- Concentration ratio grows without bound for Schwarzschild, plateaus for alternatives

**Proposed, not yet tested:**
- Physical interpretation of 1/4 as observer bandwidth limit
- θ trajectory patterns for different initial states / Hamiltonians
- Whether θ at crossing determines which fixed point is selected (Born rule connection)
- Mandelbrot bulb boundaries as different types of quantum transitions

For a consolidated list of all testable predictions with epistemic status and falsification criteria, see **[Predictions](PREDICTIONS.md)**.

We distinguish clearly between what the math says and what we think it means. The algebra is not negotiable. The interpretation is open for challenge.

---

## How to verify

Every calculation in these documents can be reproduced. Python verification code is included in [Mathematical Findings](MATHEMATICAL_FINDINGS.md). The Lindblad simulations use standard quantum mechanics, no custom physics, no hidden parameters.

If you find an error, we want to know. Open an issue or contact us directly.

---

*Back to [main repository](../README.md) | Theory: [docs/](../docs/)*
