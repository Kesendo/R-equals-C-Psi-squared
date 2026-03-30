# Experiments: Open Quantum Systems Under Dephasing

<!-- Keywords: open quantum system experiments, Lindblad dephasing simulation,
palindromic Liouvillian spectrum, quantum decoherence channel capacity,
CΨ quarter boundary, quantum state transfer spin chain, IBM quantum hardware
validation, dephasing noise information channel, quantum MIMO channel,
R=CPsi2 framework experiments, palindromic spectral symmetry verification -->

**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

This folder contains all computational experiments for the R = CΨ² project.
The project studies the spectral structure of open quantum systems under
dephasing noise and discovers that the Liouvillian eigenvalue spectrum is
**palindromically paired** (every decay rate d has a partner at 2Σγ − d).
This symmetry, verified for 54,118 eigenvalues with zero exceptions, has
consequences for quantum state transfer, decoherence thresholds, and
information channels.

All experiments are reproducible using Python, QuTiP, and NumPy.
Scripts are in `simulations/`, results in `simulations/results/`.

---

## Headline Results

Four results that a new reader should see first:

### 1. One-line formula beats 18 years of optimization (March 24, 2026)
Concentrate all dephasing on one edge qubit, protect the rest.
gamma_edge = N*gamma_base - (N-1)*epsilon. This trivially simple rule
outperforms every published dephasing optimization by two orders of
magnitude. C#-validated from N=5 (360×) through N=15 (63.5×) vs V-shape.
ENAQT literature achieves 2-3×. Information grows quadratically with
chain length (SumMI ~ 0.0053*N²) instead of decaying exponentially.
First hardware test: selective DD 2-3× on ibm_torino.

-> **[Resonant Return: from SVD to formula](RESONANT_RETURN.md)**
-> [Signal Analysis: Quadratic Scaling](SIGNAL_ANALYSIS_SCALING.md)
-> [IBM Hardware: Selective DD](IBM_SACRIFICE_ZONE.md) (Tier 2, single run)

### 2. The dephasing channel (15.5 bits capacity)
The spatial profile of dephasing rates across a qubit chain is not just
noise. It is a readable information channel with 15.5 bits of theoretical
capacity at 1% measurement noise. An external agent encoding information
in the γ profile can be decoded from internal quantum observables with
100% accuracy.

→ **[Dephasing Noise as Information Channel (γ as Signal)](GAMMA_AS_SIGNAL.md)**
→ **[Practical γ Control (+124% MI)](GAMMA_CONTROL.md)**

### 3. The CΨ = 1/4 boundary (IBM hardware validated)
The product CΨ = Tr(ρ²) × L₁/(d−1) has a critical boundary at exactly
1/4, determined by the discriminant of the self-referential purity
recursion R = C(Ψ+R)². All standard quantum channels cross this boundary.
Validated on IBM Torino at 1.9% deviation.

→ **[IBM Hardware Validation](IBM_RUN3_PALINDROME.md)**
→ **[Crossing Taxonomy (Type A/B/C)](CROSSING_TAXONOMY.md)**
→ **[Boundary Navigation (θ compass)](BOUNDARY_NAVIGATION.md)**

### 4. The palindromic spectrum (proven, N=2 through N=8)
The Liouvillian eigenvalue spectrum under local Z-dephasing is exactly
palindromic for Heisenberg/XXZ systems on any graph. The conjugation
operator Π swaps populations (immune sector) with coherences (decaying
sector), creating a time-reversal symmetry in the rescaled frame.

→ **[Π as Time Reversal](PI_AS_TIME_REVERSAL.md)**
→ **[Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md)**
→ **[Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md)**

---

## Full Experiment Index

### The Dephasing Channel (γ as readable signal)

| Experiment | Key finding |
|-----------|------------|
| [γ as Signal](GAMMA_AS_SIGNAL.md) | The bidirectional bridge: γ profiles are 100% classifiable, 15.5 bits capacity, 5 independent SVD channels, 21.5× optimization |
| [γ Control](GAMMA_CONTROL.md) | V-shape dephasing profile gives +124% MI, DD on receiver +81%, AC modulation falsified, time-resolved decoder works |
| [Relay Protocol](RELAY_PROTOCOL.md) | Staged transfer with time-dependent γ: +83% end-to-end mutual information |
| [Scaling Curve](SCALING_CURVE.md) | MI vs chain length (N=3 to N=11), hierarchy falsification, push vs pull principle |
| [Resonant Return](RESONANT_RETURN.md) | SVD-optimal profiles (10×), sacrifice-zone formula (360× at N=5 through 63.5× at N=15), frequency pulsing falsified. The formula: gamma_edge = N*gamma_base - (N-1)*epsilon |
| [Signal Analysis: Scaling](SIGNAL_ANALYSIS_SCALING.md) | Sacrifice-zone formula scaling N=2-15. Quadratic growth (SumMI ~ 0.0053*N²), constant brake (-0.020), two converging channels |
| [Temporal Sacrifice](TEMPORAL_SACRIFICE.md) | Fold catastrophe at CΨ = ¼ observed: endpoint MI peaks at exact crossing. With Bell+bath: CΨ oscillates around ¼ (81 crossings at J=5.0). MI pulses at each crossing. Damped: each cycle deposits irreversible reality |

### The CΨ = 1/4 Boundary (decoherence threshold)

| Experiment | Key finding |
|-----------|------------|
| [Crossing Taxonomy](CROSSING_TAXONOMY.md) | Three observer types: Type A (pure-Ψ, K=0.072), Type B (mixed, K=0.039), Type C (never crosses). K-invariance from Lindblad scaling |
| [Boundary Navigation](BOUNDARY_NAVIGATION.md) | θ = arctan(√(4CΨ−1)) as compass to the 1/4 transition. Triangulation: WHERE (1/4), HOW FAR (θ), HOW LONG (t_coh) |
| [Subsystem Crossing](SUBSYSTEM_CROSSING.md) | Crossing is local to entangled pairs, not a whole-system property |
| [N-Scaling Barrier](N_SCALING_BARRIER.md) | Full-system CΨ drops below 1/4 at large N due to Hilbert space dimension, but subsystem pairs still cross |
| [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) | Product states can cross 1/4 upward (entanglement generation) |
| [Noise Robustness](NOISE_ROBUSTNESS.md) | Type A/B/C taxonomy is identical under σ_x, σ_y, σ_z dephasing |
| [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) | The CΨ ≤ 1/4 bound as attractor of the self-referential map |
| [Observer Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md) | Different observers see different crossing times for the same system |
| [Observer Dependent Visibility](OBSERVER_DEPENDENT_VISIBILITY.md) | γ-dependent visibility windows |
| [Mandelbrot Connection](MANDELBROT_CONNECTION.md) | CΨ ↔ c maps the 1/4 boundary to the Mandelbrot cardioid cusp |

### Palindromic Spectral Structure

| Experiment | Key finding |
|-----------|------------|
| [Π as Time Reversal](PI_AS_TIME_REVERSAL.md) | Π maps populations (past) to coherences (future). Standing wave = interference of forward and backward modes |
| [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md) | ZZZ is universal node (classical), XX/YY are antinodes (quantum). Bell rings, GHZ is silent |
| [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) | Palindrome extends beyond Heisenberg to XY, Ising, XXZ, DM interaction. Two Π families (P1/P4) |
| [XOR Space](XOR_SPACE.md) | Where information lives in the palindrome: GHZ vs W states, Pauli weight correlation |
| [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md) | Why depolarizing noise breaks the palindrome: 1:3 per-site split vs 2:2 for Z-dephasing |
| [V-Effect](V_EFFECT_PALINDROME.md) | What happens when the palindrome breaks: immune extremes, 3× more frequencies |
| [N→∞ Palindrome](N_INFINITY_PALINDROME.md) | Thermodynamic limit: Gaussian rate density, past/future boundary blurs |
| [Error Correction](ERROR_CORRECTION_PALINDROME.md) | Palindromic protection hierarchy, optimal state (90% slow-mode), Π as Z₄ operator |
| [Cavity Modes Formula](CAVITY_MODES_FORMULA.md) | At Sigma_gamma=0: Stat(N) = Sum_J m(J,N)*(2J+1)^2 (Clebsch-Gordan). Exact for chain, lower bound for symmetric topologies. Star has N-1 harmonic frequencies, chain has rich irrational spectrum |
| [IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) | Cavity modes meet IBM hardware: sacrifice zone protects the slowest oscillating modes at 2.81x (vs 1.97x measured). Same 43 frequencies, different damping. 100% palindromic under 26x asymmetric noise |

### Star Topology and Mediator Physics

| Experiment | Key finding |
|-----------|------------|
| [Star Topology](STAR_TOPOLOGY_OBSERVERS.md) | Three conditions for observer-observer connection through mediator S |
| [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md) | CΨ windows live on 3D manifold (98% variance in 3 PCs). Two modes: glide and switch. Phase map of 4 independent mechanisms |
| [What's Inside the Windows](WHATS_INSIDE_THE_WINDOWS.md) | Skeleton + rotation decomposition of CΨ visibility windows |
| [Quantum Sonar](QUANTUM_SONAR.md) | Passive detection of hidden observers through spectral shifts. IBM investigation: qubit detuning dominates |
| [Theta-Palindrome-Echo](THETA_PALINDROME_ECHO.md) | θ connects to the channel (r=0.87 with fidelity), not to the echo |
| [Orphaned Results](ORPHANED_RESULTS.md) | Topology as gatekeeper, antiferromagnet crossing, echo characterization |
| [Optimal QST Encoding](OPTIMAL_QST_ENCODING.md) | Negative result: standard encoding already near-optimal |

### IBM Quantum Hardware

| Experiment | Key finding |
|-----------|------------|
| **[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md)** | **All IBM data combined: r* threshold at precision 0.000014, fold one-way, sacrifice MI gradient, 12 permanent crossers (24,073 records, 133 qubits, 181 days)** |
| [IBM Run 3: Palindrome Validation](IBM_RUN3_PALINDROME.md) | CΨ = 1/4 crossing confirmed at 1.9% deviation on IBM Torino (Eagle r3, 127 qubits) |
| [IBM Sacrifice-Zone](IBM_SACRIFICE_ZONE.md) | Selective DD beats uniform DD by 2-3.2× at all 5 time points on ibm_torino. First hardware test of spatial noise engineering (Tier 2, single run, caveats apply) |
| [IBM Quantum Tomography](IBM_QUANTUM_TOMOGRAPHY.md) | Single-qubit state tomography runs on IBM hardware |
| [Fixed Point Shadow](FIXED_POINT_SHADOW.md) | Shadow investigation, IBM skeleton analysis |

### Benchmarks and Comparisons

| Experiment | Key finding |
|-----------|------------|
| [QST Bridge](QST_BRIDGE.md) | Connecting to 20 years of quantum state transfer literature |
| [Localizable Entanglement](LOCALIZABLE_ENTANGLEMENT_BENCHMARK.md) | LE vs CΨ comparison: three-layer separation (CoA/LE/CΨ) |

### Additional Experiments

| Experiment | Key finding |
|-----------|------------|
| [Standing Wave (Two Observers)](STANDING_WAVE_TWO_OBSERVERS.md) | Standing wave pattern with two observer perspectives |
| [Bridge Fingerprints](BRIDGE_FINGERPRINTS.md) | Different states produce different bridge signatures |
| [Bridge Closure](BRIDGE_CLOSURE.md) | J=0 bridge is dead (no-signalling holds exactly) |
| [No-Signalling Boundary](NO_SIGNALLING_BOUNDARY.md) | CΨ drops when B measures, but A's reduced state is unchanged |
| [Coherence Density](COHERENCE_DENSITY.md) | Coherence per qubit analysis |
| [Signal Processing View](SIGNAL_PROCESSING_VIEW.md) | Two-sector frequency structure: f(c+) and f(c−) as matched filter bank |
| [When Ψ Matters](WHEN_PSI_MATTERS.md) | AND-gate justification for CΨ product |
| [Minimum Crossing Energy](MINIMUM_CROSSING_ENERGY.md) | Energy requirements for boundary crossing |
| [Simulation Evidence](SIMULATION_EVIDENCE.md) | Comprehensive simulation results |
| [Residual Analysis](RESIDUAL_ANALYSIS.md) | Post-crossing residual coherence analysis |
| [Mathematical Findings](MATHEMATICAL_FINDINGS.md) | Collected mathematical results |
| [Algebraic Exploration](ALGEBRAIC_EXPLORATION.md) | Algebraic structure exploration |
| [Born Rule Mirror](BORN_RULE_MIRROR.md) | Connection to the Born rule |
| [Decoherence Relativity](DECOHERENCE_RELATIVITY.md) | Observer-dependent decoherence rates |
| [Gravitational Invariance](GRAVITATIONAL_INVARIANCE.md) | Gravitational analogy exploration |
| [Observer Gravity Bridge](OBSERVER_GRAVITY_BRIDGE.md) | γ gradient as gravitational analogue |
| [Operator Feedback](OPERATOR_FEEDBACK.md) | Operator feedback dynamics |
| [Universal Quantum Lifetime](UNIVERSAL_QUANTUM_LIFETIME.md) | Universal lifetime scaling |
| [Why the Sum](WHY_THE_SUM.md) | Why Σγ appears in the palindromic sum |
| [Dyad Experiment](DYAD_EXPERIMENT.md) | Two-observer dynamics |
| [QKD Eavesdropping Forensics](QKD_EAVESDROPPING_FORENSICS.md) | Application to quantum key distribution |

---

## How to Read This

**If you are new to the project:** Start with the three headline results
above, then read the [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
for the core theorem.

**If you are a quantum information researcher:** The CΨ = 1/4 boundary
section and the IBM validation are the most directly relevant. The
palindromic structure connects to Liouvillian symmetry classification
(cf. Haga et al. 2023, Buca and Prosen 2012).

**If you are a signal processing engineer:** Start with
[γ as Signal](GAMMA_AS_SIGNAL.md) and its signal engineering perspective.
The system is a quantum MIMO channel with palindromic matched filters.

**If you want to reproduce results:** Every experiment links to its
simulation script. All use QuTiP + NumPy. Typical runtime: seconds to
minutes on a standard laptop.

**If you want a guided reading path:** See the [Reading Guide](../docs/READING_GUIDE.md),
which organizes the experiments into three stories: the proof, the
application, and the ontology.
