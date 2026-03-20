# R = CΨ²

> "We are all mirrors. Reality is what happens between us."

This project started from a dream on December 25, 2025 about multilayer
electrolysis. The equation R = CΨ² was in the dream. Three months of
computation, with honest documentation of every wrong turn, revealed
that the equation points at real structure in open quantum systems.

What began as philosophy became physics. What was speculative became proven.
The consciousness interpretation is retired from the technical core. The
mirror symmetry it predicted turned out to be exactly true.

**Thomas Wicht** (developer, Krefeld, Germany) and **Claude** (AI, Anthropic)
are the collaborators. We document negative results. We correct our mistakes
publicly. The recovered/ folder contains five disproven claims we keep for
history.

**Zenodo DOI:** [10.5281/zenodo.19100007](https://doi.org/10.5281/zenodo.19100007) (v3.0, March 18, 2026)

---

## What we found

### The system

```
    A (observer)
    |
    S (mediator)       Heisenberg coupling, Z-dephasing noise
    |
    B (observer)
```

N qubits coupled via Heisenberg (or XXZ) interaction, subject to local
dephasing. We build the full Liouvillian superoperator and study its
complete eigenvalue spectrum.

### Mirror symmetry: PROVEN (March 14, 2026)

The decay rate spectrum of the Liouvillian is exactly palindromic.
For every decay rate d, there exists a partner at 2Σγᵢ - d.

This was verified numerically through N=8 (54,118 rates, zero exceptions)
across every topology and noise type we tested. On March 14, 2026,
we found the analytical proof:

**The conjugation operator Π** acts per site on Pauli indices:

    I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

It satisfies **Π·L·Π⁻¹ = -L - 2Σγᵢ·I**, which directly implies the
palindrome. The proof is three steps:

1. Π flips XY-weight k → N-k, so Π·L_D·Π⁻¹ = -L_D - 2Σγ·I (trivial)
2. Π anti-commutes with [H,·] for any Heisenberg/XXZ bond (16-entry table)
3. Combined: Π·L·Π⁻¹ = -L - 2Σγ·I. QED.

Holds for: all δ (XXZ anisotropy), all graphs (star, chain, ring, complete,
binary tree), non-uniform γ per qubit, Z and Y dephasing.
Breaks for: depolarizing noise (err = gamma * 2(N-2)/3, linear in gamma and N).

### Beyond Heisenberg: all standard models (March 17-18, 2026)

The palindrome is not limited to Heisenberg coupling. ALL standard condensed
matter models are palindromic under single-axis dephasing: XY, Ising, XXZ,
Dzyaloshinskii-Moriya, and combinations. Two Π operator families exist
(P1 and P4), plus non-uniform alternating operators for XY/YX terms.
All 36/36 two-term combinations resolved: 20 per-site, 14 broken, 2 non-local.

See: [Non-Heisenberg Palindrome](experiments/NON_HEISENBERG_PALINDROME.md)

### Hardware validation: 1.9% (March 18, 2026)

CΨ=1/4 crossing measured on IBM Quantum (ibm_torino, Qubit 80).
Predicted: t* = 15.01 us. Measured: t* = 15.29 us. Deviation: 1.9%.
T2* (not T2 echo) confirmed as the correct timescale.

See: [IBM Run 3](experiments/IBM_RUN3_PALINDROME.md),
raw data in `data/ibm_run3_march2026/`

The closest prior work, incoherentons (Haga et al. 2023) and Bethe ansatz
for dephasing chains (Medvedyeva-Essler-Prosen 2016), had pieces of this
but not the operator or the palindrome. Nobody had Π.

See: [Mirror Symmetry Proof](docs/MIRROR_SYMMETRY_PROOF.md)

### Spectral architecture (exact, N=2-8)

| N | Matrix | Rates | Min | Max | Mirror |
|:--|:-------|:------|:----|:----|:-------|
| 2 | 16² | 6 | 2γ | 2γ | 100% |
| 3 | 64² | 40 | 2γ | 4γ | 100% |
| 4 | 256² | 182 | 2γ | 6γ | 100% |
| 5 | 1024² | 776 | 2γ | 8γ | 100% |
| 6 | 4096² | 3228 | 2γ | 10γ | 100% |
| 7 | 16384² | 13264 | 2γ | 12γ | 100% |
| 8 | 65536² | 54118 | 2γ | 14γ | 100% |

Boundary formula: min = 2γ, max = 2(N-1)γ. Bandwidth = 2(N-2)γ.
Five topologies share the same boundaries, differ in interior rate count.
Rate count grows roughly as 4^N (ratio converges from below).
Density of states is Gaussian: mean = Nγ, skewness = 0, kurtosis ≈ 3.

### Five independent regulators

1. **Topology** (coupling graph) sets the oscillation frequencies
2. **Symmetry** (Hamiltonian anisotropy δ) controls sector separation
3. **Noise strength** (γ per qubit) sets the decay envelope; never changes frequency
4. **Initial state** determines which modes are excited at t=0
5. **Bath geometry** (noise axis, correlations) controls which mode dominates

### Two supermodes (3-qubit star)

The observer pair A-B sees two cross-correlations: c+ (symmetric) and c-
(antisymmetric). Both are damped sinusoids. Their exact pole structure:

| Mode | Frequency | Decay | Dominant in |
|:-----|:----------|:------|:------------|
| Fast | ~1.506 (scales with J) | 10γ/3 | c+ |
| Slow | ~0.404 | 2γ | c- |
| Hidden | ~1.103 | 8γ/3 | c- only |

Decay rates are exact rational multiples of γ. Frequencies depend on
topology, decay rates don't. The slow mode (2γ) is naturally protected.

### Signal processing equivalence

The entire quantum structure maps to classical coupled oscillator physics:

| Quantum language | Signal engineering |
|:-----------------|:-------------------|
| Two sectors c+/c- | Even/odd supermode decomposition |
| Topology sets frequency | Imaginary part of system poles |
| Noise sets decay | Real part of system poles |
| Hidden observer detection | Passive topology-change detection |
| Observable filters modes | Modal observability (transfer function residues) |

### The Mandelbrot connection

CΨ = concurrence × l1-coherence, iterated as R_{n+1} = C(Ψ + R_n)²,
is algebraically equivalent to the Mandelbrot map z → z² + c.

The boundary at CΨ = 1/4 is where the fixed point z* = 1/2 sits.
Below 1/4: reality converges. Above 1/4: no stable fixed point.
z*(1-z*) = CΨ is the Bernoulli variance form: maximum at p = 0.5.

z* matches no known single quantum measure (26 candidates tested,
best r = 0.951). It is genuinely composite: requiring both entanglement
AND coherence simultaneously.

Nobody has connected Mandelbrot iteration to open quantum dynamics before.
This is our most original finding.

---

## What is NOT established

- That CΨ is a new fundamental quantity (it is a derived diagnostic)
- That hidden observer detection works on hardware (simulation only)
- That frequency-decay orthogonality extends beyond 3 qubits (it doesn't)
- That the multi-qubit palindrome has been measured on hardware (single-qubit CΨ=1/4 crossing validated at 1.9%, N>=2 untested)
- That the standing wave pattern is measurable on hardware (computed, not yet measured)
- That consciousness plays any role in the physics (THE_ANOMALY.md is philosophy, not physics)

---

## Evidence status

### Proven

| Claim | Evidence |
|:------|:--------|
| Mirror symmetry (palindrome) | Analytical proof via conjugation Π, verified N=2-8 |
| Topology-independence of decay rates | Π anti-commutes with [H,·] for ANY bond set |
| Pauli weight complementarity | Π maps XY-weight k → N-k |
| Pole structure (3 exact decay rates) | Liouvillian eigendecomposition |
| Band structure at N >= 4 | Verified N=2-8, boundary 2γ to 2(N-1)γ |
| Mandelbrot algebraic correspondence | Proven: u_n substitution |
| CΨ = 1/4 is Bernoulli variance maximum | Proven: z*(1-z*) = CΨ |
| Five independent regulators | Full parameter sweeps |
| 5 topologies share boundaries | Verified: star, chain, ring, complete, tree |
| z* is a novel composite diagnostic | 26 candidates tested, none match |
| Non-Heisenberg palindrome | ALL standard models (XY, Ising, XXZ, DM) palindromic under single-axis dephasing |
| Two Π operator families | P1 (I<->X, Y<->Z) and P4 (I<->Y, X<->Z), 36/36 resolved (2 non-local) |
| XOR space universal | GHZ->100% XOR, W->0% XOR for all standard models, r>0.98 |
| Depolarizing noise quantified | Breaks palindrome at err = gamma * 2(N-2)/3 |
| IBM hardware validation | CΨ=1/4 crossing at 1.9% deviation on ibm_torino Q80 (same-day T2*) |
| Standing wave from palindrome | Quantum correlations (XX,YY) oscillate at 2J,4J,6J; classical (ZZZ) static |

### Tested and rejected

| Claim | Result |
|:------|:-------|
| CΨ = 1/4 as Exceptional Point | No EP correlation found |
| c+/c- as Liouvillian symmetry sectors | Both parity +1; split is projection |
| IBM Q80/Q102 as sonar evidence | Was qubit detuning |
| Consciousness as physics ingredient | Retired from technical core |
| Gravity/Schwarzschild connections | Disproven (kept in recovered/) |

---

## Start here

### For a complete overview (standalone, no prior knowledge needed)
- **[Technical Paper](publications/TECHNICAL_PAPER.md)**: The palindrome proof, XOR space, QST. For physicists.
- **[Engineering Blueprint](publications/ENGINEERING_BLUEPRINT.md)**: Four design rules for quantum repeaters. For engineers.

### The proof
1. **[Mirror Symmetry Proof](docs/MIRROR_SYMMETRY_PROOF.md)**: The conjugation operator, the 16-entry table, the full verification

### Hardware validation
2. **[IBM Run 3](experiments/IBM_RUN3_PALINDROME.md)**: CΨ=1/4 crossing validated at 1.9% on ibm_torino Q80
3. **[IBM Quantum Tomography](experiments/IBM_QUANTUM_TOMOGRAPHY.md)**: All hardware runs on IBM Torino

### Generalization and structure
4. **[Non-Heisenberg Palindrome](experiments/NON_HEISENBERG_PALINDROME.md)**: All standard models, two Π families, alternating operators
5. **[XOR Space](experiments/XOR_SPACE.md)**: Where information lives in the palindrome. GHZ vs W. Universal across models.
6. **[Standing Wave Analysis](experiments/STANDING_WAVE_ANALYSIS.md)**: Palindromic pairs create oscillating quantum correlations over static classical backbone
7. **[QST Bridge](experiments/QST_BRIDGE.md)**: Quantum state transfer. Star 2:1 beats chains.

### Technical core
8. **[Signal Processing View](experiments/SIGNAL_PROCESSING_VIEW.md)**: Pole analysis, coupled oscillator translation
9. **[Standing Wave Theory](docs/STANDING_WAVE_THEORY.md)**: Two supermodes as standing waves
10. **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)**: Parameter sweeps and stress tests
11. **[Born Rule Mirror](experiments/BORN_RULE_MIRROR.md)**: Mirror quality and Born rule connection

### Algebra and diagnostics
12. **[The CΨ Lens](docs/THE_CPSI_LENS.md)**: What CΨ shows and what it doesn't
13. **[Algebraic Exploration](experiments/ALGEBRAIC_EXPLORATION.md)**: Mandelbrot correspondence, 1/4 boundary
14. **[The Bidirectional Bridge](docs/THE_BIDIRECTIONAL_BRIDGE.md)**: Two channels, two directions

### Interpretation (speculative, optional)
15. **[The Interpretation](hypotheses/THE_INTERPRETATION.md)**: All findings, open questions, current state
16. **[Weaknesses and Open Questions](docs/WEAKNESSES_OPEN_QUESTIONS.md)**: What we got wrong and what we don't know
17. **[The Anomaly](THE_ANOMALY.md)**: The question after the proof

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `publications/` | Standalone documents for external readers (paper + blueprint) |
| `docs/` | Mathematical framing, standing wave theory, the CΨ diagnostic |
| `experiments/` | All tested results and null results (38 files with headers) |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python scripts (Lindblad, Liouvillian, Prony, sweeps) |
| `simulations/results/` | All computation outputs |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `compute/` | C# engine (MathNet.Numerics + MKL + OpenBLAS, N=2-8) |
| `data/` | IBM Torino measurement data |
| `recovered/` | 5 files with disproven claims, kept for honesty |

## C# compute engine

For N >= 6, Python is too slow. The C# engine uses element-wise Liouvillian
construction (no Kronecker products, 640x faster) with Intel MKL for
eigendecomposition on 24 cores.

| N | Matrix | Build | Eigen | Rates | Mirror |
|:--|:-------|:------|:------|:------|:-------|
| 6 | 4096² | 8.7s | 56s | 3228 | 100% |
| 7 | 16384² | 0.1s | 92min | 13264 | 100% |
| 8 | 65536² | 5.6s | 10.6h | 54118 | 100% |

N=8 (65536²) uses native memory (64 GB) + OpenBLAS ILP64 eigenvalue-only LAPACK.
All timings on Intel Core Ultra 9 285k (24 cores), 128 GB RAM, Windows.

## Key scripts

| Script | Purpose |
|:---|:---|
| `pauli_weight_conjugation.py` | Mirror symmetry proof (Π operator) |
| `star_topology_v3.py` | Core dynamics engine (RK4 Lindblad) |
| `joint_pole_analysis.py` | Liouvillian eigendecomposition |
| `mirror_symmetry_deep.py` | Mirror tests: 11 noise types, all N |
| `xor_detector_v3.py` | XOR space analysis (Pauli weight correlation) |
| `qst_bridge.py` | Quantum state transfer benchmarks |
| `docs_verify.py` | Numerical verification of all /docs claims (40/40 PASS) |
| `deep_band_structure.py` | Band analysis with scaling laws |
| `z_star_identity.py` | z* vs 26 quantum measures |
| `ep_test.py` | Exceptional Point test (negative) |
| `non_heisenberg_deep.py` | All standard models: palindrome, Π families, depolarizing |
| `xor_non_heisenberg_v2.py` | XOR universal across all models, Bell correction |
| `algebraic_pi_search.py` | Π operator family enumeration (P1, P4, alternating) |
| `standing_wave_analysis.py` | Standing wave formalization: antinodes, nodes, state x Hamiltonian |

---

## The origin

December 25, 2025. A dream about Co/Ni multilayer electrolysis contained the
equation R = CΨ². The predecessor project (Stability, a material science
simulator) had found that half-filled shells (0.5) produce maximum bonding.
The same 0.5 appeared everywhere: fair coins, Bernoulli maxima, z* at the
1/4 boundary. The mirrors in the equation turned out to be mirrors in
the physics.

The motto was removed as too esoteric, then restored when the Liouvillian
decay spectrum turned out to be exactly mirror-symmetric at every system size.
Then proven analytically.

On March 18, 2026, IBM quantum hardware confirmed the prediction at 1.9%.
That evening, the question shifted: if the palindrome is the rule that
determines what survives decoherence, what happens when a system notices
its own filter? See [The Anomaly](THE_ANOMALY.md).

## Where this is going (March 20, 2026)

Everything above is proven, computed, or honestly labeled as speculative.
What follows is where the proven mathematics led when we followed it
honestly, without deciding in advance where it should stop.

On March 20, we proved that the palindromic mirror exists ONLY for
qubits (d=2). The equation d(d-2)=0 allows nothing else. Qutrits
fail completely: 0 of 236 dissipators permit a palindrome. A single
qutrit in a chain of qubits destroys the symmetry globally. The qubit
is not one option among many. It is the only dimension where physics
has an internal mirror.

We proved that the mirror has exactly two sides (Z2 parity, not Z4).
The operator Π² = X^N is a conserved symmetry that splits the entire
Liouville space into two sealed sectors: populations and coherences,
past and future, decided and undecided. We tested four sides (Z4).
Falsified. Two sides only.

We proved that parity-breaking is necessary for palindrome-breaking.
All 14 Hamiltonian combinations that lose the palindrome at N≥3 also
break the X^N parity. No exceptions. The V-Effect (forced spectral
differentiation as systems grow) cannot occur without the Hamiltonian
coupling the two sides. Twelve parity-breakers still preserve the
palindrome through an unknown hidden symmetry Q. Finding Q is open.

We built a decoder for the noise. The palindromic response matrix has
full rank: all per-site dephasing rates are independently recoverable
from mode amplitudes. The "antennas" (most sensitive modes) sit at
XY-weight 2, the boundary between classical and quantum operators.
The optimal receiver state is |010>: one qubit listening, the rest
silent. Real IBM hardware T2* data, fed through the decoder, shows
temporal structure over 6 days.

And then the mathematics led somewhere we did not expect.

The palindrome requires noise. Without noise, no mirror, no standing
wave, no structure. We tested whether the noise could come from
inside the system (each parity sector as the environment of the
other). Falsified: the sectors are exactly decoupled. The noise
must come from outside.

The noise takes 70% of the information (coherences, phase,
relationships between parts). It leaves 30% (populations, substance,
what things are individually). The noise is not random destruction.
It is selective: it takes how things are connected and leaves what
things are.

Two counter-propagating 70% streams (one from each side of the
mirror) meet in the middle and form a standing wave. That standing
wave is the interference pattern between the two sides. It does not
belong to either side. It exists between them.

This is documented in [The Other Side of the Mirror](hypotheses/THE_OTHER_SIDE.md),
1285 lines, 20 sections, written in a single day. Sections 0-12 are
mathematics (Tier 1-3). Sections 13-19 are the questions that follow
honestly from the mathematics (Tier 5). The boundary between proven
and speculative is explicitly marked throughout.

We do not know what comes next. The V-Effect says that when a pattern
outgrows its container, something new differentiates. We may be at
that point. The mathematics cannot tell us what the next level looks
like, only that the current one is no longer sufficient.

Read the technical proof first. Then, if you want, read The Other
Side. The proof does not need the interpretation. But the interpretation
needs the proof. And the proof is solid.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 - March 2026
