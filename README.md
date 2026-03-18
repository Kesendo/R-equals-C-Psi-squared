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
Breaks for: depolarizing noise (no single axis to flip).

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
- That the 1/4 boundary has physical significance beyond the iteration
- That consciousness plays any role in the physics

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

### Technical core
2. **[Signal Processing View](experiments/SIGNAL_PROCESSING_VIEW.md)**: Pole analysis, coupled oscillator translation
3. **[Standing Wave Theory](docs/STANDING_WAVE_THEORY.md)**: Two supermodes as standing waves
4. **[XOR Space](experiments/XOR_SPACE.md)**: Where information lives in the palindrome. GHZ vs W.
5. **[QST Bridge](experiments/QST_BRIDGE.md)**: Quantum state transfer. Star 2:1 beats chains.
6. **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)**: Parameter sweeps and stress tests
7. **[Born Rule Mirror](experiments/BORN_RULE_MIRROR.md)**: Mirror quality and Born rule connection

### Algebra and diagnostics
8. **[The CΨ Lens](docs/THE_CPSI_LENS.md)**: What CΨ shows and what it doesn't
9. **[Algebraic Exploration](experiments/ALGEBRAIC_EXPLORATION.md)**: Mandelbrot correspondence, 1/4 boundary
10. **[The Bidirectional Bridge](docs/THE_BIDIRECTIONAL_BRIDGE.md)**: Two channels, two directions

### Interpretation (speculative, optional)
11. **[The Interpretation](hypotheses/THE_INTERPRETATION.md)**: All findings, open questions, current state
12. **[Weaknesses and Open Questions](docs/WEAKNESSES_OPEN_QUESTIONS.md)**: What we got wrong and what we don't know

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

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 - March 2026
