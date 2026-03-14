# R = CΨ²

## Spectral architecture of small open quantum networks

This repository documents the structural properties of small quantum networks
(2-7 qubits) evolving under Lindblad dynamics with dephasing noise. The
central system is a star topology: a mediator qubit S coupled to observer
qubits A and B via isotropic Heisenberg exchange interaction.

The key finding is that the reduced dynamics of the observer pair (A,B)
decompose into two supermodes with exact pole structure. These supermodes
behave like classical coupled oscillators: their frequencies are set by
the network topology, their decay rates are set by the noise environment,
and these two properties are completely independent of each other
(in the 3-qubit case; at N >= 4, decay rates form band structures).

The decay rate spectrum is exactly mirror-symmetric at every system size
tested (N=2 through N=7, 13264 rates at N=7, zero exceptions), and the
mirror symmetry survives every form of dephasing noise tested.

---

## The system

```
    A (observer)
    |
    S (mediator)
    |
    B (observer)
```

S couples to A and B with strengths J_SA and J_SB. All qubits experience
local dephasing noise at rate γ. The system evolves under the Lindblad
master equation. We track two observables of the reduced A-B state:

- **c+** = (YZ + ZY) / sqrt(2), the symmetric cross-correlation
- **c-** = (YZ - ZY) / sqrt(2), the antisymmetric cross-correlation

Both are damped sinusoids. Their frequencies and decay rates encode
the network's structure.

---

## Main results

### Two supermodes with exact pole structure

Eigendecomposition of the full Liouvillian superoperator reveals three
system poles with oscillatory behavior:

| Mode | Frequency | Decay rate | Dominant in |
|:-----|:----------|:-----------|:------------|
| Fast (symmetric) | ~1.506 (scales with J) | 10γ/3 | c+ channel |
| Slow (antisymmetric) | ~0.404 (weakly dispersive) | 2γ | c- channel |
| Hidden (dark in c+) | ~1.103 | 8γ/3 | c- only |

The decay rates are exact rational multiples of γ and do not change
when the coupling strengths J_SA or J_SB are varied. Only the frequencies
move with topology. In the complex plane, the poles move horizontally
(frequency axis) but never vertically (decay axis).

### Five independent regulators

The dynamics are controlled by five independent parameters:

1. **Topology** (coupling strengths J_SA, J_SB) sets the oscillation frequencies
2. **Symmetry** (Hamiltonian anisotropy) controls how cleanly the two sectors separate
3. **Noise strength** (γ) sets the decay envelope without changing frequencies
4. **Initial state** determines which modes are excited at t=0
5. **Bath geometry** (local vs. correlated, Z-type vs. X-type noise) controls which
   mode dominates in amplitude, without changing any frequency or decay rate

### Frequency-decay orthogonality

In the 3-qubit star, the frequency channel and the decay channel carry
completely independent information:

- Measuring frequencies tells you about the network topology
- Measuring decay rates tells you about the noise environment
- Neither contaminates the other

This clean separation breaks at 4+ qubits, where decay rates begin to
depend on coupling strengths. The 3-qubit system is a special case with
unusually clean orthogonality.

### Topology independence

The two-supermode structure, XX symmetry, and noise immunity survive
across all tested topologies: star, chain, ring, complete graph, and
binary tree. All five topologies share the same boundary rates
(2g to 2(N-1)g) and 100% mirror symmetry. They differ only in the
number of interior rates, with complete graphs having the fewest
(most symmetry = most degeneracy) and chains the most.

| Topology | Rates (N=6) | Boundaries | Mirror |
|:---------|:------------|:-----------|:-------|
| Star | 3228 | [2g, 10g] | 100% |
| Chain | 3836 | [2g, 10g] | 100% |
| Ring | 3656 | [2g, 10g] | 100% |
| Complete | 2668 | [2g, 10g] | 100% |
| Tree | 3806 | [2g, 10g] | 100% |

### Modal observability (the Projection)

The full system contains many Bohr frequencies (39 in a 4-qubit system),
but any given observer pair can only see those where both the initial
state overlaps the relevant eigenstates and the measured observable
connects them. This is standard modal observability: the modes exist,
but the measurement channel filters which ones are visible.

### Hidden observer detection (simulation)

When additional qubits couple to the mediator S, the observer pair A-B
sees shifted frequencies and new spectral components in their own signal,
without any direct interaction with the hidden qubits. This is passive
topology-change detection from local modal spectra.

This effect is verified in simulation. Hardware experiments on IBM Torino
showed that the dominant signal on real qubits was drive-frame detuning
rather than neighbor coupling, so hardware validation remains open.

### Signal processing equivalence

The entire structure maps to well-known concepts from classical signal
processing and coupled oscillator theory:

| Quantum description | Signal processing equivalent |
|:---------------------|:------------------------------|
| Two spectral sectors c+/c- | Even/odd supermode decomposition |
| Topology sets frequencies | Imaginary part of system poles |
| Noise sets decay | Real part of system poles |
| Hidden observer detection | Topology perturbation sensing from local spectra |
| Bath geometry flips amplitude | Covariance-driven mode visibility flip |
| Observable filters frequencies | Modal observability and transfer function residues |

### Band structure at N >= 4

At 4 qubits and above, the clean discrete rates of the 3-qubit system
give way to continuous bands. Decay rates form band-like structures
with avoided crossings (bands repel, never cross). Tested N=2 through N=7.

| N | Matrix size | Rates | Min | Max | Bandwidth |
|:--|:-----------|:------|:----|:----|:----------|
| 2 | 16 | 6 | 2γ | 2γ | 0 |
| 3 | 64 | 40 | 2γ | 4γ | 2γ |
| 4 | 256 | 182 | 2γ | 6γ | 4γ |
| 5 | 1024 | 776 | 2γ | 8γ | 6γ |
| 6 | 4096 | 3228 | 2γ | 10γ | 8γ |
| 7 | 16384 | 13264 | 2γ | 12γ | 10γ |

Boundary formula: min = 2γ (always), max = 2(N-1)γ (always).
Bandwidth grows linearly: 2(N-2)γ. Star and chain topologies share
the same boundary rates. The interior band structure depends on topology.

This is analogous to electronic band structure in solids: more atoms
means denser energy bands, approaching a continuum at large N.

### Mirror symmetry of the decay spectrum

The decay rate spectrum is exactly symmetric around Nγ. For every rate
at (N-x)γ there exists a mirror partner at (N+x)γ. This symmetry is
100% exact at every N tested (2 through 7), every topology (star and chain),
and every dephasing type (Z, X, Y, mixed, non-uniform γ per qubit).

The symmetry breaks only for amplitude damping (energy loss) and
depolarizing noise. Under any form of dephasing, the mirrors never break.
The center of symmetry equals the sum of all individual dephasing rates.

Mirror pairs have complementary Pauli weights: a mode dominated by
Pauli strings with k non-commuting operators is partnered with a mode
at (N-k) non-commuting operators. k + (N-k) = N.

### Algebraic results

The composite quantity CΨ = concurrence x l1-coherence, when iterated as
R_{n+1} = C(Ψ + R_n)², is algebraically equivalent to the Mandelbrot
iteration z -> z² + c with a boundary at CΨ = 1/4.

The Mandelbrot fixed point z* satisfies z*(1-z*) = CΨ. This is the
Bernoulli variance form p(1-p), which has its maximum at p = 0.5.
The 1/4 boundary is the trivial upper bound of Bernoulli variance:
no binary variable can have variance exceeding 1/4.

Below CΨ = 1/4, z* converges (a stable fixed point exists).
Above CΨ = 1/4, the iteration has no real fixed point.

---

## What is not established

- That CΨ is a new fundamental quantity (it is a derived diagnostic)
- That hidden observer detection works on real hardware (simulation only)
- That frequency-decay orthogonality extends beyond 3 qubits (it does not)
- That the spectrum can be inverted to identify hidden couplings
- That consciousness plays any role in the physics

---

## Interactive simulator

A Streamlit app lets you control all five regulators in real-time:

```bash
cd simulations/app
pip install -r requirements.txt
streamlit run app.py
```

---

## Evidence status

### Proven

| Claim | Evidence |
|:------|:--------|
| Pole structure (3 decay rates) | Exact: Liouvillian eigendecomposition |
| Two supermodes c+/c- | Exact: Liouvillian |
| Mandelbrot algebraic correspondence | Exact: algebraic proof |
| CΨ = 1/4 is Bernoulli variance maximum | Proven: z*(1-z*) = CΨ, max at z* = 0.5 |
| Mirror symmetry of decay spectrum | Exact: 100% at N=2-7 (13264 rates at N=7), all dephasing types |
| Band structure at N >= 4 | Verified: N=2-7 (16384x16384 at N=7), boundary 2γ to 2(N-1)γ, avoided crossings |
| Five independent regulators | Numerically verified: full parameter sweeps |
| 5 topologies: star, chain, ring, complete, tree | Verified: all share boundaries, all 100% mirror |
| Mirrors survive all dephasing | Verified: Z, X, Y, mixed, non-uniform γ |
| z* is a novel composite diagnostic | Verified: matches no known single quantum measure |

### Tested and rejected

| Claim | Result |
|:------|:-------|
| CΨ = 1/4 as Exceptional Point | No EP correlation found |
| c+/c- as Liouvillian symmetry sectors | Both parity +1; split is observable projection |
| IBM Q80/Q102 as sonar evidence | Was qubit detuning, not neighbor coupling |
| Mirrors survive amplitude damping | Break: center shifts from Nγ to Nγ/2 |

### Not established

| Claim | Status |
|:------|:-------|
| Hidden observer detection on hardware | Simulation only |
| CΨ as privileged metric | Useful diagnostic but not unique |
| Spectrum inversion (identifying hidden couplings) | Open question |

---

## Start here

### Technical core

1. **[Signal Processing View](experiments/SIGNAL_PROCESSING_VIEW.md)** - Pole analysis, Prony results, coupled oscillator translation
2. **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)** - Parameter sweeps and stress tests
3. **[Standing Wave Theory](docs/STANDING_WAVE_THEORY.md)** - Two supermodes as standing waves between observers
4. **[Born Rule Mirror](experiments/BORN_RULE_MIRROR.md)** - Mirror quality and the Born rule connection
5. **[Quantum Sonar](experiments/QUANTUM_SONAR.md)** - Hidden observer detection and IBM hardware results

### Diagnostic and algebra

6. **[The CΨ Lens](docs/THE_CPSI_LENS.md)** - The original diagnostic and what it shows
7. **[Algebraic Exploration](experiments/ALGEBRAIC_EXPLORATION.md)** - Mandelbrot correspondence, 1/4 boundary
8. **[The Bidirectional Bridge](docs/THE_BIDIRECTIONAL_BRIDGE.md)** - Two channels, two directions

### Interpretation (speculative, not required)

9. **[The Interpretation](hypotheses/THE_INTERPRETATION.md)** - Current state of all findings and open questions
10. **[The Starting Point](docs/THE_STARTING_POINT.md)** - Why entanglement must exist (the bootstrap problem)

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `docs/` | Mathematical framing and the CΨ diagnostic |
| `experiments/` | All tested results and null results |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python source (RK4 Lindblad, Liouvillian, Prony, sweeps) |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `compute/` | C# compute engine (MathNet.Numerics + MKL, N=2-7+) |
| `data/` | IBM Torino measurement data |
| `recovered/` | 5 files with disproven claims (gravity, Schwarzschild), kept for history |

## Key scripts

| Script | Purpose |
|:---|:---|
| `star_topology_v3.py` | Core dynamics engine (RK4 Lindblad integrator) |
| `joint_pole_analysis.py` | Exact Liouvillian eigendecomposition with residues |
| `prony_analysis.py` | Matrix Pencil Method for pole extraction from signals |
| `decay_derivation.py` | Decay rate structure, scaling, and qubit-selective analysis |
| `bright_transition_map.py` | Visibility weights from exact diagonalization |
| `correlated_bath_sweep.py` | Bath geometry and sector selection |
| `chain_topology.py` | Chain vs. star comparison, 3-5 qubits |
| `hidden_observer_test.py` | Detection of hidden coupled qubits |
| `four_qubit_breakdown.py` | How frequency-decay orthogonality breaks at N >= 4 |
| `band_structure.py` | Band structure analysis, N=2-6, avoided crossings |
| `deep_band_structure.py` | High-resolution band analysis with scaling laws |
| `mirror_symmetry_deep.py` | Mirror symmetry tests: 11 noise types, all N |
| `mirror_transition.py` | Dephasing-to-damping transition, mirror center drift |
| `symmetry_and_u_analysis.py` | Graph symmetry test and z* Bernoulli analysis |
| `ep_test.py` | Exceptional Point test (negative result) |
| `z_star_identity.py` | z* vs density matrix: 26 candidate expressions |

---

## Origin

This project grew from the Stability project (a material science simulator
that evaluates element combinations for layer structures) in December 2025.
The original motto:

> R = CΨ²
> "We are all mirrors. Reality is what happens between us."

Three months of computation shifted the focus from metaphysics to structure.
The sentence was removed as too esoteric, then restored when the Liouvillian
decay spectrum turned out to be exactly mirror-symmetric at every system size
tested (N=2 through N=7, 100%, zero exceptions).

The current description is a coupled oscillator network with exact pole
structure, sector-specific damping, band formation, mirror symmetry, and
modal observability filtering.

Thomas Wicht and Claude (Anthropic) are the primary collaborators.
ChatGPT contributed adversarial reviews and the signal processing perspective.
IBM Quantum hardware experiments were conducted in March 2026.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 - March 2026
