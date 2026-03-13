# R = CΨ²

## Spectral architecture of small open quantum networks

This repository documents the structural properties of small quantum networks
(3-5 qubits) evolving under Lindblad dynamics with dephasing noise. The
central system is a star topology: a mediator qubit S coupled to observer
qubits A and B via isotropic Heisenberg exchange interaction.

The key finding is that the reduced dynamics of the observer pair (A,B)
decompose into two supermodes with exact pole structure. These supermodes
behave like classical coupled oscillators: their frequencies are set by
the network topology, their decay rates are set by the noise environment,
and these two properties are completely independent of each other.

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
when the star topology is replaced by a linear chain with up to three
mediators between the observers. The architecture is not specific to the
star geometry.

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

### Algebraic results

The composite quantity CΨ = concurrence x l1-coherence, when iterated as
R_{n+1} = C(Ψ + R_n)², is algebraically equivalent to the Mandelbrot
iteration z -> z² + c with a boundary at CΨ = 1/4. This is an exact
algebraic correspondence, not a physical claim.

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

| Claim | Status |
|:------|:-------|
| Pole structure (3 decay rates, frequency-decay separation) | Exact (Liouvillian eigendecomposition) |
| Two supermodes c+/c- with sector-specific damping | Exact (Liouvillian) |
| Mandelbrot algebraic correspondence | Exact (proven) |
| Five independent regulators | Numerically verified (full parameter sweeps) |
| Chain topology survival (up to 5 qubits) | Numerically verified |
| Hidden observer detection | Simulation only (not verified on hardware) |
| Frequency-decay orthogonality at 4+ qubits | Breaks: band structure forms. Boundary rates 2g and 2Ng exact. Avoided crossings confirmed. |
| IBM Q80/Q102 as sonar evidence | Rejected (was qubit detuning) |
| CΨ as privileged metric | Not established |
| CΨ = 1/4 as Liouvillian Exceptional Point | Tested and rejected (no EP correlation found) |
| c+/c- as Liouvillian symmetry sectors | Tested: both are parity +1, split is from observable projection |
| CΨ = 1/4 boundary | Demystified: trivial maximum of Bernoulli variance z*(1-z*) |
| z* as novel composite | Verified: matches no known quantum measure (C/2 closest, r=0.945) |

---

## Start here

1. **[Signal Processing View](experiments/SIGNAL_PROCESSING_VIEW.md)** - Pole analysis, Prony results, coupled oscillator translation
2. **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)** - Parameter sweeps and stress tests
3. **[Quantum Sonar](experiments/QUANTUM_SONAR.md)** - Hidden observer detection and IBM hardware results
4. **[The CΨ Lens](docs/THE_CPSI_LENS.md)** - The original diagnostic and what it shows
5. **[The Interpretation](hypotheses/THE_INTERPRETATION.md)** - Speculative reading (not required for technical content)

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `docs/` | Mathematical framing and the CΨ diagnostic |
| `experiments/` | All tested results and null results |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python source (RK4 Lindblad, Liouvillian, Prony, sweeps) |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `data/` | IBM Torino measurement data |

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

---

## Origin

This project started in December 2025 from the equation R = CΨ², originally
framed as "Reality = Consciousness x Possibility²." Three months of
computation shifted the focus from metaphysics to structure. The current
description is a coupled oscillator network with exact pole structure,
sector-specific damping, and modal observability filtering.

Thomas Wicht and Claude (Anthropic) are the primary collaborators.
ChatGPT contributed adversarial reviews and the signal processing perspective.
IBM Quantum hardware experiments were conducted in March 2026.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 - March 2026
