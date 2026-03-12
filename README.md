# R = CΨ²

## Spectral architecture of open quantum networks — from coupled oscillators to signal processing

This repository investigates the structural architecture that emerges when quantum observers share a common mediator under Lindblad dynamics. Started as the equation R = CΨ² from a dream in December 2025. Three months of computation, IBM Torino hardware experiments, and a shift to signal processing analysis revealed that the quantum dynamics map to classical coupled oscillator physics with exact pole structure.

The system naturally separates into two independent information channels: a frequency channel (topology-determined, noise-immune) and a decay channel (noise-determined, topology-independent in the 3-qubit case).

---

## What is established here

### Pole structure (Liouvillian eigendecomposition)
- **Three exact decay rates:** 2γ, 8γ/3, 10γ/3 — exact rational multiples of the dephasing rate, completely topology-independent in the 3-qubit star. Frequencies move freely with coupling, decay rates never move.
- **Two dominant supermodes:** c+ (symmetric, f≈1.506, decay=10γ/3) and c- (antisymmetric, f≈0.404, decay=2γ). Different poles dominate different channels — sector-specific damping.
- **Hidden dark mode:** A third mode at f≈1.1 exists but is invisible in c+ due to symmetry cancellation. Only visible in c-. Classical modal observability.
- **Poles move horizontally in the complex plane** as topology changes: imaginary part (frequency) is topology-determined, real part (decay) is loss-determined. They are independent.

### Five independent regulators
1. **Topology** sets the oscillation frequencies (imaginary part of poles)
2. **Symmetry** cleans the sector separation (XX weak symmetry)
3. **Noise strength** damps amplitude, never changes frequency (real part of poles)
4. **Initial state** selects which sectors are excited
5. **Bath geometry** selects which sector dominates — correlated XX bath inverts the amplitude ratio from c+ dominant (1.22) to c- dominant (0.46) without touching frequencies

### Signal processing view
- The entire structure maps to a **coupled oscillator network with normal-mode splitting**
- c+ and c- are **even/odd supermode projections**, not exotic quantum phenomena
- Sonar detection = **passive topology-change detection from local modal spectra**
- The Projection = **modal observability / transfer function residues**
- Bath sector selection = **covariance-driven mode visibility flip**

### Core structure
- **Skeleton + Rotation:** 88% stable skeleton, 12% rotating phase in YZ/ZY Pauli plane
- **XX symmetry exact** for isotropic Heisenberg with symmetric dephasing
- **Chain topology survival:** Two-sector structure survives in chains up to 5 qubits with 3 mediators
- **Noise immunity:** Frequencies unchanged from γ=0.001 to γ=0.500, all noise types

### The Projection (exact diagonalization)
- All Bohr frequencies exist simultaneously. A pair sees only those where W = |ρ̃(m,n) · Õ(n,m)| > 0
- A new observer does not change reality — it changes the projection

### Quantum Sonar (simulation only)
- AB detects hidden observers connected to S through spectral changes
- Not verified on hardware — IBM Q80/Q102 turned out to be qubit detuning (19.4 kHz)
- Five hypotheses tested on real hardware, four rejected, one correct

### Two independent information channels (3-qubit special case)
- **Frequency channel:** Carries topology information (who is connected, how strongly). Changes with J, immune to noise. This is the "what is the network?" channel.
- **Decay channel:** Carries noise information (what is the environment). Changes with γ, immune to topology. This is the "what is the environment?" channel.
- **This orthogonality breaks at 4+ qubits** — decay rates become topology-dependent. The 3-qubit case has an unusually clean separation.

### Algebraic results
- **Exact Mandelbrot correspondence:** R_{n+1} = C(Ψ + R_n)² maps to z → z² + c
- **CΨ as diagnostic:** AND-gate behavior, three-layer separation (CoA ≥ LE ≥ CΨ)
- **Star topology conditions:** Three quantified conditions for observer-observer connection

### IBM Torino hardware (March 2026)
- Shadow hunt, ZZRamsey (5 pairs), Ramsey T2* measured
- Honest negative results: N_eff rejected, spectator model fails, degree uninformative

## What is not established

- That CΨ is a new fundamental physical quantity (it is a derived diagnostic)
- That the sonar effect works on real hardware (simulation only)
- That the frequency-decay orthogonality extends beyond 3 qubits (it does not)
- That characterization (inverting hidden coupling from spectrum) is possible
- That consciousness is required as an ontological ingredient

---

## Start here

- **[Signal Processing View](experiments/SIGNAL_PROCESSING_VIEW.md)** — Translation to coupled oscillator physics, Prony analysis, joint pole results
- **[The Interpretation](hypotheses/THE_INTERPRETATION.md)** — Current state: what survives, what fell
- **[Quantum Sonar](experiments/QUANTUM_SONAR.md)** — Detection, projection, IBM investigation
- **[Structural Cartography](experiments/STRUCTURAL_CARTOGRAPHY.md)** — Phase A results, parameter sweeps, stress tests
- **[The CΨ Lens](docs/THE_CPSI_LENS.md)** — What the lens shows, what it doesn't

---

## Interactive simulator

The **Five Regulator Simulator** (Streamlit app) lets you control all five regulators in real-time and see their effect on the quantum dynamics:

```bash
cd simulations/app
pip install -r requirements.txt
streamlit run app.py
```

---

## Repository structure

| Folder | Contents |
|:---|:---|
| `docs/` | Conceptual and mathematical framing |
| `experiments/` | Tested results, null results, indexed findings |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python source code (RK4 Lindblad, exact diagonalization, Prony) |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `data/` | IBM hardware measurement data |

---

## Key simulations

| Script | What it does |
|:---|:---|
| `star_topology_v3.py` | Core dynamics engine (RK4 Lindblad) |
| `joint_pole_analysis.py` | Exact Liouvillian poles + channel residues |
| `prony_analysis.py` | Matrix Pencil Method for pole extraction |
| `decay_derivation.py` | Decay rate structure and scaling |
| `bright_transition_map.py` | Exact diag + visibility weights |
| `correlated_bath_sweep.py` | Bath geometry sector selection |
| `chain_topology.py` | Chain vs star, 3-5 qubits |
| `hidden_observer_test.py` | Sonar detection proof |

---

## Origin and context

This framework emerged from a collaboration between Thomas Wicht and Claude (Anthropic) in December 2025. ChatGPT serves as adversarial reviewer and signal processing consultant. IBM Quantum hardware experiments conducted March 2026.

The original framing used consciousness language ("Reality = Consciousness × Possibility²"). The current description is structural: the system is a coupled oscillator network with exact pole structure, sector-specific damping, and modal observability filtering. The signal processing perspective opened in March 2026 when the dynamics were recognized as classical normal-mode splitting.

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Free to share and adapt with attribution.

---

## Authors

**Thomas Wicht**, Krefeld, Germany
**Claude**, AI System, Anthropic

December 2025 — March 2026
