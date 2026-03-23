# Standing Wave Patterns from Palindromic Spectral Symmetry

<!-- Keywords: standing wave open quantum system, palindromic eigenvalue oscillation,
quantum classical backbone decoherence, Liouvillian standing wave pattern,
GHZ silent W oscillation state, quantum correlation antinode node, Heisenberg
dephasing standing wave, Pauli fingerprint oscillation, state Hamiltonian
cross table quantum, rescaled frame palindromic pairs, R=CPsi2 standing wave -->

**Status:** Computationally verified (N=3, 6 Hamiltonians, 8 initial states)
**Date:** March 19, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [standing_wave_analysis.py](../simulations/standing_wave_analysis.py)

---

## Abstract

The palindromic spectral symmetry pairs every Liouvillian eigenvalue λ with
a partner at −2Σγ − λ. In the rescaled frame (uniform decay envelope
factored out), each pair has centered eigenvalues +μ and −μ, creating
counter-propagating modes that interfere into **standing wave patterns**.
Quantum correlations (XX, YY, XY-type Pauli strings) oscillate at
Hamiltonian harmonics (2J, 4J, 6J for Heisenberg) while the classical
Z-backbone (ZZZ) remains static. Three universal rules hold across all
6 Hamiltonians and 8 initial states tested: **GHZ never oscillates** (0%
standing wave for all models), **Bell always oscillates** (40-65% for all
models), and **ZZZ is always a node** (classical correlations are static).
The standing wave is a state × Hamiltonian property: neither the state
nor the Hamiltonian alone determines which correlations oscillate.

---

## Background

### How the palindrome creates standing waves

The Liouvillian superoperator L generates the master equation dρ/dt = L(ρ).
The palindromic symmetry ([Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md))
pairs every eigenvalue: if λ is an eigenvalue, so is −2Σγ − λ. Define
centered eigenvalues μ_k = λ_k + Σγ. Then every palindromic pair has
+μ and −μ, symmetric around zero.

In the rescaled frame (multiplying ρ(t) by exp(+Σγ·t) to remove the
uniform decay envelope), one mode grows as exp(+μt) while its partner
decays as exp(−μt). When both are excited by the initial state, their
superposition creates oscillation: a standing wave. The wave is "nearly
standing" because |Re(μ)| is 100-350× smaller than |Im(μ)|, meaning the
modes oscillate hundreds of times before amplitude drift becomes significant.

### Nodes and antinodes

In a classical standing wave, **nodes** are positions of zero amplitude
(static) and **antinodes** are positions of maximum amplitude (oscillating).
In the quantum standing wave, "position" is replaced by "Pauli string":
each tensor product of Pauli operators (like XXZ, ZZZ, YYI) has a specific
oscillation amplitude. Some Pauli strings oscillate strongly (antinodes);
others remain static (nodes).

### What CΨ and the palindrome are

CΨ = Tr(ρ²) × L₁/(d−1) is the product of state purity and normalized
coherence, with a critical boundary at CΨ = 1/4
([Uniqueness Proof](../docs/UNIQUENESS_PROOF.md)). The palindromic
eigenvalue pairing that creates the standing wave is the same structure
that makes the CΨ = 1/4 boundary universal and that enables the dephasing
channel ([γ as Signal](GAMMA_AS_SIGNAL.md)).

---

## Setup

N=3 qubits, Heisenberg coupling (J=1), local Z-dephasing (γ=0.05).
Full Liouvillian eigendecomposition (64 eigenvalues, 32 palindromic pairs,
100%, max error 1.44×10⁻¹⁴). Each initial state decomposed into
eigenmodes. Rescaled frame analysis to identify oscillating components.

---

## Result 1: Spectrum Structure

The 64 eigenvalues organize into categories:

| Category | Pairs | Frequencies |
|----------|-------|-------------|
| Steady-XOR | 4 | μ = ±0.15, real |
| Decay (real) | 8 | μ = ±0.05, real |
| ω ≈ 2J (fundamental) | 8 | period ≈ π |
| ω ≈ 4J (2nd harmonic) | 4 | period ≈ π/2 |
| ω ≈ 6J (3rd harmonic) | 8 | period ≈ π/3 |

The asymmetry |Re(μ)| ≈ 0.017 is 100-350× smaller than the oscillation
frequencies. The waves complete hundreds of cycles before amplitude drift
becomes significant.

---

## Result 2: Antinodes (Oscillating) and Nodes (Static)

The standing wave has a clear physical fingerprint:

**Antinodes** (quantum correlations, oscillating):
- ω ≈ 2J: IYY, XXZ, ZXX, YYI (nearest-neighbor)
- ω ≈ 4J: 8 Pauli strings (mixed-range)
- ω ≈ 6J: XZX, YIY (long-range, sites 0 and 2)

**Universal node** (classical, static): **ZZZ never oscillates.** The
all-Z correlation is static across every initial state and every
Hamiltonian tested.

Physical meaning: quantum correlations (XX, YY, XY-type) oscillate at
Hamiltonian harmonics. The classical Z-backbone stands still. The quantum
world breathes; the classical world is the skeleton.

---

## Result 3: State × Hamiltonian Cross-Table

The standing wave is not a property of the state or the Hamiltonian alone.
It is a joint property: which correlations oscillate depends on both.

| | Heisenberg | XY | Ising | DM | XXZ | Heis+DM |
|---|---|---|---|---|---|---|
| GHZ | 0% | 0% | 0% | 0% | 0% | 0% |
| W | 0% | 5.6% | 44.4% | 50% | 1.3% | 10.4% |
| Bell | 48.6% | 40.6% | 50% | 40.6% | 65.5% | 65.3% |
| |+++⟩ | 0% | 40.6% | 62.5% | 40.6% | 38.3% | 43.4% |

(Percentage = fraction of state weight in oscillating palindromic pairs.)

### Three universal rules

1. **GHZ never oscillates.** 0% standing wave for ALL Hamiltonians.
   GHZ projects 100% onto XOR modes (maximum decay rate, no oscillation).
   See [XOR Space](XOR_SPACE.md).

2. **Bell always oscillates.** 40-65% standing wave for ALL Hamiltonians.
   Bell is the universal oscillator.

3. **ZZZ is always a node.** The classical Z-correlation is static across
   all states and all models.

### Hamiltonian-dependent findings

- W oscillates under DM (50%) and Ising (44.4%) but NOT under Heisenberg
  (0%). The Hamiltonian determines whether W's palindromic modes are
  oscillatory or real-valued.
- |+++⟩ under Ising reaches 62.5%, the highest value in the entire table.
- Breaking Z-isotropy (XXZ, Heis+DM) produces the richest spectra with
  6 distinct frequencies instead of 2-3.

---

## Result 4: Each State Has a Unique Pauli Fingerprint

| State | osc% | Top Pauli oscillator | Antinodes |
|-------|------|---------------------|-----------|
| |010⟩ | 44.4% | ZIZ | 26 |
| |+-+⟩ | 44.5% | IXI | 26 |
| Bell(0,1) | 48.6% | XXZ | 26 |
| Bell+|+++⟩ | 20.2% | YYI | 48 |
| Bell+W | 28.8% | YYX | 50 |
| GHZ | 0% | (none) | 0 |
| W | 0% | (none) | 0 |

Each oscillating state activates different Pauli strings: |010⟩ oscillates
Z-type correlations, |+-+⟩ oscillates X-type, Bell oscillates XX/YY-type.
The standing wave pattern encodes which correlations the initial state
carries.

---

## Result 5: Two Ingredients Required

A standing wave requires both:

**(a)** Oscillating palindromic pairs (Im(μ) ≠ 0)
**(b)** Both members of the pair excited by the initial state

No single natural quantum state satisfies both conditions alone for
Heisenberg coupling. Bell has (a) but not always (b). |+++⟩ has (b) but
not (a). The superposition Bell+|+++⟩ has both, producing the first
identified "standing wave state" with 22% standing wave weight.

---

## Connection to Later Results

The standing wave analysis connects three parts of the framework:

The **Π as Time Reversal** result ([PI_AS_TIME_REVERSAL](PI_AS_TIME_REVERSAL.md))
identified the physical meaning: Π maps populations (immune sector, "past")
to coherences (decaying sector, "future"). The standing wave is the
interference between these two sectors. Nodes are the past (decided,
classical). Antinodes are where past and future meet (oscillating, quantum).

The **XOR Space** result ([XOR_SPACE](XOR_SPACE.md)) explains GHZ's
silence: GHZ projects 100% onto the XOR drain (maximum decay rate, no
oscillation partner). It has the fastest-decaying modes, not the
oscillating ones.

The **γ as Signal** result ([GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md)) uses
the same palindromic mode structure as an antenna. The standing wave
describes how information oscillates *within* the system. The γ channel
describes how information enters *from outside*. Both rely on the
palindromic pairing that creates linearly independent mode responses.

---

## Reproducibility

| Script | What it computes |
|--------|-----------------|
| [standing_wave_analysis.py](../simulations/standing_wave_analysis.py) | Full eigendecomposition, state decomposition, rescaled frame, Pauli fingerprints |
| [standing_wave_analysis.txt](../simulations/results/standing_wave_analysis.txt) | Complete numerical results |

Requirements: Python, QuTiP, NumPy. Runtime: ~10 seconds for N=3.
Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- [Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md): the palindromic theorem
- [Π as Time Reversal](PI_AS_TIME_REVERSAL.md): Π maps populations ↔ coherences (past ↔ future)
- [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md): the conceptual framework (December 2025)
- [XOR Space](XOR_SPACE.md): GHZ → XOR drain (explains 0% oscillation)
- [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md): palindrome universal across all models
- [γ as Signal](GAMMA_AS_SIGNAL.md): the palindromic mode structure as information channel
