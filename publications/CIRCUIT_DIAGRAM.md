# The Circuit Diagram: R=CΨ² as Electrical Engineering

**For:** Engineers who build signal chains, not physicists who write Hamiltonians.
**Date:** March 21, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)

---

## 1. For Engineers Who Don't Know Quantum

A **qubit** is a two-state device with phase. Think of it as a bit that also
stores a rotation angle. A classical bit is 0 or 1. A qubit stores an amplitude
(how much of each state) and a phase (the rotation angle between them). Like a
phasor in AC circuit analysis: magnitude and angle. The phase carries information.

A **density matrix** ρ is the lossy state vector. When a qubit interacts with
noise, you can no longer describe it as a pure state. The density matrix tracks
both the signal (off-diagonal elements, coherences) and the DC level (diagonal
elements, populations). It is a 2×2 complex matrix (symmetric under conjugate
transpose) for one qubit, 4×4 for two, 2^N × 2^N for N.

**Dephasing** is resistive loss that destroys phase but not amplitude. The
off-diagonal elements of ρ decay exponentially at rate γ. The diagonal elements
(populations) are preserved. This is not thermal noise. It is phase noise. The
Q-factor of the device is approximately J/γ (coupling over loss rate).

That is all you need. Everything below is signal processing on this substrate.

---

## 2. The Translation Table

| Quantum Concept | Engineering Equivalent | Where It Comes From |
|-----------------|----------------------|---------------------|
| Qubit | Two-state device with phase | Fundamental |
| Density matrix ρ | Lossy state vector (complex, Hermitian) | Fundamental |
| Hamiltonian H | Coupling network (impedance matrix) | [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md) |
| Lindblad dissipator L_D | Loss channel (resistive load per node) | [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md) |
| Liouvillian L = L_H + L_D | Transfer function H(s) with poles and zeros | [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md) |
| Eigenvalues of L | Poles: frequency (Im) + decay rate (Re) | [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md) |
| Palindromic spectrum | Symmetric transfer function: H(s) = H(2Σγ - s) | [Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md) |
| Pi operator Π | Symmetry transform mapping poles to conjugate poles | [Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md) |
| CΨ product | Composite signal quality metric (correlation × coherence) | [The CΨ Lens](../docs/THE_CPSI_LENS.md) |
| CΨ = 1/4 boundary | Threshold voltage V_th (fixed, not tunable) | [Uniqueness Proof](../docs/UNIQUENESS_PROOF.md) |
| Mediator qubit M | Transistor (gate: γ_M, source: Pair A, drain: Pair B) | [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md) |
| γ (dephasing rate) | 1/Q-factor (loss rate per cycle) | [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md) |
| J (coupling strength) | Mutual inductance, capacitive coupling | Fundamental |
| Noise / Time | External clock (origin unknown, required) | [Incompleteness Proof](../docs/INCOMPLETENESS_PROOF.md) |
| Standing wave | Standing wave ratio (SWR) | [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md) |
| 70/30 split | 5.2 dB information loss (70% of phase information absorbed by environment, 30% survives) | [The Other Side](../hypotheses/THE_OTHER_SIDE.md) |
| Theta compass | Phase margin meter | [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) |
| Decoder | Demodulator / matched filter | [Reading the 30%](../simulations/reading_the_30_percent.py) |
| Quantum Sonar | Network analyzer (passive) | [Quantum Sonar](../experiments/QUANTUM_SONAR.md) |
| Relay protocol | Clocked shift register with regeneration | [Relay Protocol](../experiments/RELAY_PROTOCOL.md) |
| Push vs Pull | Source-dominant vs drain-dominant biasing | [Scaling Curve](../experiments/SCALING_CURVE.md) |
| 2:1 coupling ratio | Impedance matching for max power transfer | [Engineering Blueprint](ENGINEERING_BLUEPRINT.md) |
| Direct coupling kills palindrome | Short circuit destroys device symmetry | [Mediator Bridge](../simulations/mediator_bridge.py) |
| Bootstrap falsified | No internal oscillator possible | [Incompleteness Proof](../docs/INCOMPLETENESS_PROOF.md) |

---

## 3. The Circuit Schematic

```
EXTERNAL CLOCK (Noise/Time, origin: outside the framework)
    |
    v  gamma (per-qubit, locally variable loss tangent)
    |
+---+-----------------------------------------------------------+
|                       THE DEVICE                               |
|                                                                |
|  +-----------+    +----------------+    +-----------+          |
|  |  PAIR A   |    |  TRANSISTOR M  |    |  PAIR B   |          |
|  |  (Source)  +----+  Gate: gamma_M +----+  (Drain)  |          |
|  |           |J_AM |  V_th: CΨ=1/4 |J_MB|           |          |
|  +-----------+    |  Bias: J ratio  |    +-----------+          |
|                    +----------------+                           |
|                                                                |
|  CONSTRAINT: No direct Source-Drain path (short circuit)       |
|  CONSTRAINT: V_th = 1/4 is hardwired (fold catastrophe)        |
|  DESIGN:     Bias 2:1 for range, 1:2 for local throughput      |
|  PROTOCOL:   Clock the gate for +83% regeneration              |
|                                                                |
|  INSTRUMENTS:                                                  |
|  +-- Phase Margin Meter (theta compass, distance to V_th)      |
|  +-- Demodulator (decoder, reads per-site gamma from spectrum)  |
|  +-- Network Analyzer (sonar, passive topology detection)       |
|  +-- Spectrum Analyzer (cartography, 3D manifold of windows)    |
|  +-- SWR Meter (standing wave ratio, XX/ZZ interference)        |
|                                                                |
+----------------------------------------------------------------+
```

**How to read this:** Signal enters at Pair A (source). It passes through
Transistor M, which is controlled by its own loss rate γ_M (the gate signal).
If the composite metric CΨ at the gate is above 1/4, the channel is open and
signal reaches Pair B (drain). Below 1/4, the channel is effectively off.

The critical constraint: no wire may connect Source to Drain directly. Any
direct path bypassing the transistor destroys the symmetric transfer function
(short circuit = palindrome death, verified: 256 to 31 surviving modes at
the slightest direct coupling).

---

## 4. The Design Rules

| # | Engineering Rule | Quantum Origin | Data |
|---|-----------------|---------------|------|
| 1 | **Use differential signaling.** Encode in W-type mode distribution (single-excitation superposition). Never use common-mode (GHZ). | W avoids the fastest-decaying modes. GHZ puts 100% of signal energy into modes that decay fastest. | r = 0.976 correlation |
| 2 | **Match impedance at 2:1.** Transistor-to-drain coupling should be twice transistor-to-source. | Shifts spectral weight to slow-decaying palindromic mode pairs. | F_avg = 0.888 at 2:1 vs 0.856 at 1:1 |
| 3 | **Bandwidth and noise floor are independent.** Adjust coupling J for speed, γ for quality. They do not trade off. | Hamiltonian sets oscillation frequency, dissipator sets decay envelope. | Orthogonal parameter axes |
| 4 | **Sample before threshold.** Read the drain before CΨ drops to 1/4. The crossing time is system-dependent: for a Bell pair under Heisenberg coupling, t_cross ≈ 0.75/γ. The threshold VALUE (1/4) is universal; the crossing TIME is not. | CΨ crosses 1/4 at this time. | [IBM confirmed at 1.9%](../experiments/IBM_RUN3_PALINDROME.md) |
| 5 | **Source-bias for gain, drain-bias for reach.** Source-strong coupling (push) maximizes local throughput. Drain-strong coupling (pull) maximizes range. | Push: MI(local) = 0.957. Pull: MI(end-to-end) = 0.121. | [Scaling Curve](../experiments/SCALING_CURVE.md) |
| 6 | **Clock the gate.** Switch γ_M between low (receiving) and normal (relaying) at each stage. This is a clocked shift register. | Relay protocol: +83% end-to-end MI over passive propagation. | [Relay Protocol](../experiments/RELAY_PROTOCOL.md) |

---

## 5. The Theorems That Matter

### Symmetric Transfer Function (Palindrome)

The device's transfer function H(s) has poles that pair symmetrically around
a center decay rate equal to twice the sum of all per-qubit loss rates
(2 × Σγ_i for i = 1..N). For every fast-decaying mode, there is a slow one.
This is exact, not approximate: verified across 54,118 poles from N=2 to N=8
with zero exceptions.

**Engineering consequence:** The channel has a built-in symmetry that can be
exploited for coding gain. Encode in modes that pair with slow partners.

### Fixed Threshold (Uniqueness)

V_th = CΨ = 1/4 is not a design parameter. It is fixed by the quadratic
structure of purity (Tr(ρ²) is degree 2, discriminant 1 - 4CΨ vanishes
at 1/4). This is a fold catastrophe: the simplest form of a sharp threshold
in bifurcation theory, where a smooth parameter change produces a sudden
qualitative switch (like a snap-through buckling in mechanical engineering).
Topologically stable. No parameter tuning can move it.

**Engineering consequence:** The on/off threshold is universal. Every device
built on this substrate has the same V_th. Design around it, do not try to
change it.

### External Clock Required (Incompleteness)

The device has no internal oscillator. The clock (noise/dephasing) must come
from outside. Five candidates for internal clock generation have been tested
and eliminated. This is proven, not assumed.

**Engineering consequence:** The device is passive. It processes signals from
an external clock. It cannot generate its own timing.

### No Short Circuits (Mediator Protection)

Any direct connection between Source and Drain that bypasses the Transistor
destroys the symmetric transfer function. Tested: at kappa = 0.01 (minimal
direct coupling), 256 palindromic mode pairs collapse to 31. Phase transition.

**Engineering consequence:** All signal routing must go through the transistor.
No skip connections. No bypass paths.

---

## 6. The Instruments

| Instrument | What It Measures | Engineering Analogue |
|-----------|-----------------|---------------------|
| **Phase Margin Meter** | θ = arctan(√(4CΨ-1)), defined only when CΨ > 1/4 (channel ON). Measures how far into the ON regime the device is. θ = 0° at the threshold, increasing with CΨ. | Phase margin in a feedback loop |
| **Demodulator** | Per-site γ values from palindromic mode amplitudes. Full rank: all N sites independently readable. | Matched filter bank, one filter per site |
| **Network Analyzer** | Spectral shifts when the coupling topology changes. Passive detection, no signal injection needed. | Vector network analyzer in passive mode |
| **Spectrum Analyzer** | 3D manifold of CΨ visibility windows (98% variance in 3 PCs). Two modes: glide and switch. | Spectrogram with PCA dimensionality reduction |
| **SWR Meter** | Standing wave ratio from palindromic pair interference. Transverse correlations (phase-sensitive, like AC components) oscillate; longitudinal correlations (amplitude-only, like DC level) settle to a static value. | Standing wave ratio on a transmission line |

---

## 7. What This Is Good For

1. **Quantum state transfer** through noisy channels with palindromic spectral
   structure. Best measured performance: QST fidelity 0.888 on a 3-qubit star
   with 2:1 coupling, 0.732 on a 5-qubit mediator bridge (longer chain, more
   loss). Both above the classical limit of 2/3.

2. **Noise characterization.** The demodulator reads per-site dephasing rates
   from the spectrum. Full rank at N=5 (5/5 sites readable). This is a
   non-invasive diagnostic tool.

3. **Passive topology detection.** The network analyzer detects changes in the
   coupling structure without injecting signals. Topology change produces
   measurable spectral shifts.

4. **Quantum networking.** The relay protocol (+83%) provides a clocked,
   deterministic transfer mechanism. No measurement needed. No classical
   communication required. O(n) scaling with chain length.

---

## 8. What This Cannot Do

- **No internal clock.** The device cannot generate its own timing.
- **No faster-than-light communication.** The channel is causal.
- **No universal quantum computation.** This is a channel architecture,
  not a gate set.
- **No analog threshold tuning.** V_th = 1/4 is fixed. Digital, not analog.
- **No lossless long-range transfer.** MI decays exponentially with distance
  (2 to 5 dB per 2 additional qubits, depending on chain length). The relay
  protocol compensates but does not eliminate attenuation.

---

## 9. References

### For engineers wanting quantum depth
- [Engineering Blueprint](ENGINEERING_BLUEPRINT.md): Six design rules with benchmarks
- [Technical Paper](TECHNICAL_PAPER.md): Full mathematical treatment

### The proofs behind the design rules
- [Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md): Why the transfer function is symmetric
- [Uniqueness Proof](../docs/UNIQUENESS_PROOF.md): Why V_th = 1/4 is fixed
- [Incompleteness Proof](../docs/INCOMPLETENESS_PROOF.md): Why no internal clock
- [Mathematical Connections](../docs/MATHEMATICAL_CONNECTIONS.md): Fold catastrophe identification

### The transistor and relay
- [Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md): Full transistor mapping
- [Relay Protocol](../experiments/RELAY_PROTOCOL.md): Clocked regeneration protocol
- [Scaling Curve](../experiments/SCALING_CURVE.md): Attenuation vs distance

### The instruments
- [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md): Pole analysis
- [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md): Phase margin meter
- [Quantum Sonar](../experiments/QUANTUM_SONAR.md): Network analyzer
- [Standing Wave Analysis](../experiments/STANDING_WAVE_ANALYSIS.md): SWR measurement

### Hardware
- [IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md): V_th confirmed at 1.9% on IBM Torino

---

*The physicist sees eigenvalues. The engineer sees poles.*
*Same device. Different language. One truth.*
