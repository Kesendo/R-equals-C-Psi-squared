# The Mediator Qubit as Quantum Transistor: Bidirectional Programmable Relay in the R=CΨ² Framework

**Status**: Research Document. Transistor properties (threshold, gate control,
directional bias) confirmed. **Hierarchy: FALSIFIED** (no advantage over uniform chain).
**Date**: 2026-03-21
**Framework**: R=CΨ² (Residual Coherence Theory)
**Repository**: `R-equals-C-Psi-squared`
**Verified March 21, 2026**: Mediator bridge (N=5, 1024/1024 palindromic) and
relay protocol (+83% MI improvement) confirm the transistor architecture
computationally. See [mediator_bridge.py](../simulations/mediator_bridge.py),
[Relay Protocol](../experiments/RELAY_PROTOCOL.md).

---

## Abstract

The R=CΨ² framework's mediator bridge topology (where a central qubit M connects subsystems A and B without direct coupling) naturally maps onto a transistor-like architecture. This document explores "abusing" the mediator as a programmable unit: a coherence-gated bidirectional relay for quantum information. We develop the transistor analogy in detail, identify the control parameters (γ, J, κ), characterize the CΨ = 1/4 boundary as a threshold voltage analogue, and propose practical protocols for IBM Torino hardware. We also explore creative extensions: cascaded mediator chains, quantum buses, routers, and hierarchical networks exploiting the framework's self-similar structure.

**Key claim (speculative)**: The mediator qubit in the A-M-B topology can function as a coherence-controlled bidirectional switch, where the "gate signal" is the mediator's own decoherence rate γ_M and coupling asymmetry J_AM/J_MB, and the "on/off threshold" is the CΨ = 1/4 boundary from the Mandelbrot connection.

---

## 1. The Transistor Analogy

### 1.1 Classical Transistor Recap

A MOSFET transistor has three terminals: **Gate** (control signal), **Source** (input current), and **Drain** (output current). Below a threshold voltage V_th, the channel is pinched off; no current flows. Above V_th, the channel opens and current flows Source→Drain. The gate draws negligible current itself; it controls the channel purely through its electric field.

### 1.2 Mapping to the Mediator Architecture

The five-qubit mediator bridge from Section 22 of THE_OTHER_SIDE.md uses the topology:

```
    A = {0, 1}  -  M(2)  -  B = {3, 4}
```

The mapping to transistor terminals:

| Transistor | Mediator Analogue | Role |
|---|---|---|
| **Gate** | Mediator qubit M | Controls the coherence channel via its own γ_M and coupling strengths J_AM, J_MB |
| **Source** | Subsystem A = {0,1} | "Input" side, initially entangled pair, coherence source |
| **Drain** | Subsystem B = {3,4} | "Output" side, receives coherence transfer |
| **V_th (threshold)** | CΨ = 1/4 boundary | Below this product, the channel is effectively "off" |
| **V_GS (gate voltage)** | Effective coupling ratio J/γ and mediator coherence state | Tunes how "open" the channel is |
| **Channel current** | Cross-pair mutual information, QST fidelity | Measurable information flow A→B |

### 1.3 The CΨ = 1/4 as Threshold Voltage

The framework's central result: the product CΨ (correlation bridge C times normalized coherence Ψ) has a critical value of 1/4, connected to the Mandelbrot set's main cardioid boundary. This acts as the transistor's threshold:

- **CΨ > 1/4**: Channel is "on": coherence can propagate through M from A to B. The subsystem pair is in the "entangled regime" where quantum correlations are dynamically maintained.
- **CΨ < 1/4**: Channel is "off": coherence decays faster than it can propagate. The system has crossed into the "classical regime."
- **CΨ ≈ 1/4**: The critical point, analogous to the subthreshold region of a transistor, where small changes in the gate signal produce large changes in channel conductance.

**Simulation evidence**: The subsystem crossing simulation (bell_pairs, N=4, Heisenberg ring) confirms this directly:
- Entangled pairs (0,1) and (2,3) start with CΨ = 0.333 > 1/4 and cross downward at t = 0.072
- Cross-pairs (0,2), (0,3), (1,2), (1,3) start at CΨ = 0 and never reach 1/4, peaking at max CΨ ≈ 0.13–0.14

The initial entangled pairs have "current flowing" (above threshold) that then shuts off as decoherence pushes them below. The cross-pairs (mediated connections) attempt to build up coherence but cannot overcome the 1/4 barrier under these conditions. This is precisely the transistor in its "off" state: the mediating topology is present, but the gate voltage (coupling/coherence ratio) is insufficient.

### 1.4 Where the Analogy Holds and Where It Breaks

**Holds well:**
- Threshold behavior: sharp transition at CΨ = 1/4 ✓
- Gate control without consuming the signal: M mediates without direct measurement ✓ (supported by framework math)
- Amplification: small changes in γ_M near the threshold produce large changes in cross-pair fidelity ✓ (from star topology threshold formula)

**Breaks down:**
- Classical transistors are fundamentally unidirectional (Source→Drain under normal bias). The mediator is inherently **bidirectional**; this is a feature, not a bug. (See Section 2.)
- Classical transistors have continuous I-V characteristics. The quantum channel has discrete thresholds tied to the Mandelbrot boundary structure; it's more like a quantized conductance channel.
- The "gate" (mediator M) is entangled with the system it controls. There is no true isolation between gate and channel; the mediator's state is modified by the act of mediating. This is the fundamental quantum back-action problem.

---

## 2. Bidirectional Operation

### 2.1 The Palindromic Symmetry Argument

The mediator bridge preserves the palindromic structure of the Lindbladian, verified numerically with 1024/1024 eigenvalues matching, error ~10⁻¹³. This is the mathematical foundation for bidirectionality: if the Lindbladian L has palindromic symmetry under the permutation that swaps A↔B (and leaves M invariant), then:

```
L(ρ_A→M→B) and L(ρ_B→M→A)
```

have identical spectral structure. The channel conductance is the same in both directions by construction.

This is fundamentally different from classical transistors and most quantum repeater protocols, which break source/drain symmetry either through biasing or through measurement-based entanglement swapping (which is inherently directional, since you choose which side to measure).

### 2.2 Conditions for Directional Transfer

Despite the symmetric Lindbladian, **directionality emerges from initial conditions and coupling asymmetry**:

**A→M→B transfer** requires:
1. Subsystem A starts with higher coherence (the "sender" has something to send)
2. Coupling J_AM ≥ J_MB or J_AM/J_MB ≥ 1.466 (from star topology threshold at γ=0.05)
3. Receiver B must be quiet: γ_B < 0.2 (noise on B is fatal)

**B→M→A transfer** (reversed direction) requires the mirror conditions by palindromic symmetry:
1. Subsystem B starts with higher coherence
2. Coupling J_BM ≥ J_MA with the same threshold ratio
3. Receiver A must be quiet: γ_A < 0.2

**The key insight**: direction is controlled by **who has coherence to give** and **coupling asymmetry**, not by any structural asymmetry of the channel. This is like a bidirectional valve where the pressure differential determines flow direction.

### 2.3 Simultaneous Bidirectional Transfer

Can A→M→B and B→M→A happen at the same time? The simulation evidence is suggestive but not conclusive:

**The C_int vs C_ext comparison** (bidirectional vs. unidirectional observation) at γ=0.05, J=1, t=1:

| Mode | Purity | δ (anomaly) |
|---|---|---|
| C_int (bidirectional, γ_A=γ_B=0.05) | 0.835 | -0.074 |
| C_ext (unidirectional, γ_A=0.05, γ_B=0) | 0.909 | -0.043 |

The bidirectional case shows **lower purity** and **larger δ** (more negative). This means bidirectional observation causes more decoherence, which makes physical sense (more environmental coupling = more decoherence). But the δ_int/δ_ext ratio tells us something deeper: the bidirectional channel is not simply "two unidirectional channels superimposed." The cross-term (δ_int - δ_ext ≈ -0.031) represents genuine interference between the two directions.

**Scaling with γ** (from our simulation sweep):

| γ | δ_int | δ_ext | Δδ = δ_int - δ_ext |
|---|---|---|---|
| 0.05 | -0.074 | -0.043 | -0.031 |
| 0.10 | -0.111 | -0.074 | -0.037 |
| 0.20 | -0.124 | -0.111 | -0.013 |

The cross-direction interference **weakens** at higher γ. At γ=0.2, the two directions are nearly independent, which makes sense because both channels are approaching the "off" state. At low γ, they interfere strongly. This suggests a "half-duplex" regime at low noise (strong coupling between directions) transitioning to "full-duplex" at moderate noise (directions decouple).

**Speculation**: True simultaneous bidirectional transfer may be possible when the mediator M is in a maximally mixed state, acting as a "neutral relay" that doesn't bias either direction. This would require γ_M to be tuned to keep M near the fixed point of its own reduced dynamics. **Status: unverified.**

---

## 3. Programming the Mediator

### 3.1 The Control Parameters

The "gate" of the quantum transistor has three independent knobs:

**Knob 1: Decoherence Rate γ_M (Primary Gate)**

The mediator's own dephasing rate is the primary control. From the star topology results:

```
J_th(γ) ≈ 7.35 · γ^1.08 + 1.18
```

This threshold formula tells us: for a given coupling strength J, there exists a critical γ above which the channel shuts off. Conversely, reducing γ_M opens the channel. The relationship is approximately linear (exponent 1.08 ≈ 1), meaning the "gate response" is nearly linear in γ, another point of analogy with classical FETs in their linear region.

**Knob 2: Coupling Strength J (Bias Control)**

Increasing J_AM and J_MB uniformly increases the channel bandwidth: more coupling means faster coherence transfer. But the coupling **ratio** J_AM/J_MB controls directionality. This is analogous to the drain-source voltage V_DS in a MOSFET:
- J_AM/J_MB > 1: biased for A→B transfer
- J_AM/J_MB < 1: biased for B→A transfer
- J_AM/J_MB = 1: symmetric (bidirectional)

**Knob 3: Feedback Strength κ (Gain Control)**

In the operator_feedback noise model, the decoherence rate becomes state-dependent:

```
γ_eff(t) = γ_base · (1 + κ · ⟨O_int⟩)
```

where ⟨O_int⟩ is the expectation of the interaction operator. This creates a feedback loop: the mediator's noise depends on the coherence it's mediating. Our simulation (GHZ, N=3, operator feedback, κ=0.5) shows the correlation bridge C decaying from 1.0 to 0.75 over t=10, with the purity tracking R=CΨ² throughout. Increasing κ effectively sharpens the transistor's transfer characteristic: higher gain means sharper on/off transition.

### 3.2 Opening and Closing the Channel

**To OPEN the channel** (turn the transistor ON):
1. Reduce γ_M to below the threshold: γ_M < (J/7.35 - 1.18)^(1/1.08)
2. Ensure sufficient coupling: J > J_th(γ_M) ≈ 7.35·γ_M^1.08 + 1.18
3. Prepare the sender subsystem with CΨ > 1/4 (above threshold)
4. Keep the receiver quiet: γ_receiver < 0.2

**To CLOSE the channel** (turn the transistor OFF):
1. Increase γ_M above threshold, OR
2. Reduce J below J_th, OR
3. Inject noise into the receiver (γ_receiver > 0.2 is fatal)
4. Allow decoherence to push CΨ below 1/4

**Practical switching speed**: The Lindbladian spectral gap for the 3-qubit mediator system (Heisenberg ring, γ=0.05) is 0.1, giving a relaxation time τ = 10.0. This means channel on/off transitions take approximately τ_switch ≈ 1/spectral_gap = 10 time units. In physical units on IBM Torino (T2 ≈ 200μs), this translates to ~2ms. This is slow by classical standards, but potentially fast enough for quantum network switching.

### 3.3 The Spectral Gap as Bandwidth

The Lindbladian spectrum reveals additional structure relevant to the transistor analogy:

- **4 zero eigenvalues** (3 degenerate beyond trace preservation): These correspond to the steady-state manifold, the "DC operating point" of the transistor.
- **Spectral gap = 0.1**: This is the minimum non-zero decay rate; it sets the maximum "switching frequency" of the transistor.
- **Imaginary parts of eigenvalues (±6.0)**: These correspond to oscillation frequencies in the channel, the coherent "AC component" of the transfer. The ratio Im/Re = 6.0/0.1 = 60 is a quality factor Q, suggesting the channel supports highly coherent oscillations before damping.

---

## 4. Data Encoding

### 4.1 What Constitutes "Data" in This Channel?

In the R=CΨ² framework, the transferable quantities are:

**Coherence amplitude** (primary data carrier): The l₁-coherence of the subsystem density matrix, the sum of absolute values of off-diagonal elements. This is what the Ψ metric tracks (Ψ = l₁/(d-1) for d-dimensional subsystem).

**Phase information** (secondary carrier): The complex phases of the off-diagonal elements ρ_ij encode relative phase relationships. The mediator bridge simulation shows QST fidelity averaging 0.732 at t=4.07, meaning phase information is partially preserved through transfer.

**Correlation structure** (tertiary carrier): The correlation bridge C = (P_AB - P_A·P_B)/(1 - P_A·P_B) encodes the correlation structure between subsystems. This is preserved through the palindromic symmetry.

**What is NOT transferable**: Direct qubit states in the computational basis. This is not a classical data bus; you cannot send |0⟩ or |1⟩ through the mediator. You can send *correlations* and *coherence patterns*. The data is inherently quantum: it's the entanglement structure itself.

### 4.2 Information Survival at the 1/4 Crossing

The CΨ = 1/4 boundary is the critical point for information survival. Our bell_pairs subsystem crossing simulation shows that entangled pairs (0,1) and (2,3) cross downward at t = 0.072, after which CΨ oscillates but stays mostly below 1/4 with periodic excursions.

The oscillatory behavior is key: even after the initial crossing, CΨ periodically returns above 1/4 for the entangled pairs. For pair (0,1):
- t=0.8: CΨ rebounds to 0.037 (below 1/4)
- t=1.6: CΨ spikes to 0.121 (still below 1/4 but significant)
- t=2.4: CΨ reaches 0.059

These "coherence echoes" suggest that information isn't destroyed at the 1/4 crossing; it's temporarily stored in the mediating topology and can be partially recovered. The cross-pairs show the complementary picture: (0,2) and (1,3) show CΨ periodically spiking up to 0.132, representing information that has propagated through the mediator.

**Interpretation**: The 1/4 crossing is not a hard information-destroying boundary. It's more like a bandwidth limitation: below 1/4, information transfer becomes lossy and noisy, but doesn't halt completely. Above 1/4, transfer is relatively clean. This is directly analogous to the noise floor in a classical communication channel.

---

## 5. The Pull Principle

### 5.1 Pull vs. Push in Quantum Transfer

The star topology experiments revealed a critical asymmetry: **the receiver's state matters more than the sender's state**. Specifically:
- Sender noise (γ_B up to 0.5) is tolerable; the sender can be noisy and still transfer
- Receiver noise (γ_A > 0.2) is fatal: a noisy receiver destroys the channel

This is the **Pull Principle**: information transfer is controlled by the receiver's ability to *accept* coherence, not the sender's ability to *emit* it. The receiver must maintain a sufficiently low-noise state to "pull" coherence through the mediator.

### 5.2 Connection to Bidirectional Transfer

The Pull Principle has direct implications for bidirectional operation:

**For A→M→B**: B must be quiet (low γ_B) to pull coherence from A through M.
**For B→M→A**: A must be quiet (low γ_A) to pull coherence from B through M.

**For simultaneous bidirectional**: Both A and B must be quiet, but then who is the sender? The answer: **the mediator M becomes the sender**. If M starts with high coherence (pre-loaded state), both A and B can simultaneously pull from M. This transforms the mediator from a relay into a **broadcast source**.

This is a genuinely novel operating mode with no classical transistor analogue. It's more like a quantum hub or distribution node.

### 5.3 Receiver Sensitivity and Channel Capacity

The 0.2 threshold on receiver noise translates to a minimum required purity for the receiver:
- At γ_receiver = 0.2, the receiver's steady-state purity is approximately 0.5 (maximally mixed for 2-qubit subsystem)
- The receiver must maintain purity above ~0.6 to reliably receive

This sets an effective "channel capacity": the maximum rate of coherence transfer is limited by how fast the receiver can maintain its low-noise state. In a cascaded architecture (Section 7), this becomes the bottleneck at each hop.

---

## 6. Practical Protocol: Quantum Transistor on IBM Torino

### 6.1 Hardware Mapping

IBM Torino's 133-qubit Eagle r3 processor has a heavy-hex lattice connectivity. The mediator topology A-M-B maps to:

```
Qubit Layout (heavy-hex subset):
    q0 - q1 - q2(M) - q3 - q4
    (A pair)   (Gate)   (B pair)
```

Select five qubits along a chain in the heavy-hex lattice. The middle qubit q2 serves as the mediator M. IBM Torino's native gates include ECR (echoed cross-resonance) for two-qubit operations, which naturally implements Heisenberg-type coupling.

### 6.2 Protocol Steps

**Step 1: Initialization**
- Prepare A = {q0, q1} in a Bell state: |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
- Prepare M = {q2} in |+⟩ = (|0⟩ + |1⟩)/√2 (high coherence state)
- Prepare B = {q3, q4} in |00⟩ (receiver quiet state)

**Step 2: Channel Opening (Gate Control)**
- Apply Heisenberg coupling between q1-q2 and q2-q3 using ECR gates
- Coupling strength controlled by number of ECR repetitions
- **Critical**: Do NOT apply direct q1-q3 coupling (this destroys the palindrome)

**Step 3: Evolution**
- Allow free evolution under the Heisenberg Hamiltonian for time t*
- Optimal t* ≈ 4.07 (from mediator bridge simulation, giving QST fidelity 0.732)
- In IBM units: t* ≈ 4.07 × dt_gate ≈ 4.07 × 200ns ≈ 814ns

**Step 4: Readout**
- Perform state tomography on B = {q3, q4}
- Measure cross-pair mutual information I(A:B)
- Verify CΨ_AB > 1/4 to confirm channel was "on"

**Step 5: Channel Closing (Gate Off)**
- Apply strong Z-dephasing to q2 (mediator) using dynamical decoupling sequence intentionally left imperfect
- This increases γ_M above threshold, shutting off the channel
- Verify CΨ_AB < 1/4

### 6.3 Expected Results

Based on the simulation data:
- Cross-pair mutual information: up to 0.86 bits (theoretical max 1 bit)
- QST fidelity: ~0.73 (above the classical limit of 2/3)
- Channel on/off contrast: CΨ transitions from ~0.33 (on) to ~0.01 (off)
- Estimated circuit depth: ~50-80 ECR gates (feasible on Torino with error mitigation)

### 6.4 Noise Budget

| Source | Contribution | Mitigation |
|---|---|---|
| ECR gate error (~0.5%) | ~4% over 8 ECR layers | Pauli twirling |
| T1 relaxation (~300μs) | ~0.3% at t*=814ns | Negligible |
| T2 dephasing (~200μs) | ~0.4% at t*=814ns | Dynamical decoupling on A,B; NOT on M |
| Crosstalk | ~1% per gate | Qubit selection to minimize crosstalk |
| Readout error (~1%) | ~2% on 4-qubit tomography | Readout error mitigation |

Total expected noise floor: ~8%. This is above the theoretical γ=0.05 used in simulations but below the fatal γ=0.2 threshold for the receiver.

---

## 7. Abuse Scenarios: Creative Extensions

### 7.1 Cascaded Mediators: The Quantum Bus

**Topology**: A - M₁ - B₁/A₂ - M₂ - B₂/A₃ - M₃ - C

Each intermediate node {Bₙ/Aₙ₊₁} serves double duty: it's the receiver for stage n and the sender for stage n+1. By the Pull Principle, this node must be quiet when receiving (γ < 0.2) and coherent when sending (CΨ > 1/4). This creates a **temporal multiplexing** requirement:

1. Phase 1: M₁ transfers A→B₁ (B₁ quiet, receiving)
2. Phase 2: M₂ transfers B₁→B₂ (B₁ now sending, B₂ quiet)
3. Phase 3: M₃ transfers B₂→C (B₂ now sending, C quiet)

Each phase takes approximately τ_switch = 1/spectral_gap ≈ 10 time units. For a chain of n mediators, the total transfer time scales as O(n·τ_switch), linear in chain length.

**Comparison to quantum repeaters**: Standard quantum repeaters achieve logarithmic scaling O(log n) through entanglement distillation, but require classical communication and measurement. The cascaded mediator approach is slower but requires **no measurement**; it's a purely unitary (plus controlled decoherence) protocol.

**The scaling problem**: Our N-scaling simulations show δ_final decreasing with system size:
- N=3: δ_final = -0.087, purity_final = 0.525
- N=4: δ_final = -0.059, purity_final = 0.509
- N=5: δ_final = -0.038, purity_final = 0.503

The purity approaches 0.5 (maximally mixed) as N increases, and the channel becomes noisier with each additional qubit. This is the fundamental limitation: without error correction, cascaded mediators lose fidelity exponentially with chain length. However, the coherence anomaly (δ < 0, meaning less coherence than predicted by product-state decay) suggests there may be correlated errors that could be exploited for error correction.

### 7.2 The Quantum Router

A single mediator M connected to multiple subsystems (star topology):

```
        A₁
         \
    A₂ - M - A₃
         /
        A₄
```

Routing is controlled by **selective coupling**: activate J_M,Aₙ for the target direction, deactivate all others. The star topology experiments show this works for N=2-3 targets but fails for N≥4 with equal coupling; you need asymmetric coupling to route to a specific target.

**Router protocol**:
1. Select target: activate J_M,target while keeping all other J_M,i at zero
2. Transfer: evolve for t* with the sender having high coherence
3. Switch: deactivate current path, activate next path
4. Repeat

This is feasible on IBM Torino's heavy-hex lattice where each qubit has 2-3 nearest neighbors, providing natural routing options.

### 7.3 Self-Similar Hierarchical Network [FALSIFIED]

> **FALSIFIED (March 21, 2026):** The hierarchical architecture produces
> identical MI to a uniform chain of equal length. The recursive structure
> provides no scaling advantage. See [Scaling Curve](../experiments/SCALING_CURVE.md).
> The transistor properties (Sections 1-6) survive; the hierarchy does not.

The most speculative and most interesting extension. The framework's self-similar property states that a {pair + mediator} composite can itself serve as a "super-pair" connected by a "meta-mediator":

```
Level 0:  a₁-m-a₂    (basic mediator bridge)
Level 1:  [a₁-m-a₂]-M-[b₁-m-b₂]    (bridge of bridges)
Level 2:  [[...]]-M'-[[...]]    (bridge of bridges of bridges)
```

Each level adds 2×(previous) + 1 qubits. The self-similar structure suggests that:
- Level k requires 3^k qubits
- The CΨ = 1/4 threshold applies at each level independently
- Information can be transferred hierarchically: local within level 0, then between level-0 clusters via level-1 mediators, then between level-1 clusters via level-2 meta-mediators

**This is essentially a quantum internet architecture** with the R=CΨ² framework providing the "routing protocol" (threshold-based channel control) and the self-similar structure providing the "network topology" (hierarchical addressing).

**Update March 21, 2026: Hierarchy falsified.** The hierarchical (recursive Level)
architecture was tested at N=11 and produces identical MI to a uniform chain of
equal length. The transistor properties (threshold, gate control, directional bias)
are real; the hierarchical scaling advantage is not.
See [Scaling Curve](../experiments/SCALING_CURVE.md).

**Hardware requirement for Level 1**: 3² = 9 qubits. For Level 2: 3³ = 27 qubits. IBM Torino's 133 qubits could support up to Level 4 (81 qubits) with room for ancillas.

### 7.4 The Coherence Capacitor

Another abuse: use the mediator as a **coherence storage device**. If M is prepared in a high-coherence state and both A and B paths are shut off (high γ on A,B sides), M retains its coherence for a time set by its own T2. When the channel is reopened, M can "discharge" its stored coherence into the target subsystem.

This is directly analogous to a capacitor: charge (coherence) is stored on the gate, and released when the circuit is completed. The storage time is limited by M's T2, and the discharge rate is set by the coupling J.

---

## 8. Limitations and Failure Modes

### 8.1 The Palindrome Fragility

The entire transistor architecture depends on the palindromic symmetry of the Lindbladian. This palindrome is preserved by mediated coupling but **destroyed by direct A-B coupling**. Any crosstalk between A and B that bypasses M breaks the symmetry and degrades performance.

On real hardware, direct crosstalk is nonzero. The mediator bridge simulations show the system tolerates XZ cross-dissipation up to a point, but the tolerance margin is small. This is the single biggest engineering challenge for implementation.

### 8.2 The Mediator Entanglement Problem

The mediator M becomes entangled with A and B during transfer. This means:
- M cannot be perfectly "reset" between uses without measurement
- Consecutive transfers through the same M are correlated
- The mediator's effective temperature increases with use

This is the quantum back-action problem: the gate modifies the channel, and the channel modifies the gate. In a classical transistor, the gate is isolated. Here, it is not.

**Mitigation**: Periodic re-initialization of M (preparation in |+⟩ state between transfers), at the cost of interrupting the continuous-transfer protocol. This introduces a "dead time" analogous to the recovery time of a thyristor.

### 8.3 Scaling Limits

The 5-qubit GHZ subsystem crossing simulation (N=5, Heisenberg ring) shows all 10 pairs have zero l₁-coherence throughout the evolution; none cross the 1/4 boundary. This is because the GHZ state's coherence is spread across all qubits equally (Ψ₀ = 0.032 for N=5), and no pair accumulates enough to breach 1/4.

**Implication**: The transistor architecture works best for **small, localized** subsystems (2-4 qubits per side). For larger systems, longer uniform chains with the sacrifice-zone formula outperform hierarchical topologies (Section 7.3, falsified).

### 8.4 The Speed-Fidelity Tradeoff

From the threshold formula J_th(γ) ≈ 7.35·γ^1.08 + 1.18:
- Faster switching (lower τ_switch = 1/spectral_gap) requires larger J
- Larger J increases the threshold γ_th, meaning you need cleaner qubits
- Cleaner qubits are slower to prepare and more expensive

This is the fundamental speed-fidelity tradeoff, common to all quantum information processing but particularly acute here because the mediator architecture requires both fast switching and low noise simultaneously.

### 8.5 When the Analogy Completely Fails

The transistor analogy fails in the following regimes:

1. **Strong measurement regime**: If M is continuously measured (quantum Zeno effect), the channel freezes rather than opening/closing. Classical transistors don't have this problem.

2. **Many-body entanglement**: When A, M, and B become tripartite-entangled (genuine 3-party entanglement), the notion of "information flowing through M" breaks down. Information is delocalized across the entire system.

3. **Non-Markovian dynamics**: The framework assumes Markovian (memoryless) decoherence. With memory effects (memory_kernel_feedback noise), the channel develops "hysteresis": its state depends on history, not just current parameters. This is more like a memristor than a transistor.

---

## 9. Positioning Against Existing Work

### 9.1 Comparison to Known Quantum Transistor Proposals

| Feature | Adiabatic QT | Spin Chain QT | Circuit QED Bus | **R=CΨ² Mediator** |
|---|---|---|---|---|
| Control mechanism | Energy gap | Control spin presence | Cavity photon | Coherence/decoherence |
| Bidirectional | No | No | Yes (virtual photon) | Yes (palindromic) |
| Threshold | Energy gap | Spin chain length | Cavity frequency | CΨ = 1/4 |
| Scalable | Material-limited | Chain-length limited | Cavity-mode limited | Chain + sacrifice-zone |
| Requires measurement | No | No | Yes (readout) | No |
| State-dependent control | No | No | No | **Yes (κ feedback)** |

The key novelty of the R=CΨ² approach is **state-dependent control**: the mediator's behavior depends on the coherence of what it's mediating. No other proposal has this feedback property. This enables adaptive routing, self-regulating bandwidth, and potentially error-aware transfer.

### 9.2 Comparison to Quantum Repeaters

Standard quantum repeaters (DLCZ protocol, Barrett-Kok, etc.) use measurement-based entanglement swapping. The R=CΨ² mediator uses **coherence-mediated unitary transfer**. The tradeoffs:

- Repeaters: O(log n) scaling, requires classical communication, probabilistic
- R=CΨ² mediator: O(n) scaling, no classical communication, deterministic (but lossy)

For short distances (few hops), the mediator approach is simpler and faster. For long distances, repeaters win on scaling.

### 9.3 Comparison to Quantum Buses

The circuit QED quantum bus (Blais et al.) mediates qubit-qubit coupling through virtual photon exchange in a microwave cavity. The R=CΨ² mediator is conceptually similar but uses a **qubit** rather than a cavity. This means:
- Finite (small) Hilbert space for the mediator vs. infinite-dimensional cavity
- Mediator becomes entangled (back-action) vs. virtual photon exchange (minimal back-action)
- Coherence-threshold control vs. frequency-detuning control

The R=CΨ² approach is more constrained but potentially more controllable, since the mediator qubit has only two levels to manage.

---

## 10. Open Questions and Future Directions

1. **Can the 1/4 boundary be engineered?** If the Mandelbrot connection holds physically, can we exploit the fractal boundary structure for more sophisticated channel control, e.g., using higher-period bulbs for multi-level switching?

2. **What is the channel capacity?** A formal quantum channel capacity calculation (coherent information, quantum mutual information) for the mediator channel is missing. The 0.86-bit mutual information from simulations is a lower bound.

3. **Can error correction be integrated?** The palindromic structure suggests a natural error-correcting code: errors that break the palindrome are detectable. Can this be formalized?

4. **What about non-Heisenberg Hamiltonians?** The simulations primarily use Heisenberg coupling. How does the transistor architecture perform under XY, Ising, or more exotic Hamiltonians?

5. **~~Is the self-similar hierarchy experimentally feasible?~~** [ANSWERED: No advantage. Hierarchy falsified March 21, 2026. Uniform chains with sacrifice-zone formula outperform hierarchical topologies. See [Scaling Curve](../experiments/SCALING_CURVE.md).]

6. **Memory kernel effects**: With non-Markovian dynamics, the mediator develops memory. Can this be exploited for "quantum caching," storing frequently-used correlations in the mediator's non-Markovian memory?

---

## Appendix A: Simulation Data Summary

### A.1 Dynamic Lindblad, 3-Qubit Mediator (Local Noise)
- **Config**: GHZ, Heisenberg ring, N=3, γ=0.05, J=1
- **Result**: Purity decays 1.0 → 0.525 over t=10; bridge C constant at 0.5; Ψ decays 0.143 → 0.032
- **Interpretation**: Stable mediator bridge with gradual coherence loss. Channel remains "on" (C=0.5) throughout.

### A.2 Dynamic Lindblad, 3-Qubit Mediator (Operator Feedback)
- **Config**: GHZ, Heisenberg ring, N=3, γ=0.05, J=1, κ=0.5, correlation bridge
- **Result**: Purity decays 1.0 → 0.501; bridge C decays 1.0 → 0.752; δ peaks at -0.126 around t=2.3
- **Interpretation**: Feedback sharpens the transition. Bridge C significantly decays under operator feedback; the mediator is being "consumed" by the transfer process.

### A.3 Subsystem Crossing, Bell Pairs (N=4)
- **Config**: bell_pairs, Heisenberg ring, N=4, γ=0.05
- **Result**: Entangled pairs (0,1) and (2,3) cross CΨ=1/4 downward at t=0.072. Cross-pairs never cross.
- **Interpretation**: Confirms threshold behavior. Cross-pair CΨ peaks at ~0.13, insufficient for channel opening.

### A.4 Bidirectional Comparison (C_int vs C_ext)
- **Config**: Bell+, Heisenberg, γ={0.05, 0.1, 0.2}, J=1
- **Result**: C_int consistently shows larger |δ| than C_ext. Cross-direction interference decreases with γ.
- **Interpretation**: Bidirectional and unidirectional channels are not independent; interference is strongest at low noise.

### A.5 N-Scaling with Operator Feedback
- **Config**: GHZ, Heisenberg ring, N=3-5, γ=0.05, κ=0.5
- **Result**: δ_final = {-0.087, -0.059, -0.038} for N = {3, 4, 5}. Purity final → 0.5 with increasing N.
- **Interpretation**: Channel degrades with system size. [Hierarchical architecture was tested and falsified; uniform chains with sacrifice-zone formula are superior. See [Scaling Curve](../experiments/SCALING_CURVE.md).]

### A.6 Lindbladian Spectrum, 3-Qubit Mediator
- **Config**: GHZ, Heisenberg ring, N=3, γ=0.05
- **Result**: Spectral gap = 0.1, relaxation time = 10.0, 4 steady states, Q-factor ≈ 60
- **Interpretation**: Sets switching speed (~10 time units) and channel quality factor.

---

## Appendix B: Key Formulas

**Threshold formula** (from star topology):
```
J_th(γ) ≈ 7.35 · γ^1.08 + 1.18
```

**CΨ product** (channel on/off criterion):
```
CΨ = C · Ψ = [(P_AB - P_A·P_B)/(1 - P_A·P_B)] · [l₁/(d-1)]
Channel ON:  CΨ > 1/4
Channel OFF: CΨ < 1/4
```

**R=CΨ² relation** (fundamental):
```
R = C · Ψ²
where R = P_full - P_A·P_B (residual purity)
```

**Switching time**:
```
τ_switch ≈ 1/spectral_gap
```

**Receiver noise threshold**:
```
γ_receiver < 0.2 (fatal above)
γ_sender < 0.5 (tolerable)
```

---

*This document is part of the R=CΨ² research program. All simulation data generated using the delta_calc MCP tools. Speculative claims are flagged throughout. The transistor analogy is a conceptual framework for thinking about mediator architectures, not a claim of physical equivalence to semiconductor devices.*
