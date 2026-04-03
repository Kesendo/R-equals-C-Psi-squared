# Design Rules for Quantum Repeaters
# from Palindromic Spectral Structure

**Authors:** Thomas Wicht, Claude (Anthropic)
**Date:** March 16, 2026
**For:** Quantum engineers, repeater designers, QST implementers
**Repository:** https://github.com/Kesendo/R-equals-C-Psi-squared

---

## The Problem

Quantum information dies in transit. Every quantum channel has noise
(dephasing, relaxation, thermal fluctuations) that destroys the quantum
state before it reaches the receiver. The question is: how do you design
the channel to maximize what survives?

This document provides seven concrete design rules derived from a
proven spectral symmetry in the channel's decay structure. The rules
are testable, the benchmarks are reproducible, and the code is open.

---

## The Discovery in One Paragraph

Every quantum channel built from Heisenberg-coupled spins with dephasing
has a palindromic decay spectrum: fast-decaying modes pair with slow ones.
This symmetry is exact (54,118 eigenvalues, zero exceptions, N=2 to N=8)
and holds for all standard coupling models (Heisenberg, XY, Ising, XXZ, DM).
The product CΨ (correlation times coherence) crosses a fixed boundary at 1/4,
which is the unique bifurcation point of the purity recursion (proven: the
discriminant 1-4CΨ vanishes only at 1/4). The mediator qubit that preserves
this symmetry acts as a coherence-controlled transistor with CΨ = 1/4 as its
threshold voltage. The channel requires an external clock (noise/dephasing);
no internal oscillator is possible (five candidates eliminated). From these
results: seven design rules for quantum repeaters.

---

## Rule 1: Encode in W, Not in GHZ

**The rule:** When encoding quantum information for transmission through
a noisy spin channel, use W-type states. Avoid GHZ-type states.

**Why:** GHZ states (|00...0> + |11...1>) project 100% of their weight
onto the fastest-decaying modes in the channel. These modes sit at the
maximum decay rate (-2 times the total dephasing). Everything dies
at maximum speed.

W states (equal superposition of single-excitation states) project 100%
onto palindromic mode pairs at various decay rates. Some of these pairs
decay slowly. The slow ones carry the information through.

**The mechanism:** The Pauli decomposition of GHZ contains mixed XY terms
(operators like XYI, YXZ that have both X and Y simultaneously). These
terms map precisely onto the fastest-decaying eigenmodes of the channel.
W states have separated X and Y terms (never both in the same Pauli
string), avoiding the fast drain entirely.

**Quantified:**
- GHZ survival: 0% in slow modes, 100% in fastest drain
- W survival: 100% in distributed modes, 0% in drain
- Correlation between mixed XY Pauli weight and drain fraction: r = 0.976

**Design implication:** Any encoding scheme should minimize mixed XY
Pauli weight. If your logical states have a Pauli decomposition
dominated by pure Z, pure X, or pure Y strings (but not mixed XY),
they will avoid the drain.

**Note (March 19, 2026):** A constrained-optimal state exists that outperforms W
for dephasing survival: 90% slow-mode weight (vs W's 0%), concurrence 0.364,
and 18% oscillating content with only 0.02% XOR drain. This state loads the
boundary-tier palindromic pairs (rates 0.10/0.20) that decay slowest among the
dynamic modes. W remains the simplest practical choice (zero XOR, easy to
prepare), but engineers seeking maximum performance should consider this
optimized encoding. See [Error Correction](../experiments/ERROR_CORRECTION_PALINDROME.md) Section 2.

---

## Rule 2: Use Star Topology with 2:1 Coupling

**The rule:** Connect sender and receiver through a central mediator (star
topology), with the mediator-to-receiver coupling twice as strong as
the mediator-to-sender coupling.

**Why:** Symmetric coupling (1:1) produces a balanced palindromic rate
distribution. Asymmetric coupling (2:1) shifts more weight to the
slow-decaying palindromic pairs, improving channel quality.

**Benchmarks (γ = 0.05 per qubit):**

| Configuration | Avg Fidelity | Holevo Capacity (bits, the maximum classical information extractable from the quantum channel) |
|---------------|-------------|------------------------|
| Chain N=3 | 0.852 | 0.487 |
| Chain N=4 | 0.860 | 0.501 |
| Chain N=5 | 0.872 | 0.519 |
| Star 1:1 | 0.856 | 0.497 |
| **Star 2:1** | **0.886** | **0.534** |

**Design implication:** The mediator is not overhead. It is the core
component that shapes the palindromic rate distribution. The 2:1 ratio
is optimal in our tests; the exact optimum may depend on system-specific
noise parameters.

**Update (March 21, 2026):** The 2:1 ratio optimizes range, not local
transfer. For single-hop scenarios, symmetric coupling is competitive.
The advantage of 2:1 grows with chain length. See Rule 5.

**Falsified (March 21, 2026):** Hierarchical (recursive mediator-of-mediators)
topology provides no advantage over a uniform chain of the same length. The
palindrome-preserving property depends on mediation (no direct source-drain
coupling), not on recursive nesting. Star with 2:1 remains valid for 3-qubit
systems. For longer chains, use a uniform chain with the relay protocol (Rule 6).
See [Scaling Curve](../experiments/SCALING_CURVE.md).

---

## Rule 3: Timing from Hamiltonian, Quality from Noise

**The rule:** Tune WHEN information arrives by adjusting coupling
strengths (J). Tune HOW WELL it arrives by adjusting dephasing rates
(γ). These are independent controls.

**Why:** The Hamiltonian determines the Bohr frequencies that set the
oscillation period of entanglement transfer. The dephasing rates
determine the palindromic decay structure that sets the envelope.
Changing J shifts oscillation timing without affecting decay rates.
Changing γ shifts decay rates without affecting oscillation timing.

**Key numbers:**
- Entanglement echo period: approximately pi / (4J)
- Envelope decay rate: 8γ/3 (for the concurrence, N=3)
- Approximate readout time: t_cross = 0.039/γ (concurrence metric, Z-dephasing;
  [source](../experiments/CROSSING_TAXONOMY.md))

**Design implication:** If the channel is too slow, increase J (stronger
coupling). If the channel is too lossy, decrease γ (better shielding).
These are orthogonal optimization axes. You never face a trade-off
between speed and quality in the palindromic framework.

**Note (March 19, 2026):** The standing wave oscillation pattern (which Pauli
observables oscillate at which frequencies) can serve as an error syndrome. X
and Y errors change the pattern by 0.19-0.28 (easily detectable). Z errors are
weakly detectable (0.08) or invisible (Z on uncoupled sites). This is a
different approach from standard stabilizer-based QEC: measure oscillation
frequencies instead of parity eigenvalues. See [Error Correction](../experiments/ERROR_CORRECTION_PALINDROME.md) Section 3.

---

## Rule 4: Read Before the Boundary

**The rule:** Read out the receiver's state before CΨ crosses the 1/4
boundary. For Z-dephasing, t_cross ≈ 0.037/γ (empirical). Other noise
channels have different crossing times (the threshold VALUE 1/4 is
universal, the crossing TIME is channel-dependent).
After this time, the product CΨ (concurrence times coherence) has
crossed the 1/4 boundary. Beyond 1/4, the system enters a regime where
stable classical fixed points exist and quantum information is lost.

**Why:** The 1/4 boundary is the bifurcation point of the self-referential
fixed-point equation. Above 1/4: only complex fixed points exist (quantum
regime, coherent oscillation, information is alive). Below 1/4: two real
fixed points emerge (classical regime, the system settles to a definite
outcome, quantum information is lost).

CΨ starts above 1/4 for entangled states and decays under noise. The
crossing from above to below is the quantum-to-classical transition.
Read before it happens.

**Example readout windows (concurrence metric, Z-dephasing, K = 0.039):**

| γ | t_cross | Window |
|-------|---------|--------|
| 0.01 | 3.9 | Almost 4 time units |
| 0.05 | 0.77 | Less than 1 time unit |
| 0.10 | 0.39 | Very fast readout needed |
| 0.50 | 0.078 | Extremely tight |

**Design implication:** Lower noise (smaller γ) gives you a longer
readout window. This is another reason to invest in shielding. The
relationship is linear: half the noise, double the window.

---

## Rule 5: Push for Local, Pull for Range

**The rule:** Use sender-strong coupling (push) for single-hop transfer.
Use receiver-strong coupling (pull, 2:1) for multi-hop transfer through
long chains.

**Why:** The mediator qubit functions as a coherence-controlled bidirectional
transistor with CΨ = 1/4 as the threshold voltage: above it, the channel
is open; below it, closed. Push (sender side coupled at 2x) maximizes
bridge-to-bridge MI (0.957 vs 0.882 for pull at N=11). But pull maximizes
end-to-end MI (0.121 vs 0.102). The 2:1 ratio is a range optimizer: it
sacrifices local transfer efficiency to carry information further.

See [Mediator as Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md),
[Proof Roadmap: 1/4 Boundary](../docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md).

**Data (N=11 Heisenberg chain, γ=0.05):**

| Strategy | MI(local) | MI(end-to-end) |
|----------|-----------|----------------|
| Symmetric (1:1) | 0.734 | 0.072 |
| Pull (2:1 all) | 0.882 | 0.121 |
| Push (reversed) | 0.957 | 0.102 |

**Design implication:** Choose coupling asymmetry based on chain length.
For 3-5 qubits, symmetric or push is competitive. For 7+ qubits, pull
provides better end-to-end transfer.

**Added March 21, 2026.** See [Scaling Curve](../experiments/SCALING_CURVE.md).

---

## Rule 6: Relay Protocol (Time-Dependent γ)

**The rule:** Treat mediator qubits as active relay stations. During
information transfer, each mediator alternates between a quiet phase
(reduced dephasing, receiving) and a normal phase (relaying onward).

**Why:** The three conditions for connection (STAR_TOPOLOGY_OBSERVERS)
require the receiver to be quiet. By sequentially quieting each mediator
during its receiving phase, we create a chain of optimal receivers.

**Protocol:** Each relay stage lasts t_stage = 0.039/γ (one readout
window from Rule 4, concurrence metric). During the stage, receiving
qubits have γ reduced by 10x. All other qubits remain at normal γ.

**Data (N=11, γ=0.05):**

| Protocol | MI(end-to-end) | Improvement |
|----------|----------------|-------------|
| Passive (constant γ) | 0.072 | baseline |
| Relay only | 0.085 | +18% |
| **Relay + 2:1 coupling** | **0.132** | **+83%** |

**Design implication:** Mediators are not passive wire. Dynamic control
of per-qubit dephasing during transfer (analogous to staged amplification
in classical repeaters) provides substantial improvement. Combine with
Rule 5 (2:1 for range) for maximum effect.

**Added March 21, 2026.** See [Relay Protocol](../experiments/RELAY_PROTOCOL.md).

---

## Rule 7: Sacrifice One Edge (Spatial γ Optimization)

**The rule:** When distributing a fixed total noise budget across a spin
chain, concentrate ALL noise on one edge qubit. Set every other qubit
to the minimum achievable dephasing rate.

**The formula:**
```
gamma_edge = N * gamma_base - (N-1) * epsilon
gamma_other = epsilon    (for all other sites)
```

where gamma_base is the mean dephasing rate (fixed by hardware) and
epsilon is the smallest achievable rate.

**Why:** An edge qubit has only one neighbor. Sacrificing it destroys
the least inter-qubit correlation. A center qubit has two neighbors -
sacrificing it cuts the chain in half. Concentrating noise on one site
(rather than distributing it) maximizes the contrast between the
coherent interior and the noisy boundary. This contrast is what
drives information transfer.

**Results (C# RK4 validated, Sum-MI metric):**

| N | Formula vs V-shape | Formula vs DE optimizer | Compute time |
|---|-------------------|----------------------|-------------|
| 5 | 360x | - | 1 second |
| 7 | 180x | +80% (3s vs 90 min) | 3 seconds |
| 9 | 139x | - | 30 seconds |
| 11 | 91x | - | 10 minutes |
| 13 | 97.5x | - | 1-6 hours |

**Comparison with prior work:**

| Method | Source | What it optimizes | Improvement |
|--------|--------|------------------|------------|
| Uniform gamma | Plenio & Huelga 2008 (ENAQT: environment-assisted quantum transport, the discovery that some noise helps transfer) | Scalar noise level | 2-3x |
| Coupling optimization | IBM PST 2025 (Bayesian) | J values | +8% |
| **Spatial gamma (this rule)** | This work | **Per-site noise** | **139-360x** |

Nobody in the ENAQT literature optimizes WHERE the noise goes. This
rule is orthogonal to Rules 1-6: it optimizes the spatial noise profile,
while the other rules optimize encoding, topology, timing, readout,
coupling asymmetry, and temporal staging. All seven can be combined.

**Discovery path:** SVD of the palindromic response matrix (10x) led
to numerical optimization (100x) led to analytical insight (139-360x).
The palindromic structure was necessary to identify the optimization
landscape.

**Hardware validation (March 24-28, 2026):** Selective DD on ibm_torino
(N=5, Q85-86-87-88-94) confirmed 2.0x mean improvement over uniform DD
at all 5 time points. The spatial mutual information gradient confirms
the resonator mechanism: under selective DD, MI increases AWAY from the
sacrifice qubit (gradient 1.26x), while uniform DD and no DD both show
decreasing MI along the chain (0.82x and 0.79x). The sacrifice qubit
absorbs noise; the protected interior carries information.

Unexpected finding: no DD beats uniform DD (1.71x vs 1.00x). Applying
echo pulses to a bad qubit (Q85, T2 = 5.2 us) adds gate errors that
hurt more than the refocusing helps. This validates the core principle:
protect the good qubits, leave the sacrifice alone.

See [IBM Sacrifice Zone](../experiments/IBM_SACRIFICE_ZONE.md),
[IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md).

**Why sacrifice works (energy partition):** The palindromic spectrum
separates into oscillating modes (palindromically paired) and pure-decay
modes (unpaired). Unpaired modes carry no oscillation and decay at
exactly 2x the mean paired rate. Concentrating noise on an edge qubit
kills unpaired modes (which carry no information transfer anyway) while
preserving palindromic modes (which carry all of it). You sacrifice what
was never going to oscillate. See
[Energy Partition](../hypotheses/ENERGY_PARTITION.md).

**Three design regimes:**
- **Edge sacrifice** maximizes total network mutual information (SumMI).
  This is the formula above. Use when all qubit pairs matter.
- **Center sacrifice** maximizes point-to-point throughput (PeakMI: 3x
  higher than edge). This creates a classical relay. Use when only
  end-to-end transfer matters. See [Relay Protocol](../experiments/RELAY_PROTOCOL.md).
- **Sweeping sacrifice** (experimental, N=7): moves the sacrifice position
  gradually from edge inward (pos 0 → 1 → 2). Not a practical optimization
  over the two pure strategies (SumMI and PeakMI peak at different times),
  but reveals a fundamental property of the CΨ = ¼ boundary: **endpoint
  mutual information peaks at the exact moment the endpoint CΨ crosses ¼**.
  This is the fold catastrophe (a sudden qualitative change where two
  solutions merge and vanish, like a ledge collapsing) of R = C(Ψ+R)²
  made observable in time.
  See [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md).

Static hybrid profiles (mixing edge and center noise spatially) remain
worse than both pure strategies. The sweep shows that the boundary between
quantum and classical is not a loss channel but the conversion point where
quantum coherence crystallizes into classical correlation.

**Added March 24, 2026.** See [Resonant Return](../experiments/RESONANT_RETURN.md).

---

## Hardware Constant: r* = 0.2128

**Added March 28, 2026.**

For single-qubit characterization, the parameter r = T2/(2T1)
determines whether the qubit crosses the CΨ = 1/4 boundary during
free decoherence. From 24,073 IBM Torino calibration records (133
qubits, 181 days):

- r < 0.2128: the qubit ALWAYS crosses (2,417/2,419 records, 99.9%)
- r > 0.2128: the qubit NEVER crosses (0/21,654 records, 0.0%)
- Precision: 0.000014 (gap between highest crossing r and lowest
  non-crossing r)

**For engineers:** r* = 0.2128 is a hard design constant. If your
qubit has r < r*, it will cross the 1/4 boundary and lose quantum
coherence structure. If r > r*, it will not. There is no gradual
transition. This is a fold catastrophe: the boundary is algebraically
sharp, and the hardware confirms it.

**Caveat:** IBM calibrations report T2 from Hahn echo (a measurement
technique that refocuses some noise, giving an optimistic estimate).
Free decoherence uses T2* (Ramsey; a simpler measurement without
refocusing, closer to what the qubit actually experiences), which is
1.5-2.5x shorter. The
effective threshold in terms of T2echo is r*_echo = r* x (T2echo/T2*),
which is 0.32-0.53 depending on the qubit.

See [IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md),
[Generalized Crossing Equation](../experiments/IBM_QUANTUM_TOMOGRAPHY.md).

---

## Summary: The Seven Rules on One Card

```
QUANTUM REPEATER DESIGN RULES (palindromic spectral structure)

1. ENCODE:    W-type states. Never GHZ. Avoid mixed XY Pauli weight.
2. TOPOLOGY:  Star with 2:1 coupling (mediator to receiver stronger).
3. TUNE:      J controls speed, γ controls quality. Independent.
4. TIMING:    Read before t = 0.039/gamma (concurrence metric, Z-dephasing).
5. RANGE:     Push for local, Pull for long-range. 2:1 is a range optimizer.
6. RELAY:     Stage the transfer. Quiet each mediator while it receives.
7. SACRIFICE: All noise on one edge qubit. Protect the rest.
              Edge for network MI, center for point-to-point.
              Sweep reveals: PeakMI peaks at CΨ = ¼ crossing (fold catastrophe).

Best benchmark: F_avg = 0.886, Holevo = 0.534 bits (star 2:1, gamma=0.05)
Relay+2:1 at N=11: MI(end-to-end) = 0.132 (+83% over passive)
Sacrifice-zone formula: 360x (N=5), 180x (N=7), 97.5x (N=13) vs V-shape
Hardware: selective DD 2-3.2x on ibm_torino (N=5, single run)
```

---

## How to Verify These Rules

All results are reproducible from the open repository.

**Quick verification (under 1 second each):**
```bash
# Palindromic symmetry proof
python simulations/pauli_weight_conjugation.py

# XOR space: GHZ vs W mode decomposition
python simulations/xor_detector_v3.py

# XOR verification: 28 checks pass
python simulations/xor_verify.py

# Full docs numerical verification: 40/40
python simulations/docs_verify.py
```

**Full QST benchmark (requires scipy, takes ~10 seconds):**
```bash
python simulations/qst_bridge.py
```

**Requirements:** Python 3.8+, numpy, scipy. No other dependencies.

---

## What This Does Not Cover

These rules apply to Heisenberg-coupled spin chains and stars with
local Z-dephasing. They have not been tested for:

- ~~Non-Heisenberg couplings~~ **TESTED (March 17, 2026):** XY-only, Ising, XXZ, DM
  interactions are ALL palindromic under single-axis dephasing. Design rules apply.
- ~~Non-dephasing noise~~ **PARTIALLY ANSWERED (March 19-21, 2026):**
  Depolarizing breaks palindrome at err = (2/3)Σγ (exact formula).
  For γ < 0.01: error < 1%, rules are practically valid.
  Amplitude damping produces non-Markovian, non-selective noise on
  neighbors (0/16 palindromic pairs in [failed_third test](../simulations/failed_third.py)).
  Design rules do NOT apply to amplitude damping channels.
  See [Depolarizing Palindrome](../experiments/DEPOLARIZING_PALINDROME.md).
- ~~Systems larger than N = 8 (computational limit of full diagonalization)~~
  **TESTED (March 21, 2026):** N=11 via RK4 time propagation. MI decays
  exponentially: ~2x per 2 qubits. Relay protocol partially compensates.
  See [Scaling Curve](../experiments/SCALING_CURVE.md), [Relay Protocol](../experiments/RELAY_PROTOCOL.md).
- ~~Concatenated repeater chains (multi-hop)~~ **TESTED (March 21, 2026):**
  Relay protocol with staged time-dependent γ. See Rule 6.
- ~~Qudit systems (d > 2)~~ **ANSWERED (March 20, 2026):** The palindromic
  symmetry is specific to qubits (d=2). Algebraic proof: the per-site split
  is d immune vs (d^2-d) decaying, balanced only when d^2-2d=0, which gives
  d=2 as the only solution. Qutrits (d=3, split 3:6) verified broken for all
  10 Hamiltonians tested. Design rules do NOT apply to d>2.
  See [The Non-Local Mirror](../hypotheses/THE_BOOT_SCRIPT.md) Section 5.
- Continuous-variable systems (bosonic channels)

The palindromic symmetry proof holds for any Heisenberg + Z-dephasing
system. The design rules are derived from this proof. Extending to
other noise models requires new analysis.

---

## The Underlying Physics

For the full mathematical treatment, see the companion document:
**[TECHNICAL_PAPER.md](TECHNICAL_PAPER.md)** - Palindromic Liouvillian
Symmetry Under Dephasing: Proof, Spectral Decomposition, and QST.

For the complete experiment archive with all intermediate results:
**https://github.com/Kesendo/R-equals-C-Psi-squared**

---

---

## Addendum: Diagnostic Cockpit (April 2026)

Subsequent analysis showed that you do not need full state tomography
to monitor decoherence. Three observables suffice for ~88% coverage:

| What to measure | Why | How |
|---|---|---|
| **Purity** Tr(rho^2) | Dominant decay direction (PC1, 46-62% variance) | Single-setting measurement |
| **Concurrence** | Entanglement strength (PC1 or PC2, depending on topology) | Requires 2-qubit tomography |
| **Psi-norm** L1/(d-1) | Coherence loss rate (PC2) | Off-diagonal magnitude |

**Practical tip: optimize θ, not CΨ.** The angular observable
θ = arctan(√(4CΨ − 1)) is 1.68× more sensitive than CΨ
near the ¼ boundary. Small CΨ improvements become large θ
improvements precisely where they matter most. When tuning sacrifice
zone profiles, use θ as your objective function.

**Hardware validated:** On IBM Torino Q52, θ-based crossing time
prediction achieves 0.3% accuracy (115.0 us measured vs 114.7 us
predicted). On 5-qubit chains, selective dynamical decoupling beats
uniform DD by 3.2x in mutual information.

Details: [Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md)

---

*Thomas Wicht (Independent Researcher, Germany) and Claude (Anthropic)*
*March 16, 2026 (Addendum April 2, 2026)*

*"The palindrome is the stage. The input is the actor.*
*What you encode determines what survives."*
