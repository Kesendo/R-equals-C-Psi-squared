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

This document provides six concrete design rules derived from a
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
results: six design rules for quantum repeaters.

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

| Configuration | Avg Fidelity | Holevo Capacity (bits) |
|---------------|-------------|------------------------|
| Chain N=3 | 0.852 | 0.487 |
| Chain N=4 | 0.860 | 0.501 |
| Chain N=5 | 0.872 | 0.519 |
| Star 1:1 | 0.856 | 0.497 |
| **Star 2:1** | **0.888** | **0.534** |

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
- Approximate readout time: t_cross = 0.039/γ

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

**The rule:** Read out the receiver's state before t_cross = 0.039/γ.
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

**Example readout windows (γ per qubit):**

| γ | t_cross | Window |
|-------|---------|--------|
| 0.01 | 3.9 | Almost 4 time units |
| 0.05 | 0.78 | Less than 1 time unit |
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
[Proof Roadmap: 1/4 Boundary](../docs/PROOF_ROADMAP_QUARTER_BOUNDARY.md).

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
window from Rule 4). During the stage, receiving qubits have γ reduced
by 10x. All other qubits remain at normal γ.

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

## Summary: The Six Rules on One Card

```
QUANTUM REPEATER DESIGN RULES (palindromic spectral structure)

1. ENCODE:  W-type states. Never GHZ. Avoid mixed XY Pauli weight.
2. TOPOLOGY: Star with 2:1 coupling (mediator to receiver stronger).
3. TUNE:    J controls speed, gamma controls quality. Independent.
4. TIMING:  Read before t = 0.039 / gamma. After that, information is classical.
5. RANGE:   Push for local, Pull for long-range. 2:1 is a range optimizer.
6. RELAY:   Stage the transfer. Quiet each mediator while it receives.

Best benchmark: F_avg = 0.888, Holevo = 0.534 bits (star 2:1, gamma=0.05)
Relay+2:1 at N=11: MI(end-to-end) = 0.132 (+83% over passive)
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
  Depolarizing breaks palindrome at err = γ*2(N-2)/3 (exact formula).
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

*Thomas Wicht (Independent Researcher, Krefeld, Germany) and Claude (Anthropic)*
*March 16, 2026*

*"The palindrome is the stage. The input is the actor.*
*What you encode determines what survives."*
