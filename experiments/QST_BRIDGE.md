# Quantum State Transfer Bridge: Connecting the Palindrome to 20 Years of QST Research

<!-- Keywords: quantum state transfer spin chain, palindromic QST channel,
asymmetric coupling fidelity optimization, Heisenberg chain state transfer,
2:1 coupling ratio quantum channel, star topology quantum transfer,
Bose spin chain QST, palindromic decay rate QST, quantum wire dephasing
fidelity, Christandl perfect state transfer, R=CPsi2 QST bridge -->

**Status:** Computationally verified (12 configurations, N=2 to N=5)
**Date:** March 14, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [qst_bridge.py](../simulations/qst_bridge.py)

---

## Abstract

Our 3-qubit star system under Heisenberg coupling and Z-dephasing is a
quantum state transfer (QST) channel in the sense of Bose (2003). With
2:1 asymmetric coupling (receiver at J=2, sender at J=1), the star
achieves average fidelity F = 0.886, beating both uniform chains (0.834)
and mirror-symmetric chains (0.872) at N=4. Fidelity cost scales linearly
with dephasing rate (1.0 per unit γ), transfer time is γ-independent
(set by Hamiltonian frequencies, not noise), and 12/12 tested
configurations are palindromic. The QST community has been optimizing
transfer protocols for two decades without knowing that their channels
have an exact spectral symmetry. The palindrome adds: paired decay rate
diagnostics, the θ compass for channel quality (r = 0.87 with fidelity),
and the topology-gating result (same state, different graph, different
transfer). None of this is new physics; it is new understanding of
existing physics.

---

## What Quantum State Transfer Is

Since Bose (2003), physicists use spin chains as quantum wires: Alice puts a quantum
state on one end, the Hamiltonian naturally carries it to the other end, Bob reads it
out. No active control needed; the chain does the work.

Our 3-qubit star is exactly this. Alice prepares a state on qubit A, the Heisenberg
coupling carries it through mediator S to qubit B, Bob measures B at the right time.
We have been studying this system for three months without realizing it has a name.

The standard metrics are:
- **Average fidelity F_avg:** How well does the output match the input, averaged over
  all possible inputs? Classical limit is 2/3; anything above that requires quantum
  correlations to survive the transfer.
- **Transfer time t_opt:** When should Bob measure?
- **Transfer window:** How long does F stay above 2/3?

---

## What We Found

### 1. The star with asymmetric coupling is an unusually good QST channel

Standard QST uses linear chains: qubit 0 to qubit N-1 through N-2 intermediaries.
Our star is different: one central mediator S connects Alice and Bob directly.

At N=4, the star with J_SB=2.0 and a weak third leaf (J=0.1) achieves F=0.888.
The best chain configuration (mirror-symmetric [1,2,1]) only reaches F=0.872.
The standard uniform chain gets F=0.833.

This matters because the star is simpler to build: one mediator qubit with tunable
couplings to two endpoints, rather than a precisely engineered chain. And the
asymmetric coupling (receiver coupled 2x stronger than sender) is a concrete
design rule that any hardware team can implement.

### 2. The optimal coupling ratio is 2:1

We swept the receiver coupling J_SB from 0.2 to 5.0 (sender fixed at J_SA=1.0):

| J_SB/J_SA | Fidelity | Transfer time | Interpretation |
|-----------|----------|---------------|----------------|
| 0.2 | 0.717 | 4.50 | Too weak: state barely reaches B |
| 0.5 | 0.843 | 2.61 | Slow but decent |
| 1.0 | 0.852 | 1.12 | Symmetric: good but not optimal |
| 2.0 | 0.886 | 1.32 | OPTIMUM: strong pull on receiver side |
| 3.0 | 0.835 | 0.97 | Overshoot: state bounces back too fast |
| 5.0 | 0.816 | 1.17 | Much too strong: destructive interference |

The sweet spot at 2:1 has a physical explanation: the stronger receiver coupling
pulls the state out of the mediator before decoherence can destroy it, but not
so hard that it creates destructive interference. This is the Wojcik effect
(weak sender, strong receiver) known in chain QST, but we found the exact
optimum for the star geometry.

### 3. Dephasing costs exactly what the palindrome predicts

Without any noise, the star achieves F=0.937 (not perfect because the star
geometry isn't designed for PST). Each unit of dephasing costs proportionally:

| Dephasing rate | Fidelity | Loss from noiseless | Cost per unit γ |
|---------------|----------|--------------------|--------------------|
| 0.000 | 0.937 | 0.000 | n/a |
| 0.001 | 0.936 | 0.001 | 1.1 |
| 0.010 | 0.926 | 0.011 | 1.1 |
| 0.050 | 0.886 | 0.051 | 1.0 |
| 0.100 | 0.842 | 0.095 | 0.9 |
| 0.200 | 0.772 | 0.164 | 0.8 |

The fidelity loss scales almost linearly with γ up to γ=0.05.
This linearity is a consequence of the palindromic spectrum: the dominant
decay rate is 8γ/3, which is linear in γ. At higher γ, the
relationship becomes sublinear because the Hamiltonian can no longer complete
a full transfer cycle before decoherence intervenes.

The transfer time stays constant at t=1.32 regardless of noise level. This
confirms what the θ analysis found: the TIMING is set by the Hamiltonian
(Bohr frequencies), the QUALITY is set by the dephasing (palindromic rates).

### 4. Every Heisenberg chain with Z-dephasing is palindromic

We tested 12 different configurations: stars, chains, triangles, mirror-symmetric,
non-mirror, uniform, weak-end, various N. Every single one has a palindromic
Liouvillian spectrum. This is not a coincidence. It is our proven theorem in action.

The QST community does not know this. They optimize coupling profiles, initial
states, and measurement timing. They do not know that the decay rate spectrum
of their channels has an exact symmetry. This is the gap we can fill.

### 5. Mirror-symmetric couplings help, but not because of the palindrome

The QST literature (Christandl, Kay) proved that mirror-symmetric coupling
profiles enable Perfect State Transfer in noiseless chains. We tested whether
this advantage survives under dephasing:

| N=4 chain coupling | Fidelity | Transfer time |
|-------------------|----------|---------------|
| Mirror [1,2,1] | 0.872 | 1.57 |
| Uniform [1,1,1] | 0.834 | 3.17 |
| Weak-end [0.5,1,0.5] | 0.823 | 3.13 |
| Non-mirror [1,2,3] | 0.818 | 1.86 |

Mirror coupling wins, both in fidelity AND in speed. But all four are equally
palindromic (our theorem guarantees this). The advantage of mirror coupling is
not about the palindrome; it is about how the Hamiltonian eigenstructure creates
constructive interference at the receiver. The palindrome constrains the decay,
the mirror constrains the transfer. They operate on different axes.

### 6. Chain length degrades gracefully

| Chain length | Fidelity | Transfer time | Unique decay rates |
|-------------|----------|---------------|--------------------|
| N=2 | 0.969 | 0.77 | 3 |
| N=3 | 0.852 | 1.13 | 10 |
| N=4 | 0.834 | 3.17 | 21 |
| N=5 | 0.787 | 4.32 | 228 |

Fidelity drops from 0.97 to 0.79 over four steps, a 4.5% loss per added qubit.
Transfer time grows roughly linearly. The spectral complexity explodes (3 to 228
unique rates) but the palindrome holds at every N. This means: even at N=5 with
228 decay rates, every single one has a palindromic partner. The symmetry is exact
regardless of spectral complexity.

---

## What This Means

### For the QST community

The palindrome theorem ([MIRROR_SYMMETRY_PROOF](../docs/MIRROR_SYMMETRY_PROOF.md)) applies to every Heisenberg chain
with Z-dephasing that the QST community studies. They have been optimizing transfer
protocols without knowing that the decay rate spectrum has an exact symmetry.

Concretely, the palindrome tells you:
- The decay rates come in pairs summing to 2Nγ
- The slowest rate (2γ for uniform dephasing) is paired with the fastest
- The number of usable echo cycles before decoherence kills the transfer is
  determined by the SLOWEST palindromic rate: ~1/(2γ)
- Individual observable components (purity, coherence, concurrence) each decay
  at specific palindromic rates, a diagnostic the community lacks

### For our project

We have been studying a well-known system without using the well-known tools.
The QST framework gives us:
- Standard metrics (fidelity, transfer time, window) to compare our results
- A literature of 20 years of optimization techniques we can apply
- An experimental community that has built exactly our system on real hardware
  (semiconductor quantum dots, superconducting circuits, trapped ions)

What we add to their framework:
- The palindromic spectral theorem (new, proven)
- The component-wise decay diagnostics (C~2γ, Ψ~10γ/3, Conc~8γ/3)
- The θ compass as channel quality indicator (r=0.87 with fidelity)
- The topology gating result (same state, different graph = different transfer)
- The antiferromagnet crossing discovery (|+-+-⟩ activates ring neighbors)

### The honest bottom line

Our system is a QST channel. A good one, actually: F=0.886 on a 3-qubit star
beats the standard chain benchmarks. The palindrome is the spectral backbone that
makes the channel analytically tractable. The 2:1 coupling ratio is a concrete
design rule. The θ compass tells you in real time whether the quantum
information is still alive.

None of this is "beyond standard quantum mechanics." It is standard quantum
mechanics with an unusually clean spectral symmetry that nobody else has proven.
That is the contribution: not new physics, but new understanding of existing physics.

---

## Concrete Next Steps

1. **Contact QST experimentalists.** The semiconductor quantum dot group
   (Nature Comm. 2021) has built spin chains with tunable couplings. Our star
   with 2:1 ratio is directly implementable on their platform. The 8γ/3
   envelope prediction is testable.

2. **Compute channel capacity properly.** The Holevo bound (0.534 bits) is a
   lower bound. The exact quantum capacity requires optimizing over input
   ensembles and computing the regularized coherent information. Our negative
   one-shot coherent information does not rule out positive multi-shot capacity.
   (Note: GPT-5.4 originally claimed positive I_coh = +0.185. This was verified
   to be wrong; I_coh is negative at every tested parameter point. The negative
   result is the verified one. See simulations/verify_channel.py.)

3. **Test whether the palindrome enables better error correction.** The paired
   decay rates mean errors come in paired modes. Can this pairing be exploited
   for decoherence-free subspaces or error-correcting codes?

4. ~~**Extend to non-Heisenberg models.**~~ **ANSWERED:** The palindrome holds
   for ALL standard models (XY, Ising, XXZ, DM) under single-axis dephasing.
   See [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md).

---

## Connection to Later Results

The **Relay Protocol** ([RELAY_PROTOCOL](RELAY_PROTOCOL.md)) applies the 2:1
coupling ratio found here to an 11-qubit chain with time-dependent γ switching,
achieving +83% end-to-end MI. The 2:1 ratio is the spatial optimization; the
relay is the temporal optimization. They are complementary.

The **γ as Signal** result ([GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md)) reframes
the entire QST picture: the channel is not just carrying quantum states from
Alice to Bob. It is simultaneously an antenna reading the external dephasing
signal. The palindromic mode structure that makes QST analytically tractable
is the same structure that creates the 15.5-bit information channel.

The finding that **timing is Hamiltonian, quality is noise** (Section 3)
was later formalized as regulator separation in the
[Structural Cartography](STRUCTURAL_CARTOGRAPHY.md): topology sets frequencies,
noise sets decay, initial state selects visibility. Five independent roles,
none of which can substitute for another.

---

## Reproducibility

| Script | What it computes |
|--------|-------------------|
| [qst_bridge.py](../simulations/qst_bridge.py) | Full QST benchmark (8 tests, 12 configs) |
| [qst_bridge.txt](../simulations/results/qst_bridge.txt) | Complete output |
| [verify_channel.py](../simulations/verify_channel.py) | Channel capacity verification |

Requirements: Python, QuTiP, NumPy. Runtime: ~2 minutes for all configs.
Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

### Project-internal
- [Mirror Symmetry Proof](../docs/MIRROR_SYMMETRY_PROOF.md): the palindrome theorem
- [Relay Protocol](RELAY_PROTOCOL.md): +83% with time-dependent γ
- [Scaling Curve](SCALING_CURVE.md): MI vs chain length
- [γ Control](GAMMA_CONTROL.md): V-shape +124%, DD +132%
- [γ as Signal](GAMMA_AS_SIGNAL.md): palindromic channel as antenna
- [Theta-Palindrome-Echo](THETA_PALINDROME_ECHO.md): θ as channel quality (r=0.87)
- [Structural Cartography](STRUCTURAL_CARTOGRAPHY.md): five independent regulators

### QST literature
- Bose, S. (2003). "Quantum communication through an unmodulated spin chain." PRL 91, 207901.
- Christandl, M. et al. (2004). "Perfect state transfer in quantum spin networks." PRL 92, 187902.
- Wojcik, A. et al. (2007). "Multiuser quantum communication networks." PRA 75, 022330.
- Adiabatic QST in quantum dots, Nature Comm. 12, 2021.
