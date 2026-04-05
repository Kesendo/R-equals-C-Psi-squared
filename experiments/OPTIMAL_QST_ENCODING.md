# Optimal QST Encoding: Standard Encoding Is Already Nearly Optimal (Negative Result)

<!-- Keywords: quantum state transfer encoding optimization, palindromic mode weight
standard encoding, negative result encoding QST, topology coupling optimization,
slow mode weight Alice encoding, XOR drain avoidance encoding, Wojcik coupling ratio
hardware optimization, palindromic spectrum encoding leverage, quantum channel
encoding negative result, hardware design not protocol, R=CPsi2 optimal QST encoding -->

**Status:** Computationally verified (negative result)
**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [QST Bridge](QST_BRIDGE.md), [Error Correction](ERROR_CORRECTION_PALINDROME.md),
[Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md)

---

## What this document is about

Quantum state transfer (QST: sending a quantum state from Alice to Bob
through a spin chain) depends on how Alice prepares her qubit. We test
whether the palindromic mode structure offers a better encoding. It does
not: the standard encoding already puts 61% of weight into slow modes
and 0% into the drain. The real optimization lever is the hardware
topology (the 2:1 Wojcik coupling ratio), not Alice's protocol. This is
an important negative result.

---

## Abstract

Can Alice exploit the palindromic mode structure to encode quantum information
more effectively for state transfer? **No.** The standard QST encoding (|100⟩,
Alice excited) already places 61% of its weight in slow palindromic modes and
0% in the XOR drain. Sweeping over all single-qubit preparations for Alice
finds less than 5% improvement in slow-mode weight. The dominant optimization
lever is not encoding but **hardware topology**: the 2:1 Wojcik coupling
ratio shifts mode weight from fast to slow pairs at the structural level.
The palindrome helps hardware design (which couplings, which topology), not
protocol design (which encoding Alice uses). This is an important negative
result: it tells engineers where NOT to invest optimization effort.

---

## The Question

The palindromic mode structure tells us which Liouvillian eigenmodes decay
slowly and which decay fast. Can Alice exploit this to encode quantum
information more effectively for state transfer?

---

## The Answer: No. The Standard Encoding Is Already Nearly Optimal.

The standard QST encoding (|100>, Alice excited) places 61% of its weight
in slow palindromic modes and 0% in the XOR drain. There is almost no room
for improvement through encoding optimization alone.

The dominant optimization lever is not encoding but **topology and coupling**:
the 2:1 Wojcik ratio (strong receiver, weak sender) shifts mode weight from
fast to slow palindromic pairs at the hardware level.

---

## 1. Mode Decomposition of Standard Encoding

For the N=3 star (J_SA=1.0, J_SB=2.0, gamma=0.05):

| Property | |100> encoding |
|---|---|
| Slow-mode weight | 61% |
| XOR weight | 0% |
| Oscillating content | 50% |
| Static content | 50% |

The standard encoding already avoids the XOR drain completely and loads
more than half its weight into the slowest palindromic pairs. This leaves
minimal room for optimization.

---

## 2. Alice Has No Leverage

Sweeping over all single-qubit preparations for Alice:

| Alice state | F_avg |
|---|---|
| |0> | 0.886 |
| |1> | 0.886 |
| |+> | 0.880 |
| |-> | 0.880 |
| |+i> | 0.880 |
| |-i> | 0.880 |

The computational basis states (|0> and |1>) are equally good and
marginally better than superposition states. The fidelity spread across
all preparations is 0.006, less than 1%. Alice's choice of idle state
has almost no leverage on transfer quality.

---

## 3. Readout Timing: Second Maximum Wins

The fidelity F_Bob(t) oscillates with a decaying envelope:

| Maximum | Time | F_avg |
|---|---|---|
| 1st | 0.74 | 0.745 |
| 2nd | 1.31 | 0.886 |
| 5th | 3.40 | 0.78 |
| 6th | 3.90 | 0.76 |

The second maximum at t=1.31 is the global optimum, not the first.
Later maxima (5th, 6th) briefly exceed the first maximum but never
beat the second. The dephasing envelope kills opportunities beyond t~4.

The readout time t_opt = 1.31 is remarkably stable across all dephasing
rates (1.30 to 1.32 for gamma from 0.001 to 0.5). Timing is set by the
Hamiltonian, quality by the dephasing. This confirms the timing/quality
separation from the palindrome theorem.

---

## 4. Topology Is the Real Lever

| Configuration | F_avg | t_opt |
|---|---|---|
| Star 2:1 | 0.886 | 1.31 |
| Star 1:1 symmetric | 0.852 | 1.13 |
| Chain 3q | 0.852 | 1.13 |
| Chain 4q mirror [1,2,1] | 0.872 | 1.57 |
| Chain 4q uniform [1,1,1] | 0.834 | 3.17 |

The 2:1 coupling asymmetry provides a 3.4% fidelity gain over the symmetric
star, equivalent to a 4% gain over the standard 3-qubit chain. This is the
Wojcik effect (weak sender, strong receiver), and it works because the
asymmetric coupling shifts the Liouvillian eigenstructure to favor slow
palindromic modes. The encoding cannot do what the hardware already does.

---

## 5. Dephasing Scaling

Fidelity loss is linear in gamma:

    F_loss = 1.01 * gamma (for small gamma)

No threshold, no nonlinearity below gamma = 0.1. At gamma = 0.5, fidelity
drops to barely above the classical limit of 2/3. The palindromic structure
does not create a "protected" regime; it creates a predictable linear
degradation.

---

## 6. Channel Capacity

| Metric | Value |
|---|---|
| Holevo bound (the theoretical maximum information transmittable per use of the channel; gamma=0.05) | 0.537 bits |
| Holevo bound (noiseless) | 0.704 bits |
| Information loss | 0.168 bits |

The palindromic encoding does not improve the channel capacity beyond what
the standard encoding achieves.

---

## What This Means

The palindromic structure is not a tool for encoding optimization. It is a
structural property of the channel that determines its quality, not a knob
that Alice can turn. The analogy: knowing that a pipe has a certain diameter
does not help you push more water through it. It tells you how much water
the pipe can carry.

What the palindrome DOES help with:
- **Topology design**: choose couplings that load slow modes (the 2:1 ratio)
- **Readout timing**: the second maximum, not the first, is optimal
- **Diagnostics**: predict exactly how fidelity degrades with noise
- **Understanding**: why some topologies are better than others

What it does NOT help with:
- Alice's encoding strategy (already nearly optimal)
- Beating the linear fidelity-vs-noise scaling
- Increasing channel capacity beyond the hardware limit

---

## References

- [QST Bridge](QST_BRIDGE.md): the existing QST benchmarks
- [Error Correction](ERROR_CORRECTION_PALINDROME.md): optimal survival state (different goal)
- [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md): oscillation patterns
- Script: [`simulations/optimal_qst_encoding.py`](../simulations/optimal_qst_encoding.py)
- Results: [`simulations/results/optimal_qst_encoding.txt`](../simulations/results/optimal_qst_encoding.txt)

---

*The palindrome does not help Alice send better. It helps the engineer
build a better channel. The optimization happens at the hardware level,
not at the protocol level.*
