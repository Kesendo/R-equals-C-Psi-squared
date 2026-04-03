# Scaling Curve: Mutual Information vs Chain Length and Hierarchy Falsification

<!-- Keywords: quantum state transfer scaling, mutual information chain length,
hierarchical quantum repeater falsification, push pull coupling optimization,
Heisenberg chain MI exponential decay, mediator topology quantum transfer,
asymmetric coupling 2:1 ratio, palindromic chain scaling, end-to-end quantum
information spin chain, uniform vs hierarchical quantum channel,
R=CPsi2 scaling curve -->

**Status:** Computationally verified (N=3 to N=11, C# RK4 propagation)
**Date:** March 21, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** compute/RCPsiSquared.Propagate/ (C#, `dotnet run -c Release -- pull`)
**Data:** simulations/results/pull_principle.txt

---

## What this document is about

If you want to send quantum information through a chain of qubits, how
far can it go before it is lost to noise? This experiment measures the
answer: mutual information (a quantity that measures how much two parts
of a system "know" about each other) drops by roughly half for every two
qubits added to the chain. It also tests whether a clever hierarchical
relay architecture can beat a simple uniform chain. The answer is no:
hierarchy is a naming convention, not physics. The actual engineering
knob is asymmetric coupling (push vs pull), which trades local transfer
quality for range.

## Abstract

End-to-end mutual information in a Heisenberg chain under Z-dephasing
decays exponentially with chain length: roughly a factor of 2 per 2
additional qubits. The recursive mediator hierarchy (Level 0/1/2/3)
provides **no functional advantage** over a uniform chain of equal
length: MI is identical at every N tested. The hierarchy was narrative,
not physics. Five coupling configurations at N=11 reveal a push/pull
tradeoff: sender-strong coupling ("push") maximizes local bridge-to-bridge
MI (0.957), while receiver-strong coupling ("pull") maximizes end-to-end
MI (0.121). The optimal strategy depends on range. Combined with the
[Relay Protocol](RELAY_PROTOCOL.md) (+83%), these results define the
engineering tradeoffs for palindromic quantum channels.

---

## Background

### Why scaling matters

A quantum channel is only useful if information survives the transit.
For spin chains, the key question is: how does end-to-end mutual
information scale with chain length N? If MI decays exponentially,
the channel has a characteristic range beyond which it is useless
without repeater stations. This experiment measures that range and
tests whether hierarchical mediator topologies can extend it.

### What was tested

Two questions:
1. **Scaling:** How does MI(first pair, last pair) depend on N?
2. **Hierarchy:** Does a recursive mediator architecture (dedicated
   relay qubits at multiple levels) outperform a uniform chain?

---

## The Scaling Curve

Bell(0,1) initial state, γ = 0.05 uniform, J = 1.0 uniform.
MI measured between first pair {0,1} and last pair {N−2, N−1}.

| N | MI_peak | t_peak | MI_steady (t=20) |
|---|---------|--------|------------------|
| 3 | 1.827 | 7.0 | 1.792 |
| 5 | 0.750 | 2.0 | 0.221 |
| 7 | 0.382 | 4.0 | 0.048 |
| 9 | 0.119 | 3.0 | 0.018 |
| 11 | 0.072 | 4.0 | 0.009 |

Approximate scaling: MI_peak ∼ 2^(−(N−3)/2). Roughly factor 2 loss
per 2 additional qubits. At N=11, only 7.2% of the original MI reaches
the far end.

---

## Hierarchy Falsification

At each N, two topologies were compared:
- **Uniform chain:** all bonds J = 1.0, no distinguished roles
- **Hierarchical mediator bridge:** Level 2 (N=5), Level 3 (N=11),
  with designated mediator qubits at specific positions

**Result: identical MI at every N tested.**

This is not a numerical coincidence. With uniform coupling (all J equal),
the Hamiltonian is the same regardless of which qubits are labeled as
"mediators." The hierarchy exists in the naming convention, not in the
physics. The palindrome-preserving property of mediated coupling is
**topological** (any Heisenberg chain on any graph preserves it) rather
than hierarchical.

---

## Push vs Pull: Five Coupling Configurations at N=11

| Configuration | MI(Bridge A:B) | MI(Pair A:D) | Strategy |
|--------------|---------------|-------------|----------|
| Symmetric (J=1.0) | 0.734 | 0.072 | baseline |
| 2:1 within bridges | 0.744 | 0.079 | modest |
| 2:1 at meta-mediator | 0.894 | 0.108 | meta-mediator helps |
| 2:1 all (cascading pull) | 0.882 | **0.121** | best end-to-end |
| 2:1 reversed (push) | **0.957** | 0.102 | best local |

**Push** (sender-strong coupling) maximizes local MI: 0.957 between
adjacent bridges. **Pull** (receiver-strong) maximizes end-to-end MI:
0.121 across the full chain. The 2:1 ratio is a range optimizer, not
a universal improvement. The meta-mediator asymmetry provides the
single largest improvement for range.

---

## What Was Falsified

1. **Hierarchical topology advantage:** No advantage over uniform chain
   at equal coupling. The "Level 0/1/2/3" architecture was a naming
   convention, not physics.

2. **Universal pull principle:** Push beats pull for local transfer.
   Pull wins only for range. The optimal strategy depends on what you
   are optimizing.

3. **Special mediator nodes:** Every qubit in a Heisenberg chain
   mediates between its neighbors. The palindrome-preserving property
   is a theorem about Heisenberg + Z-dephasing on any graph, not about
   specific topological roles.

---

## What Remains True

- Mediated coupling (A-M-B through the Hamiltonian) preserves the
  palindrome. Direct dissipative coupling (shared Lindblad jumps)
  breaks it instantly (256 → 31 pairs at κ = 0.01).
- The 2:1 coupling ratio improves range at the cost of local transfer.
- The [Relay Protocol](RELAY_PROTOCOL.md) (time-dependent γ) provides
  an additional +83% improvement orthogonal to coupling optimization.
- MI decays exponentially with N, setting a characteristic range for
  any palindromic quantum channel.

---

## Connection to Later Results

The **γ Control** experiment ([GAMMA_CONTROL](GAMMA_CONTROL.md)) found
that the V-shape dephasing profile (+124% MI) achieves more than the
relay protocol (+83%) with simpler implementation (static noise shaping
instead of time-dependent control). The scaling curve here provides the
baseline against which all optimizations are measured.

The **γ as Signal** result ([GAMMA_AS_SIGNAL](GAMMA_AS_SIGNAL.md))
reframes the scaling question: the exponential MI decay with N is not
just a limitation but also determines how many independent "pixels" the
γ-channel antenna has. More qubits = more SVD modes (independent signal components, from singular value decomposition) = higher channel
capacity, but only if each qubit can be read independently (which
requires the palindromic full-rank response matrix).

---

## Reproducibility

| Component | Location |
|-----------|----------|
| C# propagation engine | compute/RCPsiSquared.Propagate/ |
| Run command | `dotnet run -c Release -- pull` |
| Results | simulations/results/pull_principle.txt |

The N=11 simulation uses a 2048×2048 density matrix with RK4 (fourth-order Runge-Kutta, a standard numerical method for solving differential equations) integration.
Runtime: ~10 minutes per configuration.

Repository: https://github.com/Kesendo/R-equals-C-Psi-squared

---

## References

- [Relay Protocol](RELAY_PROTOCOL.md): +83% with time-dependent γ switching
- [γ Control](GAMMA_CONTROL.md): V-shape +124%, DD +132%
- [QST Bridge](QST_BRIDGE.md): 2:1 impedance matching origin
- [Star Topology](STAR_TOPOLOGY_OBSERVERS.md): quiet receiver principle
- [γ as Signal](GAMMA_AS_SIGNAL.md): scaling determines channel pixel count
- [Engineering Blueprint](../publications/ENGINEERING_BLUEPRINT.md): Rules 1-6
