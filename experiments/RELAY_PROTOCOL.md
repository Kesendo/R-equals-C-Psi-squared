# Relay Protocol: Staged Transfer with Time-Dependent Dephasing

**Tier:** 2 (computationally verified)
**Date:** March 21, 2026
**Script:** compute/RCPsiSquared.Propagate/ (C#, `dotnet run -c Release -- pull`)
**Data:** simulations/results/pull_principle.txt

---

## Summary

A quantum state transfer protocol that treats mediator qubits as active
relay stations rather than passive wire. Each mediator alternates between
a quiet phase (low dephasing, receiving information) and a normal phase
(relaying onward). Combined with 2:1 asymmetric coupling, this improves
end-to-end mutual information by 83% over passive propagation.

---

## Setup

11-qubit Heisenberg chain with local Z-dephasing (gamma = 0.05 per qubit).
Topology: linear chain, 10 bonds, all J = 1.0 (or 2:1 asymmetric).

Initial state: Bell pair on qubits 0-1, rest in |0>.

The chain is conceptually divided into relay segments:
```
(0-1) -> m1(2) -> (3-4) -> M(5) -> (6-7) -> m2(8) -> (9-10)
Pair A    relay    Pair B   relay   Pair C    relay    Pair D
```

## The Protocol

Six relay stages, each lasting t_stage = 0.039 / gamma = 0.78 time units:

| Stage | Receiving qubits | gamma_receive | Action |
|-------|-----------------|---------------|--------|
| 1 | 2 (m1) | 0.005 | m1 receives from Pair A |
| 2 | 3, 4 (Pair B) | 0.005 | Pair B receives from m1 |
| 3 | 5 (M) | 0.005 | Meta-mediator receives from Pair B |
| 4 | 6, 7 (Pair C) | 0.005 | Pair C receives from M |
| 5 | 8 (m2) | 0.005 | m2 receives from Pair C |
| 6 | 9, 10 (Pair D) | 0.005 | Pair D receives (final) |

During each stage, the receiving qubits have dephasing reduced by 10x
(gamma = 0.005 instead of 0.05). All other qubits remain at gamma = 0.05.
Total protocol time: 6 x 0.78 = 4.68 time units.

Between stages, the propagator is reconstructed with updated dephasing
rates. The density matrix evolves continuously via RK4 integration
(dt = 0.05).

## Results

| Protocol | MI(BridgeA:BridgeB) | MI(PairA:PairD) | Improvement |
|----------|---------------------|-----------------|-------------|
| Passive (constant gamma) | 0.734 | 0.072 | baseline |
| Relay only | 0.759 | 0.085 | +18% |
| **Relay + 2:1 coupling** | 0.723 | **0.132** | **+83%** |

The relay protocol alone provides a modest improvement in end-to-end MI.
Combined with 2:1 asymmetric coupling (receiver-side of each mediator
coupled at J=2.0, sender-side at J=1.0), the improvement is substantial.

## Connection to Existing Results

The staging schedule derives from Rule 4 of the Engineering Blueprint:
readout before t = 0.039 / gamma. Each relay stage is one readout window.

The quiet-phase principle derives from the three conditions for observer
connection (STAR_TOPOLOGY_OBSERVERS.md): the receiver must be quiet
(low gamma) for information to arrive. The relay protocol applies this
condition sequentially to each mediator.

The 2:1 coupling derives from the QST Bridge finding (QST_BRIDGE.md):
mediator-to-receiver coupling should be twice mediator-to-sender.

## What This Does Not Cover

- Tested only for Z-dephasing noise on a Heisenberg chain
- Tested only at N=11 (single chain length)
- Only one staging schedule tested (equal-length stages at 0.039/gamma)
- Optimal stage timing, gamma_quiet value, and number of stages not explored
- No comparison with existing quantum repeater protocols (e.g., entanglement
  swapping, purification-based repeaters)

## Implementation

The protocol requires time-dependent control of per-qubit dephasing rates.
In a physical implementation, this corresponds to dynamically varying the
shielding or isolation of individual qubits during the transfer protocol.

The C# implementation creates a new LindbladPropagator for each stage
(different gamma array) and propagates the density matrix through each
stage sequentially. The density matrix is passed between stages without
reinitialization.
