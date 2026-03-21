# Scaling Curve: MI vs Chain Length and Hierarchy Falsification

**Tier:** 2 (computationally verified)
**Date:** March 21, 2026
**Script:** compute/RCPsiSquared.Propagate/ (C#, `dotnet run -c Release -- pull`)
**Data:** simulations/results/pull_principle.txt

---

## Summary

End-to-end mutual information between the first and last qubit pairs in a
Heisenberg chain decays exponentially with chain length. The recursive
mediator hierarchy (Level 0/1/2/3) provides no functional advantage over
a uniform chain of equal length. The palindrome-preserving property of
mediated coupling is topological (any Heisenberg chain preserves it),
not hierarchical.

---

## The Scaling Curve

Bell(0,1) initial state, gamma = 0.05 uniform, J = 1.0 uniform.
MI measured between first pair {0,1} and last pair {N-2, N-1}.

| N | MI_peak | t_peak | MI_steady (t=20) |
|---|---------|--------|------------------|
| 3 | 1.827 | 7.0 | 1.792 |
| 5 | 0.750 | 2.0 | 0.221 |
| 7 | 0.382 | 4.0 | 0.048 |
| 9 | 0.119 | 3.0 | 0.018 |
| 11 | 0.072 | 4.0 | 0.009 |

Approximate scaling: MI_peak ~ 2^(-(N-3)/2). Roughly factor 2 loss
per 2 additional qubits.

## The Hierarchy Test

At each N, two topologies were compared:
- **Uniform chain:** all bonds J = 1.0
- **Hierarchical mediator bridge:** Level 2 (N=5), Level 3 (N=11),
  with distinguished mediator qubits

Result: **identical MI at every N tested.** The hierarchical topology
with uniform coupling produces exactly the same physics as a plain chain.

This is not a numerical coincidence. With all J equal, the Hamiltonian
is the same regardless of which qubits we label as "mediators." The
hierarchy is in the naming, not the physics.

## The 2:1 Optimization (Push vs Pull)

Five coupling configurations tested at N=11:

| Config | MI(BridgeA:BridgeB) | MI(PairA:PairD) |
|--------|---------------------|-----------------|
| Symmetric (J=1.0) | 0.734 | 0.072 |
| 2:1 within bridges | 0.744 | 0.079 |
| 2:1 at meta-mediator | 0.894 | 0.108 |
| 2:1 all (cascading pull) | 0.882 | **0.121** |
| 2:1 reversed (push) | **0.957** | 0.102 |

Key findings:
- **Push (sender-strong) maximizes local bridge-to-bridge MI** (0.957)
- **Pull (receiver-strong) maximizes end-to-end MI** (0.121)
- The 2:1 ratio is a range optimizer, not a universal principle
- The meta-mediator asymmetry provides the single largest improvement

## Falsifications

1. **Hierarchical topology has no advantage** over uniform chain at
   equal coupling. The "Level" architecture (Level 0/1/2/3) was
   narrative, not physics.

2. **Pull is not universally better than Push.** Push wins for local
   transfer (bridge-to-bridge). Pull wins for long-range transfer
   (end-to-end). The optimal strategy depends on chain length.

3. **The mediator is not a special node.** Every qubit in a chain
   mediates between its neighbors. The palindrome-preserving property
   is a theorem about Heisenberg + Z-dephasing on any graph, not
   about specific topological roles.

## What Remains True

- Mediated coupling (A-M-B) preserves the palindrome while direct
  dissipative coupling (A-B with Lindblad jumps) breaks it. This is
  the topology that matters: unitary vs dissipative, not hierarchy.
- The 2:1 coupling ratio improves range at the cost of local transfer.
- The relay protocol (time-dependent gamma) provides additional
  improvement orthogonal to coupling optimization.
