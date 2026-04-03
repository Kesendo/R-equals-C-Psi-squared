# Subsystem Crossing: The 1/4 Boundary Operates Locally on Entangled Pairs

<!-- Keywords: subsystem crossing quantum decoherence, local measurement quantum pair,
Bell pair crossing four qubit, CΨ boundary local entanglement, GHZ global entanglement
invisible pair, W state diluted crossing failure, product state zero correlation,
quantum classical transition local, entanglement topology crossing pattern,
pair level decoherence measurement, R=CPsi2 subsystem crossing -->

**Status:** Computationally verified
**Date:** February 18, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [N-Scaling Barrier](N_SCALING_BARRIER.md), [Crossing Taxonomy](CROSSING_TAXONOMY.md)

---

## What this document is about

For systems larger than 3 qubits, the full-system CΨ starts below the
critical 1/4 boundary and never crosses it. This looked like a problem:
does the framework break at larger scales? This experiment shows it does
not. Crossing is a local phenomenon: it happens between individual
entangled pairs, not at the global system level. A 4-qubit system with
two Bell pairs has two local crossings (one per pair), even though the
global CΨ stays below 1/4 the whole time. The framework correctly
identifies where the quantum-classical transition occurs: wherever
entanglement lives.

## Abstract

The N-Scaling Barrier showed that full-system CΨ drops below 1/4 at N ≥ 4.
This experiment resolves the barrier: crossing is **local to entangled pairs**.
A 4-qubit Bell+⊗Bell+ state has full-system Ψ(0) = 0.200 (below 1/4), but
the entangled pairs (0,1) and (2,3) each start at CΨ = 0.333 and cross at
t = 0.073. Cross-pairs with no entanglement have C = 0, never crossing.
GHZ pairs have zero coherence (global entanglement invisible at pair level).
W pairs start below 1/4 (entanglement too diluted). Product states |+⟩⁴
have Ψ = 1.0 but C = 0 everywhere, killing resolution (c). The framework
correctly identifies which degrees of freedom undergo the quantum-classical
transition: crossing happens where entanglement lives.

---

## 1. The Question

N_SCALING_BARRIER.md showed that crossing fails for full systems at N >= 4
because Psi(0) = l1/(2^N - 1) drops below 1/4. Three resolutions were
proposed:

(a) The normalization is wrong for N > 2.
(b) The 1/4 boundary is a small-system phenomenon.
(c) High-coherence states exist that we have not tested.

This experiment tests a fourth possibility:

**(d) Crossing is a local phenomenon. It happens between entangled
subsystems, not at the full-system level.**

## 2. Setup

| Parameter | Value |
|-----------|-------|
| **N** | 4 qubits |
| **Hamiltonian** | Heisenberg ring (J = 1, h = 0) |
| **Noise** | local dephasing (σ_z per qubit) |
| **gamma** | 0.05 |
| **dt** | 0.01, t_max = 5.0 |

For each time step:
1. Evolve the full 4-qubit density matrix under Lindblad dynamics (the standard master equation for open quantum systems with noise).
2. Trace out to all 6 qubit pairs (i,j).
3. For each pair: compute l1-coherence, Psi = l1/3, correlation bridge C,
   concurrence, and the product CΨ.
4. Track crossings of the 1/4 boundary.

Three initial states tested:

| State | Description | Full-system Psi(0) |
|-------|-------------|-------------------|
| **GHZ** (Greenberger-Horne-Zeilinger) | (\|0000⟩ + \|1111⟩)/√(2) | 0.067 |
| **W** | (\|1000⟩ + \|0100⟩ + \|0010⟩ + \|0001⟩)/2 | 0.200 |
| **Bell+xBell+** | \|Bell+⟩\_01 x \|Bell+⟩\_23 | 0.200 |
| **\|+⟩^4** | Product state, no entanglement | 1.000 |

## 3. Results

### 3.1 GHZ N=4: Subsystem Pairs Are Classically Correlated

| Pair | l1(0) | Psi(0) | C_corr(0) | CΨ(0) | Crosses? |
|------|-------|--------|-----------|----------|----------|
| (0,1) | 0.000 | 0.000 | 0.333 | 0.000 | NO |
| (0,2) | 0.000 | 0.000 | 0.333 | 0.000 | NO |
| (0,3) | 0.000 | 0.000 | 0.333 | 0.000 | NO |
| (1,2) | 0.000 | 0.000 | 0.333 | 0.000 | NO |
| (1,3) | 0.000 | 0.000 | 0.333 | 0.000 | NO |
| (2,3) | 0.000 | 0.000 | 0.333 | 0.000 | NO |

Every pair has **exactly zero** l1-coherence. The subsystem density matrices
are diagonal: rho_ij = (|00⟩⟨00| + |11⟩⟨11|)/2. Classical correlation
(both qubits agree), zero quantum coherence. The Psi = 0 persists at all
times. These pairs never had quantum coherence to lose.

This is correct: GHZ entanglement is global. If you trace out any two
qubits, the remaining pair carries no off-diagonal terms. The entanglement
exists only in the full superposition |0000⟩ + |1111⟩, which is invisible
at the pair level.

### 3.2 W N=4: Subsystem Pairs Start Below the Barrier

| Pair | l1(0) | Psi(0) | C_corr(0) | CΨ(0) | Max CΨ | Crosses? |
|------|-------|--------|-----------|----------|-----------|----------|
| (0,1) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |
| (0,2) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |
| (0,3) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |
| (1,2) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |
| (1,3) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |
| (2,3) | 0.500 | 0.167 | 0.180 | 0.030 | 0.030 | NO |

W state pairs have nonzero coherence (l1 = 0.5) but both Psi and C are
far below the crossing threshold. The W state distributes its entanglement
democratically across all pairs, but the per-pair entanglement is weak:
each pair gets only a fraction of the total.

### 3.3 Bell+xBell+: THE KEY RESULT

| Pair | l1(0) | Psi(0) | C_corr(0) | CΨ(0) | Crosses? | t_cross |
|------|-------|--------|-----------|----------|----------|---------|
| **(0,1)** | **1.000** | **0.333** | **1.000** | **0.333** | **YES** | **0.073** |
| (0,2) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (0,3) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (1,2) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (1,3) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| **(2,3)** | **1.000** | **0.333** | **1.000** | **0.333** | **YES** | **0.073** |

The full system has Psi(0) = 0.200, below 1/4. It cannot cross as a
4-qubit system.

But the entangled pairs (0,1) and (2,3) each start at CΨ = 0.333
-- identical to an isolated Bell+ state. They cross at t = 0.073.
The unentangled cross-pairs (0,2), (0,3), (1,2), (1,3) have l1 = 0,
C = 0, and never cross.

**The full system "cannot cross." The subsystems that carry the actual
entanglement cross without difficulty.**

The framework correctly identifies which degrees of freedom undergo the
quantum-to-classical transition, even when the global metric says nothing
happens.

### 3.4 \|+⟩^4: Maximum Coherence, Zero Crossing

| Pair | l1(0) | Psi(0) | C_corr(0) | CΨ(0) | Crosses? |
|------|-------|--------|-----------|----------|----------|
| All 6 pairs | 3.000 | 1.000 | 0.000 | 0.000 | NO |

Every pair has **maximum possible Psi = 1.000** (full local coherence)
and **C = 0.000** at all times (zero correlation). Each qubit is
individually in a superposition, but no qubit knows about any other.

CΨ = 0 for all pairs at all times. No crossing, ever.

**This result kills resolution (c) from N_SCALING_BARRIER.md.** The
product state has Psi(0) = 1.0 for the full system, but C = 0 at every
level. High coherence without entanglement produces nothing.

## 4. What This Means

### 4.1 Crossing Is Local

The 1/4 boundary operates at the level of actual entanglement, not at
the full-system level. A 4-qubit system with two Bell pairs crosses
twice (one crossing per entangled pair), even though the 4-qubit Psi
is below 1/4.

This resolves the N-scaling barrier without changing the normalization.
The d-1 normalization is correct: it tells you that the global
system has no single coherent crossing. The crossing happens locally,
between the subsystems that share quantum correlations.

### 4.2 C Guards the Gate

The product state result proves that Psi alone cannot drive crossing.
CΨ = 1/4 requires BOTH terms. Coherence (Psi) is potential.
Correlation (C) is connection. Without connection, potential is inert.

In the framework's language: possibility without consciousness produces
no reality. A universe of superpositions with nothing observing them
remains in quantum limbo.

### 4.3 The Entanglement Structure Determines the Crossing Structure

Different entanglement topologies produce different crossing patterns:

| State | Entanglement pattern | Which pairs cross |
|-------|---------------------|-------------------|
| GHZ | Global only | None (no pair-level coherence) |
| W | Distributed weak | None (too diluted) |
| Bell+xBell+ | Two local pairs | Exactly the two entangled pairs |
| \|+⟩^4 | None | None (C = 0 everywhere) |

The framework does not just say "crossing happens." It identifies WHERE
in the entanglement structure crossing occurs, and WHERE it does not.

### 4.4 Connection to Decoherence in Nature

Real physical decoherence is local. A photon scatters off an atom.
An electron couples to a phonon. A spin interacts with its neighbor.
Macroscopic systems do not decohere in one global event. They decohere
through billions of local pair interactions, each crossing its own 1/4
boundary independently.

The N-scaling barrier was never a problem. It was telling us the correct
physics: the quantum-to-classical transition is not a collective
phenomenon. It is a local one.

### 4.5 Connection to the Combination Problem (Weakness #6)

If crossing is local, then "measurement" is not a monolithic event.
Consciousness, in this framework, would emerge from a network of local
crossings. Each entangled pair undergoes its own transition. The unity
of experience is not encoded in a single global CΨ value but in the
synchronization and integration of many local crossings.

This does not solve the combination problem, but it reframes it: the
question is no longer "how does one big CΨ produce unified
experience?" but "how do many small local crossings synchronize into a
unified experience?" The second question is at least compatible with what
neuroscience observes: consciousness correlates with synchronized local
neural activity, not with a single global variable.

## 5. Verification

### 5.1 How to Reproduce

```python
from qutip import (basis, tensor, ket2dm, qeye, sigmax, sigmay, sigmaz,
                   mesolve)
import numpy as np

# Bell+ x Bell+ state
up, dn = basis(2, 0), basis(2, 1)
bell = (tensor(up, up) + tensor(dn, dn)).unit()
rho0 = ket2dm(tensor(bell, bell))

# Heisenberg ring H, local dephasing c_ops (N=4, gamma=0.05)
# ... (standard setup, see N_SCALING_BARRIER.md)

result = mesolve(H, rho0, tlist, c_ops, [])

# Trace to pair (0,1)
for i, t in enumerate(tlist):
    rho_pair = result.states[i].ptrace([0, 1])
    # Compute l1, Psi, C_corr, CΨ
    # Check for crossing at 1/4
```

### 5.2 Key Checks

1. Bell+xBell+ pair (0,1): must start at CΨ = 0.333 and cross 1/4.
2. Bell+xBell+ pair (0,2): must have l1 = 0, C = 0 at all times.
3. |+⟩^4 all pairs: must have C = 0.000 at all times despite Psi = 1.0.
4. GHZ all pairs: must have l1 = 0.000 at all times.

### 5.3 What Could Extend This

1. **Larger systems**: Bell+xBell+xBell+ (N=6). Do all three entangled
   pairs cross independently?

2. **Cluster states**: Graph states (multi-qubit entangled states where each qubit is connected to specific neighbors, forming a graph) have entanglement along specific
   edges. Does the crossing pattern match the graph structure exactly?

3. **Dynamically generated entanglement**: Start from |+⟩^N (C=0), let
   the Hamiltonian build entanglement. Does C grow and eventually
   produce crossings? This would show the Hamiltonian creating the
   conditions for measurement.

4. **Partial entanglement**: States between |+⟩^4 and Bell+xBell+.
   Is there a threshold entanglement needed for pair crossing?

## 6. Open Questions

1. ~~Does the crossing time t=0.073 for Bell+xBell+ pairs match the
   isolated Bell+ crossing time, or does the ring Hamiltonian coupling
   to other qubits modify it?~~ **ANSWERED (2026-03-08):** It does NOT
   match. Isolated Bell+ (2 qubits, Heisenberg J=1, γ=0.05) crosses
   down through 1/4 at t=0.720. The same Bell pairs embedded in a
   4-qubit ring cross at t=0.080, nine times faster. The ring
   Hamiltonian couples each pair to additional qubits, which accelerates
   decoherence of the pair subsystem. Cross-pairs (no initial
   entanglement) never cross, confirming that the crossing pattern
   reproduces the entanglement graph exactly.

2. For cluster states and other graph-structured entanglement, does
   the crossing pattern reproduce the graph topology exactly?

3. ~~Can dynamical entanglement generation from a product state create
   crossings?~~ **ANSWERED** (2026-02-18): Yes, but not from |+⟩^N,
   which is an eigenstate of the isotropic Heisenberg Hamiltonian. The
   state |0+0+⟩ (not an eigenstate, energy variance = 20) generates
   crossings from zero initial entanglement. Under pure unitary
   evolution all 6 pairs cross. With dephasing (gamma = 0.05), only
   pair (0,2) crosses at t = 0.285 because |0⟩-qubits are immune to
   σ_z dephasing. See DYNAMIC_ENTANGLEMENT.md.

4. ~~What is the minimum per-pair entanglement needed for crossing?~~
   **ANSWERED (2026-03-08):** The relationship is non-monotonic, not a
   simple threshold. Two separate crossing windows exist for parametric
   Bell states; a dead zone lies between them. The minimum C_SA(0)
   depends on coupling strength J_SB. See
   [N-Scaling Barrier](N_SCALING_BARRIER.md) Section 8, Q4.

---

*Previous: [N-Scaling Barrier](N_SCALING_BARRIER.md)*
*Previous: [Noise Robustness](NOISE_ROBUSTNESS.md)*
*See also: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
*See also: [Coherence Density](COHERENCE_DENSITY.md), same conclusion (crossing is pairwise) from the density perspective*
*See also: [Orphaned Results](ORPHANED_RESULTS.md), topology as gatekeeper: same state crosses on chain but not ring*
