<!-- QUARTER-CURRENT -->
# Subsystem quarter readouts in named finite models

Current reading: the stored N=4 calculations compare a finite readout for
specified subsystems, states, and generators. The scalar quarter does not
decide Wootters concurrence or identify a measurement. For the analytic
Bell+ local-Z-dephasing control, `C=f>0`: concurrence remains positive at the scalar
quarter and at every finite later time, including scalar values below `1/4`.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the February narrative and tables below are retained as
the route by which the finite subsystem question was found; their boundary
labels are not current claims.

# Subsystem Crossing: finite pair-level readouts in N=4 examples

<!-- Keywords: subsystem crossing quantum decoherence, local measurement quantum pair,
Bell pair crossing four qubit, CΨ boundary local entanglement, GHZ global entanglement
invisible pair, W state diluted crossing failure, product state zero correlation,
quantum classical transition local, entanglement topology crossing pattern,
pair level decoherence measurement, R=CPsi2 subsystem crossing -->

**Status:** Finite catalogue; boundary interpretation withdrawn
**Date:** February 18, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [N-Scaling Barrier](N_SCALING_BARRIER.md), [Crossing Taxonomy](CROSSING_TAXONOMY.md)

---

## What this document is about

For the specified N=4 examples, a full-system scalar and pair-reduced scalars
can sit on different sides of 1/4. A product of two Bell pairs supplies two
pair-level crossings while its chosen full-system readout stays below 1/4.
This is a useful warning about subsystem dependence; it does not prove that
crossing is universally local or locate a quantum/classical transition.

## Abstract

The earlier N-scaling note observed a size-dependent full-system normalization.
This experiment compares that book with **pair-reduced finite readouts**.
A 4-qubit Bell+⊗Bell+ state has full-system Ψ(0) = 0.200 (below 1/4), but
the entangled pairs (0,1) and (2,3) each start at CΨ = 0.333 and cross at
t = 0.080. Cross-pairs with no entanglement have C = 0, never crossing.
GHZ pairs have zero coherence (global entanglement invisible at pair level).
W pairs start below 1/4 (entanglement too diluted). Product states |+⟩⁴
have Ψ = 1.0 but C = 0 in this construction. These rows distinguish several
state families; they do not make the scalar crossing an entanglement theorem.

---

## 1. The Question

N_SCALING_BARRIER.md showed that crossing fails for full systems at N >= 4
because Psi(0) = l1/(2^N - 1) drops below 1/4. Three resolutions were
proposed:

(a) The normalization is wrong for N > 2.
(b) The 1/4 boundary is a small-system phenomenon.
(c) High-coherence states exist that we have not tested.

This experiment tests a fourth possibility:

**(d) In these examples, do pair-reduced readouts cross when the chosen
full-system readout does not?**

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
| (0,1) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |
| (0,2) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |
| (0,3) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |
| (1,2) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |
| (1,3) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |
| (2,3) | 0.500 | 0.167 | 0.417 | 0.083 | 0.083 | NO |

(C_corr is the connected-correlator bridge, 0.417 for a W pair; the
concurrence is 0.5, and CΨ = concurrence·Ψ = 0.083. The February run
read the same pair with the pairwise bridge
(P_AB − P_A·P_B)/(1 − P_A·P_B) = 7/39 ≈ 0.18, which gives
CΨ = 7/234 ≈ 0.030
([delta_calc_pairwise_bridge.py](../simulations/delta_calc_pairwise_bridge.py));
the verdict holds in every one of these books.) W state pairs have
nonzero coherence (l1 = 0.5), but the product CΨ = 0.083 stays far below
the crossing threshold. The W
state distributes its entanglement democratically across all pairs, but
the per-pair entanglement is weak: each pair gets only a fraction of the
total.

### 3.3 Bell+xBell+: THE KEY RESULT

| Pair | l1(0) | Psi(0) | C_corr(0) | CΨ(0) | Crosses? | t_cross |
|------|-------|--------|-----------|----------|----------|---------|
| **(0,1)** | **1.000** | **0.333** | **1.000** | **0.333** | **YES** | **0.080** |
| (0,2) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (0,3) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (1,2) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| (1,3) | 0.000 | 0.000 | 0.000 | 0.000 | NO | n/a |
| **(2,3)** | **1.000** | **0.333** | **1.000** | **0.333** | **YES** | **0.080** |

The full system has Psi(0) = 0.200, below 1/4. It cannot cross as a
4-qubit system.

But the entangled pairs (0,1) and (2,3) each start at CΨ = 0.333,
identical to an isolated Bell+ state. They cross at t = 0.080 in the
concurrence book, as the reproduction script below and Open Question 1's
own nine-times ratio both give; the pairwise bridge the February run was
read with crosses the same pair at t = 0.073 under exact propagation
([delta_calc_pairwise_bridge.py](../simulations/delta_calc_pairwise_bridge.py)).
The unentangled cross-pairs (0,2), (0,3), (1,2), (1,3) have l1 = 0,
C = 0, and never cross.

**In this finite book, the selected full-system trace does not cross while the
two initially entangled pair readouts do.**

That contrast is a statement about which readout was computed. It does not say
that those degrees of freedom undergo a physical regime transition.

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

### 4.1 Pair and full-system readouts differ

The Bell+⊗Bell+ example reports two pair-readout crossings while the chosen
full-system scalar remains below 1/4. The result is finite and depends on the
state, reduction, normalization, and generator listed above.

It shows how changing the subsystem changes the scalar book without changing
the normalization. It does not prove that every entangled subsystem crosses:
the W-family rows and the analytic Bell+ concurrence control already separate
entanglement from this threshold.

### 4.2 C Guards the Gate

The product-state row is the direct arithmetic statement C=0, hence CΨ=0.
Calling Ψ "potential" and C "connection" is the historical interpretation,
not an additional theorem.

**Interpretive invitation:** the older notebook called this "possibility
without consciousness." No observer or consciousness model occurs in the
calculation.

### 4.3 The finite state catalogue gives different crossing patterns

Different entanglement topologies produce different crossing patterns:

| State | Entanglement pattern | Which pairs cross |
|-------|---------------------|-------------------|
| GHZ | Global only | None (no pair-level coherence) |
| W | Distributed weak | None (too diluted) |
| Bell+xBell+ | Two local pairs | Exactly the two entangled pairs |
| \|+⟩^4 | None | None (C = 0 everywhere) |

Within these four rows, the pair-reduced readout differs by state family. This
catalogue is not an all-state equivalence between entanglement and crossing.

### 4.4 Why the local picture invited a physical story

Physical decoherence can arise through local couplings, but those examples do
not imply that each interaction crosses this readout or that the finite N=4
calculation models a macroscopic environment.

The useful current conclusion is methodological: state a subsystem and a
normalization before comparing scalar values. Whether classicality is local or
collective is outside this computation.

### 4.5 Connection to the Combination Problem (Weakness #6)

**Historical interpretation:** the original note asked whether a network of
local crossings might offer language for measurement or consciousness. That
question remains an invitation; the finite pair catalogue neither supplies the
network dynamics nor shows that any pair undergoes such a transition.

This does not solve the combination problem, but it reframes it: the
question is no longer "how does one big CΨ produce unified
experience?" but "how do many small local crossings synchronize into a
unified experience?" The second question is at least compatible with what
neuroscience observes: consciousness correlates with synchronized local
neural activity, not with a single global variable.

## 5. Verification

### 5.1 How to Reproduce

The four result tables were first computed in February 2026 with QuTiP
`mesolve`, as sketched below, and read with the pairwise bridge
(P_AB − P_A·P_B)/(1 − P_A·P_B)
([delta_calc_pairwise_bridge.py](../simulations/delta_calc_pairwise_bridge.py)).
The committed reproduction in the concurrence book is
[`simulations/subsystem_crossing_pairs.py`](../simulations/subsystem_crossing_pairs.py),
which rebuilds all four tables (and the isolated-Bell+ baseline of Open
Question 1) from scratch. Do not confuse it with
[`simulations/subsystem_crossing.py`](../simulations/subsystem_crossing.py),
which is the *proof's* verifier (a different setup: Bell+ ⊗ |0...0⟩ and
random CPTP maps). Sketch of the original setup:

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
   produce crossings in the named scalar readout? This is a finite dynamics
   question, not a measurement criterion.

4. **Partial entanglement**: States between |+⟩^4 and Bell+xBell+.
   Is there a threshold entanglement needed for pair crossing?

## 6. Open Questions

1. ~~Does the crossing time t=0.080 for Bell+xBell+ pairs match the
   isolated Bell+ crossing time, or does the ring Hamiltonian coupling
   to other qubits modify it?~~ **ANSWERED (2026-03-08):** It does NOT
   match. Isolated Bell+ (2 qubits, Heisenberg J=1, γ=0.05) crosses
   down through 1/4 at t=0.720. The same Bell pairs embedded in a
   4-qubit ring cross at t=0.080, nine times faster. The ring
   Hamiltonian couples each pair to additional qubits and changes the
   subsystem trajectory. In this one state, ring, channel, and sampled time
   window, cross-pairs with zero initial concurrence did not cross the
   selected readout level. That finite negative row does not establish a
   graph-reconstruction law. (A note on the isolated
   baseline, because three numbers circulate for it: 0.720 is the
   crossing of THIS document's CΨ = concurrence·Ψ = f²/3, with
   f = e^(−4γt). The F25 purity-bridge book CΨ = f(1+f²)/6 crosses at
   0.747 = K_fold/γ, and the constant-bridge book CΨ = f/3 (the
   taxonomy's Type-A "correlation bridge", C held at 1) crosses at
   1.438, which is [Crossing Taxonomy](CROSSING_TAXONOMY.md)'s 1.437.
   One trajectory, three C-readings; see
   [The CΨ Lens](../docs/THE_CPSI_LENS.md) for the reading conventions.)

2. For cluster states and other graph-structured entanglement, does
   the crossing pattern reproduce the graph topology exactly?

3. ~~Can dynamical entanglement generation from a product state create
   crossings?~~ **ANSWERED** (2026-02-18): Yes, but not from |+⟩^N,
   which is an eigenstate of the isotropic Heisenberg Hamiltonian. The
   state |0+0+⟩ (not an eigenstate, energy variance = 20) generates
   crossings from zero initial entanglement. Under pure unitary
   evolution all 6 pairs cross. With dephasing (gamma = 0.05), only
   pair (0,2) crosses at t = 0.285 because |0⟩-qubits are immune to
   σ_z dephasing. See DYNAMIC_ENTANGLEMENT.md. *(The ring-(0,2)
   crossing numbers in that answer are read in the pairwise bridge, in
   which they reproduce under exact propagation; a σ_x run keeps (0,2) as
   the only crossing pair, so the immunity reason does not hold. Under the
   canonical pair-CΨ book that pair does not cross, and at this γ the
   upward crossing sits on the chain and in |+-+-⟩/|0+0-⟩ on the ring; see
   the reproduction note in DYNAMIC_ENTANGLEMENT.md.)*

4. ~~What is the minimum per-pair entanglement needed for crossing?~~
   **ANSWERED (2026-03-08):** The relationship is non-monotonic, not a
   simple threshold. Two separate crossing windows exist for parametric
   Bell states; a dead zone lies between them. The minimum C_SA(0)
   depends on coupling strength J_SB. See
   [N-Scaling Barrier](N_SCALING_BARRIER.md) Section 8, Q4.

---

## Where this went (successors)

- [Subsystem Crossing Proof](../docs/proofs/PROOF_SUBSYSTEM_CROSSING.md):
  this experiment's observation, formalized (March 2026). Its general
  "any primitive CPTP map" version was scope-retracted 2026-06-22 (Case
  C: a primitive channel with an entangled fixed point, CΨ = 0.2935,
  never crosses); what survives, proven, is exactly the scope this
  experiment measured: physical, computational-basis-aligned noise.
- [F28 in the formula registry](../docs/ANALYTICAL_FORMULAS.md#f28): the
  fixed-point absorber theorem carrying the same scoped verdict.

*Previous: [N-Scaling Barrier](N_SCALING_BARRIER.md)*
*Previous: [Noise Robustness](NOISE_ROBUSTNESS.md)*
*See also: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
*See also: [Coherence Density](COHERENCE_DENSITY.md), same conclusion (crossing is pairwise) from the density perspective*
*See also: [Orphaned Results](ORPHANED_RESULTS.md), topology as gatekeeper: at γ = 0.05 the same state crosses on the chain but not on the ring (concurrence book)*
