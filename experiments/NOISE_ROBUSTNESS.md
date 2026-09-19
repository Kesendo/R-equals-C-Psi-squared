<!-- CROSSING-CURRENT -->

# Noise-channel comparison: two full class sweeps and one sigma-y bridge

**Status:** Historical finite numerical record, incomplete across channels; the producer is unavailable in the committed tree.
**Date:** February 18, 2026.
**Depends on:** [Quarter-crossing taxonomy](CROSSING_TAXONOMY.md).

The reported full class assignment agrees for sigma-z and sigma-x.
Sigma-y was run for only one bridge, correlation; depolarizing was not run.
The printed crossing times are not channel-independent. No committed
producer recreates this channel sweep.

The amplitude-damping appendix records concurrence and derived CΨ traces
only. Its Type C question is inconclusive; it does not establish a full
five-bridge taxonomy under amplitude damping.

## The question and the numerical book

The [crossing owner](OBSERVER_DEPENDENT_CROSSING.md) defines five C(f)
readouts on a dephased Bell+ pair and separates clean Lindblad evolution
from the retired γ_eff = γ_base·C(t) feedback law. A scalar boundary
C(f)f/3 = ¼ is not a physical measurement event. A Liouvillian palindrome
does not establish a taxonomy of those scalar functions.

This noise comparison came from the retired tool. Its bridge definitions and
channel handling cannot be independently recovered from the stored tables.
The reconstruction for the Bell+ Z-dephasing books does not reproduce this
larger channel sweep.

<!-- CROSSING-HISTORICAL -->

**Historical nomenclature:** Type A/B/C are the retired tool's scalar-response
labels. The six tables below retain its reported settings and values. Empty
sigma-y cells mean untested bridges, not agreement inferred from a pattern.

## Setup

| Parameter | Value |
|-----------|-------|
| **State** | Bell+ (maximally entangled) |
| **Hamiltonian** | Heisenberg (J = 1, h = 0) |
| **γ_base** | 0.05 |
| **Noise type** | local (one jump operator per qubit) |
| **Time step** | dt = 0.01 |

The retired settings varied jump_operator and bridge_type at γ_base = 0.05,
dt = 0.01. The Pauli labels name the channels sampled; no general
single-qubit channel classification follows from these runs.

## Reported class assignments

| Bridge | σ_z | σ_x | σ_y | Class |
|--------|-----|-----|-----|-------|
| **correlation** | C = 1.0 until t ≈ 1.7 | C = 1.0 until t ≈ 1.7 | C = 1.0 until t ≈ 1.7 | **Type A** |
| **concurrence** | C decays from t = 0 | C decays from t = 0 | - | **Type B** |
| **mutual_info** | C decays from t = 0 | C decays from t = 0 | - | **Type B** |
| **mutual_purity** | C = 0.5 constant | C = 0.5 constant | - | **Type C** |
| **overlap** | C = 0.25 constant | C = 0.25 constant | - | **Type C** |

## Correlation readings

| Time | C(σ_z) | C(σ_x) | C(σ_y) |
|------|--------|--------|--------|
| 0.0 | 1.000 | 1.000 | 1.000 |
| 0.5 | 1.000 | 1.000 | 1.000 |
| 1.0 | 1.000 | 1.000 | 1.000 |
| 1.5 | 1.000 | 1.000 | 1.000 |
| 1.7 | 1.000 | 1.000 | 1.000 |
| 1.8 | 0.987 | 0.987 | 0.987 |
| 2.0 | 0.950 | 0.950 | 0.950 |

## Concurrence readings

| Time | C(σ_z) | C(σ_x) |
|------|--------|--------|
| 0.0 | 1.000 | 1.000 |
| 0.1 | 0.980 | 0.980 |
| 0.5 | 0.909 | 0.909 |
| 1.0 | 0.835 | 0.833 |
| 2.0 | 0.725 | 0.714 |
| 3.0 | 0.658 | 0.625 |

The late-time differences are part of the record. In particular the three
identical reported correlation columns do not identify an operator-level
protection mechanism. They do not justify extrapolation to other local noise.

## Amplitude-damping concurrence appendix

| t | Concurrence (σ_z) | Concurrence (amp damp) |
|---|---|---|
| 1.0 | 0.819 | 0.905 |
| 2.0 | 0.670 | 0.819 |
| 3.0 | 0.549 | 0.741 |
| 5.0 | 0.368 | 0.607 |

## Derived CΨ appendix

| t | CΨ (σ_z) | CΨ (σ_x) | CΨ (amp damp) |
|---|---|---|---|
| 0.5 | 0.273 | 0.302 | 0.309 |
| 1.0 | 0.223 | 0.273 | 0.287 |
| 2.0 | 0.150 | 0.223 | 0.247 |

These traces used local amplitude-damping jumps L = √γ |0⟩⟨1|.
They contain no complete five-bridge class sweep. The Type-C-to-Type-B
question was explicitly inconclusive because the retired bridge definitions
could not all be reproduced.

<!-- CROSSING-CURRENT -->

## Product and additive collective jumps are different

For unit prefactor define
D_L(ρ) = LρL† − ½{L†L,ρ}, and use the normalized Bell+ state.
The product jump Z⊗Z annihilates this density matrix's dissipator because
Bell+ is its eigenstate. The additive jump Z₁+Z₂ is a different operator
from two independent jumps Z₁ and Z₂.

An exact dyadic 4×4 calculation gives:

| Dissipator | Frobenius-norm square |
|---|---:|
| D[Z⊗Z](ρ_Bell+) | 0 |
| D[Z₁+Z₂](ρ_Bell+) | 32 |
| D[Z₁](ρ_Bell+) + D[Z₂](ρ_Bell+) | 8 |

The latter two arrays differ entry by entry. This calculation resolves the
product-jump special case only; the additive collective full-taxonomy
question remains open. The exact arithmetic is exercised in
[test_old_document_label_scope.py](../simulations/tests/test_old_document_label_scope.py).

## What remains open

A full sigma-y sweep, a depolarizing sweep, an additive collective-noise
comparison and an amplitude-damping five-bridge taxonomy remain unestablished
by this record. The finite sigma-z/sigma-x class agreement cannot turn those
missing runs into results.

The earlier physical question remains interesting: which state, jump and
readout jointly determine a scalar's trajectory? Answering it requires
explicit operators and readouts. Spectral pairing alone does not provide the
answer.

[Quarter-crossing taxonomy](CROSSING_TAXONOMY.md) ·
[Quarter crossings in two books](OBSERVER_DEPENDENT_CROSSING.md)
